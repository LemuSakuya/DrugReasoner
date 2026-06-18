"""
End-to-end analysis pipeline for the innovation project.

This module turns a natural-language drug analysis query into:
1. a normalized task object,
2. local and optional external evidence,
3. a signed relation subgraph,
4. a layered explanation report,
5. exportable JSON/Markdown artifacts.
"""

from __future__ import annotations

import json
import os
import re
import time
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional, Tuple

from config import Config
from search_agent import (
    DrugResult,
    InteractionResult,
    ProteinResult,
    QueryParser,
    SearchAgent,
)


@dataclass
class TaskObject:
    """Structured task shared by retrieval, prediction, and explanation steps."""

    raw_query: str
    intent: str
    entities: List[str] = field(default_factory=list)
    task_type: str = "general"
    slots: Dict[str, Any] = field(default_factory=dict)
    requires_prediction: bool = False
    requires_external_search: bool = False


@dataclass
class EvidenceRecord:
    """Traceable evidence item used by report generation."""

    source: str
    record_id: str
    statement: str
    relation_type: str = "entity"
    polarity: int = 0
    confidence: float = 0.5
    url: str = ""


@dataclass
class SignedEdge:
    """Signed subgraph edge."""

    source: str
    target: str
    relation_type: str
    polarity: int = 0
    evidence_ids: List[str] = field(default_factory=list)
    weight: float = 0.5


@dataclass
class AnalysisReport:
    """Layered report expected by the project application."""

    task: TaskObject
    answer: str
    reasoning: List[str]
    evidence: List[EvidenceRecord]
    signed_subgraph: List[SignedEdge]
    generated_content: str
    exports: Dict[str, str] = field(default_factory=dict)
    created_at: str = field(default_factory=lambda: time.strftime("%Y-%m-%d %H:%M:%S"))


POSITIVE_TERMS = (
    "协同", "激活", "增强", "治疗", "适应症", "结合", "可相互作用", "agonist",
    "activate", "increase", "enhance", "treat", "bind",
)
NEGATIVE_TERMS = (
    "拮抗", "抑制", "禁忌", "副作用", "风险", "毒性", "无明显", "antagonist",
    "inhibit", "decrease", "adverse", "toxicity", "contraindication",
)


def _infer_polarity(text: str) -> int:
    value = str(text).lower()
    if any(term.lower() in value for term in NEGATIVE_TERMS):
        return -1
    if any(term.lower() in value for term in POSITIVE_TERMS):
        return 1
    return 0


def _task_type_from_intent(intent: str) -> str:
    return {
        "drug_drug": "DDI",
        "drug_protein": "DTA/DTI",
        "protein_protein": "PPI",
        "protein_info": "protein_info",
        "drug_info": "drug_info",
        "general_pair": "relation",
    }.get(intent, "general")


def parse_task(query: str, force_external: Optional[bool] = None) -> TaskObject:
    """Parse natural language into a stable task object."""
    intent, params = QueryParser.parse(query)
    entities = [
        params.get("entity1", ""),
        params.get("entity2", ""),
    ]
    entities = [e.strip() for e in entities if e and str(e).strip()]
    text = query.lower()
    external_requested = any(k in query for k in ("文献", "外部", "PubMed", "证据", "溯源"))
    if force_external is None:
        force_external = os.getenv("ENABLE_EXTERNAL_SEARCH", "").strip() in ("1", "true", "yes")
    return TaskObject(
        raw_query=query,
        intent=intent,
        entities=entities,
        task_type=_task_type_from_intent(intent),
        slots={
            "entity1": params.get("entity1", ""),
            "entity2": params.get("entity2", ""),
            "risk_focus": any(k in query for k in ("风险", "慎用", "禁忌", "副作用")),
            "affinity_focus": any(k in text for k in ("dta", "affinity", "binding")) or "亲和力" in query,
        },
        requires_prediction=intent in ("drug_drug", "drug_protein", "protein_protein", "general_pair"),
        requires_external_search=bool(force_external or external_requested),
    )


def _resolve_pair(agent: SearchAgent, entity1: str, entity2: str) -> Tuple[str, str]:
    """Resolve fuzzy pair names to canonical local names when possible."""
    left, right = entity1, entity2
    left_drugs = agent.search_drug(entity1, top_k=1)
    right_drugs = agent.search_drug(entity2, top_k=1)
    left_prots = agent.search_protein(entity1, top_k=1)
    right_prots = agent.search_protein(entity2, top_k=1)

    if left_drugs and right_drugs:
        return left_drugs[0].name, right_drugs[0].name
    if left_drugs and right_prots:
        return left_drugs[0].name, right_prots[0].name
    if left_prots and right_drugs:
        return left_prots[0].name, right_drugs[0].name
    if left_prots and right_prots:
        return left_prots[0].name, right_prots[0].name
    return left, right


def collect_local_evidence(agent: SearchAgent, task: TaskObject, top_k: int = 5) -> List[EvidenceRecord]:
    """Collect database/model evidence from the local SearchAgent."""
    evidence: List[EvidenceRecord] = []
    response = agent.search(task.raw_query, top_k=top_k)

    for idx, result in enumerate(response.results[:top_k], 1):
        record_id = f"local-{idx}"
        if isinstance(result, InteractionResult):
            statement = f"{result.entity1} 与 {result.entity2} 的 {result.interaction_type} 结果：{result.result}"
            evidence.append(EvidenceRecord(
                source=result.source or "local",
                record_id=record_id,
                statement=statement,
                relation_type=result.interaction_type,
                polarity=_infer_polarity(result.result),
                confidence=result.confidence if result.confidence is not None else 0.75,
            ))
        elif isinstance(result, DrugResult):
            details = []
            if result.side_effects:
                details.append("副作用：" + "、".join(result.side_effects[:5]))
            if result.chemical_features:
                details.append("化学特征：" + "、".join(result.chemical_features[:5]))
            statement = f"药物匹配：{result.name} (ID: {result.drug_id})"
            if details:
                statement += "；" + "；".join(details)
            evidence.append(EvidenceRecord(
                source="local_db",
                record_id=record_id,
                statement=statement,
                relation_type="drug_info",
                polarity=_infer_polarity(statement),
                confidence=max(0.5, min(1.0, result.match_score)),
            ))
        elif isinstance(result, ProteinResult):
            aliases = "、".join(result.aliases[:5]) if result.aliases else "无"
            statement = f"蛋白匹配：{result.name}；基因名：{result.gene_name or '未知'}；别名：{aliases}"
            evidence.append(EvidenceRecord(
                source="local_db",
                record_id=record_id,
                statement=statement,
                relation_type="protein_info",
                polarity=0,
                confidence=max(0.5, min(1.0, result.match_score)),
            ))

    if len(task.entities) >= 2 and not any(e.relation_type in ("DDI", "DTI", "DTA", "PPI") for e in evidence):
        left, right = _resolve_pair(agent, task.entities[0], task.entities[1])
        result = agent.search_interaction(left, right)
        if result:
            evidence.append(EvidenceRecord(
                source=result.source or "local",
                record_id=f"local-{len(evidence) + 1}",
                statement=f"{result.entity1} 与 {result.entity2} 的 {result.interaction_type} 结果：{result.result}",
                relation_type=result.interaction_type,
                polarity=_infer_polarity(result.result),
                confidence=result.confidence if result.confidence is not None else 0.72,
            ))

    if not evidence and response.summary:
        evidence.append(EvidenceRecord(
            source="local_db",
            record_id="local-summary",
            statement=response.summary,
            relation_type=task.task_type,
            polarity=_infer_polarity(response.summary),
            confidence=0.45,
        ))
    return evidence


def fetch_pubmed_evidence(query: str, retmax: int = 3, timeout: int = 8) -> List[EvidenceRecord]:
    """Fetch lightweight PubMed evidence through NCBI E-utilities."""
    term = urllib.parse.quote(query)
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    search_url = f"{base}/esearch.fcgi?db=pubmed&retmode=json&retmax={retmax}&term={term}"
    try:
        with urllib.request.urlopen(search_url, timeout=timeout) as resp:
            payload = json.loads(resp.read().decode("utf-8"))
        ids = payload.get("esearchresult", {}).get("idlist", [])
        if not ids:
            return []
        summary_url = f"{base}/esummary.fcgi?db=pubmed&retmode=json&id={','.join(ids)}"
        with urllib.request.urlopen(summary_url, timeout=timeout) as resp:
            summary_payload = json.loads(resp.read().decode("utf-8"))
    except Exception as exc:
        return [EvidenceRecord(
            source="pubmed",
            record_id="pubmed-error",
            statement=f"外部文献检索失败：{exc}",
            relation_type="external_search",
            polarity=0,
            confidence=0.1,
        )]

    records: List[EvidenceRecord] = []
    result_map = summary_payload.get("result", {})
    for pmid in ids:
        item = result_map.get(pmid, {})
        title = re.sub(r"\s+", " ", item.get("title", "")).strip()
        journal = item.get("fulljournalname", "") or item.get("source", "")
        year = str(item.get("pubdate", ""))[:4]
        statement = title
        if journal:
            statement += f" ({journal}, {year})"
        records.append(EvidenceRecord(
            source="pubmed",
            record_id=f"PMID:{pmid}",
            statement=statement,
            relation_type="literature",
            polarity=_infer_polarity(statement),
            confidence=0.55,
            url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
        ))
    return records


def build_signed_subgraph(task: TaskObject, evidence: List[EvidenceRecord]) -> List[SignedEdge]:
    """Convert evidence records into a small signed subgraph."""
    edges: List[SignedEdge] = []
    entities = task.entities or [task.raw_query]
    if len(entities) >= 2:
        left, right = entities[0], entities[1]
    else:
        left, right = entities[0], task.task_type

    for ev in evidence:
        relation = ev.relation_type
        if relation == "drug_info" and "副作用：" in ev.statement:
            target = "side_effect"
        elif relation == "drug_info" and "化学特征：" in ev.statement:
            target = "chemical_feature"
        elif relation == "protein_info":
            target = "protein_alias"
        else:
            target = right
        edges.append(SignedEdge(
            source=left,
            target=target,
            relation_type=relation,
            polarity=ev.polarity,
            evidence_ids=[ev.record_id],
            weight=ev.confidence,
        ))
    return edges


def _polarity_text(polarity: int) -> str:
    if polarity > 0:
        return "正向/协同/激活倾向"
    if polarity < 0:
        return "负向/拮抗/抑制或风险倾向"
    return "中性或证据不足"


def generate_report(task: TaskObject, evidence: List[EvidenceRecord], subgraph: List[SignedEdge]) -> AnalysisReport:
    """Generate a deterministic layered explanation report."""
    if not evidence:
        answer = "本地知识库未命中可用证据，建议开启外部文献检索或补充数据源。"
    else:
        best = max(evidence, key=lambda e: e.confidence)
        answer = f"{task.task_type} 分析结论：{best.statement}（{_polarity_text(best.polarity)}）。"

    reasoning = [
        f"任务解析为 {task.intent}，标准任务类型为 {task.task_type}。",
        f"共聚合 {len(evidence)} 条证据，其中本地证据 {sum(1 for e in evidence if e.source != 'pubmed')} 条，外部文献证据 {sum(1 for e in evidence if e.source == 'pubmed')} 条。",
        f"符号子图包含 {len(subgraph)} 条边，可用于后续 Signed GNN 或图增强 LLM 输入。",
    ]
    if any(e.source == "prediction" for e in evidence):
        reasoning.append("部分关系来自本地预测模型，应在报告中标注为预测结果而非已知数据库事实。")
    if any(e.record_id == "pubmed-error" for e in evidence):
        reasoning.append("外部检索失败，本次报告仅使用已有证据或错误记录。")

    evidence_lines = []
    for ev in evidence[:8]:
        suffix = f" URL: {ev.url}" if ev.url else ""
        evidence_lines.append(
            f"- [{ev.record_id}] {ev.source} | {ev.relation_type} | {ev.statement}{suffix}"
        )
    graph_lines = [
        f"- {edge.source} -> {edge.target} ({edge.relation_type}, {_polarity_text(edge.polarity)}, w={edge.weight:.2f})"
        for edge in subgraph[:8]
    ]
    generated = "\n".join([
        "# 药物关系分析报告",
        "",
        "## 结构化结论",
        answer,
        "",
        "## 推理过程",
        *[f"- {line}" for line in reasoning],
        "",
        "## 证据来源",
        *(evidence_lines or ["- 暂无证据"]),
        "",
        "## 符号子图",
        *(graph_lines or ["- 暂无符号边"]),
    ])
    return AnalysisReport(
        task=task,
        answer=answer,
        reasoning=reasoning,
        evidence=evidence,
        signed_subgraph=subgraph,
        generated_content=generated,
    )


def export_report(report: AnalysisReport, export_dir: Optional[str] = None) -> Dict[str, str]:
    """Export report JSON and Markdown files."""
    export_dir = export_dir or Config.EXPORTS_DATA_DIR
    json_dir = os.path.join(export_dir, "json")
    md_dir = Config.REPORTS_DIR
    os.makedirs(json_dir, exist_ok=True)
    os.makedirs(md_dir, exist_ok=True)

    stamp = time.strftime("%Y%m%d_%H%M%S")
    safe_type = re.sub(r"[^A-Za-z0-9_-]+", "_", report.task.task_type)[:40] or "analysis"
    base_name = f"{stamp}_{safe_type}"
    json_path = os.path.join(json_dir, base_name + ".json")
    md_path = os.path.join(md_dir, base_name + ".md")
    report.exports = {"json": json_path, "markdown": md_path}

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(asdict(report), f, ensure_ascii=False, indent=2)
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(report.generated_content)
        f.write("\n")
    return report.exports


def analyze_query(
    agent: SearchAgent,
    query: str,
    fetch_external: Optional[bool] = None,
    export: bool = True,
) -> AnalysisReport:
    """Run the full natural-language to evidence/report pipeline."""
    task = parse_task(query, force_external=fetch_external)
    evidence = collect_local_evidence(agent, task)
    local_is_weak = not evidence or all(ev.confidence < 0.5 for ev in evidence)
    if task.requires_external_search or local_is_weak:
        evidence.extend(fetch_pubmed_evidence(query))
    subgraph = build_signed_subgraph(task, evidence)
    report = generate_report(task, evidence, subgraph)
    if export:
        export_report(report)
    return report


def format_report_for_chat(report: AnalysisReport) -> str:
    """Compact report string for the GUI assistant."""
    lines = [
        report.answer,
        "",
        "推理摘要：",
        *[f"- {item}" for item in report.reasoning],
    ]
    if report.exports:
        lines.extend([
            "",
            f"已导出 JSON：{report.exports.get('json')}",
            f"已生成 Markdown 报告：{report.exports.get('markdown')}",
        ])
    return "\n".join(lines)
