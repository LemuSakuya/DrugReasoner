"""
Search Agent for DeepBindDTA
=============================
Unified search module supporting drug, protein, and interaction search
with fuzzy matching, alias resolution, and natural language query parsing.

Usage:
    # 1. Standalone
    from search_agent import SearchAgent
    agent = SearchAgent(datareader=datareader)
    results = agent.search("阿司匹林 副作用")

    # 2. As LangChain tools
    from search_agent import build_search_tools
    tools = build_search_tools(agent)

    # 3. Integrated with app.py PROJECT_TOOL_REGISTRY
    from search_agent import register_to_project
    register_to_project()
"""

import re
import difflib
import os
from typing import Optional, List, Dict, Tuple, Any, Callable
from dataclasses import dataclass, field

# ── Data Classes ────────────────────────────────────────────────────────────

@dataclass
class DrugResult:
    """Structured drug search result."""
    name: str
    drug_id: str = ""
    match_type: str = "exact"       # exact | fuzzy | alias | feature | side_effect
    match_score: float = 1.0
    match_detail: str = ""
    smiles: Optional[str] = None
    side_effects: List[str] = field(default_factory=list)
    chemical_features: List[str] = field(default_factory=list)


@dataclass
class ProteinResult:
    """Structured protein search result."""
    name: str
    gene_name: str = ""
    match_type: str = "exact"
    match_score: float = 1.0
    match_detail: str = ""
    aliases: List[str] = field(default_factory=list)


@dataclass
class InteractionResult:
    """Structured interaction search result."""
    entity1: str
    entity2: str
    interaction_type: str          # DDI | DTI | PPI
    result: str
    confidence: Optional[float] = None
    source: str = "database"


@dataclass
class SearchResponse:
    """Top-level search response with parsed intent and ranked results."""
    query: str
    parsed_intent: str             # drug_info | drug_drug | drug_protein | protein_info | protein_protein | general
    results: List[Any] = field(default_factory=list)
    suggestions: List[str] = field(default_factory=list)
    summary: str = ""


# ── Natural Language Query Parser ───────────────────────────────────────────

class QueryParser:
    """Parse Chinese/English natural language queries into structured intents."""

    INTENT_PATTERNS = [
        # Drug-drug interaction
        ("drug_drug", [
            r'(.+?)(?:和|与|vs\.?|VS\.?|versus)(.+?)(?:的)?(?:药物)?(?:相互)?(?:作用|关系|影响|反应)',
            r'(.+?)(?:和|与)(.+?)(?:能否|可以|能)(?:一起|同时|联合)(?:使用|服用|用药)',
            r'比较(.+?)(?:和|与)(.+)',
            r'(.+?)\s+(?:and|vs\.?|versus)\s+(.+?)\s+(?:interaction|ddi)',
        ]),
        # Drug-protein / DTA
        ("drug_protein", [
            r'(.+?)(?:对|靶向|作用于|结合)(.+?)(?:的)?(?:亲和力|结合力|活性|作用|效果)',
            r'(.+?)(?:和|与)(.+?)(?:靶标|蛋白|受体)(?:的)?(?:亲和力|结合|作用)',
            r'预测(.+?)(?:和|与)(.+?)(?:的)?(?:DTA|亲和力|结合)',
            r'(.+?)\s+(?:targets?|binds?\s+to|affinity\s+(?:for|to))\s+(.+)',
        ]),
        # Drug info
        ("drug_info", [
            r'(?:查询|查找|搜索|查看)(.+?)(?:的)?(?:药物)?(?:信息|性质|详情|资料|数据)',
            r'(.+?)(?:的)?(?:副作用|化学特征|SMILES|分子式|分子结构)',
            r'(.+?)(?:是什么|是什么药|是何种药物)',
            r'(?:what\s+is|tell\s+me\s+about|info\s+(?:on|about))\s+(.+)',
        ]),
        # Protein info
        ("protein_info", [
            r'(?:查询|查找|搜索|查看)(.+?)(?:的)?(?:蛋白|靶标|受体)(?:信息|性质|详情)',
            r'(.+?)(?:基因|蛋白序列|氨基酸序列)',
            r'(?:what\s+is|tell\s+me\s+about)\s+(.+?)\s+(?:protein|target|receptor)',
        ]),
        # Protein-protein interaction
        ("protein_protein", [
            r'(.+?)(?:和|与)(.+?)(?:的)?(?:蛋白)?(?:相互)?(?:作用|关系)',
        ]),
    ]

    @classmethod
    def parse(cls, query: str) -> Tuple[str, Dict[str, Any]]:
        """Parse a natural language query into (intent, extracted_params)."""
        query = query.strip()
        params: Dict[str, Any] = {}

        for intent, patterns in cls.INTENT_PATTERNS:
            for pattern in patterns:
                m = re.search(pattern, query, re.IGNORECASE)
                if m:
                    groups = [g.strip() for g in m.groups() if g and g.strip()]
                    if intent in ("drug_drug", "drug_protein", "protein_protein"):
                        if len(groups) >= 2:
                            params["entity1"] = groups[0]
                            params["entity2"] = groups[1]
                        elif len(groups) == 1:
                            params["entity1"] = groups[0]
                    elif intent in ("drug_info", "protein_info"):
                        params["entity1"] = groups[0] if groups else query
                    return intent, params

        # Fallback: extract drug/protein names by common separators
        for sep in [r'\s+和\s+', r'\s+与\s+', r'\s+vs\.?\s+', r'\s+and\s+', r',', r'，']:
            parts = re.split(sep, query, maxsplit=1)
            if len(parts) == 2 and len(parts[0]) > 1 and len(parts[1]) > 1:
                params["entity1"] = parts[0].strip()
                params["entity2"] = parts[1].strip()
                return "general_pair", params

        params["entity1"] = query
        return "general", params


# ── Main Search Agent ───────────────────────────────────────────────────────

class SearchAgent:
    """
    Unified search agent for drug discovery data.

    Provides search across drugs, proteins, and their interactions,
    with fuzzy matching, alias resolution, and NL query support.
    """

    def __init__(
        self,
        datareader: Callable[[str], Any],
        extra_dependencies: Optional[Dict[str, Any]] = None
    ):
        """
        Args:
            datareader: Function that takes a table name and returns a numpy array (from DB).
            extra_dependencies: Optional dict of extra functions/state needed by the agent.
        """
        self._datareader = datareader
        self._cache: Dict[str, Any] = {}
        self._extra = extra_dependencies or {}

    # ── Cache helpers ───────────────────────────────────────────────────

    def _cached(self, key: str, loader: Callable[[], Any]) -> Any:
        if key not in self._cache:
            self._cache[key] = loader()
        return self._cache[key]

    def _load_drug_names(self) -> Tuple[List[str], List[str]]:
        """Return (drug_names, drug_ids)."""
        data = self._datareader("NamesWithID")
        return data[:, 0].tolist(), data[:, 1].tolist()

    def _load_protein_names(self) -> List[str]:
        """Return list of protein names from ProteinNameID."""
        data = self._datareader("ProteinNameID")
        return data[:, 0].tolist()

    def _load_protein_alias_map(self) -> Dict[str, str]:
        """Build lowercase alias -> canonical protein name mapping."""
        try:
            data = self._datareader("ProteinNameID")
        except Exception:
            return {}
        alias_map: Dict[str, str] = {}
        for row in data:
            row_vals = [str(v).strip() for v in row if v is not None and str(v).strip()]
            if not row_vals:
                continue
            canonical = row_vals[0]
            for v in row_vals:
                alias_map[v.lower()] = canonical
        return alias_map

    def _load_offside_names(self) -> List[str]:
        data = self._datareader("OffsideName")
        return [str(row[0]) for row in data]

    def _load_feature_names(self) -> List[str]:
        try:
            import csv
            csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    'pubchem_fingerprints.csv')
            with open(csv_path, 'r') as f:
                reader = csv.reader(f)
                return [row[1] for row in reader]
        except Exception:
            return []

    # ── Fuzzy matching ──────────────────────────────────────────────────

    @staticmethod
    def _fuzzy_filter(
        query: str,
        candidates: List[str],
        top_k: int = 10,
        cutoff: float = 0.4
    ) -> List[Tuple[str, float]]:
        """Return top-k fuzzy matches with scores."""
        if not candidates:
            return []
        q = query.strip().lower()
        # Exact match first
        exact = [c for c in candidates if c.lower() == q]
        if exact:
            return [(exact[0], 1.0)]
        # Prefix match
        prefix = [(c, 0.95) for c in candidates if c.lower().startswith(q)]
        if prefix:
            return prefix[:top_k]
        # Contains match
        contains = [(c, 0.85) for c in candidates if q in c.lower()]
        # difflib fallback
        close = difflib.get_close_matches(query, candidates, n=top_k, cutoff=cutoff)
        scored = []
        for c in close:
            score = difflib.SequenceMatcher(None, q, c.lower()).ratio()
            scored.append((c, score))
        # Merge: contains matches score higher than pure difflib
        seen = {c for c, _ in scored}
        for c, s in contains:
            if c not in seen:
                scored.append((c, s))
        scored.sort(key=lambda x: -x[1])
        return scored[:top_k]

    # ── Drug Search ─────────────────────────────────────────────────────

    def search_drug(
        self,
        query: str,
        top_k: int = 10,
        include_side_effects: bool = False,
        include_features: bool = False
    ) -> List[DrugResult]:
        """
        Search for drugs by name, ID, side effect, or chemical feature.

        Args:
            query: Drug name, ID, or property keyword.
            top_k: Max results to return.
            include_side_effects: Whether to include side effect details.
            include_features: Whether to include chemical feature details.

        Returns:
            List of DrugResult sorted by match score.
        """
        drug_names, drug_ids = self._load_drug_names()
        results: List[DrugResult] = []

        # 1. Exact name match
        q = query.strip()
        if q in drug_names:
            idx = drug_names.index(q)
            r = DrugResult(
                name=q, drug_id=str(drug_ids[idx]),
                match_type="exact", match_score=1.0
            )
            self._enrich_drug_result(r, idx, include_side_effects, include_features)
            results.append(r)

        # 2. Exact ID match
        id_matches = [(name, i) for i, (name, did) in enumerate(zip(drug_names, drug_ids))
                       if str(did) == q]
        for name, idx in id_matches:
            if not any(r.name == name for r in results):
                r = DrugResult(
                    name=name, drug_id=str(drug_ids[idx]),
                    match_type="exact", match_score=1.0, match_detail=f"ID: {q}"
                )
                self._enrich_drug_result(r, idx, include_side_effects, include_features)
                results.append(r)

        # 3. Fuzzy name match
        fuzzy = self._fuzzy_filter(q, drug_names, top_k=top_k, cutoff=0.3)
        for name, score in fuzzy:
            if not any(r.name == name for r in results):
                idx = drug_names.index(name)
                r = DrugResult(
                    name=name, drug_id=str(drug_ids[idx]),
                    match_type="fuzzy", match_score=score,
                    match_detail=f"similarity={score:.2f}"
                )
                self._enrich_drug_result(r, idx, include_side_effects, include_features)
                results.append(r)

        # 4. Search by side effect
        if len(results) < top_k:
            try:
                offside_names = self._load_offside_names()
                offside_matches = self._fuzzy_filter(q, offside_names, top_k=5, cutoff=0.3)
                if offside_matches:
                    offsides_matrix = self._datareader("Offsides").reshape((1, -1))[0]
                    n_drugs = len(drug_names)
                    n_offsides = len(offside_names)
                    offsides_matrix = offsides_matrix[:n_drugs * n_offsides].reshape((n_drugs, n_offsides))
                    for side_name, side_score in offside_matches:
                        side_idx = offside_names.index(side_name)
                        for i in range(n_drugs):
                            if offsides_matrix[i][side_idx] == 1:
                                if not any(r.name == drug_names[i] for r in results):
                                    r = DrugResult(
                                        name=drug_names[i], drug_id=str(drug_ids[i]),
                                        match_type="side_effect", match_score=side_score,
                                        match_detail=f"side_effect: {side_name}"
                                    )
                                    self._enrich_drug_result(r, i, include_side_effects, include_features)
                                    results.append(r)
                                    if len(results) >= top_k:
                                        break
                        if len(results) >= top_k:
                            break
            except Exception:
                pass

        # 5. Search by chemical feature name
        if len(results) < top_k:
            try:
                feat_names = self._load_feature_names()
                feat_matches = self._fuzzy_filter(q, feat_names, top_k=5, cutoff=0.3)
                if feat_matches:
                    feat_matrix = self._datareader("drug_881feat")
                    for feat_name, feat_score in feat_matches:
                        if feat_name not in feat_names:
                            continue
                        feat_idx = feat_names.index(feat_name)
                        if feat_idx >= feat_matrix.shape[1]:
                            continue
                        for i in range(min(len(drug_names), feat_matrix.shape[0])):
                            if feat_matrix[i][feat_idx] == 1:
                                if not any(r.name == drug_names[i] for r in results):
                                    r = DrugResult(
                                        name=drug_names[i], drug_id=str(drug_ids[i]),
                                        match_type="feature", match_score=feat_score,
                                        match_detail=f"feature: {feat_name}"
                                    )
                                    self._enrich_drug_result(r, i, include_side_effects, include_features)
                                    results.append(r)
                                    if len(results) >= top_k:
                                        break
                        if len(results) >= top_k:
                            break
            except Exception:
                pass

        results.sort(key=lambda r: -r.match_score)
        return results[:top_k]

    def _enrich_drug_result(
        self, result: DrugResult, idx: int,
        include_side_effects: bool, include_features: bool
    ):
        """Optionally fill in side effects and chemical features."""
        if include_side_effects:
            try:
                offside_names = self._datareader("OffsideName")
                offsides = self._datareader("Offsides").reshape((1, -1))[0]
                n_drugs = len(self._load_drug_names()[0])
                n_offsides = len(offside_names)
                offsides = offsides[:n_drugs * n_offsides].reshape((n_drugs, n_offsides))
                result.side_effects = [
                    str(offside_names[i][0])
                    for i in range(n_offsides) if offsides[idx][i] == 1
                ][:20]
            except Exception:
                pass
        if include_features:
            try:
                feat_names = self._load_feature_names()
                feat_matrix = self._datareader("drug_881feat")
                result.chemical_features = [
                    feat_names[i] for i in range(min(len(feat_names), feat_matrix.shape[1]))
                    if feat_matrix[idx][i] == 1
                ][:20]
            except Exception:
                pass

    # ── Protein Search ──────────────────────────────────────────────────

    def search_protein(self, query: str, top_k: int = 10) -> List[ProteinResult]:
        """
        Search for proteins by name, gene name, or alias.

        Args:
            query: Protein name, gene symbol, or alias keyword.
            top_k: Max results to return.

        Returns:
            List of ProteinResult sorted by match score.
        """
        results: List[ProteinResult] = []

        try:
            prot_data = self._datareader("ProteinNameID")
        except Exception:
            return results

        # Build alias map and candidate list
        alias_map = self._load_protein_alias_map()
        canonical_names: List[str] = []
        gene_map: Dict[str, str] = {}            # canonical -> gene
        all_aliases: Dict[str, List[str]] = {}   # canonical -> all aliases

        for row in prot_data:
            row_vals = [str(v).strip() for v in row if v is not None and str(v).strip()]
            if not row_vals:
                continue
            canonical = row_vals[0]
            canonical_names.append(canonical)
            gene = row_vals[1] if len(row_vals) > 1 else ""
            gene_map[canonical] = gene
            all_aliases[canonical] = row_vals

        q = query.strip()

        # 1. Exact alias match
        key = q.lower()
        if key in alias_map:
            canonical = alias_map[key]
            results.append(ProteinResult(
                name=canonical, gene_name=gene_map.get(canonical, ""),
                match_type="exact", match_score=1.0,
                aliases=all_aliases.get(canonical, [])
            ))

        # 2. Exact canonical match
        if q in canonical_names and not any(r.name == q for r in results):
            results.append(ProteinResult(
                name=q, gene_name=gene_map.get(q, ""),
                match_type="exact", match_score=1.0,
                aliases=all_aliases.get(q, [])
            ))

        # 3. Fuzzy match on all aliases
        all_candidates = list(set(
            v for row in prot_data
            for v in [str(x).strip() for x in row if x is not None and str(x).strip()]
        ))
        fuzzy = self._fuzzy_filter(q, all_candidates, top_k=top_k, cutoff=0.3)
        for match_name, score in fuzzy:
            canonical = alias_map.get(match_name.lower(), match_name)
            if not any(r.name == canonical for r in results):
                results.append(ProteinResult(
                    name=canonical, gene_name=gene_map.get(canonical, ""),
                    match_type="fuzzy", match_score=score,
                    match_detail=f"matched via: {match_name}",
                    aliases=all_aliases.get(canonical, [])
                ))

        results.sort(key=lambda r: -r.match_score)
        return results[:top_k]

    # ── Interaction Search ──────────────────────────────────────────────

    def search_interaction(
        self,
        entity1: str,
        entity2: str,
        interaction_type: str = "auto"
    ) -> Optional[InteractionResult]:
        """
        Search for known or predicted interaction between two entities.

        Args:
            entity1: First entity name (drug or protein).
            entity2: Second entity name (drug or protein).
            interaction_type: 'DDI', 'DTI', 'PPI', or 'auto' (auto-detect).

        Returns:
            InteractionResult if found, else None.
        """
        drug_names, drug_ids = self._load_drug_names()
        e1_is_drug = entity1 in drug_names
        e2_is_drug = entity2 in drug_names

        prot_names = self._load_protein_names()
        e1_is_prot = entity1 in prot_names
        e2_is_prot = entity2 in prot_names

        # Auto-detect interaction type
        if interaction_type == "auto":
            if e1_is_drug and e2_is_drug:
                interaction_type = "DDI"
            elif e1_is_prot and e2_is_prot:
                interaction_type = "PPI"
            elif (e1_is_drug and e2_is_prot) or (e1_is_prot and e2_is_drug):
                interaction_type = "DTI"
            else:
                # Try fuzzy resolution
                fuzzy_d1 = self._fuzzy_filter(entity1, drug_names, top_k=1, cutoff=0.5)
                fuzzy_d2 = self._fuzzy_filter(entity2, drug_names, top_k=1, cutoff=0.5)
                fuzzy_p1 = self._fuzzy_filter(entity1, prot_names, top_k=1, cutoff=0.5)
                fuzzy_p2 = self._fuzzy_filter(entity2, prot_names, top_k=1, cutoff=0.5)
                d1_score = fuzzy_d1[0][1] if fuzzy_d1 else 0
                d2_score = fuzzy_d2[0][1] if fuzzy_d2 else 0
                p1_score = fuzzy_p1[0][1] if fuzzy_p1 else 0
                p2_score = fuzzy_p2[0][1] if fuzzy_p2 else 0
                if d1_score > 0.5 and d2_score > 0.5:
                    interaction_type = "DDI"
                elif p1_score > 0.5 and p2_score > 0.5:
                    interaction_type = "PPI"
                elif (d1_score > 0.5 and p2_score > 0.5) or (p1_score > 0.5 and d2_score > 0.5):
                    interaction_type = "DTI"
                else:
                    return None

        try:
            if interaction_type == "DDI":
                if entity1 not in drug_names or entity2 not in drug_names:
                    return None
                idx1, idx2 = drug_names.index(entity1), drug_names.index(entity2)
                result = self._predict_ddi(idx1, idx2)
                return InteractionResult(
                    entity1=entity1, entity2=entity2,
                    interaction_type="DDI", result=result[0],
                    source="prediction" if "预测" in str(result) else "database"
                )

            elif interaction_type == "DTI":
                drug, prot = (entity1, entity2) if e1_is_drug else (entity2, entity1)
                if drug not in drug_names or prot not in prot_names:
                    return None
                d_idx = drug_names.index(drug)
                p_idx = prot_names.index(prot)
                action = self._lookup_dta(d_idx, p_idx)
                return InteractionResult(
                    entity1=drug, entity2=prot,
                    interaction_type="DTI", result=action
                )

            elif interaction_type == "PPI":
                if entity1 not in prot_names or entity2 not in prot_names:
                    return None
                idx1, idx2 = prot_names.index(entity1), prot_names.index(entity2)
                ppi_result = self._lookup_ppi(idx1, idx2)
                return InteractionResult(
                    entity1=entity1, entity2=entity2,
                    interaction_type="PPI", result=ppi_result[0]
                )
        except Exception:
            return None

        return None

    def _predict_ddi(self, idx1: int, idx2: int) -> List[str]:
        """Replicate app.py prediction() logic."""
        import torch
        import joblib
        try:
            drug_drug_sign = self._datareader("drug_drug_sign")
            temp = drug_drug_sign[:, 0:2].tolist()
            drug_ids = self._load_drug_names()[1]
            drug_label = [drug_ids[idx1], drug_ids[idx2]]
            ans = []
            if drug_label in temp:
                sign = drug_drug_sign[temp.index(drug_label)][2]
                ans.append('拮抗作用' if sign == 1 else '协同作用')
            else:
                f = torch.Tensor(joblib.load("feature.pkl"))
                weight = torch.Tensor(joblib.load("weight.pkl"))
                bias = torch.Tensor(joblib.load("bias.pkl"))
                edge = [idx1, idx2]
                value = torch.cat([f[edge[0]], f[edge[1]]], dim=1)
                value = torch.nn.functional.linear(value, weight.T, bias)
                pred = torch.nn.functional.log_softmax(value)
                _, topk_indices = torch.topk(pred, k=1)
                for i in topk_indices:
                    ans.append('拮抗作用' if i == 1 else '协同作用')
            return ans
        except Exception as e:
            return [f"预测失败: {e}"]

    def _lookup_dta(self, drug_idx: int, prot_idx: int) -> str:
        """Look up drug-target interaction type from DTA table."""
        try:
            DTA = self._datareader("DTA").reshape((1, -1))[0][:2089464].reshape((1443, 1448))
            dta_num = DTA[drug_idx][prot_idx]
            action_table = self._datareader("drug_protein_action")
            return str(action_table[:, 1][dta_num - 1])
        except Exception as e:
            return f"查询失败: {e}"

    def _lookup_ppi(self, idx1: int, idx2: int) -> List[str]:
        """Look up protein-protein interaction."""
        try:
            PPI = self._datareader("PPI").reshape((1, -1))[0][:1243225].reshape((1115, 1115))
            sign = int(float(PPI[idx1][idx2]))
            return ['可相互作用' if sign == 1 else '无明显互作用']
        except Exception as e:
            return [f"查询失败: {e}"]

    # ── High-level search entry ─────────────────────────────────────────

    def search(self, query: str, top_k: int = 10) -> SearchResponse:
        """
        Main search entry point. Parses natural language, dispatches to
        appropriate sub-searches, and returns a unified response.

        Args:
            query: Natural language search query (Chinese or English).
            top_k: Max results per category.

        Returns:
            SearchResponse with parsed intent and ranked results.
        """
        intent, params = QueryParser.parse(query)

        if intent == "drug_info":
            drug_query = params.get("entity1", query)
            drugs = self.search_drug(drug_query, top_k=top_k,
                                     include_side_effects=True, include_features=True)
            suggestions = self._build_suggestions(drugs, "drug")
            return SearchResponse(
                query=query, parsed_intent=intent, results=drugs,
                suggestions=suggestions,
                summary=self._format_drug_summary(drugs)
            )

        elif intent == "protein_info":
            prot_query = params.get("entity1", query)
            prots = self.search_protein(prot_query, top_k=top_k)
            suggestions = self._build_suggestions(prots, "protein")
            return SearchResponse(
                query=query, parsed_intent=intent, results=prots,
                suggestions=suggestions,
                summary=self._format_protein_summary(prots)
            )

        elif intent in ("drug_drug", "drug_protein", "protein_protein", "general_pair"):
            e1 = params.get("entity1", "")
            e2 = params.get("entity2", "")
            if e1 and e2:
                interaction = self.search_interaction(e1, e2)
                if interaction:
                    return SearchResponse(
                        query=query, parsed_intent=intent,
                        results=[interaction],
                        summary=f"{interaction.entity1} ↔ {interaction.entity2}: {interaction.result}"
                    )
            # Fallback: search each entity individually
            d1 = self.search_drug(e1, top_k=3) if e1 else []
            d2 = self.search_drug(e2, top_k=3) if e2 else []
            p1 = self.search_protein(e1, top_k=3) if e1 else []
            p2 = self.search_protein(e2, top_k=3) if e2 else []
            all_results = list(d1 + d2 + p1 + p2)
            return SearchResponse(
                query=query, parsed_intent=intent, results=all_results,
                summary=f"找到 {len(all_results)} 个可能匹配"
            )

        else:
            # General: try everything
            drugs = self.search_drug(query, top_k=5,
                                     include_side_effects=True, include_features=True)
            prots = self.search_protein(query, top_k=5)
            all_results = list(drugs) + list(prots)
            return SearchResponse(
                query=query, parsed_intent="general", results=all_results,
                suggestions=self._build_suggestions(all_results, "any"),
                summary=self._format_general_summary(drugs, prots)
            )

    # ── Formatting helpers ──────────────────────────────────────────────

    def _build_suggestions(self, results: List[Any], kind: str) -> List[str]:
        suggestions = []
        for r in results[:5]:
            if isinstance(r, DrugResult):
                suggestions.append(f"查看 {r.name} 的详细药物信息")
                suggestions.append(f"查看 {r.name} 的药物关系图")
                if r.side_effects:
                    suggestions.append(f"查看 {r.name} 的副作用详情")
            elif isinstance(r, ProteinResult):
                suggestions.append(f"查看 {r.name} 的蛋白质网络图")
                suggestions.append(f"查看 {r.name} 的蛋白互作用关系")
        return suggestions[:8]

    def _format_drug_summary(self, drugs: List[DrugResult]) -> str:
        if not drugs:
            return "未找到匹配的药物。"
        lines = [f"找到 {len(drugs)} 个药物匹配："]
        for r in drugs[:5]:
            parts = [f"• {r.name} (ID: {r.drug_id})"]
            if r.match_type != "exact":
                parts.append(f" [{r.match_type}: {r.match_detail}]")
            if r.side_effects:
                parts.append(f" | 副作用({len(r.side_effects)}项)")
            if r.chemical_features:
                parts.append(f" | 化学特征({len(r.chemical_features)}项)")
            lines.append("".join(parts))
        return "\n".join(lines)

    def _format_protein_summary(self, prots: List[ProteinResult]) -> str:
        if not prots:
            return "未找到匹配的蛋白质。"
        lines = [f"找到 {len(prots)} 个蛋白质匹配："]
        for r in prots[:5]:
            parts = [f"• {r.name}"]
            if r.gene_name:
                parts.append(f" (基因: {r.gene_name})")
            if r.aliases:
                parts.append(f" 别名: {', '.join(r.aliases[:4])}")
            if r.match_type != "exact":
                parts.append(f" [{r.match_type}]")
            lines.append("".join(parts))
        return "\n".join(lines)

    def _format_general_summary(
        self, drugs: List[DrugResult], prots: List[ProteinResult]
    ) -> str:
        parts = []
        if drugs:
            parts.append(f"药物: {len(drugs)} 个匹配")
        if prots:
            parts.append(f"蛋白质: {len(prots)} 个匹配")
        return "，".join(parts) if parts else "未找到匹配结果。"


# ── LangChain Tool Builder ──────────────────────────────────────────────────

def build_search_tools(agent: SearchAgent) -> List[Any]:
    """
    Build LangChain-compatible tool functions from a SearchAgent instance.

    Returns a list of @tool-decorated functions that can be passed to
    create_agent(system_prompt=..., tools=[...]).

    Usage:
        from langchain.agents import create_agent
        from langchain_core.tools import tool

        agent = SearchAgent(datareader=my_datareader)
        search_tools = build_search_tools(agent)
        lc_agent = create_agent(model=llm, tools=search_tools, system_prompt=...)
    """
    try:
        from langchain_core.tools import tool as lc_tool
    except ImportError:
        # Fallback: return plain functions
        return [
            lambda q: agent.search(q),
            lambda q: agent.search_drug(q),
            lambda q: agent.search_protein(q),
            lambda e1, e2: agent.search_interaction(e1, e2),
        ]

    @lc_tool
    def search_all(query: str) -> str:
        """
        统一搜索入口。搜索药物、蛋白质或相互作用。
        输入自然语言查询（中文或英文），返回结构化搜索结果。
        例如：'阿司匹林的副作用'、'吉非替尼和EGFR的结合亲和力'、'查找布洛芬的药物信息'
        """
        try:
            response = agent.search(query)
            return response.summary
        except Exception as e:
            return f"搜索失败: {e}"

    @lc_tool
    def search_drug_tool(query: str) -> str:
        """
        搜索药物信息。输入药物名称、ID、副作用关键词或化学特征名称。
        返回匹配的药物列表，包括名称、ID、副作用和化学特征。
        例如：'阿司匹林'、'头痛'（查找有头痛副作用的药物）
        """
        try:
            drugs = agent.search_drug(query, include_side_effects=True, include_features=True)
            return agent._format_drug_summary(drugs)
        except Exception as e:
            return f"药物搜索失败: {e}"

    @lc_tool
    def search_protein_tool(query: str) -> str:
        """
        搜索蛋白质/靶标信息。输入蛋白质名称、基因名或别名。
        返回匹配的蛋白质列表，包括名称、基因名和别名。
        例如：'EGFR'、'BRCA1'、'epidermal growth factor receptor'
        """
        try:
            prots = agent.search_protein(query)
            return agent._format_protein_summary(prots)
        except Exception as e:
            return f"蛋白质搜索失败: {e}"

    @lc_tool
    def search_interaction_tool(entity1: str, entity2: str) -> str:
        """
        查询两个实体之间的相互作用。自动检测 DDI/DTA/PPI 类型。
        输入两个药物名称或蛋白质名称，返回相互作用结果。
        例如：search_interaction_tool('阿司匹林', '布洛芬')
        """
        try:
            result = agent.search_interaction(entity1, entity2)
            if result:
                return f"{result.entity1} ↔ {result.entity2} ({result.interaction_type}): {result.result}"
            return f"未找到 {entity1} 与 {entity2} 之间的已知相互作用。"
        except Exception as e:
            return f"相互作用查询失败: {e}"

    return [search_all, search_drug_tool, search_protein_tool, search_interaction_tool]


# ── Integration with app.py's PROJECT_TOOL_REGISTRY ─────────────────────────

def register_to_project(agent: Optional[SearchAgent] = None) -> Dict[str, Callable]:
    """
    Register search agent functions into app.py's PROJECT_TOOL_REGISTRY format.

    Usage in app.py:
        from search_agent import SearchAgent, register_to_project
        agent = SearchAgent(datareader=datareader)
        PROJECT_TOOL_REGISTRY.update(register_to_project(agent))

    Returns a dict compatible with PROJECT_TOOL_REGISTRY schema.
    """
    _agent = agent  # captured by closures below

    def _search_all(payload: str) -> str:
        if _agent is None:
            return "搜索代理未初始化。"
        return _agent.search(payload).summary

    def _search_drug_func(payload: str) -> str:
        if _agent is None:
            return "搜索代理未初始化。"
        drugs = _agent.search_drug(payload, include_side_effects=True, include_features=True)
        return _agent._format_drug_summary(drugs)

    def _search_protein_func(payload: str) -> str:
        if _agent is None:
            return "搜索代理未初始化。"
        prots = _agent.search_protein(payload)
        return _agent._format_protein_summary(prots)

    def _search_interaction_func(payload: str) -> str:
        if _agent is None:
            return "搜索代理未初始化。"
        parts = [p.strip() for p in re.split(r'[,，\s|;；]+', payload) if p.strip()]
        if len(parts) < 2:
            return "需要两个实体名称，例如：阿司匹林, 布洛芬"
        result = _agent.search_interaction(parts[0], parts[1])
        if result:
            return f"{result.entity1} ↔ {result.entity2} ({result.interaction_type}): {result.result}"
        return f"未找到 {parts[0]} 与 {parts[1]} 之间的已知相互作用。"

    return {
        "search": _search_all,
        "search_drug": _search_drug_func,
        "search_protein": _search_protein_func,
        "search_interaction": _search_interaction_func,
    }
