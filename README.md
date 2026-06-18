# 药研智析 DrugReasoner

<p align="right">
  <b>中文</b> | <a href="README.en.md">English</a> | <a href="README.zh.md">完整中文说明</a>
</p>

## 项目简介

**药研智析 DrugReasoner** 是面向大学生创新训练项目“融合语言理解与符号关系推理的药物数据分析系统及其实现”的软件原型。系统融合本地药物数据库、符号关系图、LLMDTA 深度学习模型与 LLM Agent，用于完成药物关系查询、智能预测、自然语言交互和可解释报告生成。

核心能力：

- **自然语言解析**：将用户问题解析为药物、靶标、任务类型和证据需求等结构化参数。
- **符号关系推理**：围绕 DDI、DTI、DTA、PPI 等关系构建带极性的有向关系结果。
- **药物-靶标亲和力预测**：基于 Mol2Vec + ESM-2 双编码器和双线性注意力机制调用 LLMDTA 模型。
- **药物-药物相互作用分析**：查询或预测协同/拮抗关系，并提供关系图展示。
- **智能助手与报告导出**：通过 LangChain/HTTP LLM 接口调用本地工具，导出 Markdown 与 JSON 证据报告。

## 系统架构

```text
┌────────────────────────────────────────────────────────┐
│        GUI 主程序 app.py / 智能助手 / 报告导出         │
│ 药物查询 │ DDI 分析 │ DTA 预测 │ DTA/DDI 对比 │ Agent │
└───────────────────────────┬────────────────────────────┘
                            │
          ┌─────────────────┼─────────────────┐
          ▼                 ▼                 ▼
      MySQL 数据库      SearchAgent      LLMDTA 模型
  drug_discovery     统一检索与别名解析  Mol2Vec + ESM-2
          │                 │                 │
          ▼                 ▼                 ▼
  DDI/DTI/PPI 数据   AnalysisPipeline   预测结果 CSV
  符号关系知识底座    证据包/符号子图     可视化与报告
```

## 核心文件

| 文件 | 说明 |
| --- | --- |
| `app.py` | Tkinter GUI 主程序 |
| `analysis_pipeline.py` | 自然语言解析、证据聚合、符号子图和报告导出 |
| `search_agent.py` | 药物、蛋白、相互作用统一检索 Agent |
| `llmdta.py` | LLMDTA 神经网络模型定义 |
| `attention_blocks.py` | 双线性注意力与 Transformer 模块 |
| `dataset.py` | PyTorch 数据集与预测数据加载 |
| `model_config.py` | 模型超参数配置 |
| `utils.py` | Mol2Vec / ESM-2 特征提取与缓存 |
| `data_extractor.py` | 药物/蛋白预测输入提取工具 |
| `config.py` | 全局配置、路径、数据库和软件名称 |
| `pred.py` | 独立批量预测脚本 |
| `setup.ps1` | 一键环境配置与启动脚本 |
| `docs/PROJECT_TODO.md` | 项目配置、验收和后续研发清单 |

## 快速开始

### 环境要求

- Windows 10/11
- Python 3.9+
- Docker Desktop，用于启动本地 MySQL
- Git LFS，用于拉取大模型和数据文件

### 一键部署

```powershell
.\setup.ps1
```

首次部署并导入数据库：

```powershell
.\setup.ps1 -InitDB
```

如果已有外部 MySQL：

```powershell
.\setup.ps1 -SkipDocker
```

仅配置环境不启动 GUI：

```powershell
.\setup.ps1 -SetupOnly
```

### LLM 配置

复制环境模板：

```powershell
Copy-Item .env.example .env
notepad .env
```

DeepSeek / OpenAI-compatible 配置示例：

```env
LLM_API_KEY=你的 DeepSeek API Key
DEEPSEEK_API_KEY=你的 DeepSeek API Key
LLM_BASE_URL=https://api.deepseek.com
LLM_MODEL=deepseek-v4-pro
LLM_THINKING=disabled
```

## 常用演示流程

1. 点击“药物信息查询”，输入 `Aspirin`。
2. 点击“药物-药物分析”，输入 `Aspirin` 与 `Warfarin`。
3. 点击“DTA预测”，输入药物名称和 `davis_prots.csv` 中的 `prot_id`。
4. 打开“智能助手”，输入：

```text
Aspirin和Warfarin的相互作用并生成证据报告
```

报告会导出到：

- `reports/`
- `data/exports/json/`

## 数据库环境变量

| 变量名 | 默认值 | 说明 |
| --- | --- | --- |
| `DB_HOST` | `localhost` | MySQL 主机地址 |
| `DB_PORT` | `3306` | MySQL 端口 |
| `DB_USER` | `root` | MySQL 用户名 |
| `DB_PASSWORD` | `12345` | MySQL 密码 |
| `DB_NAME` | `drug_discovery` | 数据库名称 |

## 更多文档

- [完整中文说明](README.zh.md)
- [English README](README.en.md)
- [项目配置与验收清单](docs/PROJECT_TODO.md)
