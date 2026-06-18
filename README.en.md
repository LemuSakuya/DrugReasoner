# DrugReasoner

<p align="right">
  <a href="README.en.md"><b>English</b></a> | <a href="README.zh.md">中文</a>
</p>

## Overview

**DrugReasoner** is a drug data analysis prototype for the innovation project “Drug Data Analysis System Integrating Language Understanding and Symbolic Relation Reasoning”. It combines a local drug database, signed relation graphs, the LLMDTA deep learning model, and an LLM Agent to support drug relation search, prediction, natural-language interaction, and explainable report generation.

Core capabilities:

- **Natural-language parsing**: map user questions into structured entities, task types, and evidence requirements.
- **Symbolic relation reasoning**: handle signed and directed DDI, DTI, DTA, and PPI relations.
- **Drug-target affinity prediction**: run the LLMDTA model with Mol2Vec + ESM-2 features and bilinear attention.
- **Drug-drug interaction analysis**: query or predict synergistic/antagonistic DDI relations and visualize graphs.
- **Agent assistant and report export**: call local tools through LangChain or an HTTP LLM endpoint, then export Markdown and JSON reports.

## Architecture

```text
┌────────────────────────────────────────────────────────┐
│        GUI app.py / Assistant / Report Export          │
│ Drug Query │ DDI │ DTA │ DTA/DDI Compare │ Agent       │
└───────────────────────────┬────────────────────────────┘
                            │
          ┌─────────────────┼─────────────────┐
          ▼                 ▼                 ▼
      MySQL Database    SearchAgent       LLMDTA Model
  drug_discovery        entity search     Mol2Vec + ESM-2
          │                 │                 │
          ▼                 ▼                 ▼
  DDI/DTI/PPI data   AnalysisPipeline   prediction CSV
  signed knowledge   evidence/subgraph  visualization/report
```

## Key Files

| File | Description |
| --- | --- |
| `app.py` | Tkinter GUI application |
| `analysis_pipeline.py` | NL parsing, evidence aggregation, signed subgraph, report export |
| `search_agent.py` | Unified drug, protein, and interaction search agent |
| `llmdta.py` | LLMDTA neural network model |
| `attention_blocks.py` | Bilinear attention and Transformer modules |
| `dataset.py` | PyTorch dataset and prediction dataloader |
| `model_config.py` | Model hyperparameters |
| `utils.py` | Mol2Vec / ESM-2 feature extraction with caching |
| `data_extractor.py` | Drug/protein prediction input extraction |
| `config.py` | Global paths, DB settings, and software metadata |
| `pred.py` | Standalone prediction script |
| `setup.ps1` | One-click setup and launcher |
| `docs/PROJECT_TODO.md` | Configuration, acceptance, and roadmap checklist |

## Quick Start

### Requirements

- Windows 10/11
- Python 3.9+
- Docker Desktop for local MySQL
- Git LFS for large model and data files

### One-click setup

```powershell
.\setup.ps1
```

First setup with database import:

```powershell
.\setup.ps1 -InitDB
```

Use an existing MySQL instance:

```powershell
.\setup.ps1 -SkipDocker
```

Prepare the environment without launching the GUI:

```powershell
.\setup.ps1 -SetupOnly
```

### LLM configuration

```powershell
Copy-Item .env.example .env
notepad .env
```

DeepSeek / OpenAI-compatible example:

```env
LLM_API_KEY=your DeepSeek API key
DEEPSEEK_API_KEY=your DeepSeek API key
LLM_BASE_URL=https://api.deepseek.com
LLM_MODEL=deepseek-v4-pro
LLM_THINKING=disabled
```

## Demo Flow

1. Open “Drug Query” and search `Aspirin`.
2. Open “Drug-Drug Analysis” and enter `Aspirin` and `Warfarin`.
3. Open “DTA Prediction” and enter a drug plus a `prot_id` from `davis_prots.csv`.
4. Open the assistant and ask:

```text
Aspirin and Warfarin interaction with an evidence report
```

Reports are exported to:

- `reports/`
- `data/exports/json/`

## Database Variables

| Variable | Default | Description |
| --- | --- | --- |
| `DB_HOST` | `localhost` | MySQL host |
| `DB_PORT` | `3306` | MySQL port |
| `DB_USER` | `root` | MySQL user |
| `DB_PASSWORD` | `12345` | MySQL password |
| `DB_NAME` | `drug_discovery` | Database name |

## Project Name

Recommended folder rename:

```text
E:\VSCode\software
```

to:

```text
E:\VSCode\DrugReasoner
```

See [docs/PROJECT_TODO.md](docs/PROJECT_TODO.md).
