# DeepBindDTA

<p align="right">
  <a href="README.en.md"><b>English</b></a> | <a href="README.zh.md">中文</a>
</p>

## 📌 Project Overview

**DeepBindDTA** is an intelligent drug relationship analysis system that integrates deep learning, knowledge graphs, and large language models (LLMs). It is designed for:

- **Drug–Target Affinity (DTA) Prediction** — dual-encoder architecture (Mol2Vec + ESM-2) with bilinear attention
- **Drug–Drug Interaction (DDI) Analysis** — signed directed graph visualization
- **Drug & Protein Query** — database-backed entity search with alias resolution
- **AI Assistant** — LangChain-powered agent with local tool integration

### 🏗️ Architecture

```
┌─────────────────────────────────────────────┐
│                  GUI (app.py)               │
│  Drug Query │ DDI Analysis │ DTA Prediction │
│             │  3D Viewer   │  AI Assistant  │
└──────────────────┬──────────────────────────┘
                   │
       ┌───────────┴───────────┐
       ▼                       ▼
  MySQL Database          LLMDTA Model
  (drug_discovery)    (llmdta.py + attention_blocks.py)
       │                       │
  DDI / DTI data        Mol2Vec + ESM-2
  Knowledge Graph         features
```

### 📁 Key Files

| File | Description |
|------|-------------|
| `app.py` | Main GUI application (Tkinter) |
| `llmdta.py` | LLMDTA neural network model |
| `attention_blocks.py` | Bilinear attention & Transformer blocks |
| `dataset.py` | PyTorch Dataset and DataLoader |
| `model_config.py` | Model hyperparameters |
| `utils.py` | Feature extraction (Mol2Vec / ESM-2) with caching |
| `data_extractor.py` | Drug/protein data extraction utilities |
| `config.py` | Global configuration (paths, DB, GUI) |
| `pred.py` | Standalone batch prediction script |
| `setup.ps1` | **One-click environment setup and launcher** |
| `.env.example` | Local database and LLM environment template |

---

## 🚀 Quick Start

### Prerequisites

- Windows 10/11
- Python 3.9 ([download](https://www.python.org/downloads/release/python-3913/))
- Docker Desktop (for MySQL) — [download](https://www.docker.com/products/docker-desktop/)
- Git LFS — `git lfs install`

### Option A — One-click Setup (Recommended)

```powershell
.\setup.ps1
```

This will automatically detect Python, create or repair `.venv`, install dependencies, start MySQL, configure the optional LLM assistant, and launch the GUI.

To persist local passwords or LLM keys, copy the environment template first:

```powershell
Copy-Item .env.example .env
notepad .env
```

For DeepSeek, use the OpenAI-compatible direct API mode:

```env
LLM_API_KEY=your DeepSeek API key
DEEPSEEK_API_KEY=your DeepSeek API key
LLM_BASE_URL=https://api.deepseek.com
LLM_MODEL=deepseek-v4-pro
LLM_THINKING=disabled
```

Common modes:

```powershell
.\setup.ps1 -SkipDocker    # Use an external MySQL instance
.\setup.ps1 -InitDB        # Import SQL data on first setup
.\setup.ps1 -PredOnly      # Run batch prediction only
.\setup.ps1 -SetupOnly     # Prepare the environment without launching
```

To use a specific Python interpreter:

```powershell
$env:DEEPBIND_PYTHON = "C:\Path\To\python.exe"
.\setup.ps1
```

### Option B — Manual Setup

```powershell
# 1. Create a virtual environment and install dependencies
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt

# 2. Start MySQL
docker compose up -d

# 3. Import the database on first setup

# 4. Launch the GUI
.venv\Scripts\python.exe app.py
```

### Database Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `DB_HOST` | `localhost` | MySQL host |
| `DB_PORT` | `3306` | MySQL port |
| `DB_USER` | `root` | MySQL user |
| `DB_PASSWORD` | `12345` | MySQL password |
| `DB_NAME` | `drug_discovery` | Database name |

> Recommendation: set your own `MYSQL_ROOT_PASSWORD` / `DB_PASSWORD` in `.env` for local development, and never commit real secrets.

---

## 🔧 Git LFS

This repository uses **Git LFS** for large files (`*.pth`, `*.pkl`, `*.pt`, `*.h5`, `*.ckpt`, `*.onnx`, `*.sql`).

```bash
# After cloning, fetch all LFS objects
git lfs install
git lfs pull
```

> **Note:** If you have an existing clone, the safest update is to reclone:
> ```bash
> git clone https://github.com/LemuSakuya/DeepBindDTA--Preview-.git
> cd DeepBindDTA--Preview-
> git lfs pull
> ```

---

## 📎 Contact

Open an issue on GitHub for questions, bugs, or LFS tracking requests.
