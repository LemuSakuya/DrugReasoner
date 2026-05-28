# DeepBindDTA

<p align="right">
  <a href="README.en.md">English</a> | <a href="README.zh.md"><b>中文</b></a>
</p>

## 📌 项目简介

**DeepBindDTA** 是一个融合深度学习、知识图谱与大语言模型（LLM）的智能药物关系分析系统，主要功能包括：

- **药物-靶标亲和力（DTA）预测** — 基于 Mol2Vec + ESM-2 双编码器 + 双线性注意力机制
- **药物-药物相互作用（DDI）分析** — 有向符号图可视化
- **药物与蛋白质查询** — 支持别名解析的数据库实体搜索
- **AI 智能助手** — 基于 LangChain Agent，集成本地工具调用

### 🏗️ 系统架构

```
┌─────────────────────────────────────────────────┐
│               GUI 主程序 (app.py)               │
│  药物查询 │ DDI 分析 │ DTA 预测 │ AI 智能助手  │
└──────────────────┬──────────────────────────────┘
                   │
       ┌───────────┴───────────┐
       ▼                       ▼
  MySQL 数据库              LLMDTA 模型
  (drug_discovery)   (llmdta.py + attention_blocks.py)
       │                       │
  DDI/DTI 数据           Mol2Vec + ESM-2
  知识图谱                   特征提取
```

### 📁 核心文件说明

| 文件 | 说明 |
|------|------|
| `app.py` | GUI 主程序（Tkinter） |
| `llmdta.py` | LLMDTA 神经网络模型定义 |
| `attention_blocks.py` | 双线性注意力 & Transformer 模块 |
| `dataset.py` | PyTorch 数据集与数据加载器 |
| `model_config.py` | 模型超参数配置 |
| `utils.py` | 特征提取工具（Mol2Vec / ESM-2，含缓存机制） |
| `data_extractor.py` | 药物/蛋白质数据提取工具 |
| `config.py` | 全局配置（路径、数据库、GUI） |
| `pred.py` | 独立批量预测脚本 |
| `setup.ps1` | **一键环境配置与启动脚本** |
| `.env.example` | 本地数据库与 LLM 环境变量模板 |

---

## 🚀 快速开始

### 环境要求

- Windows 10/11
- Python 3.9+（[下载地址](https://www.python.org/downloads/)）
- Docker Desktop（用于 MySQL）— [下载地址](https://www.docker.com/products/docker-desktop/)
- Git LFS — `git lfs install`

### 一键部署（推荐）

```powershell
.\setup.ps1
```

脚本自动完成：检测 Python → 创建虚拟环境 → 安装依赖 → 启动 MySQL → 配置 LLM → 启动 GUI。

如需保存本地密码或 LLM Key，可先复制环境模板：

```powershell
Copy-Item .env.example .env
notepad .env
```

DeepSeek 用户推荐使用 OpenAI-compatible 直连配置：

```env
LLM_API_KEY=你的 DeepSeek API Key
DEEPSEEK_API_KEY=你的 DeepSeek API Key
LLM_BASE_URL=https://api.deepseek.com
LLM_MODEL=deepseek-v4-pro
LLM_THINKING=disabled
```

**常用参数：**

```powershell
.\setup.ps1 -SkipDocker    # 已有外部 MySQL，跳过 Docker
.\setup.ps1 -InitDB        # 首次部署，导入 SQL 数据
.\setup.ps1 -PredOnly      # 仅运行预测，不启动 GUI
.\setup.ps1 -SetupOnly     # 仅配置环境，不启动应用
```

**指定 Python 路径：**

```powershell
$env:DEEPBIND_PYTHON = "C:\Path\To\python.exe"
.\setup.ps1
```

### 手动部署

```powershell
# 1. 创建虚拟环境并安装依赖
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt

# 2. 启动 MySQL
docker compose up -d

# 3. 导入数据库（首次需要）
.\setup.ps1 -InitDB

# 4. 启动
.venv\Scripts\python.exe app.py
```

### 数据库环境变量

| 变量名 | 默认值 | 说明 |
|--------|--------|------|
| `DB_HOST` | `localhost` | MySQL 主机地址 |
| `DB_PORT` | `3306` | MySQL 端口 |
| `DB_USER` | `root` | MySQL 用户名 |
| `DB_PASSWORD` | `12345` | MySQL 密码 |
| `DB_NAME` | `drug_discovery` | 数据库名称 |

> 建议：本地开发请在 `.env` 中设置自己的 `MYSQL_ROOT_PASSWORD` / `DB_PASSWORD`，不要把真实密钥提交到仓库。

---

## 🔧 Git LFS 说明

本仓库使用 **Git LFS** 管理大型文件（`*.pth`, `*.pkl`, `*.pt`, `*.h5`, `*.ckpt`, `*.onnx`, `*.sql`）。

```bash
# 克隆后拉取所有 LFS 对象
git lfs install
git lfs pull
```

> **提示：** 如果已有本地克隆，最安全的更新方式是重新克隆：
> ```bash
> git clone https://github.com/LemuSakuya/DeepBindDTA--Preview-.git
> cd DeepBindDTA--Preview-
> git lfs pull
> ```

---

## 📎 联系与反馈

如有问题、Bug 报告或 LFS 追踪需求，请在 GitHub 上提交 Issue。
