# 药研智析 DrugReasoner 项目完成清单

## 一、推荐命名

- 软件中文名：药研智析
- 软件英文名：DrugReasoner
- 推荐项目文件夹名：DrugReasoner
- 对应申请书题目：融合语言理解与符号关系推理的药物数据分析系统及其实现

命名理由：原文件夹名 `software` 过于临时，旧名 `DeepBindDTA` 又偏向药物-靶标亲和力单任务。`DrugReasoner` 能覆盖自然语言解析、符号关系推理、关系预测和报告生成四个核心模块，中文名“药研智析”更适合答辩展示和软著材料。

## 二、当前已补齐内容

- 统一软件展示名为“药研智析 DrugReasoner”。
- 在 `config.py` 中新增软件副标题和推荐文件夹名。
- 首页标题对齐“语言理解 + 符号关系推理”的大创方向。
- 智能助手标题和说明改为项目品牌名。
- 帮助窗口增加兜底逻辑，缺少帮助文件时不再崩溃。
- 新增 `help_interduction.txt`，用于 GUI 帮助按钮。
- 已有 `analysis_pipeline.py` 支持自然语言解析、证据聚合、符号子图构建、Markdown/JSON 报告导出。
- 已有 `search_agent.py` 支持药物、蛋白、DDI/DTI/PPI 的统一检索入口。

## 三、你需要在本机完成的配置

- 安装 Python 3.9 或兼容版本，并确认 `python --version` 可用。
- 安装 Docker Desktop，并保持 Docker 正常运行。
- 复制 `.env.example` 为 `.env`。
- 在 `.env` 中填写本地 MySQL 密码：
  - `MYSQL_ROOT_PASSWORD`
  - `DB_PASSWORD`
- 首次部署时运行：
  ```powershell
  .\setup.ps1 -InitDB
  ```
- 日常启动时运行：
  ```powershell
  .\setup.ps1
  ```
- 如果你已经手动启动 MySQL，运行：
  ```powershell
  .\setup.ps1 -SkipDocker
  ```
- 如果要使用智能助手，填写以下 LLM 配置之一：
  ```env
  LLM_API_KEY=你的API Key
  LLM_BASE_URL=https://api.deepseek.com
  LLM_MODEL=deepseek-v4-pro
  ```
  或者使用 LangChain Provider：
  ```env
  LLM_PROVIDER_MODEL=openai:gpt-4o-mini
  OPENAI_API_KEY=你的API Key
  ```
- 如果模型大文件来自 Git LFS，克隆后运行：
  ```powershell
  git lfs install
  git lfs pull
  ```

## 四、建议你手动验收的演示流程

- 点击“帮助”，确认帮助文档能正常打开。
- 点击“药物信息查询”，输入 `Aspirin`，确认能返回药物信息。
- 点击“药物-药物分析”，输入 `Aspirin` 和 `Warfarin`，确认能返回 DDI 关系。
- 点击“DTA预测”，选择或输入一个药物和 `davis_prots.csv` 中存在的靶标 `prot_id`，确认预测文件能生成。
- 点击“DTA/DDI对比”，分别测试 DTA 和 DDI 两侧按钮。
- 打开“智能助手”，点击“环境检测”和“API测试”。
- 在智能助手输入：
  ```text
  Aspirin和Warfarin的相互作用并生成证据报告
  ```
  确认 `reports/` 和 `data/exports/json/` 中生成新报告。

## 五、后续研发任务

- 完善实体解析 Agent：覆盖药物别名、靶标别名、疾病名、突变位点和任务类型槽位。
- 扩展符号化知识底座：补充 DrugBank、CTD、PubChem、MalaCards 等数据源的来源字段、证据句、置信分数和正负极性。
- 增加 DDA 药物-疾病关系模块，使系统覆盖 DDI、DTI、DTA、DDA 四类任务。
- 将外部文献检索从 PubMed 轻量摘要扩展到可复核证据包，保存 PMID、标题、年份、摘要和检索词。
- 加入报告引用格式，使 Markdown 报告能直接用于阶段汇报。
- 增加结构化日志和案例回放，便于答辩时复现实验过程。
- 对 GUI 做响应式调整，避免不同 DPI 下按钮文字挤压。
- 增加自动化测试脚本，至少覆盖配置读取、QueryParser、SearchAgent 和 AnalysisPipeline。
- 给预测模型增加失败提示，包括缺少模型文件、缺少 RDKit、缺少 ESM 权重、CUDA/CPU 切换失败等情况。
- 准备软著材料：软件说明书、操作手册、核心源代码前后 30 页、界面截图和模块流程图。

## 六、推荐文件夹重命名方式

当前项目路径是：

```text
E:\VSCode\software
```

建议最终改为：

```text
E:\VSCode\DrugReasoner
```

如果当前 VS Code、终端或 Python 进程仍在占用该目录，先关闭相关窗口，再在 PowerShell 中从父目录执行：

```powershell
Rename-Item -LiteralPath "E:\VSCode\software" -NewName "DrugReasoner"
```

重命名后重新用 VS Code 打开：

```text
E:\VSCode\DrugReasoner
```
