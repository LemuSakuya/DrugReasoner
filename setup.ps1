<#
.SYNOPSIS
    药研智析 DrugReasoner 一键环境配置与启动

.DESCRIPTION
    自动完成：
      1. 检测/选择 Python 解释器
      2. 创建虚拟环境并安装依赖
      3. 启动 MySQL Docker 容器
      4. 初始化数据库（可选）
      5. 配置 LLM API 环境变量（可选）
      6. 启动 GUI 主程序

.PARAMETER SkipDocker
    跳过 Docker / MySQL 启动

.PARAMETER InitDB
    强制导入 SQL 数据

.PARAMETER PredOnly
    运行预测脚本而非 GUI

.PARAMETER SetupOnly
    仅配置环境，不启动应用

.EXAMPLE
    .\setup.ps1                     # 完整部署并启动 GUI
    .\setup.ps1 -SkipDocker          # 已有外部 MySQL
    .\setup.ps1 -InitDB              # 首次部署，导入数据
    .\setup.ps1 -PredOnly            # 仅运行预测
    .\setup.ps1 -SetupOnly           # 只配环境不启动
#>

param(
    [switch]$SkipDocker,
    [switch]$InitDB,
    [switch]$PredOnly,
    [switch]$SetupOnly
)

$ErrorActionPreference = 'Stop'
$root = $PSScriptRoot
Set-Location $root

function Import-DotEnv([string]$path) {
    if (-not (Test-Path $path)) { return }
    Get-Content $path | ForEach-Object {
        $line = $_.Trim()
        if (-not $line -or $line.StartsWith('#') -or $line -notmatch '=') { return }
        $name, $value = $line.Split('=', 2)
        $name = $name.Trim()
        $value = $value.Trim().Trim('"').Trim("'")
        if ($name) {
            [Environment]::SetEnvironmentVariable($name, $value, 'Process')
        }
    }
}

Import-DotEnv (Join-Path $root '.env')

# ── 颜色输出工具 ────────────────────────────────────────────────────────────
function Write-Step([string]$msg) {
    Write-Host "`n>> $msg" -ForegroundColor Cyan
}
function Write-OK([string]$msg) { Write-Host "   OK  $msg" -ForegroundColor Green }
function Write-Warn([string]$msg) { Write-Host "   ??  $msg" -ForegroundColor Yellow }
function Write-Err([string]$msg) { Write-Host "   !!  $msg" -ForegroundColor Red }

# ── 横幅 ────────────────────────────────────────────────────────────────────
Write-Host ""
Write-Host "  ========================================" -ForegroundColor Magenta
Write-Host "   药研智析 DrugReasoner  一键环境配置" -ForegroundColor Magenta
Write-Host "  ========================================" -ForegroundColor Magenta

# ══════════════════════════════════════════════════════════════════════════════
# 步骤 1：检测 Python
# ══════════════════════════════════════════════════════════════════════════════
Write-Step "步骤 1/5：检测 Python 环境"

$venvDir = Join-Path $root '.venv'
$venvPy  = Join-Path $venvDir 'Scripts\python.exe'

function Find-Python {
    # 1. 用户通过环境变量指定
    $userPy = $env:DEEPBIND_PYTHON
    if ($userPy) {
        $userPy = $userPy.Trim('"').Trim()
        if (Test-Path $userPy) {
            Write-OK "使用 DEEPBIND_PYTHON: $userPy"
            return $userPy
        }
        Write-Warn "DEEPBIND_PYTHON 指定的路径不存在: $userPy"
    }

    # 2. 查找 Python 3.9（项目推荐版本）
    $try = @(
        @{exe='py';     args=@('-3.9')},
        @{exe='python3.9'; args=@()},
        @{exe='python3';   args=@()},
        @{exe='python';    args=@()}
    )
    foreach ($t in $try) {
        try {
            $ver = & $t.exe @($t.args) --version 2>&1
            if ($LASTEXITCODE -eq 0) {
                $path = & $t.exe @($t.args) -c "import sys; print(sys.executable)" 2>$null
                $path = ($path | Select-Object -First 1).Trim()
                if ($path -and (Test-Path $path)) {
                    Write-OK "找到 $ver -> $path"
                    return $path
                }
            }
        } catch {}
    }

    # 3. 扫描常见安装路径
    $commonPaths = @(
        "$env:LOCALAPPDATA\Programs\Python\Python39\python.exe",
        "$env:LOCALAPPDATA\Programs\Python\Python310\python.exe",
        "$env:LOCALAPPDATA\Programs\Python\Python311\python.exe",
        "$env:LOCALAPPDATA\Programs\Python\Python312\python.exe",
        "$env:LOCALAPPDATA\Programs\Python\Python313\python.exe",
        "C:\Python39\python.exe",
        "C:\Python310\python.exe",
        "C:\Program Files\Python39\python.exe",
        "C:\Program Files\Python310\python.exe"
    )
    foreach ($p in $commonPaths) {
        if (Test-Path $p) {
            Write-OK "找到: $p"
            return $p
        }
    }

    return $null
}

$pythonExe = Find-Python
if (-not $pythonExe) {
    Write-Err "未找到 Python。请安装 Python 3.9+: https://www.python.org/downloads/"
    Write-Warn "或设置环境变量 DEEPBIND_PYTHON 指向 python.exe"
    exit 1
}

$pyVer = & $pythonExe --version 2>&1
if ($pyVer -notmatch '3\.(9|1[0-9])') {
    Write-Warn "推荐 Python 3.9，当前: $pyVer（可能不兼容 RDKit）"
}

# ══════════════════════════════════════════════════════════════════════════════
# 步骤 2：创建/修复虚拟环境 & 安装依赖
# ══════════════════════════════════════════════════════════════════════════════
Write-Step "步骤 2/5：Python 虚拟环境与依赖"

$needRebuild = $true
if (Test-Path $venvPy) {
    try {
        & $venvPy -c "import sys" 2>$null
        if ($LASTEXITCODE -eq 0) {
            Write-OK "虚拟环境可用，跳过创建"
            $needRebuild = $false
        }
    } catch {}
}

if ($needRebuild) {
    if (Test-Path $venvDir) {
        Write-Warn "虚拟环境已损坏，重建中..."
        Remove-Item -Recurse -Force $venvDir -ErrorAction SilentlyContinue
    }
    Write-Host "   创建 .venv ..."
    & $pythonExe -m venv $venvDir
    if ($LASTEXITCODE -ne 0) {
        Write-Err "虚拟环境创建失败"
        exit 1
    }
    Write-OK "虚拟环境创建成功"
}

# 确保 pip 最新并安装依赖
Write-Host "   升级 pip ..."
& $venvPy -m pip install --upgrade pip setuptools wheel --quiet 2>&1 | Out-Null

Write-Host "   安装依赖 (requirements.txt) ..."
& $venvPy -m pip install -r requirements.txt --quiet
if ($LASTEXITCODE -ne 0) {
    Write-Warn "部分依赖安装失败，尝试逐个安装..."
    & $venvPy -m pip install -r requirements.txt
    if ($LASTEXITCODE -ne 0) {
        Write-Err "依赖安装失败，请检查网络或 requirements.txt"
        exit 1
    }
}
Write-OK "依赖安装完成"

# ══════════════════════════════════════════════════════════════════════════════
# 步骤 3：Docker MySQL
# ══════════════════════════════════════════════════════════════════════════════
Write-Step "步骤 3/5：MySQL（Docker）"

if ($SkipDocker) {
    Write-Warn "跳过 Docker（-SkipDocker），使用外部 MySQL"
} else {
    # 检测 docker
    try {
        docker info 2>&1 | Out-Null
        if ($LASTEXITCODE -ne 0) { throw }
    } catch {
        Write-Err "Docker 未运行。请启动 Docker Desktop 或使用 -SkipDocker"
        exit 1
    }

    $container = 'deepbinddta-mysql'
    $running = docker ps --filter "name=$container" --format '{{.Names}}' 2>$null
    if ($running -eq $container) {
        Write-OK "MySQL 容器已在运行"
    } else {
        Write-Host "   启动 MySQL 容器 ..."
        docker compose up -d 2>&1
        if ($LASTEXITCODE -ne 0) {
            Write-Err "容器启动失败，请检查 docker-compose.yml"
            exit 1
        }

        # 等待就绪
        $mysqlPwd = if ($env:MYSQL_ROOT_PASSWORD) { $env:MYSQL_ROOT_PASSWORD } else { '12345' }
        Write-Host "   等待 MySQL 就绪 ..."
        $deadline = (Get-Date).AddSeconds(90)
        $ready = $false
        while ((Get-Date) -lt $deadline) {
            docker exec $container mysqladmin ping -h 127.0.0.1 -uroot "-p$mysqlPwd" --silent 2>$null
            if ($LASTEXITCODE -eq 0) { $ready = $true; break }
            Start-Sleep -Seconds 3
        }
        if (-not $ready) {
            Write-Err "MySQL 启动超时。检查: docker logs $container"
            exit 1
        }
        Write-OK "MySQL 就绪"
    }

    # 导入 SQL
    $sqlFile = Join-Path $root 'data\drug_discovery_dump.sql'
    if ($InitDB -and (Test-Path $sqlFile)) {
        Write-Host "   导入数据库 (首次约需数分钟) ..."
        $mysqlPwd = if ($env:MYSQL_ROOT_PASSWORD) { $env:MYSQL_ROOT_PASSWORD } else { '12345' }
        $dbName   = if ($env:MYSQL_DATABASE)      { $env:MYSQL_DATABASE }      else { 'drug_discovery' }

        # 创建库
        docker exec -e "MYSQL_PWD=$mysqlPwd" $container mysql -uroot -e "CREATE DATABASE IF NOT EXISTS $dbName CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;" 2>$null

        # 复制 SQL 进容器再导入（避免 Windows 重定向问题）
        $containerPath = '/tmp/dump.sql'
        docker cp "$sqlFile" "${container}:$containerPath"
        if ($LASTEXITCODE -ne 0) {
            Write-Err "复制 SQL 文件到容器失败"
            exit 1
        }
        try {
            docker exec -e "MYSQL_PWD=$mysqlPwd" $container sh -lc "mysql -uroot $dbName < $containerPath"
            if ($LASTEXITCODE -ne 0) {
                Write-Err "SQL 导入失败"
                exit 1
            }
            Write-OK "数据库导入完成"
        } finally {
            docker exec $container sh -lc "rm -f $containerPath" 2>$null
        }
    } elseif (Test-Path $sqlFile) {
        Write-OK "SQL 文件已就绪（跳过导入，加 -InitDB 可强制导入）"
    }
}

# ══════════════════════════════════════════════════════════════════════════════
# 步骤 4：LLM 配置
# ══════════════════════════════════════════════════════════════════════════════
Write-Step "步骤 4/5：LLM API 配置"

$envFile = Join-Path $root '.env'

# 检查是否已有配置
$existingKey = $env:LLM_API_KEY
$existingModel = $env:LLM_MODEL
$existingBase  = $env:LLM_BASE_URL
$existingProvider = $env:LLM_PROVIDER_MODEL

if ($existingKey -or $existingProvider) {
    Write-OK "检测到已有 LLM 环境变量，跳过配置"
    if ($existingProvider) { Write-Host "   Provider: $existingProvider" }
    if ($existingKey)      { Write-Host "   API Key: 已设置" }
} elseif (Test-Path $envFile) {
    # 读取 .env 文件检查
    $envContent = Get-Content $envFile -Raw -ErrorAction SilentlyContinue
    if ($envContent -match 'LLM_API_KEY=\S+' -or $envContent -match 'LLM_PROVIDER_MODEL=\S+') {
        Write-OK ".env 文件已有 LLM 配置"
    } else {
        $doConfig = $true
    }
} else {
    $doConfig = $true
}

if ($doConfig) {
    Write-Host ""
    Write-Host "  智能助手功能需要 LLM API，支持两种配置方式：" -ForegroundColor Yellow
    Write-Host ""
    Write-Host "  方式 A - 直连 DeepSeek / OpenAI 兼容 API（推荐 DeepSeek 用户）：" -ForegroundColor Cyan
    Write-Host "    设置 LLM_API_KEY + LLM_BASE_URL + LLM_MODEL"
    Write-Host ""
    Write-Host "  方式 B - LangChain Provider（推荐）：" -ForegroundColor Cyan
    Write-Host "    设置 LLM_PROVIDER_MODEL，如: openai:gpt-4o / groq:llama-3.1-70b / anthropic:claude-sonnet-4-6"
    Write-Host "    对应的 *_API_KEY 也需设置（OPENAI_API_KEY / GROQ_API_KEY / ANTHROPIC_API_KEY）"
    Write-Host ""

    $choice = Read-Host "  配置 LLM？[A=直连 / B=Provider / N=跳过]"
    if ($choice -eq 'A' -or $choice -eq 'a') {
        $key   = Read-Host "  API Key"
        $base  = Read-Host "  Base URL (默认 https://api.deepseek.com)"
        $model = Read-Host "  Model (默认 deepseek-v4-pro)"
        if (-not $base)  { $base  = 'https://api.deepseek.com' }
        if (-not $model) { $model = 'deepseek-v4-pro' }

        @"
# DrugReasoner LLM 配置 (DeepSeek / OpenAI-compatible API)
LLM_API_KEY=$key
DEEPSEEK_API_KEY=$key
LLM_BASE_URL=$base
LLM_MODEL=$model
LLM_THINKING=disabled
LLM_REASONING_EFFORT=
"@ | Out-File -FilePath $envFile -Encoding utf8

        $env:LLM_API_KEY  = $key
        $env:DEEPSEEK_API_KEY = $key
        $env:LLM_BASE_URL = $base
        $env:LLM_MODEL    = $model
        Write-OK "配置已保存到 .env"
    }
    elseif ($choice -eq 'B' -or $choice -eq 'b') {
        $prov = Read-Host "  Provider Model (如 openai:gpt-4o / groq:llama-3.1-70b)"
        $pkey = Read-Host "  对应的 API Key"
        if ($prov) {
            $providerPrefix = ($prov -split ':')[0].ToUpper()
            @"
# DrugReasoner LLM 配置 (LangChain Provider)
LLM_PROVIDER_MODEL=$prov
${providerPrefix}_API_KEY=$pkey
"@ | Out-File -FilePath $envFile -Encoding utf8

            $env:LLM_PROVIDER_MODEL = $prov
            if ($providerPrefix -eq 'OPENAI')   { $env:OPENAI_API_KEY   = $pkey }
            if ($providerPrefix -eq 'GROQ')     { $env:GROQ_API_KEY     = $pkey }
            if ($providerPrefix -eq 'ANTHROPIC'){ $env:ANTHROPIC_API_KEY = $pkey }
            if ($providerPrefix -eq 'DEEPSEEK') { $env:DEEPSEEK_API_KEY  = $pkey }
            Write-OK "配置已保存到 .env"
        }
    }
    else {
        Write-Warn "跳过 LLM 配置（智能助手功能将不可用）"
    }
}

# ══════════════════════════════════════════════════════════════════════════════
# 步骤 5：启动应用
# ══════════════════════════════════════════════════════════════════════════════
Write-Step "步骤 5/5：启动"

# Matplotlib 缓存固定在项目内
$mplDir = Join-Path $root '.mplconfig'
New-Item -ItemType Directory -Force -Path $mplDir | Out-Null
$env:MPLCONFIGDIR = $mplDir

if ($SetupOnly) {
    Write-OK "环境配置完成（-SetupOnly，跳过启动）"
    Write-Host ""
    Write-Host "  启动 GUI ：.\setup.ps1 -SkipDocker" -ForegroundColor Cyan
    Write-Host "  运行预测 ：.\setup.ps1 -SkipDocker -PredOnly" -ForegroundColor Cyan
    exit 0
}

if ($PredOnly) {
    $target = Join-Path $root 'pred.py'
    Write-OK "运行预测脚本: pred.py"
} else {
    $target = Join-Path $root 'app.py'
    Write-OK "启动 GUI 主程序: app.py"
}

Write-Host ""
Write-Host "  ========================================" -ForegroundColor Magenta
Write-Host ""

& $venvPy -u $target
