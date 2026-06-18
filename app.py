import csv
import tkinter
from tkinter import *
from tkinter import messagebox
from tkinter import font as tkfont
import warnings
import joblib
import networkx as nx
import pandas as pd
import pymysql
import torch
import os
import subprocess
import sys
import glob
import json
import threading
import urllib.request
import urllib.error
import re
import difflib
import time
import ctypes
import queue
from typing import Optional, List
API_TIMEOUT = 15


def _enable_windows_dpi_awareness():
    if os.name != 'nt':
        return
    try:
        ctypes.windll.shcore.SetProcessDpiAwareness(1)
    except Exception:
        try:
            ctypes.windll.user32.SetProcessDPIAware()
        except Exception:
            pass


_enable_windows_dpi_awareness()

UI_FONT_FAMILY = 'Microsoft YaHei UI'
MONO_FONT_FAMILY = 'Consolas'
try:
    from langchain.agents import create_agent
    from langchain_core.tools import tool
    from langchain_core.language_models.llms import LLM
    LANGCHAIN_AVAILABLE = True
except Exception:
    LANGCHAIN_AVAILABLE = False
import numpy as np
import plotly.graph_objects as go
from plotly.offline import plot
from PIL import Image, ImageTk
import io


def _load_dotenv_file(path=None):
    path = path or os.path.join(os.path.dirname(os.path.abspath(__file__)), '.env')
    if not os.path.exists(path):
        return
    try:
        with open(path, 'r', encoding='utf-8') as f:
            for raw_line in f:
                line = raw_line.strip()
                if not line or line.startswith('#') or '=' not in line:
                    continue
                name, value = line.split('=', 1)
                name = name.strip()
                value = value.strip().strip('"').strip("'")
                if name and name not in os.environ:
                    os.environ[name] = value
    except Exception:
        pass


_load_dotenv_file()

from config import Config

os.chdir(Config.BASE_DIR)


# Matplotlib 在 Windows 上首次导入 pyplot 可能触发字体扫描/缓存构建，表现为启动期“卡住”。
# 这里改为按需导入：只有在需要绘图时才加载 pyplot。
plt = None
Axes3D = None
MATPLOTLIB_AVAILABLE = False


def _ensure_matplotlib() -> bool:
    global plt, Axes3D, MATPLOTLIB_AVAILABLE
    if plt is not None:
        return True
    try:
        # 优先把 matplotlib 缓存放在项目内，避免落到慢盘/网络盘导致导入耗时或卡住。
        root_dir = os.path.dirname(os.path.abspath(__file__))
        os.environ.setdefault('MPLCONFIGDIR', os.path.join(root_dir, '.mplconfig'))
        try:
            os.makedirs(os.environ['MPLCONFIGDIR'], exist_ok=True)
        except Exception:
            pass

        import matplotlib.pyplot as _plt
        from mpl_toolkits.mplot3d import Axes3D as _Axes3D

        plt = _plt
        Axes3D = _Axes3D
        MATPLOTLIB_AVAILABLE = True
        return True
    except Exception:
        MATPLOTLIB_AVAILABLE = False
        return False


# RDKit用于处理分子结构（延迟导入：避免启动阶段卡在加载二进制库）
Chem = None
AllChem = None
RDKIT_AVAILABLE = False


def _ensure_rdkit() -> bool:
    global Chem, AllChem, RDKIT_AVAILABLE
    if Chem is not None and AllChem is not None:
        return True
    try:
        from rdkit import Chem as _Chem
        from rdkit.Chem import AllChem as _AllChem
        Chem = _Chem
        AllChem = _AllChem
        RDKIT_AVAILABLE = True
        return True
    except Exception:
        RDKIT_AVAILABLE = False
        return False

import data_extractor
import search_agent
import analysis_pipeline

warnings.filterwarnings('ignore', category=UserWarning, module='pandas') #去除推荐SQLAlchemy警告
warnings.filterwarnings('ignore', category=FutureWarning, module='torch') #去除torch版本警告

software_name = Config.SOFTWARE_NAME

# 全局主窗口引用，用于创建子窗口
main_root = None
agent_window_ref = None

# 预测日志缓存
_LAST_PRED_LOG = ""

# 蛋白名称别名映射缓存
_PROTEIN_ALIAS_CACHE = {
    'prot_count': None,
    'alias_map': {},
    'alias_candidates': [],
    'proteinnameid_values': set()
}

_AUTOCOMPLETE_CACHE = {
    'drug_names': None,
    'protein_candidates': None,
    'protein_count': None,
    'davis_prot_ids': None
}


def _get_drug_name_candidates():
    if _AUTOCOMPLETE_CACHE['drug_names'] is not None:
        return _AUTOCOMPLETE_CACHE['drug_names']
    names = []
    try:
        names = datareader("NamesWithID")[:, 0].tolist()
    except Exception:
        try:
            df = pd.read_csv('EGFR-Case/drug.tsv', sep='\t')
            names = df.iloc[:, 0].dropna().astype(str).tolist()
        except Exception:
            names = []
    names = sorted(set(names))
    _AUTOCOMPLETE_CACHE['drug_names'] = names
    return names


def _build_protein_alias_cache(prot_ids):
    prot_count = len(prot_ids)
    if _PROTEIN_ALIAS_CACHE['prot_count'] == prot_count:
        return _PROTEIN_ALIAS_CACHE['alias_map'], _PROTEIN_ALIAS_CACHE['alias_candidates']

    alias_map = {}
    alias_candidates = []
    proteinnameid_values = set()
    try:
        protein_table = datareader("ProteinNameID")
        for row in protein_table:
            row_vals = [str(v).strip() for v in row if v is not None and str(v).strip()]
            if not row_vals:
                continue
            for v in row_vals:
                proteinnameid_values.add(v)
            name_val = row_vals[0] if len(row_vals) > 0 else ''
            gene_val = row_vals[1] if len(row_vals) > 1 else ''

            matched_prot = None
            for v in row_vals:
                if v in prot_ids:
                    matched_prot = v
                    break
            if not matched_prot and gene_val:
                if gene_val in prot_ids:
                    matched_prot = gene_val
                else:
                    for pid in prot_ids:
                        if pid.startswith(f"{gene_val}("):
                            matched_prot = pid
                            break
            if not matched_prot and name_val:
                for pid in prot_ids:
                    if pid.lower() == name_val.lower():
                        matched_prot = pid
                        break

            if not matched_prot:
                continue

            for v in row_vals:
                alias_map[v.lower()] = matched_prot
            alias_candidates.extend(row_vals)
    except Exception:
        alias_map = {}
        alias_candidates = []
        proteinnameid_values = set()

    _PROTEIN_ALIAS_CACHE['prot_count'] = prot_count
    _PROTEIN_ALIAS_CACHE['alias_map'] = alias_map
    _PROTEIN_ALIAS_CACHE['alias_candidates'] = alias_candidates
    _PROTEIN_ALIAS_CACHE['proteinnameid_values'] = proteinnameid_values
    return alias_map, alias_candidates


def _is_in_proteinnameid(input_name, prot_ids):
    if not input_name:
        return False
    _build_protein_alias_cache(prot_ids)
    values = _PROTEIN_ALIAS_CACHE.get('proteinnameid_values', set())
    if not values:
        return False
    key = input_name.strip().lower()
    return any(key == str(v).strip().lower() for v in values)


def _get_protein_name_candidates():
    try:
        df = pd.read_csv('davis_prots.csv', sep=',')
        prot_ids = df.iloc[:, 0].dropna().astype(str).tolist()
    except Exception:
        prot_ids = []
    prot_count = len(prot_ids)
    if (_AUTOCOMPLETE_CACHE['protein_candidates'] is not None and
            _AUTOCOMPLETE_CACHE['protein_count'] == prot_count):
        return _AUTOCOMPLETE_CACHE['protein_candidates']

    _, alias_candidates = _build_protein_alias_cache(prot_ids)
    if alias_candidates:
        candidates = sorted(set(alias_candidates))
    else:
        candidates = sorted(set(prot_ids))

    # 过滤掉 Ensembl 等“技术型ID”（例如：9606.ENSP00000363117），避免下拉框出现难以理解的条目。
    # 注意：这里只影响候选展示，不影响内部别名映射/解析逻辑。
    ensembl_like = re.compile(r'^(?:\d{3,}\.)?(?:ENSP|ENST|ENSG)\d+(?:\.\d+)?$', re.IGNORECASE)
    candidates = [c for c in candidates if not ensembl_like.match(str(c).strip())]

    _AUTOCOMPLETE_CACHE['protein_candidates'] = candidates
    _AUTOCOMPLETE_CACHE['protein_count'] = prot_count
    return candidates


def _get_davis_prot_ids():
    if _AUTOCOMPLETE_CACHE['davis_prot_ids'] is not None:
        return _AUTOCOMPLETE_CACHE['davis_prot_ids']
    try:
        df = pd.read_csv('davis_prots.csv', sep=',')
        prot_ids = df.iloc[:, 0].dropna().astype(str).tolist()
    except Exception:
        prot_ids = []
    prot_ids = sorted(set(prot_ids))
    _AUTOCOMPLETE_CACHE['davis_prot_ids'] = prot_ids
    return prot_ids


def _open_select_dialog(parent, title, candidates, on_select):
    win = tkinter.Toplevel(parent)
    win.title(title)
    win.geometry('520x420')
    win.config(bg='#F5F5F5')
    win.transient(parent)
    win.grab_set()
    win.resizable(False, False)

    header = tkinter.Label(win, text=title, font=(UI_FONT_FAMILY, 13, 'bold'), bg='#F5F5F5', fg='#333333')
    header.pack(pady=10)

    search_frame = tkinter.Frame(win, bg='#F5F5F5')
    search_frame.pack(pady=6, padx=12, fill=tkinter.X)
    search_entry = tkinter.Entry(search_frame, font=(UI_FONT_FAMILY, 11), relief=tkinter.SOLID, borderwidth=1)
    search_entry.pack(side=tkinter.LEFT, fill=tkinter.X, expand=True)

    list_frame = tkinter.Frame(win, bg='#F5F5F5')
    list_frame.pack(padx=12, pady=6, fill=tkinter.BOTH, expand=True)

    listbox = tkinter.Listbox(list_frame, height=14, exportselection=False)
    listbox.pack(side=tkinter.LEFT, fill=tkinter.BOTH, expand=True)
    scrollbar = tkinter.Scrollbar(list_frame, command=listbox.yview)
    scrollbar.pack(side=tkinter.RIGHT, fill=tkinter.Y)
    listbox.config(yscrollcommand=scrollbar.set)

    items = list(candidates or [])

    def _refresh(filter_text=''):
        listbox.delete(0, tkinter.END)
        key = filter_text.strip().lower()
        for item in items:
            if not key or key in str(item).lower():
                listbox.insert(tkinter.END, item)

    def _confirm(_=None):
        if not listbox.curselection():
            return
        value = listbox.get(listbox.curselection()[0])
        try:
            on_select(value)
        except Exception:
            pass
        win.destroy()

    search_entry.bind('<KeyRelease>', lambda e: _refresh(search_entry.get()))
    listbox.bind('<Double-Button-1>', _confirm)
    listbox.bind('<Return>', _confirm)
    _refresh('')
    search_entry.focus()


def _attach_autocomplete(entry, candidates_provider, max_items=8):
    listbox = tkinter.Listbox(entry.master, height=max_items, exportselection=False)
    listbox.place_forget()

    def _update_list(_=None):
        text = entry.get().strip()
        if not text:
            listbox.place_forget()
            return
        candidates = candidates_provider() or []
        key = text.lower()
        matches = [c for c in candidates if key in str(c).lower()]
        if not matches:
            listbox.place_forget()
            return
        listbox.delete(0, tkinter.END)
        for item in matches[:max_items]:
            listbox.insert(tkinter.END, item)
        listbox.place(in_=entry, relx=0, rely=1, relwidth=1)

    def _select(_=None):
        if not listbox.curselection():
            return
        value = listbox.get(listbox.curselection()[0])
        entry.delete(0, tkinter.END)
        entry.insert(0, value)
        listbox.place_forget()
        entry.focus_set()

    def _hide_later(_=None):
        entry.after(200, listbox.place_forget)

    entry.bind('<KeyRelease>', _update_list)
    entry.bind('<FocusOut>', _hide_later)
    listbox.bind('<<ListboxSelect>>', _select)
    listbox.bind('<Return>', _select)
    entry.bind('<Down>', lambda e: (listbox.focus_set(), listbox.selection_set(0)) if listbox.winfo_ismapped() and listbox.size() else None)

##########读取数据库##########
_db_cfg = Config.get_db_config()

# Lazy DB connection: do not connect at import time.
_sql_connection = None


def get_sql_connection():
    global _sql_connection
    if _sql_connection is not None:
        return _sql_connection
    _sql_connection = pymysql.connect(
        host=_db_cfg['host'],
        user=_db_cfg['user'],
        password=_db_cfg['password'],
        db=_db_cfg['database'],
        port=_db_cfg['port'],
        autocommit=False,
        charset=_db_cfg['charset'],
    )
    return _sql_connection


_TABLE_NAME_MAP = None


def _refresh_table_name_map(conn):
    global _TABLE_NAME_MAP
    try:
        df = pd.read_sql("SHOW TABLES", conn)
        tables = [str(x) for x in df.values.flatten().tolist()]
        _TABLE_NAME_MAP = {t.lower(): t for t in tables}
    except Exception:
        _TABLE_NAME_MAP = {}


def _resolve_table_name(requested: str, conn) -> str:
    if not re.fullmatch(r"[A-Za-z0-9_]+", str(requested)):
        raise ValueError(f"Invalid table name: {requested}")
    global _TABLE_NAME_MAP
    if _TABLE_NAME_MAP is None:
        _refresh_table_name_map(conn)
    actual = _TABLE_NAME_MAP.get(str(requested).lower())
    if actual:
        return actual
    _refresh_table_name_map(conn)
    actual = _TABLE_NAME_MAP.get(str(requested).lower())
    if actual:
        return actual
    raise RuntimeError(f"Table '{requested}' not found in schema '{_db_cfg['database']}'.")


def datareader(table):
    conn = get_sql_connection()
    actual_table = _resolve_table_name(table, conn)
    schema = _db_cfg['database']
    sql = f"SELECT * FROM `{schema}`.`{actual_table}`"
    return pd.read_sql(sql, conn).values


# Lazy SearchAgent: do not create at import time.
_search_agent = None


def get_search_agent():
    global _search_agent
    if _search_agent is not None:
        return _search_agent
    _search_agent = search_agent.SearchAgent(datareader=datareader)
    return _search_agent


#########设置中文字体，解决乱码问题###########
def setup_chinese_font():
    """设置中文字体为微软雅黑"""
    try:
        if not _ensure_matplotlib():
            return
        plt.rcParams['font.sans-serif'] = [UI_FONT_FAMILY, 'Microsoft YaHei', 'SimHei', 'DejaVu Sans']
        plt.rcParams['axes.unicode_minus'] = False
        plt.rcParams['font.family'] = 'sans-serif'
    except:
        pass


def setup_tk_fonts(root, base_size=10):
    """Set crisp Windows-friendly Tk fonts after a root window exists."""
    try:
        root.tk.call('tk', 'scaling', 1.2)
    except Exception:
        pass
    font_specs = {
        'TkDefaultFont': (UI_FONT_FAMILY, base_size),
        'TkTextFont': (UI_FONT_FAMILY, base_size),
        'TkFixedFont': (MONO_FONT_FAMILY, base_size),
        'TkMenuFont': (UI_FONT_FAMILY, base_size),
        'TkHeadingFont': (UI_FONT_FAMILY, base_size + 1, 'bold'),
        'TkCaptionFont': (UI_FONT_FAMILY, base_size),
        'TkSmallCaptionFont': (UI_FONT_FAMILY, max(8, base_size - 1)),
        'TkIconFont': (UI_FONT_FAMILY, base_size),
        'TkTooltipFont': (UI_FONT_FAMILY, base_size),
    }
    for name, spec in font_specs.items():
        try:
            tkfont.nametofont(name).configure(family=spec[0], size=spec[1])
            if len(spec) > 2:
                tkfont.nametofont(name).configure(weight=spec[2])
        except Exception:
            pass


#########统一按钮样式（简化重复代码）#########
MAIN_BTN_GRADIENT = ['#A5B4FC', '#93C5FD', '#60A5FA', '#3B82F6', '#60A5FA', '#93C5FD']
PRIMARY_BUTTON_STYLE = {
    'font': (UI_FONT_FAMILY, 12),
    'bg': MAIN_BTN_GRADIENT[0],
    'fg': 'white',
    'relief': tkinter.FLAT,
    'cursor': 'hand2',
    'activebackground': MAIN_BTN_GRADIENT[1],
    'activeforeground': 'white'
}
SECONDARY_BUTTON_STYLE = {
    'font': (UI_FONT_FAMILY, 12),
    'bg': '#95A5A6',
    'fg': 'white',
    'relief': tkinter.FLAT,
    'cursor': 'hand2',
    'activebackground': '#7F8C8D',
    'activeforeground': 'white'
}


def _create_rounded_rect(canvas, x1, y1, x2, y2, radius, **kwargs):
    radius = max(0, min(radius, (x2 - x1) / 2, (y2 - y1) / 2))
    items = []
    items.append(canvas.create_rectangle(x1 + radius, y1, x2 - radius, y2, **kwargs))
    items.append(canvas.create_rectangle(x1, y1 + radius, x2, y2 - radius, **kwargs))
    items.append(canvas.create_arc(x1, y1, x1 + 2 * radius, y1 + 2 * radius, start=90, extent=90, style=tkinter.PIESLICE, **kwargs))
    items.append(canvas.create_arc(x2 - 2 * radius, y1, x2, y1 + 2 * radius, start=0, extent=90, style=tkinter.PIESLICE, **kwargs))
    items.append(canvas.create_arc(x1, y2 - 2 * radius, x1 + 2 * radius, y2, start=180, extent=90, style=tkinter.PIESLICE, **kwargs))
    items.append(canvas.create_arc(x2 - 2 * radius, y2 - 2 * radius, x2, y2, start=270, extent=90, style=tkinter.PIESLICE, **kwargs))
    return items


class RoundedButton(tkinter.Canvas):
    def __init__(self, parent, text, command=None, font=None, bg='#3B82F6', fg='white',
                 activebackground=None, activeforeground=None, radius=10,
                 width=None, height=None, padx=12, pady=6, **kwargs):
        self._text = text
        self._command = command
        self._font = font or (UI_FONT_FAMILY, 12)
        self._bg = bg
        self._fg = fg
        self._hover_bg = activebackground or bg
        self._hover_fg = activeforeground or fg
        self._radius = radius
        self._padx = padx
        self._pady = pady
        self._enabled = True

        parent_bg = '#F5F5F5'
        try:
            parent_bg = parent.cget('bg')
        except Exception:
            pass

        text_font = tkfont.Font(font=self._font)
        text_width = text_font.measure(self._text)
        text_height = text_font.metrics('linespace')
        # 用中文字宽估算（中文约 2x 拉丁字符宽）
        char_width = text_font.measure('中')
        if width is not None:
            desired_width = int(width) * char_width + self._padx * 2
        else:
            desired_width = text_width + self._padx * 2
        if height is not None:
            desired_height = int(height) * text_height + self._pady * 2
        else:
            desired_height = text_height + self._pady * 2

        super().__init__(
            parent,
            width=desired_width,
            height=desired_height,
            highlightthickness=0,
            bg=parent_bg,
            bd=0,
            **kwargs
        )
        self.configure(cursor='hand2')
        self._text_id = None
        self._shape_ids = []
        self.bind('<Enter>', self._on_enter)
        self.bind('<Leave>', self._on_leave)
        self.bind('<Button-1>', self._on_click)
        self.bind('<Configure>', self._on_resize)
        # 延迟绘制：等 widget 被 pack/place 映射后再画，否则 winfo_width 返回 1
        self.after(1, self._draw)

    def _draw(self):
        self.delete('all')
        w = self.winfo_width()
        h = self.winfo_height()
        if w <= 1:
            w = self.winfo_reqwidth()
        if h <= 1:
            h = self.winfo_reqheight()
        w = max(1, w)
        h = max(1, h)
        self._shape_ids = _create_rounded_rect(
            self, 1, 1, w - 1, h - 1, self._radius,
            fill=self._bg, outline=self._bg
        )
        self._text_id = self.create_text(
            w / 2, h / 2,
            text=self._text,
            font=self._font,
            fill=self._fg
        )

    def _on_resize(self, _):
        self._draw()

    def _on_enter(self, _):
        if not self._enabled:
            return
        for item in self._shape_ids:
            self.itemconfigure(item, fill=self._hover_bg, outline=self._hover_bg)
        if self._text_id:
            self.itemconfigure(self._text_id, fill=self._hover_fg)

    def _on_leave(self, _):
        if not self._enabled:
            return
        for item in self._shape_ids:
            self.itemconfigure(item, fill=self._bg, outline=self._bg)
        if self._text_id:
            self.itemconfigure(self._text_id, fill=self._fg)

    def _on_click(self, _):
        if not self._enabled:
            return
        if self._command:
            self._command()

    def set_state(self, enabled=True):
        self._enabled = enabled
        if not enabled:
            disabled_bg = '#CBD5E1'
            disabled_fg = '#94A3B8'
            for item in self._shape_ids:
                self.itemconfigure(item, fill=disabled_bg, outline=disabled_bg)
            if self._text_id:
                self.itemconfigure(self._text_id, fill=disabled_fg)

    def config(self, **kwargs):
        state = kwargs.pop('state', None)
        if state is not None:
            self.set_state(state != tkinter.DISABLED)
        command = kwargs.pop('command', None)
        if command is not None:
            self._command = command
        return super().config(**kwargs)

    configure = config


def create_button(parent, text, command, style='primary', **kwargs):
    base_style = PRIMARY_BUTTON_STYLE if style == 'primary' else SECONDARY_BUTTON_STYLE
    _ = kwargs.pop('animate', style == 'primary')
    merged = {**base_style, **kwargs}
    width = merged.pop('width', None)
    height = merged.pop('height', None)
    padx = merged.pop('padx', 12)
    pady = merged.pop('pady', 6)
    font = merged.pop('font', (UI_FONT_FAMILY, 12))
    bg = merged.pop('bg', '#3B82F6')
    fg = merged.pop('fg', 'white')
    activebackground = merged.pop('activebackground', MAIN_BTN_GRADIENT[-1])
    activeforeground = merged.pop('activeforeground', fg)
    btn = RoundedButton(
        parent,
        text=text,
        command=command,
        font=font,
        bg=bg,
        fg=fg,
        activebackground=activebackground,
        activeforeground=activeforeground,
        radius=10,
        width=width,
        height=height,
        padx=padx,
        pady=pady
    )
    return btn


def _get_pred_output_path():
    default_path = 'Pred_egfr_soft.csv'
    if os.path.exists(default_path):
        return default_path
    matches = sorted(glob.glob('Pred_*_soft.csv'), key=os.path.getmtime, reverse=True)
    return matches[0] if matches else default_path


def _resolve_protein_id(input_name, prot_ids):
    if not input_name:
        return None, []
    if input_name in prot_ids:
        return input_name, []
    lower_map = {p.lower(): p for p in prot_ids}
    key = input_name.lower()
    if key in lower_map:
        return lower_map[key], []
    contains_matches = [p for p in prot_ids if key in p.lower()]
    if len(contains_matches) == 1:
        return contains_matches[0], []
    alias_map, alias_candidates = _build_protein_alias_cache(prot_ids)
    if key in alias_map:
        return alias_map[key], []

    suggest_pool = alias_candidates or prot_ids
    suggestions = difflib.get_close_matches(input_name, suggest_pool, n=5, cutoff=0.6)
    return None, suggestions


def _run_pred_update_async(on_success=None, on_error=None):
    script_path = os.path.join(Config.BASE_DIR, 'pred.py')
    venv_python = os.path.join(Config.BASE_DIR, '.venv', 'Scripts', 'python.exe')
    python_exec = venv_python if os.path.exists(venv_python) else sys.executable

    def _dispatch(callback):
        if not callback:
            return
        if main_root:
            try:
                if main_root.winfo_exists():
                    try:
                        main_root.after(0, callback)
                        return
                    except Exception:
                        pass
            except Exception:
                pass
        try:
            callback()
        except Exception:
            pass
        return

    def _worker():
        try:
            if not os.path.exists(script_path):
                raise FileNotFoundError(f"未找到预测脚本：{script_path}")
            
            # 使用 Popen 以便实时打印输出
            # Add -u to force unbuffered output from the child python process
            # Read binary to avoid TextIOWrapper buffering
            process = subprocess.Popen(
                [python_exec, '-u', script_path],
                cwd=Config.BASE_DIR,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                bufsize=0  # Unbuffered binary
            )
            
            output_lines = []
            while True:
                line_bytes = process.stdout.readline()
                if not line_bytes and process.poll() is not None:
                    break
                if line_bytes:
                    line = line_bytes.decode('utf-8', errors='replace')
                    print(line, end='', flush=True)  # 实时打印到终端
                    output_lines.append(line)
            
            process.wait()
            
            global _LAST_PRED_LOG
            _LAST_PRED_LOG = "".join(output_lines)
            
            if process.returncode != 0:
                raise RuntimeError(f"pred.py 退出码：{process.returncode}")
            
            _dispatch(on_success)
        except Exception as exc:
            if on_error:
                _dispatch(lambda: on_error(exc))

    threading.Thread(target=_worker, daemon=True).start()


def _exit_app():
    global main_root
    try:
        if main_root and main_root.winfo_exists():
            main_root.destroy()
    except Exception:
        pass
    try:
        if start and start.winfo_exists():
            start.destroy()
    except Exception:
        pass
    try:
        sys.exit(0)
    except SystemExit:
        return


def _get_llm_config():
    base_url = os.getenv('LLM_BASE_URL', 'https://api.deepseek.com').strip()
    api_key = os.getenv('LLM_API_KEY', os.getenv('DEEPSEEK_API_KEY', '')).strip()
    model = os.getenv('LLM_MODEL', 'deepseek-v4-pro').strip()
    if base_url.endswith('/'):
        base_url = base_url[:-1]
    if 'api.deepseek.com' not in base_url and not base_url.endswith('/v1'):
        base_url = base_url + '/v1'
    return base_url, api_key, model


def _get_provider_key_state():
    provider_model = os.getenv('LLM_PROVIDER_MODEL', '').strip()
    if not provider_model:
        return "", ""
    if provider_model.startswith('openai:'):
        return provider_model, os.getenv('OPENAI_API_KEY', '').strip()
    if provider_model.startswith('groq:'):
        return provider_model, os.getenv('GROQ_API_KEY', '').strip()
    if provider_model.startswith('anthropic:'):
        return provider_model, os.getenv('ANTHROPIC_API_KEY', '').strip()
    if provider_model.startswith('deepseek:'):
        return provider_model, os.getenv('DEEPSEEK_API_KEY', os.getenv('LLM_API_KEY', '')).strip()
    return provider_model, os.getenv('LLM_PROVIDER_API_KEY', '').strip()


def _post_chat_completion(base_url, api_key, model, messages):
    if not api_key:
        raise RuntimeError('未检测到LLM_API_KEY环境变量')
    url = base_url.rstrip('/') + '/chat/completions'
    payload = {
        'model': model,
        'messages': messages,
        'temperature': 0.6
    }
    if model.startswith('deepseek-v4'):
        thinking = os.getenv('LLM_THINKING', '').strip().lower()
        reasoning_effort = os.getenv('LLM_REASONING_EFFORT', '').strip()
        if thinking in {'enabled', 'disabled'}:
            payload['thinking'] = {'type': thinking}
        if reasoning_effort:
            payload['reasoning_effort'] = reasoning_effort
    data = json.dumps(payload, ensure_ascii=False).encode('utf-8')
    req = urllib.request.Request(url, data=data, method='POST')
    req.add_header('Content-Type', 'application/json')
    req.add_header('Authorization', f'Bearer {api_key}')
    opener = _build_proxy_opener()
    with opener.open(req, timeout=API_TIMEOUT) as resp:
        body = resp.read().decode('utf-8')
        result = json.loads(body)
        return result['choices'][0]['message']['content']


class _HttpChatLLM(LLM):
    base_url: str
    api_key: str
    model: str
    timeout: int = API_TIMEOUT

    @property
    def _llm_type(self) -> str:
        return "http_chat"

    @property
    def _identifying_params(self):
        return {"model": self.model, "base_url": self.base_url}

    def _call(self, prompt: str, stop: Optional[List[str]] = None, run_manager=None) -> str:
        messages = [{'role': 'user', 'content': prompt}]
        response = _post_chat_completion(self.base_url, self.api_key, self.model, messages)
        if stop:
            for token in stop:
                if token in response:
                    response = response.split(token)[0]
        return response


def _extract_drug_name(query):
    patterns = [
        r'查找(.+?)的药物性质',
        r'查询(.+?)的药物性质',
        r'查找(.+?)药物性质',
        r'查询(.+?)药物性质',
        r'查找(.+?)的药物信息',
        r'查询(.+?)的药物信息',
        r'查找(.+?)药物信息',
        r'查询(.+?)药物信息',
    ]
    for pat in patterns:
        m = re.search(pat, query, re.IGNORECASE)
        if m:
            return m.group(1).strip()
    # 兜底：去掉常见前缀
    cleaned = re.sub(r'^(帮我|请|能否|可以|我要)?(查询|查找)?', '', query).strip()
    return cleaned if cleaned else None


def _get_drug_info_local(drug_name):
    try:
        NamesWithID = datareader("NamesWithID")
        Drug_Name = NamesWithID[:, 0].tolist()
        if drug_name not in Drug_Name:
            return f"未在数据库中找到药物：{drug_name}"
        flag = Drug_Name.index(drug_name)
        drug_id = NamesWithID[:, 1][flag]
        offside = OFFSIDE(flag)[:8]
        feat = drug_881feat(flag)[:8]
        smiles = get_drug_smiles(drug_name)

        lines = [
            f"药物名称：{drug_name}",
            f"药物编号：{drug_id}",
        ]
        if smiles:
            lines.append(f"SMILES：{smiles}")
        if offside:
            lines.append("副作用（部分）： " + "，".join(offside))
        if feat:
            lines.append("化学特征（部分）： " + "，".join(feat))
        return "\n".join(lines)
    except Exception as e:
        return f"本地查询失败：{e}"


def _handle_local_request(query):
    keywords = ['药物性质', '药物信息', '副作用', '化学特征', 'SMILES',
                '药物', '蛋白', '靶标', '相互作用', '亲和力', '搜索', '查找', '查询',
                '关系', '预测', '结合']
    analysis_keywords = ['证据报告', '可解释报告', '分析报告', '符号子图', '证据包', '溯源报告']
    if any(k in query for k in analysis_keywords):
        return _run_analysis_pipeline(query)
    if any(k in query for k in keywords):
        name = _extract_drug_name(query)
        if name and name in _get_drug_name_candidates():
            return _get_drug_info_local(name)
        # Fallback to search agent for broader queries
        try:
            agent = get_search_agent()
            response = agent.search(query)
            if response.results:
                return response.summary
        except Exception:
            pass
    return None


@tool
def _langchain_local_tool(query: str) -> str:
    """本地药物数据库查询。输入药物相关问题，返回药物性质/副作用/化学特征等信息。"""
    result = _handle_local_request(query)
    return result or "本地数据库未命中该问题。"


PROJECT_TOOL_REGISTRY = {
    "drug_info": _get_drug_info_local,
    "open_drug_detail": None,
    "open_drug_drug_analysis": None,
    "open_dta_predict": None,
    "open_relation_graph": None,
    "analyze_query": None,
    # Search agent tools — registered lazily via _search_project_tool
    "search": None,
    "search_drug": None,
    "search_protein": None,
    "search_interaction": None,
}


def _search_project_tool(action: str, payload: str) -> str:
    agent = get_search_agent()
    if action == "search":
        return agent.search(payload).summary
    elif action == "search_drug":
        drugs = agent.search_drug(payload, include_side_effects=True, include_features=True)
        return agent._format_drug_summary(drugs)
    elif action == "search_protein":
        prots = agent.search_protein(payload)
        return agent._format_protein_summary(prots)
    elif action == "search_interaction":
        parts = [p.strip() for p in re.split(r'[,，\s|;；]+', payload) if p.strip()]
        if len(parts) < 2:
            return "需要两个实体名称，例如：阿司匹林, 布洛芬"
        result = agent.search_interaction(parts[0], parts[1])
        if result:
            return f"{result.entity1} ↔ {result.entity2} ({result.interaction_type}): {result.result}"
        return f"未找到 {parts[0]} 与 {parts[1]} 之间的已知相互作用。"
    return f"不支持的搜索 action：{action}"


def _run_analysis_pipeline(payload: str) -> str:
    if not payload:
        return "缺少分析问题，例如：阿司匹林和布洛芬的相互作用并生成证据报告"
    try:
        agent = get_search_agent()
        report = analysis_pipeline.analyze_query(agent, payload, export=True)
        return analysis_pipeline.format_report_for_chat(report)
    except Exception as exc:
        return f"分析流水线执行失败：{exc}"


@tool
def project_action(action: str, payload: str = "") -> str:
    """
    调用项目内置功能的统一入口。
    action: 功能名（例如：drug_info）
    payload: 对应功能的输入（例如药物名称）
    """
    func = PROJECT_TOOL_REGISTRY.get(action)
    if func is None and action.startswith("search"):
        # Lazy-dispatch search agent actions
        if not payload:
            return f"缺少 payload，请提供 action={action} 的输入参数。"
        try:
            return _search_project_tool(action, payload)
        except Exception as exc:
            return f"搜索失败：{exc}"
    if not func:
        return f"不支持的 action：{action}"
    if not payload:
        return f"缺少 payload，请提供 action={action} 的输入参数。"
    try:
        return func(payload)
    except Exception as exc:
        return f"功能调用失败：{exc}"


def _open_drug_detail_window(drug_name: str) -> str:
    if not main_root:
        return "主窗口未初始化，无法打开详情窗口。"
    try:
        main_root.after(0, lambda: drug_inquiry(drug_name))
        return f"已打开药物详情窗口：{drug_name}"
    except Exception as exc:
        return f"打开窗口失败：{exc}"


PROJECT_TOOL_REGISTRY["open_drug_detail"] = _open_drug_detail_window


def _split_pair(payload: str):
    if not payload:
        return None, None
    for sep in [",", "，", " ", "|", ";", "；"]:
        if sep in payload:
            parts = [p.strip() for p in payload.split(sep) if p.strip()]
            if len(parts) >= 2:
                return parts[0], parts[1]
    return None, None


def _open_drug_drug_window(payload: str) -> str:
    if not main_root:
        return "主窗口未初始化，无法打开分析窗口。"
    drug1, drug2 = _split_pair(payload)
    if not drug1 or not drug2:
        return "需要两个药物名称，例如：阿司匹林, 氯吡格雷"
    try:
        main_root.after(0, lambda: drug_drug(main_root))
        main_root.after(0, lambda: drug_drug_inquiry(drug1, drug2))
        return f"已打开药物-药物分析窗口：{drug1} vs {drug2}"
    except Exception as exc:
        return f"打开窗口失败：{exc}"


def _open_dta_predict_window(payload: str) -> str:
    if not main_root:
        return "主窗口未初始化，无法打开预测窗口。"
    drug_name, protein_name = _split_pair(payload)
    if not drug_name or not protein_name:
        return "需要药物与靶标名称，例如：吉非替尼, EGFR"
    try:
        main_root.after(0, lambda: dta_predict(main_root))
        main_root.after(0, lambda: drug_prot_dta_predict(drug_name, protein_name))
        return f"已打开DTA预测窗口：{drug_name} vs {protein_name}"
    except Exception as exc:
        return f"打开窗口失败：{exc}"


def _open_relation_graph(payload: str) -> str:
    if not main_root:
        return "主窗口未初始化，无法打开关系图。"
    drug_name = payload.strip() if payload else ""
    if not drug_name:
        return "请提供药物名称，例如：阿司匹林"
    try:
        # 复用查询流程，先拿到药物编号与索引
        names = datareader("NamesWithID")
        drug_names = names[:, 0].tolist()
        if drug_name not in drug_names:
            return f"未在数据库中找到药物：{drug_name}"
        idx = drug_names.index(drug_name)
        drug_id = names[:, 1][idx]
        main_root.after(0, lambda: signed_graph(drug_id, idx))
        return f"已打开关系图：{drug_name}"
    except Exception as exc:
        return f"打开关系图失败：{exc}"


PROJECT_TOOL_REGISTRY["open_drug_drug_analysis"] = _open_drug_drug_window
PROJECT_TOOL_REGISTRY["open_dta_predict"] = _open_dta_predict_window
PROJECT_TOOL_REGISTRY["open_relation_graph"] = _open_relation_graph
PROJECT_TOOL_REGISTRY["analyze_query"] = _run_analysis_pipeline


def _build_llm():
    """构建LangChain模型，优先使用init_chat_model，其次回退到HTTP LLM"""
    provider_model = os.getenv('LLM_PROVIDER_MODEL', '').strip()
    if provider_model and ':' in provider_model:
        try:
            from langchain.chat_models import init_chat_model
            return init_chat_model(provider_model)
        except Exception:
            pass
    base_url, api_key, model = _get_llm_config()
    if not api_key:
        raise RuntimeError("未检测到LLM_API_KEY环境变量")
    return _HttpChatLLM(base_url=base_url, api_key=api_key, model=model, timeout=API_TIMEOUT)


def _build_langchain_agent(system_prompt):
    """构建LangChain Agent（LangChain 1.0 create_agent风格）"""
    if not LANGCHAIN_AVAILABLE:
        raise RuntimeError("LangChain未安装，请先安装 langchain")
    llm = _build_llm()
    return create_agent(
        model=llm,
        tools=[_langchain_local_tool, project_action],
        system_prompt=system_prompt
    )


def _run_langchain_agent(messages, system_prompt):
    """调用LangChain Agent，返回最终回复内容"""
    agent = _build_langchain_agent(system_prompt)
    result = agent.invoke({"messages": messages})
    if isinstance(result, dict) and "messages" in result and result["messages"]:
        last = result["messages"][-1]
        content = getattr(last, "content", None)
        if content:
            return content
    if isinstance(result, str):
        return result
    return str(result)


def _run_simple_chat_completion(messages):
    base_url, api_key, model = _get_llm_config()
    if not api_key:
        raise RuntimeError('未检测到LLM_API_KEY环境变量')
    return _post_chat_completion(base_url, api_key, model, messages)


def _run_agent_reply(messages, system_prompt):
    if LANGCHAIN_AVAILABLE:
        try:
            return _run_langchain_agent(messages, system_prompt)
        except Exception as exc:
            # LangChain异常时回退到HTTP直连
            try:
                return _run_simple_chat_completion(messages)
            except Exception:
                raise exc
    return _run_simple_chat_completion(messages)


def _test_llm_api():
    base_url, api_key, model = _get_llm_config()
    if not api_key:
        return False, '未检测到LLM_API_KEY环境变量'
    try:
        reply = _post_chat_completion(
            base_url,
            api_key,
            model,
            [{'role': 'user', 'content': 'ping'}]
        )
        return True, f'连接成功，模型返回: {reply[:80]}'
    except Exception as e:
        return False, f'连接失败：{e}'


def _build_models_url(base_url):
    base_url = base_url.rstrip('/')
    if 'api.deepseek.com' in base_url:
        return base_url + '/models'
    if not base_url.endswith('/v1'):
        base_url = base_url + '/v1'
    return base_url + '/models'


def _test_network_connectivity():
    base_url, _, _ = _get_llm_config()
    url = _build_models_url(base_url)
    try:
        req = urllib.request.Request(url, method='GET')
        opener = _build_proxy_opener()
        with opener.open(req, timeout=API_TIMEOUT) as resp:
            return True, f'网络连通，状态码：{resp.status}'
    except urllib.error.HTTPError as e:
        if e.code in (401, 403):
            return True, f'网络连通，服务返回认证状态码：{e.code}'
        return False, f'服务返回错误状态码：{e.code}'
    except Exception as e:
        return False, f'网络不可达：{e}'


def _mask_key(key):
    if not key:
        return 'EMPTY'
    if len(key) <= 8:
        return key[0] + '***' + key[-1]
    return f"{key[:4]}***{key[-4:]}"


def _get_proxy_env():
    http_proxy = os.getenv('HTTP_PROXY', '') or os.getenv('http_proxy', '')
    https_proxy = os.getenv('HTTPS_PROXY', '') or os.getenv('https_proxy', '')
    return http_proxy, https_proxy


def _build_proxy_opener():
    http_proxy, https_proxy = _get_proxy_env()
    proxies = {}
    if http_proxy:
        proxies['http'] = http_proxy
    if https_proxy:
        proxies['https'] = https_proxy
    if proxies:
        handler = urllib.request.ProxyHandler(proxies)
        return urllib.request.build_opener(handler)
    return urllib.request.build_opener()


def open_agent_window_legacy(parent):
    win = tkinter.Toplevel(parent)
    win.title('智能助手')
    win.geometry('820x620')
    win.minsize(760, 560)
    win.config(bg='#F5F5F5')
    win.transient(parent)
    win.grab_set()
    win.resizable(True, True)

    title_label = tkinter.Label(
        win,
        text='LangChain 智能助手',
        font=(UI_FONT_FAMILY, 16, 'bold'),
        bg='#F5F5F5',
        fg='#333333'
    )
    title_label.pack(pady=12)

    input_frame = tkinter.Frame(win, bg='#F5F5F5')
    input_frame.pack(side=tkinter.BOTTOM, fill=tkinter.X, padx=15, pady=10)
    input_frame.grid_columnconfigure(0, weight=1)
    for col in range(1, 5):
        input_frame.grid_columnconfigure(col, weight=0)

    status_bar = tkinter.Frame(win, bg='#F5F5F5')
    status_bar.pack(side=tkinter.BOTTOM, fill=tkinter.X, padx=15, pady=(0, 6))

    input_hint = tkinter.Label(
        status_bar,
        text='在下方输入问题并点击发送',
        font=(UI_FONT_FAMILY, 10),
        bg='#F5F5F5',
        fg='#666666'
    )
    input_hint.pack(side=tkinter.LEFT)

    base_url, api_key, model = _get_llm_config()
    provider_model, provider_key = _get_provider_key_state()
    key_state = 'OK' if api_key else '缺失'
    provider_state = 'OK' if provider_key else '缺失'
    agent_lib = 'OK' if LANGCHAIN_AVAILABLE else '无'
    agent_mode = 'LangChain' if LANGCHAIN_AVAILABLE else 'HTTP'
    if provider_model:
        status_text = f'状态：就绪 | Provider:{provider_model} | PKEY:{provider_state} | LangChain:{agent_lib} | Mode:{agent_mode}'
    else:
        status_text = f'状态：就绪 | KEY:{key_state} | LangChain:{agent_lib} | Mode:{agent_mode}'
    status_var = tkinter.StringVar(value=status_text)
    status_label = tkinter.Label(
        status_bar,
        textvariable=status_var,
        font=(UI_FONT_FAMILY, 10),
        bg='#F5F5F5',
        fg='#4B5563'
    )
    status_label.pack(side=tkinter.RIGHT)

    # 调试：强制显示状态栏（避免被布局挤掉）
    status_bar.after(500, lambda: status_bar.pack(side=tkinter.BOTTOM, fill=tkinter.X, padx=15, pady=(0, 6)))

    chat_frame = tkinter.Frame(win, bg='#F5F5F5')
    chat_frame.pack(side=tkinter.TOP, padx=15, pady=5, fill=tkinter.BOTH, expand=True)

    chat_text = tkinter.Text(
        chat_frame,
        font=(UI_FONT_FAMILY, 11),
        bg='white',
        fg='#333333',
        relief=tkinter.SOLID,
        borderwidth=1,
        wrap=tkinter.WORD
    )
    chat_text.pack(side=tkinter.LEFT, fill=tkinter.BOTH, expand=True)
    chat_text.config(state=tkinter.DISABLED)

    scrollbar = tkinter.Scrollbar(chat_frame, command=chat_text.yview)
    scrollbar.pack(side=tkinter.RIGHT, fill=tkinter.Y)
    chat_text.config(yscrollcommand=scrollbar.set)
    input_var = tkinter.StringVar()
    input_entry = tkinter.Entry(
        input_frame,
        textvariable=input_var,
        font=(UI_FONT_FAMILY, 12),
        relief=tkinter.SOLID,
        borderwidth=1
    )
    input_entry.grid(row=0, column=0, padx=(0, 10), pady=(0, 8), ipady=6, sticky='ew')

    system_prompt = (
        '你是药物数据分析助手。可以基于用户问题提供查询建议、解释结果、'
        '以及给出可执行的下一步。对药物性质/药物信息/副作用/化学特征等问题，'
        '优先结合本地数据库信息回答。你可以在需要时调用工具：'
        'project_action(action, payload)，目前支持 action='
        'drug_info、open_drug_detail、open_drug_drug_analysis、open_dta_predict、open_relation_graph、'
        'analyze_query（自然语言解析-证据检索-符号子图-可解释报告导出）、'
        'search（统一搜索）、search_drug（药物搜索）、search_protein（蛋白质搜索）、search_interaction（相互作用搜索）。'
        '\n搜索工具说明：'
        '\n- search: 支持自然语言查询，如"查找头痛相关的药物"、"阿司匹林和布洛芬的相互作用"'
        '\n- search_drug: 按药物名/副作用/化学特征搜索药物'
        '\n- search_protein: 按蛋白名/基因名/别名搜索蛋白质'
        '\n- search_interaction: 查询两个实体（药物/蛋白）间的 DDI/DTA/PPI 相互作用'
        '\n- analyze_query: 生成结构化结论、证据来源、符号子图和Markdown/JSON报告'
    )
    messages = [{'role': 'system', 'content': system_prompt}]

    def append_message(role, content):
        chat_text.config(state=tkinter.NORMAL)
        prefix = '你：' if role == 'user' else '助手：'
        chat_text.insert(tkinter.END, f'{prefix}{content}\n\n')
        chat_text.see(tkinter.END)
        chat_text.config(state=tkinter.DISABLED)

    def append_assistant_start():
        chat_text.config(state=tkinter.NORMAL)
        chat_text.insert(tkinter.END, '助手：')
        chat_text.see(tkinter.END)
        chat_text.config(state=tkinter.DISABLED)

    def append_assistant_delta(delta):
        chat_text.config(state=tkinter.NORMAL)
        chat_text.insert(tkinter.END, delta)
        chat_text.see(tkinter.END)
        chat_text.config(state=tkinter.DISABLED)

    def append_assistant_end():
        chat_text.config(state=tkinter.NORMAL)
        chat_text.insert(tkinter.END, '\n\n')
        chat_text.see(tkinter.END)
        chat_text.config(state=tkinter.DISABLED)

    def send_message():
        try:
            text = input_var.get().strip()
            if not text:
                status_var.set('状态：请输入内容')
                return
            status_var.set('状态：已触发发送...')
            input_var.set('')
            append_message('user', text)
            messages.append({'role': 'user', 'content': text})
            send_btn.config(state=tkinter.DISABLED)
            status_var.set('状态：发送中...')

            local_result = _handle_local_request(text)
            if local_result:
                append_message('assistant', local_result)
                messages.append({'role': 'assistant', 'content': local_result})
                status_var.set('状态：完成')
                send_btn.config(state=tkinter.NORMAL)
                return
            if not LANGCHAIN_AVAILABLE and not api_key:
                append_message('assistant', 'LangChain不可用且未检测到 LLM_API_KEY，无法调用模型。')
                status_var.set('状态：缺少依赖或Key')
                send_btn.config(state=tkinter.NORMAL)
                return
            append_message('assistant', '已收到问题，正在调用模型，请稍候...')
        except Exception as e:
            status_var.set(f'状态：发送前错误 {e}')
            messagebox.showerror('发送失败', f'发送前发生错误：{e}')
            return

        def worker():
            max_retries = 2
            backoff = 1.2
            last_error = None
            for attempt in range(max_retries + 1):
                try:
                    reply = _run_agent_reply(messages, system_prompt)
                    messages.append({'role': 'assistant', 'content': reply})
                    win.after(0, lambda: append_message('assistant', reply))
                    win.after(0, lambda: status_var.set('状态：完成'))
                    win.after(0, lambda: send_btn.config(state=tkinter.NORMAL))
                    return
                except Exception as e:
                    last_error = e
                    if attempt < max_retries:
                        wait_s = backoff * (attempt + 1)
                        win.after(0, lambda: status_var.set(f'状态：失败，{wait_s:.1f}s 后重试...'))
                        threading.Event().wait(wait_s)
                        continue
            win.after(0, lambda: append_message('assistant', f'调用失败：{last_error}'))
            win.after(0, lambda: status_var.set('状态：失败，请检查API或网络'))
            win.after(0, lambda: send_btn.config(state=tkinter.NORMAL))
            return

        threading.Thread(target=worker, daemon=True).start()

    env_btn = create_button(
        input_frame,
        text='环境检测',
        command=lambda: None,
        width=8,
        height=1,
        animate=False
    )
    env_btn.grid(row=1, column=0, sticky='w', padx=(0, 8))

    proxy_btn = create_button(
        input_frame,
        text='代理检测',
        command=lambda: None,
        width=8,
        height=1,
        animate=False
    )
    proxy_btn.grid(row=1, column=1, padx=(0, 8))

    net_btn = create_button(
        input_frame,
        text='网络测试',
        command=lambda: None,
        width=8,
        height=1,
        animate=False
    )
    net_btn.grid(row=1, column=2, padx=(0, 8))

    test_btn = create_button(
        input_frame,
        text='API测试',
        command=lambda: None,
        width=8,
        height=1,
        animate=False
    )
    test_btn.grid(row=1, column=3, padx=(0, 8))

    send_btn = create_button(
        input_frame,
        text='发送',
        command=send_message,
        width=8,
        height=1
    )
    send_btn.grid(row=0, column=1, columnspan=3, pady=(0, 8), sticky='e')
    input_entry.bind('<Return>', lambda e: send_message())
    input_entry.focus()
    win.bind('<Button-1>', lambda e: input_entry.focus_set())

    def run_api_test():
        status_var.set('状态：测试中...')
        append_message('assistant', 'API测试已启动，请稍候...')
        test_btn.config(state=tkinter.DISABLED)

        done_flag = {'done': False}
        def watchdog():
            if not done_flag['done']:
                append_message('assistant', '提示：API测试已超时，请检查网络或Key')
                status_var.set('状态：API测试超时')
                test_btn.config(state=tkinter.NORMAL)
        win.after(30000, watchdog)

        def worker():
            try:
                ok, msg = _test_llm_api()
            except Exception as e:
                ok, msg = False, f'执行异常：{e}'
            def ui_update():
                done_flag['done'] = True
                if ok:
                    append_message('assistant', f'API测试成功：{msg}')
                    status_var.set('状态：API测试成功')
                else:
                    append_message('assistant', f'API测试失败：{msg}')
                    status_var.set('状态：API测试失败')
                test_btn.config(state=tkinter.NORMAL)
            try:
                win.after(0, ui_update)
            except Exception:
                pass

        threading.Thread(target=worker, daemon=True).start()

    test_btn.config(command=run_api_test)

    def run_env_check():
        base_url, api_key, model = _get_llm_config()
        provider_model, provider_key = _get_provider_key_state()
        masked = _mask_key(api_key)
        masked_provider = _mask_key(provider_key)
        lc_state = 'OK' if LANGCHAIN_AVAILABLE else '无'
        if provider_model:
            append_message(
                'assistant',
                f'环境检测：PKEY={masked_provider} | Provider={provider_model} | LangChain={lc_state}'
            )
        else:
            append_message(
                'assistant',
                f'环境检测：KEY={masked} | LangChain={lc_state} | BASE={base_url} | MODEL={model}'
            )
        status_var.set('状态：环境检测完成')

    env_btn.config(command=run_env_check)

    def run_proxy_check():
        http_proxy, https_proxy = _get_proxy_env()
        http_proxy = http_proxy or '未设置'
        https_proxy = https_proxy or '未设置'
        append_message('assistant', f'代理检测：HTTP_PROXY={http_proxy} | HTTPS_PROXY={https_proxy}')
        status_var.set('状态：代理检测完成')

    proxy_btn.config(command=run_proxy_check)

    def run_net_test():
        status_var.set('状态：网络测试中...')
        append_message('assistant', '网络测试已启动，请稍候...')
        net_btn.config(state=tkinter.DISABLED)

        done_flag = {'done': False}
        def watchdog():
            if not done_flag['done']:
                append_message('assistant', '提示：网络测试已超时，请检查网络或代理')
                status_var.set('状态：网络测试超时')
                net_btn.config(state=tkinter.NORMAL)
        win.after(30000, watchdog)

        def worker():
            ok, msg = _test_network_connectivity()
            def ui_update():
                done_flag['done'] = True
                if ok:
                    append_message('assistant', f'网络测试成功：{msg}')
                    status_var.set('状态：网络测试成功')
                else:
                    append_message('assistant', f'网络测试失败：{msg}')
                    status_var.set('状态：网络测试失败')
                net_btn.config(state=tkinter.NORMAL)
            try:
                win.after(0, ui_update)
            except Exception:
                pass

        threading.Thread(target=worker, daemon=True).start()

    net_btn.config(command=run_net_test)


def open_agent_window(parent):
    global agent_window_ref
    try:
        if agent_window_ref is not None and agent_window_ref.winfo_exists():
            agent_window_ref.lift()
            agent_window_ref.focus_force()
            return
    except Exception:
        agent_window_ref = None

    win = tkinter.Toplevel(parent)
    agent_window_ref = win
    win.title('智能助手')
    setup_tk_fonts(win, base_size=12)
    win.geometry('1040x760')
    win.minsize(920, 660)
    win.config(bg='#ECF2F8')
    win.resizable(True, True)
    try:
        win.lift()
        win.focus_force()
    except Exception:
        pass

    bg_color = '#ECF2F8'
    panel_bg = '#FFFFFF'
    soft_bg = '#F7FAFC'
    border = '#D9E2EC'
    text_color = '#1F2937'
    muted = '#64748B'
    primary = '#2563EB'
    primary_hover = '#1D4ED8'
    success = '#059669'
    danger = '#DC2626'
    ui_font = UI_FONT_FAMILY

    base_url, api_key, model = _get_llm_config()
    provider_model, provider_key = _get_provider_key_state()
    key_state = '已配置' if api_key else '缺失'
    provider_state = '已配置' if provider_key else '缺失'
    agent_lib = '可用' if LANGCHAIN_AVAILABLE else '未安装'
    agent_mode = 'LangChain + HTTP回退' if LANGCHAIN_AVAILABLE else 'HTTP直连'
    if provider_model:
        status_text = f'Provider {provider_model} | PKEY {provider_state} | LangChain {agent_lib}'
    else:
        status_text = f'{model} | KEY {key_state} | {agent_mode}'
    status_var = tkinter.StringVar(value=f'就绪 · {status_text}')
    ui_queue = queue.Queue()
    window_alive = {'value': True}

    def ui_call(func, *args, **kwargs):
        ui_queue.put((func, args, kwargs))

    def process_ui_queue():
        if not window_alive['value']:
            return
        try:
            while True:
                func, args, kwargs = ui_queue.get_nowait()
                try:
                    func(*args, **kwargs)
                except Exception:
                    pass
        except queue.Empty:
            pass
        try:
            win.after(50, process_ui_queue)
        except Exception:
            window_alive['value'] = False

    def on_close():
        global agent_window_ref
        window_alive['value'] = False
        try:
            win.destroy()
        except Exception:
            pass
        agent_window_ref = None

    win.protocol('WM_DELETE_WINDOW', on_close)
    win.bind('<Escape>', lambda _event: on_close())

    shell = tkinter.Frame(win, bg=bg_color)
    shell.pack(fill=tkinter.BOTH, expand=True)
    shell.grid_rowconfigure(1, weight=1)
    shell.grid_columnconfigure(0, weight=1)

    header = tkinter.Frame(shell, bg=panel_bg, highlightbackground=border, highlightthickness=1)
    header.grid(row=0, column=0, sticky='ew', padx=18, pady=(18, 10))
    header.grid_columnconfigure(0, weight=1)

    title_block = tkinter.Frame(header, bg=panel_bg)
    title_block.grid(row=0, column=0, sticky='w', padx=18, pady=14)
    create_button(
        title_block,
        text='返回',
        command=on_close,
        style='secondary',
        font=(ui_font, 12),
        bg='#EEF2F7',
        fg='#334155',
        activebackground='#E2E8F0',
        width=7,
        height=1,
        padx=6,
        pady=3,
        animate=False
    ).pack(anchor='w', pady=(0, 8))
    tkinter.Label(
        title_block,
        text=f'{Config.SOFTWARE_NAME} 智能助手',
        font=(ui_font, 22, 'bold'),
        bg=panel_bg,
        fg=text_color
    ).pack(anchor='w')
    tkinter.Label(
        title_block,
        text='自然语言解析、符号关系推理、证据报告与 API 诊断',
        font=(ui_font, 12),
        bg=panel_bg,
        fg=muted
    ).pack(anchor='w', pady=(3, 0))

    toolbar = tkinter.Frame(header, bg=panel_bg)
    toolbar.grid(row=0, column=1, sticky='e', padx=18, pady=14)

    body = tkinter.Frame(shell, bg=bg_color)
    body.grid(row=1, column=0, sticky='nsew', padx=18, pady=(0, 10))
    body.grid_rowconfigure(0, weight=1)
    body.grid_columnconfigure(0, weight=1)

    chat_card = tkinter.Frame(body, bg=panel_bg, highlightbackground=border, highlightthickness=1)
    chat_card.grid(row=0, column=0, sticky='nsew')
    chat_card.grid_rowconfigure(0, weight=1)
    chat_card.grid_columnconfigure(0, weight=1)

    chat_text = tkinter.Text(
        chat_card,
        font=(ui_font, 13),
        bg=panel_bg,
        fg=text_color,
        relief=tkinter.FLAT,
        borderwidth=0,
        wrap=tkinter.WORD,
        padx=18,
        pady=16,
        insertbackground=text_color,
        spacing1=4,
        spacing3=10
    )
    chat_text.grid(row=0, column=0, sticky='nsew')
    scrollbar = tkinter.Scrollbar(chat_card, command=chat_text.yview)
    scrollbar.grid(row=0, column=1, sticky='ns')
    chat_text.config(yscrollcommand=scrollbar.set)
    chat_text.tag_configure('user_name', foreground=primary, font=(ui_font, 12, 'bold'), spacing1=8)
    chat_text.tag_configure('assistant_name', foreground=success, font=(ui_font, 12, 'bold'), spacing1=8)
    chat_text.tag_configure('system_name', foreground=muted, font=(ui_font, 12, 'bold'), spacing1=8)
    chat_text.tag_configure('user_body', foreground='#111827', lmargin1=16, lmargin2=16, rmargin=70, spacing3=12)
    chat_text.tag_configure('assistant_body', foreground='#111827', lmargin1=16, lmargin2=16, rmargin=40, spacing3=12)
    chat_text.tag_configure('system_body', foreground=muted, lmargin1=16, lmargin2=16, rmargin=40, spacing3=12)
    chat_text.config(state=tkinter.DISABLED)

    composer = tkinter.Frame(shell, bg=panel_bg, highlightbackground=border, highlightthickness=1)
    composer.grid(row=2, column=0, sticky='ew', padx=18, pady=(0, 18))
    composer.grid_columnconfigure(0, weight=1)

    quick_row = tkinter.Frame(composer, bg=panel_bg)
    quick_row.grid(row=0, column=0, columnspan=2, sticky='ew', padx=14, pady=(12, 6))

    input_text = tkinter.Text(
        composer,
        height=3,
        font=(ui_font, 13),
        bg='#F8FAFC',
        fg=text_color,
        relief=tkinter.FLAT,
        borderwidth=0,
        wrap=tkinter.WORD,
        padx=12,
        pady=10,
        insertbackground=text_color
    )
    input_text.grid(row=1, column=0, sticky='ew', padx=(14, 10), pady=(0, 12))
    input_text.focus_set()

    send_btn = create_button(
        composer,
        text='发送',
        command=lambda: send_message(),
        font=(ui_font, 13, 'bold'),
        bg=primary,
        activebackground=primary_hover,
        width=7,
        height=1,
        padx=8,
        pady=5
    )
    send_btn.grid(row=1, column=1, sticky='se', padx=(0, 14), pady=(0, 12))

    status_bar = tkinter.Frame(shell, bg=bg_color)
    status_bar.grid(row=3, column=0, sticky='ew', padx=20, pady=(0, 12))
    status_bar.grid_columnconfigure(0, weight=1)
    tkinter.Label(
        status_bar,
        textvariable=status_var,
        font=(ui_font, 11),
        bg=bg_color,
        fg=muted,
        anchor='w'
    ).grid(row=0, column=0, sticky='w')
    tkinter.Label(
        status_bar,
        text='Enter 发送 · Shift+Enter 换行',
        font=(ui_font, 11),
        bg=bg_color,
        fg=muted,
        anchor='e'
    ).grid(row=0, column=1, sticky='e')

    system_prompt = (
        '你是药物数据分析助手。可以基于用户问题提供查询建议、解释结果、'
        '以及给出可执行的下一步。对药物性质/药物信息/副作用/化学特征等问题，'
        '优先结合本地数据库信息回答。你可以在需要时调用工具：'
        'project_action(action, payload)，目前支持 action='
        'drug_info、open_drug_detail、open_drug_drug_analysis、open_dta_predict、open_relation_graph、'
        'analyze_query（自然语言解析-证据检索-符号子图-可解释报告导出）、'
        'search（统一搜索）、search_drug（药物搜索）、search_protein（蛋白质搜索）、search_interaction（相互作用搜索）。'
        '\n搜索工具说明：'
        '\n- search: 支持自然语言查询，如"查找头痛相关的药物"、"阿司匹林和布洛芬的相互作用"'
        '\n- search_drug: 按药物名/副作用/化学特征搜索药物'
        '\n- search_protein: 按蛋白名/基因名/别名搜索蛋白质'
        '\n- search_interaction: 查询两个实体（药物/蛋白）间的 DDI/DTA/PPI 相互作用'
        '\n- analyze_query: 生成结构化结论、证据来源、符号子图和Markdown/JSON报告'
    )
    messages = [{'role': 'system', 'content': system_prompt}]

    def append_message(role, content):
        content = str(content).strip()
        if not content:
            return
        chat_text.config(state=tkinter.NORMAL)
        if role == 'user':
            chat_text.insert(tkinter.END, '你\n', 'user_name')
            chat_text.insert(tkinter.END, content + '\n\n', 'user_body')
        elif role == 'system':
            chat_text.insert(tkinter.END, '系统\n', 'system_name')
            chat_text.insert(tkinter.END, content + '\n\n', 'system_body')
        else:
            chat_text.insert(tkinter.END, '助手\n', 'assistant_name')
            chat_text.insert(tkinter.END, content + '\n\n', 'assistant_body')
        chat_text.see(tkinter.END)
        chat_text.config(state=tkinter.DISABLED)

    def set_busy(is_busy, text=''):
        send_btn.config(state=tkinter.DISABLED if is_busy else tkinter.NORMAL)
        input_text.config(state=tkinter.DISABLED if is_busy else tkinter.NORMAL)
        if text:
            status_var.set(text)
        if not is_busy:
            input_text.focus_set()

    def get_input():
        return input_text.get('1.0', tkinter.END).strip()

    def clear_input():
        input_text.delete('1.0', tkinter.END)

    def send_message():
        text = get_input()
        if not text:
            status_var.set('请输入内容')
            input_text.focus_set()
            return
        clear_input()
        append_message('user', text)
        messages.append({'role': 'user', 'content': text})
        set_busy(True, '正在处理...')

        local_result = _handle_local_request(text)
        if local_result:
            append_message('assistant', local_result)
            messages.append({'role': 'assistant', 'content': local_result})
            set_busy(False, '完成 · 已使用本地工具')
            return

        if not LANGCHAIN_AVAILABLE and not api_key:
            append_message('assistant', '未检测到 LLM_API_KEY，且 LangChain 不可用。请检查 .env 或点击环境检测。')
            set_busy(False, '缺少 API 配置')
            return

        append_message('system', '正在调用模型，请稍候...')

        def worker():
            max_retries = 1
            last_error = None
            for attempt in range(max_retries + 1):
                try:
                    reply = _run_agent_reply(messages, system_prompt)
                    messages.append({'role': 'assistant', 'content': reply})
                    ui_call(append_message, 'assistant', reply)
                    ui_call(set_busy, False, '完成 · 模型已响应')
                    return
                except Exception as exc:
                    last_error = exc
                    if attempt < max_retries:
                        ui_call(status_var.set, '调用失败，正在重试...')
                        threading.Event().wait(1.2)
            ui_call(append_message, 'assistant', f'调用失败：{last_error}')
            ui_call(set_busy, False, '调用失败 · 请检查 API、网络或模型名')

        threading.Thread(target=worker, daemon=True).start()

    def insert_quick(text):
        input_text.delete('1.0', tkinter.END)
        input_text.insert('1.0', text)
        input_text.focus_set()

    def bind_enter(event):
        if event.state & 0x0001:
            return None
        send_message()
        return 'break'

    input_text.bind('<Return>', bind_enter)

    examples = [
        '查询 Aspirin 的药物信息',
        'Aspirin 和 Warfarin 的相互作用',
        'EGFR 相关靶标信息',
        'Aspirin和Warfarin的相互作用并生成证据报告',
    ]
    for item in examples:
        create_button(
            quick_row,
            text=item,
            command=lambda value=item: insert_quick(value),
            style='secondary',
            font=(ui_font, 11),
            bg='#EEF2F7',
            fg='#334155',
            activebackground='#E2E8F0',
            width=None,
            height=1,
            padx=8,
            pady=3,
            animate=False
        ).pack(side=tkinter.LEFT, padx=(0, 8))

    def run_api_test():
        status_var.set('API 测试中...')
        append_message('system', 'API 测试已启动。')
        api_btn.config(state=tkinter.DISABLED)

        def worker():
            try:
                ok, msg = _test_llm_api()
            except Exception as exc:
                ok, msg = False, f'执行异常：{exc}'
            def ui_update():
                append_message('assistant', ('API测试成功：' if ok else 'API测试失败：') + msg)
                status_var.set('API 测试成功' if ok else 'API 测试失败')
                api_btn.config(state=tkinter.NORMAL)
            ui_call(ui_update)

        threading.Thread(target=worker, daemon=True).start()

    def run_env_check():
        base_url_now, api_key_now, model_now = _get_llm_config()
        provider_model_now, provider_key_now = _get_provider_key_state()
        lc_state = '可用' if LANGCHAIN_AVAILABLE else '未安装'
        if provider_model_now:
            msg = f'Provider={provider_model_now}\nPKEY={_mask_key(provider_key_now)}\nLangChain={lc_state}'
        else:
            msg = (
                f'BASE={base_url_now}\n'
                f'MODEL={model_now}\n'
                f'KEY={_mask_key(api_key_now)}\n'
                f'LangChain={lc_state}\n'
                f'HTTP_PROXY={_get_proxy_env()[0] or "未设置"}\n'
                f'HTTPS_PROXY={_get_proxy_env()[1] or "未设置"}'
            )
        append_message('assistant', '环境检测：\n' + msg)
        status_var.set('环境检测完成')

    def run_net_test():
        status_var.set('网络测试中...')
        append_message('system', '网络测试已启动。')
        net_btn.config(state=tkinter.DISABLED)

        def worker():
            ok, msg = _test_network_connectivity()
            def ui_update():
                append_message('assistant', ('网络测试成功：' if ok else '网络测试失败：') + msg)
                status_var.set('网络测试成功' if ok else '网络测试失败')
                net_btn.config(state=tkinter.NORMAL)
            ui_call(ui_update)

        threading.Thread(target=worker, daemon=True).start()

    env_btn = create_button(
        toolbar,
        text='环境检测',
        command=run_env_check,
        style='secondary',
        font=(ui_font, 12),
        bg='#EEF2F7',
        fg='#334155',
        activebackground='#E2E8F0',
        width=8,
        height=1,
        padx=6,
        pady=4,
        animate=False
    )
    env_btn.pack(side=tkinter.LEFT, padx=(0, 8))

    net_btn = create_button(
        toolbar,
        text='网络测试',
        command=run_net_test,
        style='secondary',
        font=(ui_font, 12),
        bg='#EEF2F7',
        fg='#334155',
        activebackground='#E2E8F0',
        width=8,
        height=1,
        padx=6,
        pady=4,
        animate=False
    )
    net_btn.pack(side=tkinter.LEFT, padx=(0, 8))

    api_btn = create_button(
        toolbar,
        text='API测试',
        command=run_api_test,
        font=(ui_font, 12, 'bold'),
        bg=primary,
        activebackground=primary_hover,
        width=7,
        height=1,
        padx=6,
        pady=4,
        animate=False
    )
    api_btn.pack(side=tkinter.LEFT)

    append_message(
        'assistant',
        '你好，我可以帮你查询药物/靶标信息、分析 DDI/DTA/PPI 关系，也可以生成证据报告。'
        '\n先点右上角“环境检测”或“API测试”可以确认模型配置是否可用。'
    )
    process_ui_queue()


##########四个板块的鼠标响应##########
def Mouse_Click_drug(event, drug_num):  # 关联鼠标点击事件
    drug_x, drug_y, drug_l, drug_h, drug_r = 380, 320, 200, 60, 5
    if drug_x <= event.x <= drug_x + drug_l and drug_y <= event.y <= drug_y + drug_h:  # 响应的位置
        drug_inquiry(drug_num)


def Mouse_over_drug(event, canvas_drug, r1, r2, t):  # 关联鼠标经过事件
    drug_x, drug_y, drug_l, drug_h, drug_r = 380, 320, 200, 60, 5
    if drug_x <= event.x <= drug_x + drug_l and drug_y <= event.y <= drug_y + drug_h:  # 响应的位置
        canvas_drug.itemconfigure(r1, outline='black')  # 重设外框颜色
        canvas_drug.itemconfigure(r2, outline='black')  # 重设内框颜色
        canvas_drug.itemconfigure(t, fill='black')  # 重设显示文本颜色
        canvas_drug.configure(cursor='hand2')  # 重设鼠标样式
    else:
        canvas_drug.itemconfigure(r1, outline='white')  # 重设外框颜色
        canvas_drug.itemconfigure(r2, outline='white')  # 重设内框颜色
        canvas_drug.itemconfigure(t, fill='white')  # 重设显示文本颜色
        canvas_drug.configure(cursor='arrow')  # 重设鼠标样式


def Mouse_Click_drug_drug(event, drug_num1, drug_num2):  # 关联鼠标点击事件
    drug_drug_x, drug_drug_y, drug_drug_l, drug_drug_h, drug_drug_r = 380, 360, 200, 60, 5
    if drug_drug_x <= event.x <= drug_drug_x + drug_drug_l and drug_drug_y <= event.y <= drug_drug_y + drug_drug_h:  # 响应的位置
        drug_drug_inquiry(drug_num1, drug_num2)


def Mouse_over_drug_drug(event, canvas_drug_drug, r1, r2, t):  # 关联鼠标经过事件
    drug_drug_x, drug_drug_y, drug_drug_l, drug_drug_h, drug_drug_r = 380, 360, 200, 60, 5
    if drug_drug_x <= event.x <= drug_drug_x + drug_drug_l and drug_drug_y <= event.y <= drug_drug_y + drug_drug_h:  # 响应的位置
        canvas_drug_drug.itemconfigure(r1, outline='black')  # 重设外框颜色
        canvas_drug_drug.itemconfigure(r2, outline='black')  # 重设内框颜色
        canvas_drug_drug.itemconfigure(t, fill='black')  # 重设显示文本颜色
        canvas_drug_drug.configure(cursor='hand2')  # 重设鼠标样式
    else:
        canvas_drug_drug.itemconfigure(r1, outline='white')  # 重设外框颜色
        canvas_drug_drug.itemconfigure(r2, outline='white')  # 重设内框颜色
        canvas_drug_drug.itemconfigure(t, fill='white')  # 重设显示文本颜色
        canvas_drug_drug.configure(cursor='arrow')  # 重设鼠标样式


def Mouse_Click_drug_prot(event, drug_num1, drug_num2):  # 关联鼠标点击事件
    drug_drug_x, drug_drug_y, drug_drug_l, drug_drug_h, drug_drug_r = 380, 360, 200, 60, 5
    if drug_drug_x <= event.x <= drug_drug_x + drug_drug_l and drug_drug_y <= event.y <= drug_drug_y + drug_drug_h:  # 响应的位置
        drug_prot_dta_predict(drug_num1, drug_num2)


def Mouse_over_drug_prot(event, canvas_drug_drug, r1, r2, t):  # 关联鼠标经过事件
    drug_drug_x, drug_drug_y, drug_drug_l, drug_drug_h, drug_drug_r = 380, 360, 200, 60, 5
    if drug_drug_x <= event.x <= drug_drug_x + drug_drug_l and drug_drug_y <= event.y <= drug_drug_y + drug_drug_h:  # 响应的位置
        canvas_drug_drug.itemconfigure(r1, outline='black')  # 重设外框颜色
        canvas_drug_drug.itemconfigure(r2, outline='black')  # 重设内框颜色
        canvas_drug_drug.itemconfigure(t, fill='black')  # 重设显示文本颜色
        canvas_drug_drug.configure(cursor='hand2')  # 重设鼠标样式
    else:
        canvas_drug_drug.itemconfigure(r1, outline='white')  # 重设外框颜色
        canvas_drug_drug.itemconfigure(r2, outline='white')  # 重设内框颜色
        canvas_drug_drug.itemconfigure(t, fill='white')  # 重设显示文本颜色
        canvas_drug_drug.configure(cursor='arrow')  # 重设鼠标样式


def Mouse_Click_drug_protein(event, drug_num, protein_num):  # 关联鼠标点击事件
    drug_protein_x, drug_protein_y, drug_protein_l, drug_protein_h, drug_protein_r = 380, 360, 200, 60, 5
    if drug_protein_x <= event.x <= drug_protein_x + drug_protein_l and drug_protein_y <= event.y <= drug_protein_y + drug_protein_h:  # 响应的位置
        drug_protein_inquiry(drug_num, protein_num)


def Mouse_over_drug_protein(event, canvas_drug_protein, r1, r2, t):  # 关联鼠标经过事件
    drug_protein_x, drug_protein_y, drug_protein_l, drug_protein_h, drug_protein_r = 380, 360, 200, 60, 5
    if drug_protein_x <= event.x <= drug_protein_x + drug_protein_l and drug_protein_y <= event.y <= drug_protein_y + drug_protein_h:  # 响应的位置
        canvas_drug_protein.itemconfigure(r1, outline='black')  # 重设外框颜色
        canvas_drug_protein.itemconfigure(r2, outline='black')  # 重设内框颜色
        canvas_drug_protein.itemconfigure(t, fill='black')  # 重设显示文本颜色
        canvas_drug_protein.configure(cursor='hand2')  # 重设鼠标样式
    else:
        canvas_drug_protein.itemconfigure(r1, outline='white')  # 重设外框颜色
        canvas_drug_protein.itemconfigure(r2, outline='white')  # 重设内框颜色
        canvas_drug_protein.itemconfigure(t, fill='white')  # 重设显示文本颜色
        canvas_drug_protein.configure(cursor='arrow')  # 重设鼠标样式


def Mouse_over_protein_protein(event, canvas_protein_protein, r1, r2, t):
    protein_protein_x, protein_protein_y, protein_protein_l, protein_protein_h, protein_protein_r = 380, 360, 200, 60, 5
    if protein_protein_x <= event.x <= protein_protein_x + protein_protein_l and protein_protein_y <= event.y <= protein_protein_y + protein_protein_h:  # 响应的位置
        canvas_protein_protein.itemconfigure(r1, outline='black')  # 重设外框颜色
        canvas_protein_protein.itemconfigure(r2, outline='black')  # 重设内框颜色
        canvas_protein_protein.itemconfigure(t, fill='black')  # 重设显示文本颜色
        canvas_protein_protein.configure(cursor='hand2')  # 重设鼠标样式
    else:
        canvas_protein_protein.itemconfigure(r1, outline='white')  # 重设外框颜色
        canvas_protein_protein.itemconfigure(r2, outline='white')  # 重设内框颜色
        canvas_protein_protein.itemconfigure(t, fill='white')  # 重设显示文本颜色
        canvas_protein_protein.configure(cursor='arrow')  # 重设鼠标样式


##########很多嵌套计算模块#########
def PPI(flag_Protein1, flag_Protein2):  # 查询PPI的积极消极作用
    PPI = datareader("PPI").reshape((1, -1))[0][:1243225].reshape((1115, 1115))
    # 蛋白质蛋白质互作用矩阵：(1115,1115)
    PPI_sign = int(float(PPI[flag_Protein1][flag_Protein2]))
    if PPI_sign == 1:
        return ['可相互作用']
    else:
        return ['无明显互作用']


def protein_network(protein_name, flag):
    # 输入药物编号与flag序号，返回一张图
    if not _ensure_matplotlib():
        messagebox.showerror("错误", "Matplotlib 未安装或加载失败，无法绘图。")
        return
    G = nx.Graph()
    color_list = ["green"]
    G.add_node(protein_name)
    ProteinNameID = datareader("ProteinNameID")
    # 蛋白质名称与编号：(1115,5)
    rows_name = ProteinNameID[:, 0]
    rows_PPI = datareader("PPI").reshape((1, -1))[0][0:1243225].reshape((1115, 1115))
    # 蛋白质蛋白质互作用矩阵：(1115,1115)
    PPI = rows_PPI[:, 0]
    pos, neg = [], []
    for i in range(len(PPI)):
        if PPI[i] == 1:
            pos.append(rows_name[i])
        elif PPI[i] == -1:
            neg.append(rows_name[i])
    for i in range(len(pos)):
        G.add_edge(protein_name, str(pos[i]))
        color_list.append("red")
    for i in range(len(neg)):
        G.add_edge(protein_name, str(neg[i]))
        color_list.append("blue")
    nx.draw(G, node_color=color_list, with_labels=True)
    plt.show()
    return


def search_drug_protein(flag_drug, flag_protein):  # 查询药物与蛋白质的互作用类型
    DTA = datareader("DTA").reshape((1, -1))[0][0:2089464].reshape((1443, 1448))
    # 药物靶标互作用矩阵：(1443,1448)
    DTA_num = DTA[flag_drug][flag_protein]
    drug_protein_action = datareader("drug_protein_action")
    # 药物与蛋白质相互作用：(43,2)
    rows_Drug_Protein_Action = drug_protein_action[:, 1]
    return rows_Drug_Protein_Action[DTA_num - 1]


def prediction(src, dst):  # 查询或预测药物互作用
    drug_drug_sign = datareader("drug_drug_sign")
    # 药物-药物-符号：(707867,3)
    temp = drug_drug_sign[:, 0:2].tolist()
    drug_label = [src[0], dst[0]]
    ans = []
    if drug_label in temp:
        sign = drug_drug_sign[temp.index(drug_label)][2]
        if sign == 1:
            ans.append('拮抗作用')
        else:
            ans.append('协同作用')
    else:
        # messagebox.showinfo(title='信息提示', message='数据库种未找到对应信息，正在为您预测！')
        f = torch.Tensor(joblib.load("feature.pkl"))
        weight = torch.Tensor(joblib.load("weight.pkl"))
        bias = torch.Tensor(joblib.load("bias.pkl"))
        edge = [src, dst]
        value = torch.cat([f[edge[0]], f[edge[1]]], dim=1)
        value = torch.nn.functional.linear(value, weight.T, bias)
        pred = torch.nn.functional.log_softmax(value)
        topk_out, topk_indices = torch.topk(pred, k=1)
        for i in topk_indices:
            if i == 1:
                ans.append('拮抗作用')
            else:
                ans.append('协同作用')
    return ans


def signed_graph(drug_num, flag):
    """创建交互式药物关系图（使用Plotly）"""
    try:
        G = nx.Graph()
        G.add_node(drug_num)
        NamesWithID = datareader("NamesWithID")
        rows_name = NamesWithID[:, 0]
        DDI = datareader("DDI").reshape((1, -1))[0][0:2082249].reshape((1443, 1443))
        DDI = DDI[flag]
        
        pos, neg = [], []
        for i in range(len(DDI)):
            if DDI[i] == 1:
                pos.append(rows_name[i])
            elif DDI[i] == -1:
                neg.append(rows_name[i])
        
        for i in range(len(pos)):
            G.add_edge(drug_num, str(pos[i]))
        for i in range(len(neg)):
            G.add_edge(drug_num, str(neg[i]))
        
        # 计算节点数量
        num_nodes = len(G.nodes()) - 1
        num_pos = len(pos)
        num_neg = len(neg)
        
        # 创建分层圆形布局（与原来相同）
        pos_layout = {}
        pos_layout[drug_num] = (0, 0)
        
        if num_nodes <= 20:
            radius = 3.0
            angles = np.linspace(0, 2 * np.pi, num_nodes, endpoint=False)
            for i, node in enumerate(list(G.nodes())[1:]):
                pos_layout[node] = (radius * np.cos(angles[i]), radius * np.sin(angles[i]))
        elif num_nodes <= 60:
            n1 = int(num_nodes * 0.3)
            n2 = num_nodes - n1
            radius1, radius2 = 4.0, 6.5
            angles1 = np.linspace(0, 2 * np.pi, n1, endpoint=False)
            angles2 = np.linspace(0, 2 * np.pi, n2, endpoint=False)
            nodes_list = list(G.nodes())[1:]
            for i in range(n1):
                pos_layout[nodes_list[i]] = (radius1 * np.cos(angles1[i]), radius1 * np.sin(angles1[i]))
            for i in range(n2):
                pos_layout[nodes_list[n1 + i]] = (radius2 * np.cos(angles2[i]), radius2 * np.sin(angles2[i]))
        else:
            n1, n2, n3 = int(num_nodes * 0.2), int(num_nodes * 0.3), num_nodes - int(num_nodes * 0.2) - int(num_nodes * 0.3)
            radius1, radius2, radius3 = 5.0, 8.0, 11.5
            nodes_list = list(G.nodes())[1:]
            pos_nodes = [nodes_list[i] for i in range(num_pos)]
            neg_nodes = [nodes_list[num_pos + i] for i in range(num_neg)]
            mixed_nodes = []
            max_len = max(len(pos_nodes), len(neg_nodes))
            for i in range(max_len):
                if i < len(pos_nodes):
                    mixed_nodes.append(pos_nodes[i])
                if i < len(neg_nodes):
                    mixed_nodes.append(neg_nodes[i])
            angles1 = np.linspace(0, 2 * np.pi, n1, endpoint=False)
            angles2 = np.linspace(0, 2 * np.pi, n2, endpoint=False)
            angles3 = np.linspace(0, 2 * np.pi, n3, endpoint=False)
            for i in range(min(n1, len(mixed_nodes))):
                pos_layout[mixed_nodes[i]] = (radius1 * np.cos(angles1[i]), radius1 * np.sin(angles1[i]))
            for i in range(min(n2, len(mixed_nodes) - n1)):
                pos_layout[mixed_nodes[n1 + i]] = (radius2 * np.cos(angles2[i]), radius2 * np.sin(angles2[i]))
            for i in range(min(n3, len(mixed_nodes) - n1 - n2)):
                pos_layout[mixed_nodes[n1 + n2 + i]] = (radius3 * np.cos(angles3[i]), radius3 * np.sin(angles3[i]))
        
        # 准备Plotly数据
        node_x = []
        node_y = []
        node_text = []
        node_color = []
        node_size = []
        node_info = []
        
        for node in G.nodes():
            x, y = pos_layout[node]
            node_x.append(x)
            node_y.append(y)
            node_text.append(node)
            
            if node == drug_num:
                node_color.append('green')
                node_size.append(25)
                node_info.append(f'<b>中心药物</b><br>{node}<br>药物编号: {drug_num}')
            elif node in pos:
                node_color.append('red')
                node_size.append(15)
                node_info.append(f'<b>协同作用药物</b><br>{node}<br>关系类型: 协同作用')
            else:
                node_color.append('blue')
                node_size.append(15)
                node_info.append(f'<b>拮抗作用药物</b><br>{node}<br>关系类型: 拮抗作用')
        
        # 创建边的数据
        edge_x = []
        edge_y = []
        edge_info = []
        edge_color = []
        
        for edge in G.edges():
            x0, y0 = pos_layout[edge[0]]
            x1, y1 = pos_layout[edge[1]]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
            
            if edge[0] == drug_num:
                target = edge[1]
            else:
                target = edge[0]
            
            if target in pos:
                edge_color.append('#ff6b6b')
                edge_info.append(f'{drug_num} ↔ {target}<br>关系: 协同作用')
            else:
                edge_color.append('#4dabf7')
                edge_info.append(f'{drug_num} ↔ {target}<br>关系: 拮抗作用')
        
        # 创建Plotly图形
        fig = go.Figure()
        
        # 添加边（使用scatter绘制线条）
        for i in range(0, len(edge_x), 3):
            if i + 2 < len(edge_x):
                fig.add_trace(go.Scatter(
                    x=[edge_x[i], edge_x[i+1]],
                    y=[edge_y[i], edge_y[i+1]],
                    mode='lines',
                    line=dict(width=2, color=edge_color[i//3] if i//3 < len(edge_color) else 'gray'),
                    hoverinfo='skip',
                    showlegend=False
                ))
        
        # 添加节点
        fig.add_trace(go.Scatter(
            x=node_x,
            y=node_y,
            mode='markers+text',
            marker=dict(
                size=node_size,
                color=node_color,
                line=dict(width=2, color='black')
            ),
            text=node_text,
            textposition="middle center",
            textfont=dict(size=10, family="Microsoft YaHei", color="white"),
            hovertext=node_info,
            hoverinfo='text',
            name='药物节点',
            showlegend=False
        ))
        
        # 更新布局
        fig.update_layout(
            title=dict(
                text=f'药物关系图 - {drug_num}<br><sub>红色:协同作用 | 蓝色:拮抗作用 | 绿色:中心药物</sub>',
                font=dict(size=20, family="Microsoft YaHei", color="black"),
                x=0.5
            ),
            showlegend=False,
            hovermode='closest',
            margin=dict(b=20, l=5, r=5, t=80),
            annotations=[dict(
                text="提示: 鼠标悬停查看详细信息 | 滚轮缩放 | 拖拽平移",
                showarrow=False,
                xref="paper", yref="paper",
                x=0.5, y=-0.05,
                xanchor="center", yanchor="top",
                font=dict(size=12, color="gray")
            )],
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor='white',
            width=1200,
            height=900
        )
        
        # 显示图形（在浏览器中打开）
        fig.show()
        
    except Exception as e:
        messagebox.showerror("错误", f"创建交互式关系图时出错: {str(e)}")
    return


def feat_signed_graph(drug_name, flag):
    """创建交互式特征关系图（使用Plotly）"""
    try:
        G = nx.Graph()
        G.add_node(drug_name)
        with open('pubchem_fingerprints.csv', 'r') as NamesWithID:
            reader = csv.reader(NamesWithID)
            rows_name = [row[1] for row in reader]
        drug_881feat = datareader("drug_881feat")
        Drug_881feat = drug_881feat[flag]
        feat = []
        for i in range(len(Drug_881feat)):
            if Drug_881feat[i] == 1:
                feat.append(rows_name[i])
        for i in range(len(feat)):
            G.add_edge(drug_name, str(feat[i]))
        
        # 计算节点数量
        num_nodes = len(G.nodes()) - 1
        
        # 创建分层圆形布局
        pos_layout = {}
        pos_layout[drug_name] = (0, 0)
        
        if num_nodes <= 20:
            radius = 3.0
            angles = np.linspace(0, 2 * np.pi, num_nodes, endpoint=False)
            nodes_list = list(G.nodes())[1:]
            for i, node in enumerate(nodes_list):
                pos_layout[node] = (radius * np.cos(angles[i]), radius * np.sin(angles[i]))
        elif num_nodes <= 60:
            n1 = int(num_nodes * 0.3)
            n2 = num_nodes - n1
            radius1, radius2 = 4.0, 6.5
            angles1 = np.linspace(0, 2 * np.pi, n1, endpoint=False)
            angles2 = np.linspace(0, 2 * np.pi, n2, endpoint=False)
            nodes_list = list(G.nodes())[1:]
            for i in range(n1):
                pos_layout[nodes_list[i]] = (radius1 * np.cos(angles1[i]), radius1 * np.sin(angles1[i]))
            for i in range(n2):
                pos_layout[nodes_list[n1 + i]] = (radius2 * np.cos(angles2[i]), radius2 * np.sin(angles2[i]))
        else:
            n1 = int(num_nodes * 0.2)
            n2 = int(num_nodes * 0.3)
            n3 = num_nodes - n1 - n2
            radius1, radius2, radius3 = 5.0, 8.0, 11.5
            angles1 = np.linspace(0, 2 * np.pi, n1, endpoint=False)
            angles2 = np.linspace(0, 2 * np.pi, n2, endpoint=False)
            angles3 = np.linspace(0, 2 * np.pi, n3, endpoint=False)
            nodes_list = list(G.nodes())[1:]
            for i in range(n1):
                pos_layout[nodes_list[i]] = (radius1 * np.cos(angles1[i]), radius1 * np.sin(angles1[i]))
            for i in range(n2):
                pos_layout[nodes_list[n1 + i]] = (radius2 * np.cos(angles2[i]), radius2 * np.sin(angles2[i]))
            for i in range(n3):
                pos_layout[nodes_list[n1 + n2 + i]] = (radius3 * np.cos(angles3[i]), radius3 * np.sin(angles3[i]))
        
        # 准备Plotly数据
        node_x = []
        node_y = []
        node_text = []
        node_color = []
        node_size = []
        node_info = []
        
        for node in G.nodes():
            x, y = pos_layout[node]
            node_x.append(x)
            node_y.append(y)
            node_text.append(node)
            
            if node == drug_name:
                node_color.append('blue')
                node_size.append(25)
                node_info.append(f'<b>中心药物</b><br>{node}<br>特征关系图')
            else:
                node_color.append('red')
                node_size.append(15)
                node_info.append(f'<b>化学特征</b><br>{node}<br>关联药物: {drug_name}')
        
        # 创建边的数据
        edge_x = []
        edge_y = []
        edge_color = []
        
        for edge in G.edges():
            x0, y0 = pos_layout[edge[0]]
            x1, y1 = pos_layout[edge[1]]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
            edge_color.append('#ff6b6b')
        
        # 创建Plotly图形
        fig = go.Figure()
        
        # 添加边
        for i in range(0, len(edge_x), 3):
            if i + 2 < len(edge_x):
                fig.add_trace(go.Scatter(
                    x=[edge_x[i], edge_x[i+1]],
                    y=[edge_y[i], edge_y[i+1]],
                    mode='lines',
                    line=dict(width=2, color='#ff6b6b'),
                    hoverinfo='skip',
                    showlegend=False
                ))
        
        # 添加节点
        fig.add_trace(go.Scatter(
            x=node_x,
            y=node_y,
            mode='markers+text',
            marker=dict(
                size=node_size,
                color=node_color,
                line=dict(width=2, color='black')
            ),
            text=node_text,
            textposition="middle center",
            textfont=dict(size=10, family="Microsoft YaHei", color="white"),
            hovertext=node_info,
            hoverinfo='text',
            name='特征节点',
            showlegend=False
        ))
        
        # 更新布局
        fig.update_layout(
            title=dict(
                text=f'特征关系图 - {drug_name}<br><sub>蓝色:中心药物 | 红色:相关化学特征</sub>',
                font=dict(size=20, family="Microsoft YaHei", color="black"),
                x=0.5
            ),
            showlegend=False,
            hovermode='closest',
            margin=dict(b=20, l=5, r=5, t=80),
            annotations=[dict(
                text="提示: 鼠标悬停查看详细信息 | 滚轮缩放 | 拖拽平移",
                showarrow=False,
                xref="paper", yref="paper",
                x=0.5, y=-0.05,
                xanchor="center", yanchor="top",
                font=dict(size=12, color="gray")
            )],
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor='white',
            width=1200,
            height=900
        )
        
        # 显示图形
        fig.show()
        
    except Exception as e:
        messagebox.showerror("错误", f"创建交互式特征关系图时出错: {str(e)}")
    return


def OFFSIDE(flag):  # 寻找指定flag的offside列表
    offside = []
    OffsideNames = datareader("OffsideName")
    # 副作用名：(14341,1)
    Offsides = datareader("Offsides").reshape((1, -1))[0][0:20694063].reshape((1443, 14341))
    # 药物副作用信息：(1443,14341)
    offside_index = Offsides[flag]
    for i in range(14341):
        if offside_index[i] == 0:
            offside.append(OffsideNames[i][0])
    return offside  # 返回offside嵌套列表


def drug_881feat(flag):
    drug_881feat = []
    with open('pubchem_fingerprints.csv', 'r') as feat:
        reader = csv.reader(feat)
        featname_881 = [row[1] for row in reader]  # 881种特征
    rows_881feat = datareader("drug_881feat")
    # 药物881特征：(1443,881)
    feat_index = rows_881feat[flag]
    for i in range(881):
        if feat_index[i] == 1:
            drug_881feat.append(featname_881[i])
    return drug_881feat  # 返回drug_881feat嵌套列表


##########四个模块的计算结果函数##########
def get_drug_smiles(drug_name):
    """获取药物的SMILES字符串"""
    try:
        # 方法1: 从Case/drug.tsv文件读取（优先，因为这是用户查询时生成的）
        try:
            df = pd.read_csv('Case/drug.tsv', sep='\t')
            if 'SMILES' in df.columns and 'Name' in df.columns:
                match = df[df['Name'].str.contains(drug_name, case=False, na=False)]
                if not match.empty:
                    smiles = match.iloc[0]['SMILES']
                    if pd.notna(smiles) and smiles.strip():
                        return smiles.strip()
        except Exception as e:
            print(f"从Case/drug.tsv读取失败: {e}")
        
        # 方法2: 从EGFR-Case/drug.tsv读取
        try:
            df = pd.read_csv('EGFR-Case/drug.tsv', sep='\t')
            if 'SMILES' in df.columns and 'Name' in df.columns:
                match = df[df['Name'].str.contains(drug_name, case=False, na=False)]
                if not match.empty:
                    smiles = match.iloc[0]['SMILES']
                    if pd.notna(smiles) and smiles.strip():
                        return smiles.strip()
        except Exception as e:
            print(f"从EGFR-Case/drug.tsv读取失败: {e}")
        
        # 方法3: 从kiba数据文件读取
        try:
            df = pd.read_csv('kiba/kiba_drugs.csv', sep=',')
            if 'SMILES' in df.columns:
                # 尝试匹配药物名称或ID
                for col in df.columns:
                    match = df[df[col].astype(str).str.contains(drug_name, case=False, na=False)]
                    if not match.empty and 'SMILES' in match.columns:
                        smiles = match.iloc[0]['SMILES']
                        if pd.notna(smiles) and smiles.strip():
                            return smiles.strip()
        except Exception as e:
            print(f"从kiba数据读取失败: {e}")
        
        # 方法4: 尝试从数据库的DTA表或其他表读取（如果有drug_smile字段）
        try:
            # 检查是否有包含SMILES的数据表
            conn = get_sql_connection()
            sql_query = "SHOW TABLES LIKE '%drug%'"
            tables_df = pd.read_sql(sql_query, conn)
            for table in tables_df.values.flatten():
                try:
                    # 尝试读取表的前几行，检查是否有SMILES相关字段
                    sample = pd.read_sql(f"SELECT * FROM {table} LIMIT 1", conn)
                    for col in sample.columns:
                        if 'smile' in col.lower() or 'smiles' in col.lower():
                            # 找到包含SMILES的字段，尝试查询
                            query = f"SELECT `{col}` FROM `{table}` WHERE `drug_id` LIKE %s OR `drug_name` LIKE %s LIMIT 1"
                            result = pd.read_sql(query, conn, params=[f'%{drug_name}%', f'%{drug_name}%'])
                            if not result.empty:
                                smiles = result.iloc[0][col]
                                if pd.notna(smiles) and str(smiles).strip():
                                    return str(smiles).strip()
                except:
                    continue
        except Exception as e:
            print(f"从数据库读取失败: {e}")
        
        return None
    except Exception as e:
        print(f"获取SMILES时出错: {e}")
        return None


def show_3d_molecule_structure(drug_name, flag):
    """显示3D分子结构图（基于真实SMILES）"""
    try:
        if not _ensure_rdkit():
            messagebox.showerror("错误", "RDKit未安装或加载失败，无法显示3D分子结构。\n请安装: conda install -c conda-forge rdkit")
            return
        
        setup_chinese_font()

        # 获取SMILES字符串
        smiles = get_drug_smiles(drug_name)
        
        if not smiles:
            # 如果找不到SMILES，尝试从其他数据源获取
            messagebox.showwarning("警告", f"未找到药物 '{drug_name}' 的SMILES数据。\n将尝试使用数据库中的信息。")
            # 可以添加其他获取SMILES的方法
            return
        
        # 使用RDKit解析SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            messagebox.showerror("错误", f"无法解析SMILES字符串: {smiles}")
            return
        
        # 添加氢原子
        mol = Chem.AddHs(mol)
        
        # 生成3D坐标
        try:
            # 使用ETKDG方法生成3D坐标
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
        except:
            # 如果ETKDG失败，使用基本方法
            try:
                AllChem.EmbedMolecule(mol)
                AllChem.MMFFOptimizeMolecule(mol)
            except:
                messagebox.showerror("错误", "无法生成3D坐标")
                return
        
        # 获取原子坐标
        conf = mol.GetConformer()
        num_atoms = mol.GetNumAtoms()

        # 原子坐标
        atoms_pos = np.zeros((num_atoms, 3))
        atom_symbols = []
        atom_colors_map = {
            'C': '#909090',  # 灰色
            'O': '#FF0000',  # 红色
            'N': '#3050F8',  # 蓝色
            'H': '#FFFFFF',  # 白色
            'S': '#FFFF30',  # 黄色
            'P': '#FF8000',  # 橙色
            'F': '#90E050',  # 绿色
            'Cl': '#1FF01F', # 绿色
            'Br': '#A62929', # 深红色
            'I': '#940094',  # 紫色
        }
        atom_sizes_map = {
            'C': 300, 'O': 280, 'N': 260, 'H': 180,
            'S': 320, 'P': 340, 'F': 240, 'Cl': 300,
            'Br': 360, 'I': 400
        }
        
        for i in range(num_atoms):
            pos = conf.GetAtomPosition(i)
            atoms_pos[i] = [pos.x, pos.y, pos.z]
            atom = mol.GetAtomWithIdx(i)
            atom_symbols.append(atom.GetSymbol())
        
        # 创建3D图形
        fig = plt.figure(figsize=(14, 12))
        ax = fig.add_subplot(111, projection='3d')
        
        # 计算中心点并居中显示
        center = atoms_pos.mean(axis=0)
        atoms_pos = atoms_pos - center
        
        # 绘制化学键
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            bond_type = bond.GetBondType()
            
            # 根据键类型设置颜色和宽度
            if bond_type == Chem.BondType.SINGLE:
                color = 'gray'
                width = 2
            elif bond_type == Chem.BondType.DOUBLE:
                color = 'darkgray'
                width = 3
            elif bond_type == Chem.BondType.TRIPLE:
                color = 'black'
                width = 4
            else:
                color = 'gray'
                width = 2
            
            ax.plot([atoms_pos[begin_idx, 0], atoms_pos[end_idx, 0]],
                   [atoms_pos[begin_idx, 1], atoms_pos[end_idx, 1]],
                   [atoms_pos[begin_idx, 2], atoms_pos[end_idx, 2]],
                   color=color, linewidth=width, alpha=0.8)

        # 绘制原子
        for i in range(num_atoms):
            symbol = atom_symbols[i]
            color = atom_colors_map.get(symbol, '#808080')  # 默认灰色
            size = atom_sizes_map.get(symbol, 250)
            
            ax.scatter(atoms_pos[i, 0], atoms_pos[i, 1], atoms_pos[i, 2],
                      c=color, s=size, alpha=0.9, edgecolors='black', linewidths=1.5)
            
            # 标注原子（可选，对于大分子可能太密集）
            if num_atoms <= 50:  # 只对较小的分子标注
                ax.text(atoms_pos[i, 0], atoms_pos[i, 1], atoms_pos[i, 2],
                       f' {symbol}', fontsize=8, fontweight='bold')
        
        # 设置坐标轴
        max_range = np.max(np.abs(atoms_pos)) * 1.2
        ax.set_xlim([-max_range, max_range])
        ax.set_ylim([-max_range, max_range])
        ax.set_zlim([-max_range, max_range])
        
        ax.set_xlabel('X (Å)', fontsize=12, fontfamily='Microsoft YaHei', labelpad=10)
        ax.set_ylabel('Y (Å)', fontsize=12, fontfamily='Microsoft YaHei', labelpad=10)
        ax.set_zlabel('Z (Å)', fontsize=12, fontfamily='Microsoft YaHei', labelpad=10)
        
        # 标题
        title = f'{drug_name} - 3D分子结构\nSMILES: {smiles}'
        ax.set_title(title, fontsize=16, fontfamily='Microsoft YaHei', fontweight='bold', pad=20)

        # 添加图例
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#909090', label='碳 (C)'),
            Patch(facecolor='#FF0000', label='氧 (O)'),
            Patch(facecolor='#3050F8', label='氮 (N)'),
            Patch(facecolor='#FFFFFF', label='氢 (H)'),
            Patch(facecolor='#FFFF30', label='硫 (S)'),
            Patch(facecolor='#FF8000', label='磷 (P)'),
        ]
        ax.legend(handles=legend_elements, loc='upper left', fontsize=10, prop={'family': 'Microsoft YaHei'})
        
        # 设置背景色
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    except Exception as e:
        import traceback
        error_msg = f"显示3D分子结构时出错: {str(e)}\n{traceback.format_exc()}"
        messagebox.showerror("错误", error_msg)


def show_3d_network_graph(drug_name, flag):
    """显示3D药物关系网络图"""
    try:
        setup_chinese_font()

        # 获取药物关系数据
        NamesWithID = datareader("NamesWithID")
        rows_name = NamesWithID[:, 0]
        DDI = datareader("DDI").reshape((1, -1))[0][0:2082249].reshape((1443, 1443))
        DDI = DDI[flag]

        # 收集相关药物
        related_drugs = []
        for i in range(min(50, len(DDI))):
            if DDI[i] != 0 and i != flag:
                related_drugs.append((rows_name[i], DDI[i]))

        # 创建3D图
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        # 中心节点（查询的药物）
        ax.scatter([0], [0], [0], c='green', s=300, alpha=0.9, label=drug_name)

        # 相关药物节点
        if related_drugs:
            n_related = len(related_drugs)
            angles = np.linspace(0, 2 * np.pi, n_related, endpoint=False)

            radius = 8
            x = radius * np.cos(angles)
            y = radius * np.sin(angles)
            z = np.random.uniform(-3, 3, n_related)  # 添加高度变化

            # 绘制相关节点和连接线
            for i, (drug, relation) in enumerate(related_drugs):
                color = 'red' if relation == 1 else 'blue'
                ax.scatter(x[i], y[i], z[i], c=color, s=150, alpha=0.7)

                # 添加连接线
                line_color = 'red' if relation == 1 else 'blue'
                ax.plot([0, x[i]], [0, y[i]], [0, z[i]],
                        color=line_color, alpha=0.5, linewidth=2)

        ax.set_xlabel('X Axis', fontsize=12, labelpad=10)
        ax.set_ylabel('Y Axis', fontsize=12, labelpad=10)
        ax.set_zlabel('Z Axis', fontsize=12, labelpad=10)

        title = f'{drug_name} - 3D Drug Network'
        ax.set_title(title, fontsize=14, pad=20)

        plt.tight_layout()
        plt.show()

    except Exception as e:
        messagebox.showerror("错误", f"显示3D网络图时出错: {str(e)}")


def drug_inquiry(drug_name):  # 单药物计算查询函数
    # 接下来判断数据库中是否能找到相关信息
    NamesWithID = datareader("NamesWithID")
    # 药物及药物编号：(1443,2)
    Drug_Name = NamesWithID[:, 0].tolist()
    if Drug_Name.count(drug_name) == 0:
        messagebox.showwarning(title='错误提示', message='您输入的药物编号有误，请重新输入！')
        return
    else:
        flag = Drug_Name.index(drug_name)
    # 获取信息
    Drug_ID = NamesWithID[:, 1]
    drug_id = Drug_ID[flag]  # 实际药编号
    offside = OFFSIDE(flag)  # 副作用列表
    feat = drug_881feat(flag)

    # 显示新窗口（子窗口，简化版）
    global main_root
    win_drug_inquiry = tkinter.Toplevel(main_root)
    win_drug_inquiry.title('药物信息查询结果')
    win_drug_inquiry.geometry('600x650')  # 简化后的窗口大小
    win_drug_inquiry.config(bg='#F5F5F5')
    win_drug_inquiry.transient(main_root)
    win_drug_inquiry.grab_set()
    win_drug_inquiry.resizable(False, False)

    # 标题
    title_label = tkinter.Label(win_drug_inquiry,
                               text=f'{drug_name} 查询结果',
                               font=(UI_FONT_FAMILY, 16, 'bold'),
                               bg='#F5F5F5',
                               fg='#333333')
    title_label.pack(pady=15)

    # 信息显示区域（使用Frame组织）
    info_frame = tkinter.Frame(win_drug_inquiry, bg='#F5F5F5')
    info_frame.pack(pady=10, padx=20, fill=tkinter.BOTH, expand=True)

    # 基本信息（简化显示）
    basic_frame = tkinter.LabelFrame(info_frame, text='基本信息', font=(UI_FONT_FAMILY, 12), bg='#F5F5F5', fg='#333333')
    basic_frame.pack(fill=tkinter.X, pady=5)
    
    tkinter.Label(basic_frame, text=f'药物编号：{drug_id}', font=(UI_FONT_FAMILY, 11), bg='#F5F5F5', anchor='w').pack(fill=tkinter.X, padx=10, pady=3)
    tkinter.Label(basic_frame, text=f'药物名称：{drug_name}', font=(UI_FONT_FAMILY, 11), bg='#F5F5F5', anchor='w').pack(fill=tkinter.X, padx=10, pady=3)

    # 副作用（简化显示）
    offside_frame = tkinter.LabelFrame(info_frame, text='副作用', font=(UI_FONT_FAMILY, 12), bg='#F5F5F5', fg='#333333')
    offside_frame.pack(fill=tkinter.BOTH, expand=True, pady=5)
    
    txt_offside = tkinter.Text(offside_frame,
                               font=(UI_FONT_FAMILY, 10),
                                       bg='white',
                               fg='#333333',
                               relief=tkinter.SOLID,
                               borderwidth=1,
                               wrap=tkinter.WORD,
                               height=4)
    txt_offside.pack(fill=tkinter.BOTH, expand=True, padx=5, pady=5)
    for i in offside[:6]:  # 减少显示数量
        txt_offside.insert(tkinter.END, i + '\n')
    txt_offside.config(state=tkinter.DISABLED)

    # 化学特征（简化显示）
    feat_frame = tkinter.LabelFrame(info_frame, text='化学特征', font=(UI_FONT_FAMILY, 12), bg='#F5F5F5', fg='#333333')
    feat_frame.pack(fill=tkinter.BOTH, expand=True, pady=5)
    
    txt_feat = tkinter.Text(feat_frame,
                           font=(UI_FONT_FAMILY, 10),
                                  bg='white',
                           fg='#333333',
                           relief=tkinter.SOLID,
                           borderwidth=1,
                           wrap=tkinter.WORD,
                           height=4)
    txt_feat.pack(fill=tkinter.BOTH, expand=True, padx=5, pady=5)
    for i in feat[:8]:  # 减少显示数量
        txt_feat.insert(tkinter.END, i + '\n')
    txt_feat.config(state=tkinter.DISABLED)

    # 可视化按钮区域（简化，使用网格布局）
    viz_frame = tkinter.LabelFrame(info_frame, text='可视化', font=(UI_FONT_FAMILY, 12), bg='#F5F5F5', fg='#333333')
    viz_frame.pack(fill=tkinter.X, pady=5)

    btn_frame = tkinter.Frame(viz_frame, bg='#F5F5F5')
    btn_frame.pack(pady=10)

    # 第一行按钮
    btn_row1 = tkinter.Frame(btn_frame, bg='#F5F5F5')
    btn_row1.pack(pady=3)
    for col in range(2):
        btn_row1.grid_columnconfigure(col, weight=1, uniform='viz_row1')
    
    create_button(btn_row1, text='3D分子',
                  command=lambda: show_3d_molecule_structure(drug_name, flag),
                  font=(UI_FONT_FAMILY, 10), width=12, padx=5, pady=3).grid(row=0, column=0, padx=6)
    create_button(btn_row1, text='3D网络',
                  command=lambda: show_3d_network_graph(drug_name, flag),
                  font=(UI_FONT_FAMILY, 10), width=12, padx=5, pady=3).grid(row=0, column=1, padx=6)

    # 第二行按钮
    btn_row2 = tkinter.Frame(btn_frame, bg='#F5F5F5')
    btn_row2.pack(pady=3)
    for col in range(2):
        btn_row2.grid_columnconfigure(col, weight=1, uniform='viz_row2')
    
    create_button(btn_row2, text='药物关系图',
                  command=lambda: signed_graph(drug_id, flag),
                  font=(UI_FONT_FAMILY, 10), width=12, padx=5, pady=3).grid(row=0, column=0, padx=6)
    create_button(btn_row2, text='特征关系图',
                  command=lambda: feat_signed_graph(drug_name, flag),
                  font=(UI_FONT_FAMILY, 10), width=12, padx=5, pady=3).grid(row=0, column=1, padx=6)

    # 返回按钮
    btn_return = create_button(
        win_drug_inquiry,
        text='返回',
        command=win_drug_inquiry.destroy,
        style='secondary',
        width=15,
        pady=5
    )
    btn_return.pack(pady=15)

def drug_drug_inquiry(drug1_name, drug2_name):  # 药物-药物计算查询函数
    NamesWithID = datareader("NamesWithID")
    # 药物及药物编号：(1443,2)
    Drug_Name = NamesWithID[:, 0].tolist()
    if (drug2_name not in Drug_Name) or (drug1_name not in Drug_Name):
        messagebox.showwarning(title='错误提示', message='您输入的药物名称有误，请重新输入！')
        return
    else:
        flag1 = Drug_Name.index(drug1_name)
        flag2 = Drug_Name.index(drug2_name)
    # 显示新窗口（子窗口）
    global main_root
    win_drug_drug_inquiry = tkinter.Toplevel(main_root)
    win_drug_drug_inquiry.transient(main_root)
    win_drug_drug_inquiry.grab_set()
    win_drug_drug_inquiry.title('药物-药物分析结果')
    win_drug_drug_inquiry.geometry('500x400')
    win_drug_drug_inquiry.config(bg='#F5F5F5')
    win_drug_drug_inquiry.resizable(False, False)
    # 结果显示（简化）
    result_frame = tkinter.Frame(win_drug_drug_inquiry, bg='#F5F5F5')
    result_frame.pack(pady=30, padx=20, fill=tkinter.BOTH, expand=True)
    
    # 标题
    title_label = tkinter.Label(result_frame,
                               text=f'{drug1_name} 与 {drug2_name}',
                               font=(UI_FONT_FAMILY, 14, 'bold'),
                               bg='#F5F5F5',
                               fg='#333333')
    title_label.pack(pady=10)
    
    # 预测结果（突出显示）
    result_text = prediction([flag1], [flag2])
    result_label = tkinter.Label(result_frame,
                                 text=f'互作用关系：{result_text}',
                                 font=(UI_FONT_FAMILY, 16, 'bold'),
                                 bg='#E8F4F8',
                                 fg='#2C3E50',
                                 relief=tkinter.SOLID,
                                 borderwidth=2,
                                 padx=20,
                                 pady=15)
    result_label.pack(pady=20, fill=tkinter.X)
    
    # 按钮区域（简化）
    btn_frame = tkinter.Frame(result_frame, bg='#F5F5F5')
    btn_frame.pack(pady=20)
    for col in range(2):
        btn_frame.grid_columnconfigure(col, weight=1, uniform='drug_drug_btns')
    
    create_button(
        btn_frame,
        text='查看药物1详情',
        command=lambda: drug_inquiry(drug1_name),
        font=(UI_FONT_FAMILY, 11),
        width=15,
        padx=10,
        pady=5
    ).grid(row=0, column=0, padx=8)
    
    create_button(
        btn_frame,
        text='查看药物2详情',
        command=lambda: drug_inquiry(drug2_name),
        font=(UI_FONT_FAMILY, 11),
        width=15,
        padx=10,
        pady=5
    ).grid(row=0, column=1, padx=8)


def dta_predicts(drug_name, protein_name, resolved_prot=None):
    global main_root
    win_drug_inquiry = tkinter.Toplevel(main_root)
    win_drug_inquiry.title('DTA预测结果')
    win_drug_inquiry.geometry('600x350')
    win_drug_inquiry.config(bg='#F5F5F5')
    win_drug_inquiry.transient(main_root)
    win_drug_inquiry.grab_set()
    win_drug_inquiry.resizable(False, False)

    # 标题
    title_label = tkinter.Label(win_drug_inquiry,
                               text='DTA预测结果',
                               font=(UI_FONT_FAMILY, 16, 'bold'),
                               bg='#F5F5F5',
                               fg='#333333')
    title_label.pack(pady=20)

    # 旧版三行信息布局：药物 / 靶标 / DTA
    info_frame = tkinter.Frame(win_drug_inquiry, bg='#F5F5F5')
    info_frame.pack(pady=10, padx=30, fill=tkinter.BOTH, expand=True)

    def _row(parent, label_text, value_text, row_index):
        tkinter.Label(
            parent,
            text=label_text,
            fg='#333333',
            bg='#F5F5F5',
            font=(UI_FONT_FAMILY, 13),
            anchor='w'
        ).grid(row=row_index, column=0, sticky='w', pady=10)
        val = tkinter.Label(
            parent,
            text=value_text,
            fg='#333333',
            bg='#F5F5F5',
            font=(UI_FONT_FAMILY, 13),
            anchor='w'
        )
        val.grid(row=row_index, column=1, sticky='w', pady=10)
        return val

    info_frame.grid_columnconfigure(0, weight=0)
    info_frame.grid_columnconfigure(1, weight=1)

    _row(info_frame, '药物名称：', str(drug_name), 0)
    _row(info_frame, '靶标名称：', str(protein_name), 1)
    dta_value_label = _row(info_frame, 'DTA：', '计算中...', 2)
    
    # 记录开始时间
    start_time = time.time()
    start_pred_path = _get_pred_output_path()
    try:
        start_pred_mtime = os.path.getmtime(start_pred_path) if os.path.exists(start_pred_path) else 0
    except Exception:
        start_pred_mtime = 0
    
    # 返回按钮
    btn_return = create_button(
        win_drug_inquiry,
        text='返回',
        command=win_drug_inquiry.destroy,
        style='secondary',
        width=15,
        pady=5
    )
    btn_return.pack(pady=15)

    _poll_state = {'done': False}

    def _read_pred_and_update(force=False):
        pred_path = _get_pred_output_path()
        elapsed_time = time.time() - start_time
        try:
            if not os.path.exists(pred_path):
                if force:
                    dta_value_label.config(text=f'文件未找到: {pred_path}', fg='red')
                return False

            try:
                mtime = os.path.getmtime(pred_path)
            except Exception:
                mtime = 0

            # 如果文件还没更新（仍是旧结果），先继续等待；force=True 时跳过该判断
            if (not force) and (mtime <= start_pred_mtime) and (elapsed_time < 120):
                return False

            try:
                df_pred = pd.read_csv(pred_path)
            except Exception as e:
                if force:
                    dta_value_label.config(text=f'读取结果失败: {str(e)}', fg='red')
                return False

            target_prot = resolved_prot or protein_name
            matched = pd.DataFrame()
            if not df_pred.empty and 'drug_id' in df_pred.columns and 'prot_id' in df_pred.columns:
                # 1. Direct match
                matched = df_pred[(df_pred['drug_id'] == drug_name) & (df_pred['prot_id'] == target_prot)]

                # 2. Loose match (strip whitespace)
                if matched.empty:
                    mask_drug = df_pred['drug_id'].astype(str).str.strip() == str(drug_name).strip()
                    mask_prot = df_pred['prot_id'].astype(str).str.strip() == str(target_prot).strip()
                    matched = df_pred[mask_drug & mask_prot]

                # 3. Fallback: single row
                if matched.empty and len(df_pred) == 1:
                    matched = df_pred
            elif len(df_pred) == 1:
                matched = df_pred

            if not matched.empty:
                if 'pred' in matched.columns:
                    val = matched['pred'].iloc[0]
                else:
                    val = matched.iloc[0, -1]
                try:
                    dtan = f"{float(val):.4f}"
                except Exception:
                    dtan = str(val)
            else:
                dtan = '未匹配到结果'

            dta_value_label.config(text=f'{dtan} (耗时: {elapsed_time:.1f}s)', fg='#333333')
            return True
        except Exception as e:
            if force:
                dta_value_label.config(text=f'界面错误: {str(e)}', fg='red')
            return False

    def _poll_pred_result():
        if _poll_state.get('done'):
            return
        try:
            if not win_drug_inquiry.winfo_exists():
                _poll_state['done'] = True
                return
        except Exception:
            _poll_state['done'] = True
            return

        if _read_pred_and_update(force=False):
            _poll_state['done'] = True
            return
        if time.time() - start_time > 120:
            dta_value_label.config(text='超时：未获取到预测结果', fg='red')
            _poll_state['done'] = True
            return
        win_drug_inquiry.after(1000, _poll_pred_result)

    def _load_pred_result():
        # pred.py 正常退出后，强制读取并更新一次
        _poll_state['done'] = True
        _read_pred_and_update(force=True)

    def _handle_pred_error(exc):
        _poll_state['done'] = True
        dta_value_label.config(text='预测失败，请查看环境或数据配置', fg='red')
        detail = _LAST_PRED_LOG.strip()
        if detail:
            messagebox.showwarning(title='错误提示', message=f'运行预测失败：{exc}\n\n{detail}')
        else:
            messagebox.showwarning(title='错误提示', message=f'运行预测失败：{exc}')

    _run_pred_update_async(on_success=_load_pred_result, on_error=_handle_pred_error)
    # 同时启动轮询：即使 pred.py 未及时退出，也能在结果文件更新后刷新界面
    win_drug_inquiry.after(1000, _poll_pred_result)


def drug_prot_dta_predict(drug_name, protein_name):
    NamesWithID = datareader("NamesWithID")
    input_file = 'EGFR-Case/drug.tsv'
    df = pd.read_csv(input_file, sep="\t")
    Drug_Name = df.iloc[:, 0]
    found1 = df.iloc[:, 0].str.contains(drug_name, na=False).any()
    input_file1 = 'davis_prots.csv'
    df1 = pd.read_csv(input_file1, sep=",")
    Protein_Name = df1.iloc[:, 0].tolist()
    resolved_prot, prot_suggestions = _resolve_protein_id(protein_name, Protein_Name)
    if not found1:
        messagebox.showwarning(title='错误提示', message='您输入的药物名称有误，请重新输入！')
        return
    if not resolved_prot:
        if _is_in_proteinnameid(protein_name, Protein_Name):
            tip = '该名称存在于 ProteinNameID，但不在 davis_prots.csv 中，无法进行当前模型预测。请改用 davis_prots.csv 的 prot_id（第一列）。'
        else:
            tip = '靶标名称需使用 davis_prots.csv 的 prot_id（第一列），或数据库 ProteinNameID 中的名称。'
        if prot_suggestions:
            tip += f"\n可能的候选：{', '.join(prot_suggestions)}"
        messagebox.showwarning(title='错误提示', message=f'您输入的靶标名称有误，请重新输入！\n{tip}')
        return
    data_extractor.find_drug(drug_name, input_file)
    data_extractor.find_protein(resolved_prot, input_file1)
    # 直接进入结果窗口并开始预测（不再弹出“准备预测DTA”确认窗）
    dta_predicts(drug_name, protein_name, resolved_prot)
    return


def drug_protein_inquiry(drug_name, protein_name):  # 药物-靶标计算查询函数
    NamesWithID = datareader("NamesWithID")
    # 药物及药物编号：(1443,2)
    Drug_Name = NamesWithID[:, 0].tolist()
    ProteinNameID = datareader("ProteinNameID")
    # 蛋白质名称与编号：(1115,5)
    Protein_Name = ProteinNameID[:, 0].tolist()
    if Drug_Name.count(drug_name) == 0 or Protein_Name.count(protein_name) == 0:
        messagebox.showwarning(title='错误提示', message='您输入的药物名称有误，请重新输入！')
        return
    else:
        flag_drug = Drug_Name.index(drug_name)
        flag_Protein = Protein_Name.index(protein_name)

    Drug_ID = NamesWithID[:, 1]
    drug_id = Drug_ID[flag_drug]  # 实际药编号
    # 显示新窗口（子窗口）
    global main_root
    win_drug_protein_inquiry = tkinter.Toplevel(main_root)
    win_drug_protein_inquiry.transient(main_root)
    win_drug_protein_inquiry.grab_set()
    win_drug_protein_inquiry.title(software_name + '查询结果')
    win_drug_protein_inquiry.geometry('850x480')
    win_drug_protein_inquiry.config(background="white")

    # 互作用关系预测显示
    lb_drug_protein_inquiry = tkinter.Label(win_drug_protein_inquiry,
                                            text='您所查询的药物与靶标互作用关系为：',
                                            fg='black',
                                            bg='white',
                                            font=(UI_FONT_FAMILY, 26))
    lb_drug_protein_inquiry.place(relx=0.05, rely=0.15, relwidth=0.7, relheight=0.15)
    lb_sign_graph = tkinter.Label(win_drug_protein_inquiry,
                                text=search_drug_protein(flag_drug, flag_Protein),
                                fg='red',
                                bg='white',
                                font=(UI_FONT_FAMILY, 40))
    lb_sign_graph.place(relx=0.4, rely=0.3, relwidth=0.6, relheight=0.2)
    # 查询药物的信息按钮
    btn_drug1_inquiry = create_button(win_drug_protein_inquiry,
                                    text='点击查询药物的相关信息',
                                    command=lambda: drug_inquiry(drug_name),
                                    bg='white',
                                    fg='black',
                                    font=(UI_FONT_FAMILY, 28),
                                    style='secondary')
    btn_drug1_inquiry.place(relx=0.2, rely=0.55, relwidth=0.6, relheight=0.1)
    # 查询靶标的信息按钮
    btn_drug2_inquiry = create_button(win_drug_protein_inquiry,
                                    text='点击查询靶标的关系图',
                                    command=lambda: protein_network(protein_name, flag_Protein),
                                    bg='white',
                                    fg='black',
                                    font=(UI_FONT_FAMILY, 28),
                                    style='secondary')
    btn_drug2_inquiry.place(relx=0.2, rely=0.7, relwidth=0.6, relheight=0.1)
    
    # 添加返回按钮
    btn_return = create_button(win_drug_protein_inquiry,
                               text='返回',
                               command=win_drug_protein_inquiry.destroy,
                               bg='lightgray',
                               fg='black',
                               font=(UI_FONT_FAMILY, 16),
                               style='secondary')
    btn_return.place(relx=0.4, rely=0.85, relwidth=0.2, relheight=0.08)


##########四个模块的输入界面函数##########
def drug(root):  # 单药物输入界面（简化版）
    # 弹出新窗口
    win_drug = tkinter.Toplevel(root)
    win_drug.geometry('500x300')
    win_drug.title('药物信息查询')
    win_drug.config(bg='#F5F5F5')
    win_drug.transient(root)
    win_drug.grab_set()
    win_drug.resizable(False, False)

    # 标题
    title_label = tkinter.Label(win_drug,
                               text='药物信息查询',
                               font=(UI_FONT_FAMILY, 16, 'bold'),
                               bg='#F5F5F5',
                               fg='#333333')
    title_label.pack(pady=20)

    # 输入区域
    input_frame = tkinter.Frame(win_drug, bg='#F5F5F5')
    input_frame.pack(pady=20)

    # 标签
    label = tkinter.Label(input_frame,
                         text='药物名称：',
                         font=(UI_FONT_FAMILY, 12),
                         bg='#F5F5F5',
                         fg='#333333')
    label.pack(side=tkinter.LEFT, padx=10)

    # 输入框
    inp_drug_num = tkinter.Entry(input_frame,
                                 font=(UI_FONT_FAMILY, 12),
                                 width=25,
                                 relief=tkinter.SOLID,
                                 borderwidth=1)
    inp_drug_num.pack(side=tkinter.LEFT, padx=10)
    inp_drug_num.focus()
    _attach_autocomplete(inp_drug_num, _get_drug_name_candidates)

    create_button(
        input_frame,
        text='选择',
        command=lambda: _open_select_dialog(
            win_drug,
            '选择药物',
            _get_drug_name_candidates(),
            lambda v: (inp_drug_num.delete(0, tkinter.END), inp_drug_num.insert(0, v))
        ),
        style='secondary',
        width=6,
        height=1
    ).pack(side=tkinter.LEFT, padx=6)

    # 按钮区域
    button_frame = tkinter.Frame(win_drug, bg='#F5F5F5')
    button_frame.pack(pady=30)
    for col in range(2):
        button_frame.grid_columnconfigure(col, weight=1, uniform='drug_btns')

    # 查询按钮
    btn_search = create_button(
        button_frame,
        text='查询',
        command=lambda: drug_inquiry(inp_drug_num.get()),
        width=12,
        height=1
    )
    btn_search.grid(row=0, column=0, padx=12)

    # 返回按钮
    btn_return = create_button(
        button_frame,
        text='返回',
        command=win_drug.destroy,
        style='secondary',
        width=12,
        height=1
    )
    btn_return.grid(row=0, column=1, padx=12)

    # 绑定回车键
    inp_drug_num.bind('<Return>', lambda e: drug_inquiry(inp_drug_num.get()))


def drug_drug(root):  # 药物-药物输入界面（简化版）
    # 弹出新窗口
    win_drug_drug = tkinter.Toplevel(root)
    win_drug_drug.geometry('500x350')
    win_drug_drug.title('药物-药物分析')
    win_drug_drug.config(bg='#F5F5F5')
    win_drug_drug.transient(root)
    win_drug_drug.grab_set()
    win_drug_drug.resizable(False, False)

    # 标题
    title_label = tkinter.Label(win_drug_drug,
                               text='药物-药物相互作用分析',
                               font=(UI_FONT_FAMILY, 16, 'bold'),
                               bg='#F5F5F5',
                               fg='#333333')
    title_label.pack(pady=20)

    # 输入区域1
    input_frame1 = tkinter.Frame(win_drug_drug, bg='#F5F5F5')
    input_frame1.pack(pady=15)

    label1 = tkinter.Label(input_frame1,
                          text='药物1名称：',
                          font=(UI_FONT_FAMILY, 12),
                          bg='#F5F5F5',
                          fg='#333333',
                          width=12)
    label1.pack(side=tkinter.LEFT, padx=10)

    inp1_drug_drug_num = tkinter.Entry(input_frame1,
                                      font=(UI_FONT_FAMILY, 12),
                                      width=25,
                                      relief=tkinter.SOLID,
                                      borderwidth=1)
    inp1_drug_drug_num.pack(side=tkinter.LEFT, padx=10)
    inp1_drug_drug_num.focus()
    _attach_autocomplete(inp1_drug_drug_num, _get_drug_name_candidates)
    create_button(
        input_frame1,
        text='选择',
        command=lambda: _open_select_dialog(
            win_drug_drug,
            '选择药物1',
            _get_drug_name_candidates(),
            lambda v: (inp1_drug_drug_num.delete(0, tkinter.END), inp1_drug_drug_num.insert(0, v))
        ),
        style='secondary',
        width=6,
        height=1
    ).pack(side=tkinter.LEFT, padx=6)

    # 输入区域2
    input_frame2 = tkinter.Frame(win_drug_drug, bg='#F5F5F5')
    input_frame2.pack(pady=15)

    label2 = tkinter.Label(input_frame2,
                          text='药物2名称：',
                          font=(UI_FONT_FAMILY, 12),
                          bg='#F5F5F5',
                          fg='#333333',
                          width=12)
    label2.pack(side=tkinter.LEFT, padx=10)

    inp2_drug_drug_num = tkinter.Entry(input_frame2,
                                      font=(UI_FONT_FAMILY, 12),
                                      width=25,
                                      relief=tkinter.SOLID,
                                      borderwidth=1)
    inp2_drug_drug_num.pack(side=tkinter.LEFT, padx=10)
    _attach_autocomplete(inp2_drug_drug_num, _get_drug_name_candidates)
    create_button(
        input_frame2,
        text='选择',
        command=lambda: _open_select_dialog(
            win_drug_drug,
            '选择药物2',
            _get_drug_name_candidates(),
            lambda v: (inp2_drug_drug_num.delete(0, tkinter.END), inp2_drug_drug_num.insert(0, v))
        ),
        style='secondary',
        width=6,
        height=1
    ).pack(side=tkinter.LEFT, padx=6)

    # 按钮区域
    button_frame = tkinter.Frame(win_drug_drug, bg='#F5F5F5')
    button_frame.pack(pady=30)
    for col in range(2):
        button_frame.grid_columnconfigure(col, weight=1, uniform='dta_btns')
    for col in range(2):
        button_frame.grid_columnconfigure(col, weight=1, uniform='drug_drug_btns')

    # 查询按钮
    btn_search = create_button(
        button_frame,
        text='查询',
        command=lambda: drug_drug_inquiry(inp1_drug_drug_num.get(), inp2_drug_drug_num.get()),
        width=12,
        height=1
    )
    btn_search.grid(row=0, column=0, padx=12)

    # 返回按钮
    btn_return = create_button(
        button_frame,
        text='返回',
        command=win_drug_drug.destroy,
        style='secondary',
        width=12,
        height=1
    )
    btn_return.grid(row=0, column=1, padx=12)

    # 绑定回车键
    inp1_drug_drug_num.bind('<Return>', lambda e: inp2_drug_drug_num.focus())
    inp2_drug_drug_num.bind('<Return>', lambda e: drug_drug_inquiry(inp1_drug_drug_num.get(), inp2_drug_drug_num.get()))


def dta_predict(root):  # DTA预测输入界面（简化版）
    # 弹出新窗口
    win_drug_drug = tkinter.Toplevel(root)
    win_drug_drug.geometry('500x350')
    win_drug_drug.title('DTA预测')
    win_drug_drug.config(bg='#F5F5F5')
    win_drug_drug.transient(root)
    win_drug_drug.grab_set()
    win_drug_drug.resizable(False, False)

    # 标题
    title_label = tkinter.Label(win_drug_drug,
                               text='药物-靶标亲和力预测',
                               font=(UI_FONT_FAMILY, 16, 'bold'),
                               bg='#F5F5F5',
                               fg='#333333')
    title_label.pack(pady=20)

    # 输入区域1
    input_frame1 = tkinter.Frame(win_drug_drug, bg='#F5F5F5')
    input_frame1.pack(pady=15)

    label1 = tkinter.Label(input_frame1,
                          text='药物名称：',
                          font=(UI_FONT_FAMILY, 12),
                          bg='#F5F5F5',
                          fg='#333333',
                          width=12)
    label1.pack(side=tkinter.LEFT, padx=10)

    inp1_drug_drug_num = tkinter.Entry(input_frame1,
                                      font=(UI_FONT_FAMILY, 12),
                                      width=25,
                                      relief=tkinter.SOLID,
                                      borderwidth=1)
    inp1_drug_drug_num.pack(side=tkinter.LEFT, padx=10)
    inp1_drug_drug_num.focus()
    _attach_autocomplete(inp1_drug_drug_num, _get_drug_name_candidates)
    create_button(
        input_frame1,
        text='选择',
        command=lambda: _open_select_dialog(
            win_drug_drug,
            '选择药物',
            _get_drug_name_candidates(),
            lambda v: (inp1_drug_drug_num.delete(0, tkinter.END), inp1_drug_drug_num.insert(0, v))
        ),
        style='secondary',
        width=6,
        height=1
    ).pack(side=tkinter.LEFT, padx=6)

    # 输入区域2
    input_frame2 = tkinter.Frame(win_drug_drug, bg='#F5F5F5')
    input_frame2.pack(pady=15)

    label2 = tkinter.Label(input_frame2,
                          text='靶标名称：',
                          font=(UI_FONT_FAMILY, 12),
                          bg='#F5F5F5',
                          fg='#333333',
                          width=12)
    label2.pack(side=tkinter.LEFT, padx=10)

    inp2_drug_drug_num = tkinter.Entry(input_frame2,
                                      font=(UI_FONT_FAMILY, 12),
                                      width=25,
                                      relief=tkinter.SOLID,
                                      borderwidth=1)
    inp2_drug_drug_num.pack(side=tkinter.LEFT, padx=10)
    _attach_autocomplete(inp2_drug_drug_num, _get_protein_name_candidates)
    create_button(
        input_frame2,
        text='选择',
        command=lambda: _open_select_dialog(
            win_drug_drug,
            '选择靶标（支持别名）',
            _get_protein_name_candidates(),
            lambda v: (inp2_drug_drug_num.delete(0, tkinter.END), inp2_drug_drug_num.insert(0, v))
        ),
        style='secondary',
        width=10,
        height=1
    ).pack(side=tkinter.LEFT, padx=6)

    # 按钮区域
    button_frame = tkinter.Frame(win_drug_drug, bg='#F5F5F5')
    button_frame.pack(pady=30)

    def _submit_dta_predict():
        drug_name = inp1_drug_drug_num.get().strip()
        protein_name = inp2_drug_drug_num.get().strip()
        if not drug_name or not protein_name:
            messagebox.showwarning(title='错误提示', message='请输入药物名称和靶向名称')
            return
        drug_prot_dta_predict(drug_name, protein_name)

    # 预测按钮
    btn_predict = create_button(
        button_frame,
        text='预测',
        command=_submit_dta_predict,
        width=12,
        height=1
    )
    btn_predict.grid(row=0, column=0, padx=12)

    # 返回按钮
    btn_return = create_button(
        button_frame,
        text='返回',
        command=win_drug_drug.destroy,
        style='secondary',
        width=12,
        height=1
    )
    btn_return.grid(row=0, column=1, padx=12)

    # 绑定回车键
    inp1_drug_drug_num.bind('<Return>', lambda e: inp2_drug_drug_num.focus())
    inp2_drug_drug_num.bind('<Return>', lambda e: _submit_dta_predict())


def dta_ddi_compare(root):  # DTA与DDI预测对比界面（左右对比）
    win_compare = tkinter.Toplevel(root)
    win_compare.geometry('980x520')
    win_compare.title('DTA vs DDI 预测对比')
    win_compare.config(bg='#F5F5F5')
    win_compare.transient(root)
    win_compare.grab_set()
    win_compare.resizable(False, False)

    title_label = tkinter.Label(
        win_compare,
        text='DTA / DDI 预测对比',
        font=(UI_FONT_FAMILY, 16, 'bold'),
        bg='#F5F5F5',
        fg='#333333'
    )
    title_label.pack(pady=15)

    container = tkinter.Frame(win_compare, bg='#F5F5F5')
    container.pack(fill=tkinter.BOTH, expand=True, padx=20, pady=10)
    container.grid_columnconfigure(0, weight=1, uniform='compare_cols')
    container.grid_columnconfigure(1, weight=1, uniform='compare_cols')
    container.grid_rowconfigure(0, weight=1)

    left_frame = tkinter.Frame(container, bg='white', relief=tkinter.SOLID, borderwidth=1)
    left_frame.grid(row=0, column=0, sticky='nsew', padx=(0, 10))

    right_frame = tkinter.Frame(container, bg='white', relief=tkinter.SOLID, borderwidth=1)
    right_frame.grid(row=0, column=1, sticky='nsew', padx=(10, 0))

    # 左侧：DTA预测
    tkinter.Label(
        left_frame,
        text='DTA预测',
        font=(UI_FONT_FAMILY, 14, 'bold'),
        bg='white',
        fg='#333333'
    ).pack(pady=(12, 8))

    left_input1 = tkinter.Frame(left_frame, bg='white')
    left_input1.pack(pady=6, padx=12, fill=tkinter.X)
    tkinter.Label(left_input1, text='药物名称：', font=(UI_FONT_FAMILY, 11), bg='white', fg='#333333', width=10, anchor='w').pack(side=tkinter.LEFT)
    dta_drug_entry = tkinter.Entry(left_input1, font=(UI_FONT_FAMILY, 11), relief=tkinter.SOLID, borderwidth=1)
    dta_drug_entry.pack(side=tkinter.LEFT, fill=tkinter.X, expand=True)
    _attach_autocomplete(dta_drug_entry, _get_drug_name_candidates)
    create_button(
        left_input1,
        text='选择',
        command=lambda: _open_select_dialog(
            win_compare,
            '选择药物',
            _get_drug_name_candidates(),
            lambda v: (dta_drug_entry.delete(0, tkinter.END), dta_drug_entry.insert(0, v))
        ),
        style='secondary',
        width=6,
        height=1
    ).pack(side=tkinter.LEFT, padx=6)

    left_input2 = tkinter.Frame(left_frame, bg='white')
    left_input2.pack(pady=6, padx=12, fill=tkinter.X)
    tkinter.Label(left_input2, text='靶标名称：', font=(UI_FONT_FAMILY, 11), bg='white', fg='#333333', width=10, anchor='w').pack(side=tkinter.LEFT)
    dta_prot_entry = tkinter.Entry(left_input2, font=(UI_FONT_FAMILY, 11), relief=tkinter.SOLID, borderwidth=1)
    dta_prot_entry.pack(side=tkinter.LEFT, fill=tkinter.X, expand=True)
    _attach_autocomplete(dta_prot_entry, _get_protein_name_candidates)
    create_button(
        left_input2,
        text='选择',
        command=lambda: _open_select_dialog(
            win_compare,
            '选择可预测靶标',
            _get_davis_prot_ids(),
            lambda v: (dta_prot_entry.delete(0, tkinter.END), dta_prot_entry.insert(0, v))
        ),
        style='secondary',
        width=10,
        height=1
    ).pack(side=tkinter.LEFT, padx=6)

    dta_result_label = tkinter.Label(
        left_frame,
        text='结果：-       ',
        font=(UI_FONT_FAMILY, 12, 'bold'),
        bg='#E8F4F8',
        fg='#2C3E50',
        relief=tkinter.SOLID,
        borderwidth=1,
        padx=10,
        pady=10
    )
    dta_result_label.pack(pady=15, padx=12, fill=tkinter.X)

    def _run_dta_predict():
        drug_name = dta_drug_entry.get().strip()
        protein_name = dta_prot_entry.get().strip()
        if not drug_name or not protein_name:
            messagebox.showwarning(title='错误提示', message='请填写药物名称与靶标名称')
            return
        try:
            df = pd.read_csv('EGFR-Case/drug.tsv', sep="\t")
            df1 = pd.read_csv('davis_prots.csv', sep=",")
        except Exception as exc:
            messagebox.showwarning(title='错误提示', message=f'读取数据文件失败：{exc}')
            return
        if not df.iloc[:, 0].str.contains(drug_name, na=False).any():
            messagebox.showwarning(title='错误提示', message='您输入的药物名称有误，请重新输入！')
            return
        prot_ids = df1.iloc[:, 0].tolist()
        resolved_prot, prot_suggestions = _resolve_protein_id(protein_name, prot_ids)
        if not resolved_prot:
            if _is_in_proteinnameid(protein_name, prot_ids):
                tip = '该名称存在于 ProteinNameID，但不在 davis_prots.csv 中，无法进行当前模型预测。请改用 davis_prots.csv 的 prot_id（第一列）。'
            else:
                tip = '靶标名称需使用 davis_prots.csv 的 prot_id（第一列），或数据库 ProteinNameID 中的名称。'
            if prot_suggestions:
                tip += f"\n可能的候选：{', '.join(prot_suggestions)}"
            messagebox.showwarning(title='错误提示', message=f'您输入的靶标名称有误，请重新输入！\n{tip}')
            return
        try:
            data_extractor.find_drug(drug_name, 'EGFR-Case/drug.tsv')
            data_extractor.find_protein(resolved_prot, 'davis_prots.csv')
        except Exception as exc:
            messagebox.showwarning(title='错误提示', message=f'准备数据失败：{exc}')
            return
        
        start_time = time.time()
        dta_result_label.config(text='正在更新预测（预计耗时约30秒），请稍候...')

        def _load_pred_result():
            pred_path = _get_pred_output_path()
            elapsed_time = time.time() - start_time
            try:
                if not os.path.exists(pred_path):
                    dta_result_label.config(text='预测结果文件不存在')
                    return
                df_pred = pd.read_csv(pred_path)
                if not df_pred.empty and 'drug_id' in df_pred.columns and 'prot_id' in df_pred.columns:
                    matched = df_pred[(df_pred['drug_id'] == drug_name) & (df_pred['prot_id'] == resolved_prot)]
                else:
                    matched = pd.DataFrame()
                if not matched.empty:
                    if 'pred' in matched.columns:
                        dtan = matched['pred'].iloc[0]
                    else:
                        dtan = matched.iloc[0, -1]
                elif df_pred.empty:
                    dtan = 'N/A'
                elif 'pred' in df_pred.columns:
                    dtan = df_pred['pred'].iloc[0]
                else:
                    dtan = df_pred.iloc[0, -1]
            except Exception:
                dtan = 'N/A'
            dta_result_label.config(text=f'预测DTA值：{dtan} (耗时: {elapsed_time:.1f}s)')

        def _handle_pred_error(exc):
            dta_result_label.config(text='预测失败，请查看环境或数据配置')
            detail = _LAST_PRED_LOG.strip()
            if detail:
                messagebox.showwarning(title='错误提示', message=f'运行预测失败：{exc}\n\n{detail}')
            else:
                messagebox.showwarning(title='错误提示', message=f'运行预测失败：{exc}')

        _run_pred_update_async(on_success=_load_pred_result, on_error=_handle_pred_error)

    create_button(
        left_frame,
        text='开始预测',
        command=_run_dta_predict,
        width=12,
        height=1
    ).pack(pady=6)

    # 右侧：DDI预测
    tkinter.Label(
        right_frame,
        text='DDI符号关系预测',
        font=(UI_FONT_FAMILY, 14, 'bold'),
        bg='white',
        fg='#333333'
    ).pack(pady=(12, 8))

    right_input1 = tkinter.Frame(right_frame, bg='white')
    right_input1.pack(pady=6, padx=12, fill=tkinter.X)
    tkinter.Label(right_input1, text='药物1名称：', font=(UI_FONT_FAMILY, 11), bg='white', fg='#333333', width=10, anchor='w').pack(side=tkinter.LEFT)
    ddi_drug1_entry = tkinter.Entry(right_input1, font=(UI_FONT_FAMILY, 11), relief=tkinter.SOLID, borderwidth=1)
    ddi_drug1_entry.pack(side=tkinter.LEFT, fill=tkinter.X, expand=True)
    _attach_autocomplete(ddi_drug1_entry, _get_drug_name_candidates)
    create_button(
        right_input1,
        text='选择',
        command=lambda: _open_select_dialog(
            win_compare,
            '选择药物1',
            _get_drug_name_candidates(),
            lambda v: (ddi_drug1_entry.delete(0, tkinter.END), ddi_drug1_entry.insert(0, v))
        ),
        style='secondary',
        width=6,
        height=1
    ).pack(side=tkinter.LEFT, padx=6)

    right_input2 = tkinter.Frame(right_frame, bg='white')
    right_input2.pack(pady=6, padx=12, fill=tkinter.X)
    tkinter.Label(right_input2, text='药物2名称：', font=(UI_FONT_FAMILY, 11), bg='white', fg='#333333', width=10, anchor='w').pack(side=tkinter.LEFT)
    ddi_drug2_entry = tkinter.Entry(right_input2, font=(UI_FONT_FAMILY, 11), relief=tkinter.SOLID, borderwidth=1)
    ddi_drug2_entry.pack(side=tkinter.LEFT, fill=tkinter.X, expand=True)
    _attach_autocomplete(ddi_drug2_entry, _get_drug_name_candidates)
    create_button(
        right_input2,
        text='选择',
        command=lambda: _open_select_dialog(
            win_compare,
            '选择药物2',
            _get_drug_name_candidates(),
            lambda v: (ddi_drug2_entry.delete(0, tkinter.END), ddi_drug2_entry.insert(0, v))
        ),
        style='secondary',
        width=6,
        height=1
    ).pack(side=tkinter.LEFT, padx=6)

    ddi_result_label = tkinter.Label(
        right_frame,
        text='结果：-       ',
        font=(UI_FONT_FAMILY, 12, 'bold'),
        bg='#E8F4F8',
        fg='#2C3E50',
        relief=tkinter.SOLID,
        borderwidth=1,
        padx=10,
        pady=10
    )
    ddi_result_label.pack(pady=15, padx=12, fill=tkinter.X)

    def _run_ddi_predict():
        drug1_name = ddi_drug1_entry.get().strip()
        drug2_name = ddi_drug2_entry.get().strip()
        if not drug1_name or not drug2_name:
            messagebox.showwarning(title='错误提示', message='请填写两个药物名称')
            return
        names = datareader("NamesWithID")
        drug_names = names[:, 0].tolist()
        if drug1_name not in drug_names or drug2_name not in drug_names:
            messagebox.showwarning(title='错误提示', message='您输入的药物名称有误，请重新输入！')
            return
        flag1 = drug_names.index(drug1_name)
        flag2 = drug_names.index(drug2_name)
        result = prediction([flag1], [flag2])
        result_text = result[0] if result else 'N/A'
        ddi_result_label.config(text=f'预测DDI关系：{result_text}')

    create_button(
        right_frame,
        text='开始预测',
        command=_run_ddi_predict,
        width=12,
        height=1
    ).pack(pady=6)

    btn_return = create_button(
        win_compare,
        text='返回',
        command=win_compare.destroy,
        style='secondary',
        width=12,
        height=1
    )
    btn_return.pack(pady=12)


##########响应函数##########
def Mouse_over(event, canvas_main, rectangle_1, rectangle_2, text):  # 鼠标经过响应函数
    x = [130, 380, 630]
    y = [260, 260, 260]
    l, h, r = 200, 60, 5
    if x[0] <= event.x <= x[0] + l and y[0] <= event.y <= y[0] + h:  # 响应的位置
        canvas_main.itemconfigure(rectangle_1[0], outline='black')  # 重设外框颜色
        canvas_main.itemconfigure(rectangle_2[0], outline='black')  # 重设内框颜色
        canvas_main.itemconfigure(text[0], fill='black')  # 重设显示文本颜色
        canvas_main.configure(cursor='hand2')  # 重设鼠标样式
    elif x[1] <= event.x <= x[1] + l and y[1] <= event.y <= y[1] + h:  # 响应的位置
        canvas_main.itemconfigure(rectangle_1[1], outline='black')  # 重设外框颜色
        canvas_main.itemconfigure(rectangle_2[1], outline='black')  # 重设内框颜色
        canvas_main.itemconfigure(text[1], fill='black')  # 重设显示文本颜色
        canvas_main.configure(cursor='hand2')  # 重设鼠标样式
    elif x[2] <= event.x <= x[2] + l and y[2] <= event.y <= y[2] + h:  # 响应的位置
        canvas_main.itemconfigure(rectangle_1[2], outline='black')  # 重设外框颜色
        canvas_main.itemconfigure(rectangle_2[2], outline='black')  # 重设内框颜色
        canvas_main.itemconfigure(text[2], fill='black')  # 重设显示文本颜色
        canvas_main.configure(cursor='hand2')  # 重设鼠标样式
    else:
        canvas_main.itemconfigure(rectangle_1[0], outline='white')  # 重设外框颜色
        canvas_main.itemconfigure(rectangle_2[0], outline='white')  # 重设内框颜色
        canvas_main.itemconfigure(text[0], fill='white')  # 重设显示文本颜色
        canvas_main.itemconfigure(rectangle_1[1], outline='white')  # 重设外框颜色
        canvas_main.itemconfigure(rectangle_2[1], outline='white')  # 重设内框颜色
        canvas_main.itemconfigure(text[1], fill='white')  # 重设显示文本颜色
        canvas_main.itemconfigure(rectangle_1[2], outline='white')  # 重设外框颜色
        canvas_main.itemconfigure(rectangle_2[2], outline='white')  # 重设内框颜色
        canvas_main.itemconfigure(text[2], fill='white')  # 重设显示文本颜色


def Mouse_Click(event, root):  # 点击响应函数
    x = [130, 380, 630]
    y = [260, 260, 260]
    l, h, r = 200, 60, 5

    command = [drug, drug_drug, dta_predict]
    if x[0] <= event.x <= x[0] + l and y[0] <= event.y <= y[0] + h:  # 响应的位置
        command[0](root)
    elif x[1] <= event.x <= x[1] + l and y[1] <= event.y <= y[1] + h:  # 响应的位置
        command[1](root)
    elif x[2] <= event.x <= x[2] + l and y[2] <= event.y <= y[2] + h:  # 响应的位置
        command[2](root)


def help():  # 帮助响应函数
    parent = main_root if main_root else start
    win_help = tkinter.Toplevel(parent) if parent else tkinter.Tk()
    setup_tk_fonts(win_help, base_size=12)
    win_help.geometry('500x650')
    win_help.config(background="white")
    win_help.title('帮助文档')
    if parent:
        win_help.transient(parent)
    help_path = os.path.join(Config.BASE_DIR, 'help_interduction.txt')
    try:
        with open(help_path, encoding='utf-8') as file:
            content = file.read()
    except FileNotFoundError:
        content = (
            f'{Config.SOFTWARE_NAME}\n'
            f'{Config.SOFTWARE_SUBTITLE}\n\n'
            '主要功能：\n'
            '1. 药物信息查询：检索药物编号、副作用、化学特征和分子结构。\n'
            '2. 药物-药物分析：查询或预测 DDI 的协同/拮抗关系并展示符号图。\n'
            '3. DTA 预测：调用 LLMDTA 模型预测药物-靶标亲和力。\n'
            '4. 智能助手：将自然语言问题解析为结构化任务，聚合证据并导出报告。\n\n'
            '首次使用请先运行 setup.ps1，配置 MySQL、依赖和可选的 LLM API Key。'
        )
    lb_help = tkinter.Label(win_help,
                            text='帮助使用文档',
                            fg='red',
                            bg='white',
                            font=(UI_FONT_FAMILY, 28),
                            relief=FLAT)
    lb_help.place(relx=0, rely=0, relwidth=1, relheight=0.15)
    txt_help_interduction = tkinter.Text(win_help,
                                        fg='black',
                                        bg='white',
                                        font=(UI_FONT_FAMILY, 18),
                                        relief=FLAT)
    txt_help_interduction.place(relx=0, rely=0.15, relwidth=1, relheight=0.85)
    txt_help_interduction.insert(END, content + '\n')
    # 创建一个滚动条控件，默认为垂直方向
    sbar1_help = tkinter.Scrollbar(win_help)
    # 将滚动条放置在右侧，并设置当窗口大小改变时滚动条会沿着垂直方向延展
    sbar1_help.pack(side=RIGHT, fill=Y)
    # 使用 command 关联控件的 yview方法
    sbar1_help.config(command=txt_help_interduction.yview)
    
    # 添加返回按钮
    btn_return = create_button(win_help,
                               text='返回',
                               command=win_help.destroy,
                               bg='lightgray',
                               fg='black',
                               font=(UI_FONT_FAMILY, 16),
                               style='secondary')
    btn_return.place(relx=0.4, rely=0.92, relwidth=0.2, relheight=0.06)


def ROOT():
    global bg, main_root
    root = tkinter.Tk()
    setup_tk_fonts(root, base_size=12)
    main_root = root  # 保存主窗口引用
    root.title(software_name)
    root.geometry(Config.MAIN_WINDOW_SIZE)  # 简化后的窗口大小
    root.config(bg='#F5F5F5')  # 简洁的浅灰色背景
    root.resizable(False, False)

    def return_to_start():
        root.destroy()
        try:
            start.deiconify()
        except Exception:
            pass

    root.protocol('WM_DELETE_WINDOW', return_to_start)
    root.bind('<Escape>', lambda _event: return_to_start())

    # 顶部按钮区（替代菜单栏）
    top_bar = tkinter.Frame(root, bg='#F5F5F5')
    top_bar.pack(fill=tkinter.X, padx=10, pady=10)
    create_button(top_bar, text='返回', command=return_to_start, style='secondary', width=8, font=(UI_FONT_FAMILY, 13), animate=False).pack(side=tkinter.LEFT)
    create_button(top_bar, text='帮助', command=help, width=8, font=(UI_FONT_FAMILY, 13), animate=False).pack(side=tkinter.LEFT, padx=8)
    create_button(top_bar, text='智能助手', command=lambda: open_agent_window(root), width=10, font=(UI_FONT_FAMILY, 13), animate=False).pack(side=tkinter.LEFT, padx=8)
    create_button(top_bar, text='关闭', command=_exit_app, style='secondary', width=8, font=(UI_FONT_FAMILY, 13), animate=False).pack(side=tkinter.RIGHT)

    # 标题（简化）
    title_label = tkinter.Label(root,
                               text=software_name,
                               font=(UI_FONT_FAMILY, 24, 'bold'),
                               bg='#F5F5F5',
                               fg='#333333')
    title_label.pack(pady=(32, 6))
    subtitle_label = tkinter.Label(root,
                                  text=Config.SOFTWARE_SUBTITLE,
                                  font=(UI_FONT_FAMILY, 13),
                                  bg='#F5F5F5',
                                  fg='#64748B')
    subtitle_label.pack()

    # 按钮容器
    button_frame = tkinter.Frame(root, bg='#F5F5F5')
    button_frame.pack(pady=30)
    for col in range(3):
        button_frame.grid_columnconfigure(col, weight=1, uniform='main_buttons')

    # 三个主要功能按钮（等宽等距）
    btn_drug = create_button(
        button_frame,
        text='药物信息查询',
        command=lambda: drug(root),
        font=(UI_FONT_FAMILY, 16),
        width=18,
        height=2
    )
    btn_drug.grid(row=0, column=0, padx=22, pady=10)

    btn_drug_drug = create_button(
        button_frame,
        text='药物-药物分析',
        command=lambda: drug_drug(root),
        font=(UI_FONT_FAMILY, 16),
        width=18,
        height=2
    )
    btn_drug_drug.grid(row=0, column=1, padx=22, pady=10)

    btn_dta = create_button(
        button_frame,
        text='DTA预测',
        command=lambda: dta_predict(root),
        font=(UI_FONT_FAMILY, 16),
        width=18,
        height=2
    )
    btn_dta.grid(row=0, column=2, padx=22, pady=10)

    btn_compare = create_button(
        button_frame,
        text='DTA/DDI对比',
        command=lambda: dta_ddi_compare(root),
        font=(UI_FONT_FAMILY, 16),
        width=22,
        height=2
    )
    btn_compare.grid(row=1, column=0, columnspan=3, pady=12)

    # 说明文字（简化）
    info_label = tkinter.Label(root,
                              text='请选择功能模块',
                              font=(UI_FONT_FAMILY, 14),
                              bg='#F5F5F5',
                              fg='#666666')
    info_label.pack(pady=20)

def Mouse_Click_start(event):  # 关联鼠标点击事件
    global root
    start_x, start_y, start_l, start_h, start_r = 600, 350, 230, 50, 5
    start_help_x, start_help_y, start_help_l, start_help_h, start_help_r = 600, 420, 230, 50, 5
    if start_x <= event.x <= start_x + start_l and start_y <= event.y <= start_y + start_h:  # 响应的位置
        start.withdraw()
        ROOT()
    elif start_help_x <= event.x <= start_help_x + start_help_l and start_help_y <= event.y <= start_help_y + start_help_h:  # 响应的位置
        help()


def Mouse_over_start(event, canvas_drug, r11, r12, t1, r21, r22, t2):  # 关联鼠标经过事件
    start_x, start_y, start_l, start_h, start_r = 600, 350, 230, 50, 5
    start_help_x, start_help_y, start_help_l, start_help_h, start_help_r = 600, 420, 230, 50, 5
    if start_x <= event.x <= start_x + start_l and start_y <= event.y <= start_y + start_h:  # 响应的位置
        if not _start_hover['main']:
            _start_hover['main'] = True
            animate_start_button(canvas_drug, r11, r12, t1, 'main')
        _start_hover['help'] = False
        canvas_drug.configure(cursor='hand2')  # 重设鼠标样式
    elif start_help_x <= event.x <= start_help_x + start_help_l and start_help_y <= event.y <= start_help_y + start_help_h:  # 响应的位置
        if not _start_hover['help']:
            _start_hover['help'] = True
            animate_start_button(canvas_drug, r21, r22, t2, 'help')
        _start_hover['main'] = False
        canvas_drug.configure(cursor='hand2')  # 重设鼠标样式
    else:
        _start_hover['main'] = False
        _start_hover['help'] = False
        for item in r11:
            canvas_drug.itemconfigure(item, outline=START_BTN_OUTLINE, fill=START_BTN_FILL)  # 重设外框颜色
        for item in r12:
            canvas_drug.itemconfigure(item, outline=START_BTN_OUTLINE, fill=START_BTN_FILL)  # 重设内框颜色
        canvas_drug.itemconfigure(t1, fill=START_BTN_TEXT)  # 重设显示文本颜色
        for item in r21:
            canvas_drug.itemconfigure(item, outline=START_BTN_OUTLINE, fill=START_BTN_FILL)  # 重设外框颜色
        for item in r22:
            canvas_drug.itemconfigure(item, outline=START_BTN_OUTLINE, fill=START_BTN_FILL)  # 重设内框颜色
        canvas_drug.itemconfigure(t2, fill=START_BTN_TEXT)  # 重设显示文本颜色
        canvas_drug.configure(cursor='arrow')  # 恢复默认


# 视觉风格
START_BTN_FILL = '#A5B4FC'
START_BTN_OUTLINE = '#E0E7FF'
START_BTN_TEXT = '#FFFFFF'
START_TITLE_COLOR = '#0F172A'
START_TITLE_SHADOW = '#F8FAFC'
START_BTN_GRADIENT = ['#A5B4FC', '#93C5FD', '#60A5FA', '#3B82F6', '#60A5FA', '#93C5FD']
START_BTN_OUTLINE_HOVER = '#F1F5FF'
START_BTN_TEXT_HOVER = '#FFFFFF'

_start_hover = {'main': False, 'help': False}
_start_anim_id = {'main': None, 'help': None}
_start_anim_idx = {'main': 0, 'help': 0}
start = None

def create_title_with_shadow(canvas, x, y, text, font):
    canvas.create_text(x + 2, y + 2, text=text, fill=START_TITLE_SHADOW, font=font)
    return canvas.create_text(x, y, text=text, fill=START_TITLE_COLOR, font=font)

def create_start_button(canvas, x, y, w, h, r, text):
    rect1 = _create_rounded_rect(
        canvas, x, y, x + w, y + h, r,
        width=2, outline=START_BTN_OUTLINE, fill=START_BTN_FILL
    )
    rect2 = _create_rounded_rect(
        canvas, x + r, y + r, x + w - r, y + h - r, max(2, r - 2),
        width=1.5, outline=START_BTN_OUTLINE, fill=START_BTN_FILL
    )
    text_id = canvas.create_text(x + 0.5 * w, y + 0.5 * h, text=text,
                                 font=(UI_FONT_FAMILY, 24, 'bold'), fill=START_BTN_TEXT)
    return rect1, rect2, text_id

def animate_start_button(canvas, rect1, rect2, text_id, key):
    if not _start_hover[key]:
        return
    color = START_BTN_GRADIENT[-1]
    for item in rect1:
        canvas.itemconfigure(item, fill=color, outline=START_BTN_OUTLINE_HOVER)
    for item in rect2:
        canvas.itemconfigure(item, fill=color, outline=START_BTN_OUTLINE_HOVER)
    canvas.itemconfigure(text_id, fill=START_BTN_TEXT_HOVER)


def main():
    """启动 Tkinter 入口窗口。"""
    global start

    start = Tk()  # 创建Tk控件
    setup_tk_fonts(start, base_size=12)
    start.geometry(Config.START_WINDOW_SIZE)  # 设置窗口大小及位置
    start.title(software_name)  # 设置窗口标题
    start.protocol('WM_DELETE_WINDOW', _exit_app)

    canvas_start = Canvas(start, highlightthickness=0)  # 创建Canvas控件，并设置边框厚度为0
    canvas_start.place(width=1000, height=670)  # 设置Canvas控件大小及位置
    bg_path = os.path.join(Config.PIC_DIR, '开始2.png')
    bg = PhotoImage(file=bg_path)  # 导入背景图
    canvas_start.create_image(496, 328, image=bg)  # 添加背景图片

    create_title_with_shadow(canvas_start, 730, 210, '药研智析', (UI_FONT_FAMILY, 38, 'bold'))
    create_title_with_shadow(canvas_start, 730, 260, 'DrugReasoner', (UI_FONT_FAMILY, 30, 'bold'))
    create_title_with_shadow(canvas_start, 730, 306, '语言理解 × 符号关系推理', (UI_FONT_FAMILY, 18, 'bold'))

    start_x, start_y, start_l, start_h, start_r = 585, 350, 260, 56, 6
    start_help_x, start_help_y, start_help_l, start_help_h, start_help_r = 585, 420, 260, 56, 6
    start_rectangle_1, start_rectangle_2, start_label_text = create_start_button(
        canvas_start, start_x, start_y, start_l, start_h, start_r, '点击进入应用'
    )
    start_help_rectangle_1, start_help_rectangle_2, start_help_label_text = create_start_button(
        canvas_start, start_help_x, start_help_y, start_help_l, start_help_h, start_help_r, '点击阅读帮助'
    )

    canvas_start.bind('<Button-1>', lambda event: Mouse_Click_start(event))  # 关联鼠标点击事件
    canvas_start.bind('<Motion>', lambda event: Mouse_over_start(event, canvas_start, start_rectangle_1, start_rectangle_2,
                                                                start_label_text,
                                                                start_help_rectangle_1, start_help_rectangle_2,
                                                                start_help_label_text))  # 关联鼠标经过事件

    start.mainloop()


if __name__ == '__main__':
    main()
