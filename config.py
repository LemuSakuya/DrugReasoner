"""
项目配置文件
统一管理所有配置项
"""
import os


class Config:
    """项目配置类"""
    
    # 数据库配置
    DB_HOST = os.getenv('DB_HOST', 'localhost')
    DB_USER = os.getenv('DB_USER', 'root')
    DB_PASSWORD = os.getenv('DB_PASSWORD', os.getenv('MYSQL_ROOT_PASSWORD', '12345'))
    DB_NAME = os.getenv('DB_NAME', 'drug_discovery')
    DB_PORT = int(os.getenv('DB_PORT', 3306))
    DB_CHARSET = 'utf8mb4'
    
    @classmethod
    def get_db_config(cls):
        """获取数据库配置字典"""
        return {
            'host': cls.DB_HOST,
            'user': cls.DB_USER,
            'password': cls.DB_PASSWORD,
            'database': cls.DB_NAME,
            'port': cls.DB_PORT,
            'charset': cls.DB_CHARSET
        }
    
    # 文件路径配置
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    DATA_DIR = os.path.join(BASE_DIR, 'data')
    RAW_DATA_DIR = os.path.join(DATA_DIR, 'raw')
    PROCESSED_DATA_DIR = os.path.join(DATA_DIR, 'processed')
    CACHE_DATA_DIR = os.path.join(DATA_DIR, 'cache')
    CASES_DATA_DIR = os.path.join(DATA_DIR, 'cases')
    EXPORTS_DATA_DIR = os.path.join(DATA_DIR, 'exports')

    # 原始数据目录
    RAW_DTI_DIR = os.path.join(RAW_DATA_DIR, 'dti')
    RAW_DDI_DIR = os.path.join(RAW_DATA_DIR, 'ddi')
    RAW_ENTITY_DIR = os.path.join(RAW_DATA_DIR, 'entity')
    RAW_KNOWLEDGE_DIR = os.path.join(RAW_DATA_DIR, 'knowledge')

    RAW_DAVIS_DIR = os.path.join(RAW_DTI_DIR, 'davis')
    RAW_KIBA_DIR = os.path.join(RAW_DTI_DIR, 'kiba')
    RAW_BINDINGDB_DIR = os.path.join(RAW_DTI_DIR, 'bindingdb')

    RAW_DDINTER_DIR = os.path.join(RAW_DDI_DIR, 'ddinter')
    RAW_DRUGBANK_DIR = os.path.join(RAW_DDI_DIR, 'drugbank')
    RAW_MECDDI_DIR = os.path.join(RAW_DDI_DIR, 'mecddi')
    RAW_OFFSIDES_DIR = os.path.join(RAW_DDI_DIR, 'offsides')

    RAW_DRUG_ALIAS_DIR = os.path.join(RAW_ENTITY_DIR, 'drug_alias')
    RAW_PROTEIN_ALIAS_DIR = os.path.join(RAW_ENTITY_DIR, 'protein_alias')

    RAW_DGIDB_DIR = os.path.join(RAW_KNOWLEDGE_DIR, 'dgidb')
    RAW_CTD_DIR = os.path.join(RAW_KNOWLEDGE_DIR, 'ctd')
    RAW_STITCH_DIR = os.path.join(RAW_KNOWLEDGE_DIR, 'stitch')

    # 处理后数据目录
    PROCESSED_QUERY_DIR = os.path.join(PROCESSED_DATA_DIR, 'query')
    PROCESSED_PREDICTION_DIR = os.path.join(PROCESSED_DATA_DIR, 'prediction')
    PROCESSED_EXPLANATION_DIR = os.path.join(PROCESSED_DATA_DIR, 'explanation')

    QUERY_KNOWN_DDI_DIR = os.path.join(PROCESSED_QUERY_DIR, 'known_ddi')
    QUERY_KNOWN_DTI_DIR = os.path.join(PROCESSED_QUERY_DIR, 'known_dti')
    QUERY_ENTITY_MAP_DIR = os.path.join(PROCESSED_QUERY_DIR, 'entity_map')

    PREDICTION_DTI_DIR = os.path.join(PROCESSED_PREDICTION_DIR, 'dti')
    PREDICTION_DDI_DIR = os.path.join(PROCESSED_PREDICTION_DIR, 'ddi')

    EXPLANATION_CORPORA_DIR = os.path.join(PROCESSED_EXPLANATION_DIR, 'corpora')
    EXPLANATION_TEMPLATES_DIR = os.path.join(PROCESSED_EXPLANATION_DIR, 'templates')

    # 缓存目录
    CACHE_EMBEDDINGS_DIR = os.path.join(CACHE_DATA_DIR, 'embeddings')
    CACHE_FEATURES_DIR = os.path.join(CACHE_DATA_DIR, 'features')
    CACHE_RETRIEVAL_DIR = os.path.join(CACHE_DATA_DIR, 'retrieval')

    # 演示与案例目录
    CASE_GENERAL_DIR = os.path.join(CASES_DATA_DIR, 'general')
    CASE_EGFR_DIR = os.path.join(CASES_DATA_DIR, 'egfr')
    CASE_DEMO_DIR = os.path.join(CASES_DATA_DIR, 'demo')

    # 导出目录
    EXPORTS_CSV_DIR = os.path.join(EXPORTS_DATA_DIR, 'csv')
    EXPORTS_JSON_DIR = os.path.join(EXPORTS_DATA_DIR, 'json')
    EXPORTS_LOGS_DIR = os.path.join(EXPORTS_DATA_DIR, 'logs')
    REPORTS_DIR = os.path.join(BASE_DIR, 'reports')

    MODELS_DIR = os.path.join(BASE_DIR, 'savemodel')
    OUTPUT_DIR = os.path.join(BASE_DIR, 'output')
    PIC_DIR = os.path.join(BASE_DIR, 'pic')
    CASE_DIR = os.path.join(BASE_DIR, 'Case')

    LEGACY_CASE_DIR = CASE_DIR
    # 软件信息
    SOFTWARE_NAME = '药研智析 DrugReasoner'
    SOFTWARE_SUBTITLE = '融合语言理解与符号关系推理的药物数据分析系统'
    PROJECT_FOLDER_NAME = 'DrugReasoner'
    
    # 模型配置
    DEFAULT_MODEL_VERSION = "LLMDTA_v1"
    
    # GUI配置
    MAIN_WINDOW_SIZE = '860x620'
    START_WINDOW_SIZE = '993x663'
    
    # MySQL路径（用于导入SQL）
    MYSQL_PATHS = [
        r"C:\Program Files\MySQL\MySQL Server 9.5\bin\mysql.exe",
        r"C:\Program Files (x86)\MySQL\MySQL Server 9.5\bin\mysql.exe",
        r"C:\Program Files\MySQL\MySQL Server 8.0\bin\mysql.exe",
        r"C:\Program Files (x86)\MySQL\MySQL Server 8.0\bin\mysql.exe",
    ]


# 创建必要的目录
for dir_path in [
    Config.DATA_DIR,
    Config.RAW_DATA_DIR,
    Config.PROCESSED_DATA_DIR,
    Config.CACHE_DATA_DIR,
    Config.CASES_DATA_DIR,
    Config.EXPORTS_DATA_DIR,
    Config.RAW_DTI_DIR,
    Config.RAW_DDI_DIR,
    Config.RAW_ENTITY_DIR,
    Config.RAW_KNOWLEDGE_DIR,
    Config.RAW_DAVIS_DIR,
    Config.RAW_KIBA_DIR,
    Config.RAW_BINDINGDB_DIR,
    Config.RAW_DDINTER_DIR,
    Config.RAW_DRUGBANK_DIR,
    Config.RAW_MECDDI_DIR,
    Config.RAW_OFFSIDES_DIR,
    Config.RAW_DRUG_ALIAS_DIR,
    Config.RAW_PROTEIN_ALIAS_DIR,
    Config.RAW_DGIDB_DIR,
    Config.RAW_CTD_DIR,
    Config.RAW_STITCH_DIR,
    Config.PROCESSED_QUERY_DIR,
    Config.PROCESSED_PREDICTION_DIR,
    Config.PROCESSED_EXPLANATION_DIR,
    Config.QUERY_KNOWN_DDI_DIR,
    Config.QUERY_KNOWN_DTI_DIR,
    Config.QUERY_ENTITY_MAP_DIR,
    Config.PREDICTION_DTI_DIR,
    Config.PREDICTION_DDI_DIR,
    Config.EXPLANATION_CORPORA_DIR,
    Config.EXPLANATION_TEMPLATES_DIR,
    Config.CACHE_EMBEDDINGS_DIR,
    Config.CACHE_FEATURES_DIR,
    Config.CACHE_RETRIEVAL_DIR,
    Config.CASE_GENERAL_DIR,
    Config.CASE_EGFR_DIR,
    Config.CASE_DEMO_DIR,
    Config.EXPORTS_CSV_DIR,
    Config.EXPORTS_JSON_DIR,
    Config.EXPORTS_LOGS_DIR,
    Config.REPORTS_DIR,
    Config.OUTPUT_DIR,
]:
    os.makedirs(dir_path, exist_ok=True)

