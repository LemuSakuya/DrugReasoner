import torch
from torch.utils.data import Dataset
import numpy as np
import pandas as pd

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)

_FALLBACK_WARNED = set()


def _select_cached_id(requested_id, cache_dict, cache_name):
    if not isinstance(cache_dict, dict):
        return requested_id
    vec_dict = cache_dict.get("vec_dict", {})
    if requested_id in vec_dict:
        return requested_id
    if not vec_dict:
        raise KeyError(f"{cache_name} 缓存为空，无法匹配 {requested_id}。请检查预训练特征文件是否存在且可用。")
    fallback_id = next(iter(vec_dict))
    warn_key = (cache_name, requested_id)
    if warn_key not in _FALLBACK_WARNED:
        _FALLBACK_WARNED.add(warn_key)
        print(f"[WARN] {cache_name} 缓存缺少 {requested_id}，回退使用 {fallback_id}。")
    return fallback_id


def matrix_pad(arr, max_len):
    dim = arr.shape[-1]
    length = arr.shape[0]
    if length < max_len:
        new_arr = np.zeros((max_len, dim))
        vec_mask = np.zeros((max_len))
        new_arr[:length] = arr
        vec_mask[:length] = 1
        return new_arr, vec_mask
    else:
        new_arr = arr[:max_len]
        vec_mask = np.ones((max_len))
        return new_arr, vec_mask


def my_collate_fn4pred(batch_data, device, hp, mol2vec_dict, protvec_dict, isEsm=False):
    batch_size = len(batch_data)
    drug_substruc_max = hp.substructure_max_len
    protein_max = hp.prot_max_len
    mol2vec_dim = hp.mol2vec_dim
    protvec_dim = hp.protvec_dim

    b_drug_vec = torch.zeros((batch_size, mol2vec_dim), dtype=torch.float32)
    b_prot_vec = torch.zeros((batch_size, protvec_dim), dtype=torch.float32)
    b_drug_mask = torch.zeros((batch_size, drug_substruc_max), dtype=torch.float32)
    b_prot_mask = torch.zeros((batch_size, protein_max), dtype=torch.float32)
    b_drug_mat = torch.zeros((batch_size, drug_substruc_max, mol2vec_dim), dtype=torch.float32)
    b_prot_mat = torch.zeros((batch_size, protein_max, protvec_dim), dtype=torch.float32)
    b_label = torch.zeros(batch_size)

    for i, pair in enumerate(batch_data):
        drug_id, prot_id = pair.iloc[0], pair.iloc[1]
        drug_id = str(drug_id)
        prot_id = str(prot_id)
        drug_id = _select_cached_id(drug_id, mol2vec_dict, "Mol2vec")
        prot_id = _select_cached_id(prot_id, protvec_dict, "ESM")
        drug_vec = mol2vec_dict["vec_dict"][drug_id]
        prot_vec = protvec_dict["vec_dict"][prot_id]
        drug_mat = mol2vec_dict["mat_dict"][drug_id]
        prot_mat = protvec_dict["mat_dict"][prot_id]
        drug_mat_pad, drug_mask = matrix_pad(drug_mat, drug_substruc_max)
        prot_mat_pad, prot_mask = matrix_pad(prot_mat, protein_max)

        b_drug_vec[i] = torch.from_numpy(drug_vec)
        b_prot_vec[i] = torch.from_numpy(prot_vec)
        b_drug_mat[i] = torch.from_numpy(drug_mat_pad)
        b_drug_mask[i] = torch.from_numpy(drug_mask)
        b_prot_mat[i] = torch.from_numpy(prot_mat_pad)
        b_prot_mask[i] = torch.from_numpy(prot_mask)

    b_drug_vec = b_drug_vec.to(device)
    b_prot_vec = b_prot_vec.to(device)
    b_drug_mat = b_drug_mat.to(device)
    b_drug_mask = b_drug_mask.to(device)
    b_prot_mat = b_prot_mat.to(device)
    b_prot_mask = b_prot_mask.to(device)
    return b_drug_vec, b_prot_vec, b_drug_mat, b_drug_mask, b_prot_mat, b_prot_mask


class CustomDataSet(Dataset):
    def __init__(self, dataset, hp):
        self.hp = hp
        self.dataset = dataset

    def __getitem__(self, index):
        return self.dataset.iloc[index, :]

    def __len__(self):
        return len(self.dataset)
