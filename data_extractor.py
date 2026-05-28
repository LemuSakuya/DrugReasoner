import csv
import pandas as pd


def find_drug(names, input_file):
    output_file = 'Case/drug.tsv'
    target_name = names
    found = False
    df = pd.read_csv(input_file, sep="\t")
    found = df.iloc[:, 0].str.contains(target_name, na=False).any()
    if found:
        matching_rows = df[df.iloc[:, 0].str.contains(target_name, na=False)]
        matching_rows.to_csv(output_file, sep='\t', index=False, encoding='utf-8')
        print(f"已将 {input_file} 中名字为 {target_name} 的行保存到 {output_file}")
    else:
        print(f"在 {input_file} 中未找到名字为 {target_name} 的行")


def find_protein(names, input_file):
    output_file = 'Case/prot.tsv'
    target_name = names
    df = pd.read_csv(input_file, sep=",")
    found = (df.iloc[:, 0] == target_name).any()
    if found:
        matching_rows = df[df.iloc[:, 0] == target_name]
        matching_rows.to_csv(output_file, sep='\t', index=False, encoding='utf-8')
        print(f"已将 {input_file} 中名字为 {target_name} 的行保存到 {output_file}")
    else:
        print(f"在 {input_file} 中未找到名字为 {target_name} 的行")
