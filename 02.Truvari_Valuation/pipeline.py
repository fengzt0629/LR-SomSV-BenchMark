import os
import matplotlib.pyplot as plt
import pandas as pd
from multiprocessing import Pool

from config import base_path,truvari,bgzip,tabix,work_dir,ref_name,comp_path,hg38_polish_chain,head_txt,soft_ware
from function import VCF_to_bed


# IF ref is not hg38 ,need this function to Get location in hg38
comp_path_hg38 = VCF_to_bed(comp_path,hg38_polish_chain,head_txt)



output = os.path.join(work_dir,'truvari_dir',ref_name,'different_way_dir')
if not os.path.exists(output):
    os.makedirs(output)



os.system('{bgzip} -c {comp_path}  > {comp_path}.gz'.format(comp_path=comp_path, bgzip=bgzip))
os.system('{bgzip} -c {base_path}  > {base_path}.gz'.format(base_path=base_path, bgzip=bgzip))
os.system('{tabix} -p vcf {comp_path}.gz'.format(tabix=tabix, comp_path=comp_path))
os.system('{tabix} -p vcf {base_path}.gz'.format(tabix=tabix, base_path=base_path))
os.system('{truvari} bench -b {base_path}.gz -c {comp_path}.gz -o {output}  --pctseq 0' .format(truvari=truvari, base_path=base_path, comp_path=comp_path, output=output))
#  --pctseq 0 即不考虑序列信息




###  整理召回率和精度情况汇总
import json

recall_pre = os.path.join(work_dir,'truvari_dir',"recall_precision_2.new_golden_vcf.txt") #new_   old_
mapping = 'minimap2'
call = soft_ware
method = soft_ware
with open(recall_pre, 'w') as f:
    f.write('ref\tmapping\tcall\tmethod\tTP\tFP\tFN\tprecision\trecall\tf1\n')
    summary = os.path.join(output,"summary.json")
    with open(summary, 'r') as f1:
        data = json.load(f1)
        recall = data['recall']
        precision = data['precision']
        f1 = data['f1']
        # print(recall, precision, f1)
        # with open(recall_pre, 'a') as f2:
        f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ref_name, mapping, call, method,
                                                                  data['TP-base'], data['FP'], data['FN'],
                                                                  data['precision'], data['recall'],
                                                                  data['f1']))














