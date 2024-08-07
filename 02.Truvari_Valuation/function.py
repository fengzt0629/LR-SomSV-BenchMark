

import os

import matplotlib.pyplot as plt
import pandas as pd
from multiprocessing import Pool
import re
from config import liftover


def transform(jasmine_bed,out_jasmine_bed,unmmap,hg38_polish_chain):
    os.system('{liftover} -bedPlus=3 {inputbed} {chainfile} {out_bed} {unmapped_bed}'.format(liftover=liftover, inputbed=jasmine_bed, chainfile=hg38_polish_chain, out_bed=out_jasmine_bed, unmapped_bed=unmmap))



def VCF_to_bed(vcf_DIR,hg38_polish_chain,head_txt):
    vcf = vcf_DIR
    bed = vcf.replace('.vcf', '.bed')
    f = vcf
    print(f)
    with open(f, 'r') as fin:
        records = [x.strip().split("\t") for x in fin.readlines() if not re.search('##', x)]
    vcfDf = pd.DataFrame.from_records(records[1:])
    print(records[0])
    vcfDf.columns = records[0]
    vcfDf['SVType'] = vcfDf['INFO'].apply(lambda x: x.split("SVTYPE=")[-1].split(";")[0])
    vcfDf['POS2'] = vcfDf.apply(lambda x: int(x['POS']) + 1 if (x['SVType'] in ['INS', 'BND', 'INV']) else
    x['INFO'].split("END=")[-1].split(";")[0], axis=1)
    # 根据records排列的顺序，将vcfDf中的列顺序调整为records中的顺序
    records_2 = ['#CHROM', 'POS', 'POS2'] + records[0][2:]
    vcfDf = vcfDf[records_2]
    vcfDf.sort_values(by=['#CHROM', 'POS', 'POS2'], inplace=True)
    vcfDf.to_csv(bed, sep='\t', index=False, header=False)
    jasmine_bed = bed
    out_jasmine_bed = bed.replace('.bed', '_hg38.bed')
    unmmap = bed.replace('.bed', '_hg38_unmmap.bed')
    transform(jasmine_bed, out_jasmine_bed, unmmap, hg38_polish_chain)
    out_jasmine_bed_df = pd.read_csv(out_jasmine_bed, sep='\t', header=None)
    out_jasmine_bed_df[2] = out_jasmine_bed_df[2].astype(int)
    out_jasmine_bed_df[8] = out_jasmine_bed_df.apply(lambda row: re.sub(r'END=\d+', 'END=' + str(row[2]), row[8]),axis=1)
    out_jasmine_bed_df.drop([2], axis=1, inplace=True)
    out_jasmine_bed_df = out_jasmine_bed_df.sort_values(by=[0, 1])
    out_jasmine_bed_df.to_csv(out_jasmine_bed, sep='\t', index=False, header=False)
    head_vcf = bed.replace('.bed', '_hg38.head.vcf')
    os.system('cat {head_txt} {out_jasmine_bed} > {head_vcf}'.format(head_txt=head_txt, out_jasmine_bed=out_jasmine_bed,
                                                                     head_vcf=head_vcf))
    comp_path_hg38 = head_vcf
    return comp_path_hg38






