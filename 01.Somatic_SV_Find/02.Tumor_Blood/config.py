import os

import matplotlib.pyplot as plt
import pandas as pd
from multiprocessing import Pool

samtools = "samtools" # Version 1.9
bedtools = "bedtools" # Version 2.30.0
sniffles2 = "sniffles" #Version 2.0.6
minimap2 = 'minimap2' #Version 2.17
seqkit = 'seqkit' #Version 2.4.0
tr_bed = "human_GRCh38_no_alt_analysis_set.trf.bed" #from sniffels

ref_fasta = 'hg38_fa'
sampleid = 'Sample'
ref_name = 'hg38_fa'
tumor_fastq = 'tumor_fastq'
blood_fastq = 'blood_fastq'


# work_dir = '/NAS/wg_liuxy/02.benchmark_project/04.minimap2_test_data/'

work_dir = 'work_dir'

data_dir = os.path.join(work_dir, '01.Golden.Standard')
program_dir = os.path.join(work_dir, '02.program')
ms_dir = os.path.join(data_dir, 'minimap2_sniffles')
coverage_dir = os.path.join(data_dir,'coverage')
somatic_dir = os.path.join(data_dir, "Somatic_minimap2_2.4"+ref_name)
