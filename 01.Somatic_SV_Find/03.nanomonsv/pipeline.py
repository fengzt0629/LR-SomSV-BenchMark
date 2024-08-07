


import os
import matplotlib.pyplot as plt
import pandas as pd
from multiprocessing import Pool

from config import ref_fasta,sampleid,somatic_dir,program_dir,blood_fastq,tumor_fastq,ms_dir,ref_name
from config import samtools,sniffles2,tr_bed,nanomonsv
from function import bam_tag,ms_bam



###########################################################################################################################
################################################ using  add tag bam in Merge_way
############################################################################################################################

sniffles2_sample_dir = os.path.join(ms_dir,ref_name,sampleid, 'nanomonsv')

if os.path.exists(sniffles2_sample_dir):
    os.system("rm -rf %s" % sniffles2_sample_dir)
if not os.path.exists(sniffles2_sample_dir):
    os.makedirs(sniffles2_sample_dir)

script_ms = os.path.join(ms_dir,ref_name,sampleid, '%s_bam_nanomonsv.2nd.sh' % sampleid)

bam_tag_dir_tumor = os.path.join(somatic_dir, 'tumor', sampleid,
                                 '{sampleid}_tumor_minimap2.2.26_sorted_tag.bam'.format(sampleid=sampleid))
bam_tag_dir_blood = os.path.join(somatic_dir, 'blood', sampleid,
                                 '{sampleid}_blood_minimap2.2.26_sorted_tag.bam'.format(sampleid=sampleid))

output_prefix_tumor = os.path.join(sniffles2_sample_dir,'tumor')
if not os.path.exists(output_prefix_tumor):
    os.makedirs(output_prefix_tumor)
output_prefix_tumor_1 = output_prefix_tumor + '/' + sampleid+'_tumor'

output_prefix_blood = os.path.join(sniffles2_sample_dir,'blood')
if not os.path.exists(output_prefix_blood):
    os.makedirs(output_prefix_blood)
output_prefix_blood_1 = output_prefix_blood + '/' + sampleid+'_blood'


# Write nanomosv sh

with open(script_ms, 'w') as out:
    out.write("#! /bin/bash" + '\n')
    out.write("source activate  nanomonsv" + "\n")
    out.write('''echo "$(date) 1. Start to tumor nanomonsv parse: %s" ''' % sampleid + '\n')
    out.write("nanomonsv parse {bam_file} {output_prefix_tumor_1}  ".format(
        nanomonsv=nanomonsv, bam_file = bam_tag_dir_tumor,output_prefix_tumor_1 = output_prefix_tumor_1) + '\n')
    out.write('''echo "$(date) 1. Finish to tumor nanomonsv parse: %s" ''' % sampleid + '\n')
    out.write('''echo "$(date) 2. Start to blood nanomonsv parse: %s" ''' % sampleid + '\n')
    out.write("nanomonsv parse {bam_file} {output_prefix_blood_1} ".format(
        nanomonsv=nanomonsv, bam_file = bam_tag_dir_blood,output_prefix_blood_1 = output_prefix_blood_1) + '\n')
    out.write('''echo "$(date) 2. Finish to blood nanomonsv parse: %s" ''' % sampleid + '\n')
    out.write('''echo "$(date) 3. Start to nanomonsv get : %s" ''' % sampleid + '\n')
    out.write(
        "nanomonsv get {output_prefix_tumor_1} {tumor_bam_file} {ref}  --control_prefix {output_prefix_blood_1} --control_bam {blood_bam_file} --single_bnd --use_racon --threads 200 --control_panel_prefix /NAS/wg_liuxy/02.benchmark_project/04.minimap2_test_data/01.Golden.Standard/04.minimap2_nanomonsv_2/COLO829/02.control_panel/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control/hprc_year1_data_freeze_nanopore_minimap2_2_24_merge_control  ".format(
            nanomonsv=nanomonsv, output_prefix_tumor_1 = output_prefix_tumor_1,tumor_bam_file = bam_tag_dir_tumor ,ref = ref_fasta, output_prefix_blood_1 = output_prefix_blood_1, blood_bam_file = bam_tag_dir_blood ) + '\n') # --minsvlen default 35bp --min_tumor_variant_read_num 3 --max_control_variant_read_num 0 --min_indel_size 50  --threads 50 --median_mapQ_thres 10
    out.write('''echo "$(date) 3. Finish to nanomonsv get: %s" ''' % sampleid + '\n')





## Run sh
script_ms = os.path.join(ms_dir,ref_name,sampleid, '%s_bam_nanomonsv.2nd.sh' % sampleid)
stdout = script_ms.replace(".sh", ".o")
stderr = script_ms.replace(".sh", ".e")
for std in [stdout, stderr]:
    if os.path.exists(std):
        os.system("rm %s" % std)
# pools.apply_async(ms_bam, args=(script_ms,stdout,stderr))
ms_bam(script_ms, stdout, stderr)



