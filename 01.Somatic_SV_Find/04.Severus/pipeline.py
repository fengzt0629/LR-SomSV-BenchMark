

import os
import matplotlib.pyplot as plt
import pandas as pd
from multiprocessing import Pool

from config import ref_fasta,sampleid,somatic_dir,program_dir,blood_fastq,tumor_fastq,ms_dir,ref_name
from config import samtools,sniffles2,tr_bed,nanomonsv,severus
from function import bam_tag,ms_bam



###########################################################################################################################
################################################ using  add tag bam in Merge_way
############################################################################################################################


sniffles2_sample_dir = os.path.join(ms_dir,ref_name,sampleid, 'severus')

if os.path.exists(sniffles2_sample_dir):
    os.system("rm -rf %s" % sniffles2_sample_dir)
if not os.path.exists(sniffles2_sample_dir):
    os.makedirs(sniffles2_sample_dir)
sv_vcf2_T_N = os.path.join(sniffles2_sample_dir, '%s_T_N_merge_minimap2_severus_v2.vcf' % sampleid)


script_ms = os.path.join(ms_dir,ref_name,sampleid, '%s_bam_severus.2nd.sh' % sampleid)

bam_tag_dir_tumor = os.path.join(somatic_dir, 'tumor', sampleid,
                                 '{sampleid}_tumor_minimap2.2.26_sorted_tag.bam'.format(sampleid=sampleid))
bam_tag_dir_blood = os.path.join(somatic_dir, 'blood', sampleid,
                                 '{sampleid}_blood_minimap2.2.26_sorted_tag.bam'.format(sampleid=sampleid))

cmd_1 = "source activate severus_env"
cmd = "{severus} --target-bam {tumorbam} --control-bam {bloodbam} --out-dir {outdir} -t 30 --output-read-ids ".format(
    severus=severus, tumorbam=bam_tag_dir_tumor, bloodbam=bam_tag_dir_blood, sampleid=sampleid, outdir=sniffles2_sample_dir,
    ref=ref_fasta)  # defult mapq=10 svlen=50

# Write Severus sh

with open(script_ms, 'w') as out:
    out.write("#! /bin/bash" + '\n')
    out.write(cmd_1 + "\n")
    out.write('''echo "$(date) 2. Start to severus :%s" ''' % sampleid + '\n')
    out.write(cmd + '\n')  # --minsvlen default 35bp
    out.write('''echo "$(date) 2. Finish to severus :%s:" ''' % sampleid + '\n')


## Run sh
script_ms = os.path.join(ms_dir,ref_name,sampleid, '%s_bam_severus.2nd.sh' % sampleid)
stdout = script_ms.replace(".sh", ".o")
stderr = script_ms.replace(".sh", ".e")
for std in [stdout, stderr]:
    if os.path.exists(std):
        os.system("rm %s" % std)
# pools.apply_async(ms_bam, args=(script_ms,stdout,stderr))
ms_bam(script_ms, stdout, stderr)













