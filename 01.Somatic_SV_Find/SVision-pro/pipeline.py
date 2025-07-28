import os
import matplotlib.pyplot as plt
import pandas as pd
from multiprocessing import Pool

from config import ref_fasta,sampleid,somatic_dir,program_dir,blood_fastq,\
    tumor_fastq,ms_dir,ref_name,model_path,tag
from config import samtools,sniffles2,tr_bed,nanomonsv,severus
from function import bam_tag,ms_bam



svision_sample_dir = os.path.join(ms_dir,ref_name,sampleid, 'svisionpro')


script_ms = os.path.join(ms_dir,ref_name,sampleid, '%s_bam_severus.2nd.sh' % sampleid)

bam_tag_dir_tumor = os.path.join(somatic_dir, 'tumor', sampleid,
                                 '{sampleid}_tumor_minimap2.2.26_sorted_tag.bam'.format(sampleid=sampleid))
bam_tag_dir_blood = os.path.join(somatic_dir, 'blood', sampleid,
                                 '{sampleid}_blood_minimap2.2.26_sorted_tag.bam'.format(sampleid=sampleid))

cmd_1 = "source activate svision-pro-env" # version=2.3

cmd = "SVision-pro --target_path {tumor}  --base_path {normal} " \
                  "--genome_path {Ref} --preset error-prone --min_supp 3 --model_path {model} " \
                  "--max_sv_size 1000000 --out_path {out_path} --sample_name {sampleid} --process_num 60 " \
                  "--detect_mode somatic".format(tumor=bam_tag_dir_tumor,normal=bam_tag_dir_blood,Ref=ref_fasta,model=model_path,
                        out_path=svision_sample_dir,sampleid=tag)

# Write SVision-pro sh

with open(script_ms, 'w') as out:
    out.write("#! /bin/bash" + '\n')
    out.write(cmd_1 + "\n")
    out.write('''echo "$(date) 2. Start to SVision-pro :%s" ''' % sampleid + '\n')
    out.write(cmd + '\n')  # --minsvlen default 35bp
    out.write('''echo "$(date) 2. Finish to SVision-pro :%s:" ''' % sampleid + '\n')


## Run sh
script_ms = os.path.join(ms_dir,ref_name,sampleid, '%s_bam_svisionpro.2nd.sh' % sampleid)
stdout = script_ms.replace(".sh", ".o")
stderr = script_ms.replace(".sh", ".e")
for std in [stdout, stderr]:
    if os.path.exists(std):
        os.system("rm %s" % std)
# pools.apply_async(ms_bam, args=(script_ms,stdout,stderr))
ms_bam(script_ms, stdout, stderr)




