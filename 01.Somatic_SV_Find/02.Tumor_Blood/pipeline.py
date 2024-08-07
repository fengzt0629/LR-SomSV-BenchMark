
import os
import matplotlib.pyplot as plt
import pandas as pd
from multiprocessing import Pool

from config import ref_fasta,sampleid,somatic_dir,program_dir,blood_fastq,tumor_fastq,ms_dir,ref_name
from config import samtools,sniffles2,jasmine,SURVIVOR,tr_bed
from function import bam_tag,ms_bam



# map and add tag
###########################################################################################################################
################################################ using Merge_way add tag bam
############################################################################################################################

###########################################################################################################################
###############################################       call SV and tumor-blood
###########################################################################################################################

bam_tag_dir_tumor = os.path.join(somatic_dir, 'tumor', sampleid,
                                 '{sampleid}_tumor_minimap2.2.26_sorted_tag.bam'.format(sampleid=sampleid))
bam_tag_dir_blood = os.path.join(somatic_dir, 'blood', sampleid,
                                 '{sampleid}_blood_minimap2.2.26_sorted_tag.bam'.format(sampleid=sampleid))
sniffles2_sample_dir = os.path.join(ms_dir,ref_name,sampleid, 'sniffles2')
jasmines_sample_dir = os.path.join(ms_dir, ref_name, sampleid, 'jasmine')
SURVIVOR_sample_dir = os.path.join(ms_dir, ref_name, sampleid, 'SURVIVOR')

for sample_dir in [sniffles2_sample_dir,jasmines_sample_dir,SURVIVOR_sample_dir]:
    if not os.path.exists(sample_dir):
        os.makedirs(sample_dir)

script_ms = os.path.join(ms_dir,ref_name,sampleid, '%s_bam_%s_sniffles.sh' % (sampleid, 'sniffles2_jasmine_and_SURVIVOR'))
file_list = os.path.join(jasmines_sample_dir, 'file_list.txt')
out_file = os.path.join(jasmines_sample_dir, 'jasmine.sniffles.reslut.vcf')
out_file_survivor = os.path.join(SURVIVOR_sample_dir, 'SURVIVOR.sniffles.reslut.vcf')
sv_vcf2_T_N_tumor = os.path.join(sniffles2_sample_dir, '%s_tumor_minimap2_sniffles_v2.vcf' % sampleid)
sv_vcf2_T_N_blood = os.path.join(sniffles2_sample_dir, '%s_blood_minimap2_sniffles_v2.vcf' % sampleid)

with open(file_list, 'w') as out:
    out.write(sv_vcf2_T_N_tumor + '\n')
    out.write(sv_vcf2_T_N_blood + '\n')

with open(script_ms, 'w') as out:
    out.write("#! /bin/bash" + '\n')
    out.write('''echo "$(date) 1. Start to sniffles tumor: %s" ''' % sampleid + '\n')
    cmd_tumor = "{sniffles} -i {bam_sort} -v {vcf} --minsupport 3 --cluster-binsize 100 --minsvlen 50  -t 50 --mapq 40 --min-alignment-length 1000 --output-rnames --allow-overwrite --reference {ref}".format(
        sniffles=sniffles2, bam_sort=bam_tag_dir_tumor, vcf=sv_vcf2_T_N_tumor, tr=tr_bed,
        ref=ref_fasta)  # --minsvlen default 35bp --tandem-repeats {tr}
    out.write(cmd_tumor + '\n')  # --minsvlen default 35bp
    out.write('''echo "$(date) 1. Finish to sniffles tumor: %s" ''' % sampleid + '\n')
    out.write('''echo "$(date) 2. Start to sniffles blood: %s" ''' % sampleid + '\n')
    cmd_blood = "{sniffles} -i {bam_sort} -v {vcf} --minsupport 3 --cluster-binsize 100 --minsvlen 50  -t 50 --mapq 40 --min-alignment-length 1000 --output-rnames --allow-overwrite --reference {ref}".format(
        sniffles=sniffles2, bam_sort=bam_tag_dir_blood, vcf=sv_vcf2_T_N_blood, tr=tr_bed,
        ref=ref_fasta)  # --minsvlen default 35bp --tandem-repeats {tr}
    out.write(cmd_blood + '\n')  # --minsvlen default 35bp
    out.write('''echo "$(date) 2. Finish to sniffles blood: %s" ''' % sampleid + '\n')
    # out.write("source activate jasmine_env" + "\n") ##手动进入jasmine环境
    out.write('''echo "$(date) 3. Start to jasmine: %s" ''' % sampleid + '\n')
    cmd_jasmine = '{jasmine} threads=30 min_seq_id=0.9 k_jaccard=9 ignore_strand max_dist=500 file_list={sample_vcf} out_file={out_vcf}'.format(
        jasmine=jasmine, sample_vcf=file_list, out_vcf=out_file)
    out.write(cmd_jasmine + '\n')
    out.write('''echo "$(date) 3. Finish to jasmine: %s" ''' % sampleid + '\n')
    out.write('''echo "$(date) 4. Start to SURVIVOR: %s" ''' % sampleid + '\n')
    cmd_SURVIVOR = '{SURVIVOR} merge {sample_files} 1000 2 1 1 0 30 {sample_merged_vcf}'.format(SURVIVOR=SURVIVOR,
                                                                                                sample_files=file_list,
                                                                                                sample_merged_vcf=out_file_survivor)
    out.write(cmd_SURVIVOR + '\n')
    out.write('''echo "$(date) 4. Finish to SURVIVOR: %s" ''' % sampleid + '\n')



script_ms = os.path.join(jasmines_sample_dir, '%s_bam_%s_sniffles.sh' % (sampleid, 'sniffles2_jasmine_and_SURVIVOR'))
## Run sh
stdout = script_ms.replace(".sh", ".o")
stderr = script_ms.replace(".sh", ".e")
for std in [stdout, stderr]:
    if os.path.exists(std):
        os.system("rm %s" % std)
# pools.apply_async(ms_bam, args=(script_ms,stdout,stderr))
ms_bam(script_ms, stdout, stderr)















