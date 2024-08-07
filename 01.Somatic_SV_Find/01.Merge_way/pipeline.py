
import os
import matplotlib.pyplot as plt
import pandas as pd
from multiprocessing import Pool

from config import ref_fasta,sampleid,somatic_dir,program_dir,blood_fastq,tumor_fastq,ms_dir,ref_name
from config import samtools,sniffles2,cuteSV,tr_bed
from function import bam_tag,ms_bam,sniffles2_vcftobed,cuteSV_vcftobed



# map and add tag

for tag in ['tumor','blood']:
    if tag == 'tumor':
        fq = tumor_fastq
    else:
        fq = blood_fastq
    fq_dir = 'NONE'
    script_tag = os.path.join(program_dir, "ONT_bam_add_tag_minimap2_%s.sh" % tag) # sh in our 01.Merge_way
    somatic_sample_dir = os.path.join(somatic_dir,tag,sampleid)
    if not os.path.exists(somatic_sample_dir):
        os.makedirs(somatic_sample_dir)
    new_ID = '%s_%s' % (sampleid, tag)
    stdout = os.path.join(somatic_sample_dir, '%s_add_tag.o' % new_ID)
    stderr = os.path.join(somatic_sample_dir, '%s_add_tag.e' % new_ID)
    bam_tag(script_tag,fq,new_ID,somatic_sample_dir,fq_dir,ref_fasta,stdout,stderr)



# merge and call SV

sniffles2_sample_dir = os.path.join(ms_dir,  ref_name,sampleid, 'sniffles2')
if not os.path.exists(sniffles2_sample_dir):
    os.makedirs(sniffles2_sample_dir)

sv_vcf2_T_N = os.path.join(sniffles2_sample_dir, '%s_T_N_merge_minimap2_sniffles_v2.vcf' % sampleid)
bam_tag_dir_tumor = os.path.join(somatic_dir,'tumor',sampleid,'{sampleid}_tumor_minimap2.2.26_sorted_tag.bam'.format(sampleid=sampleid))
bam_tag_dir_blood = os.path.join(somatic_dir,'blood',sampleid,'{sampleid}_blood_minimap2.2.26_sorted_tag.bam'.format(sampleid=sampleid))
bam_merge_file= os.path.join(somatic_dir, 'Merge_bam', ref_name,sampleid)

if not os.path.exists(bam_merge_file):
    os.makedirs(bam_merge_file)

bam_merge_T_N = os.path.join(bam_merge_file, '%s_minimap2_sorted_T_N_merge.bam' % sampleid)
script_ms = os.path.join(ms_dir,ref_name,sampleid, '%s_bam_merge_sniffles.sh' % sampleid)

with open(script_ms, 'w') as out:
    out.write("#! /bin/bash" + '\n')
    out.write('''echo "$(date) 1. Start to merge bam: %s" ''' % sampleid + '\n')
    out.write("{samtools} merge -@ 40 -h {bam_N} {out_bam} {bam_N} {bam_T}".format(
        samtools=samtools, bam_N=bam_tag_dir_blood, bam_T=bam_tag_dir_tumor, out_bam=bam_merge_T_N) + '\n')
    out.write('''/NAS/wg_fzt/software/samtools-1.9/samtools index -@ 25 %s \n''' % bam_merge_T_N)
    out.write('''echo "$(date) 1. Finish to merge bam: %s" ''' % sampleid + '\n')
    out.write('''echo "$(date) 2. Start to sniffles2: %s" ''' % sampleid + '\n')
    out.write(
        "{sniffles} -i {bam_sort} -v {vcf} -t 50 --minsupport 1 --mapq 10 --min-alignment-length 1000 "
        "--output-rnames --allow-overwrite --long-ins-length 100000 --reference {ref}".format(
            sniffles=sniffles2, bam_sort=bam_merge_T_N, vcf=sv_vcf2_T_N, ref= ref_fasta) + '\n') # --minsvlen default 35bp
    out.write('''echo "$(date) 2. Finish to sniffles2: %s" ''' % sampleid + '\n')





## Run sh
stdout = script_ms.replace(".sh", ".o")
stderr = script_ms.replace(".sh", ".e")
for std in [stdout, stderr]:
    if os.path.exists(std):
        os.system("rm %s" % std)
# pools.apply_async(ms_bam, args=(script_ms,stdout,stderr))
ms_bam(script_ms,stdout,stderr)



#Filter somatic SV
sv_bed2_T_N = os.path.join(sniffles2_sample_dir, '%s_T_N_merge_minimap2_sniffles_v2.bed' % sampleid)
sniffles2_vcftobed(sv_vcf2_T_N, sv_bed2_T_N, 3)






###################################################################################################################################
##############          Cute SV
################################################################################################################################

cuteSV_sample_dir = os.path.join(ms_dir, ref_name,sampleid, 'cuteSV')
if not os.path.exists(cuteSV_sample_dir):
    os.makedirs(cuteSV_sample_dir)
sv_vcf2_T_N = os.path.join(cuteSV_sample_dir, '%s_T_N_merge_minimap2_cuteSV_v2.vcf' % sampleid)
bam_merge_file= os.path.join(somatic_dir, 'Merge_bam', sampleid)
bam_merge_T_N = os.path.join(bam_merge_file, '%s_minimap2_sorted_T_N_merge.bam' % sampleid)
cuteSV_script_ms = os.path.join(ms_dir,ref_name,sampleid,'%s_bam_cuteSV.2nd.sh' % sampleid)
with open(cuteSV_script_ms, 'w') as out:
    out.write("#! /bin/bash" + '\n')
    out.write('''echo "$(date) 1. Start to CuteSV: %s" ''' % sampleid + '\n')
    out.write(
        "{cuteSV} -t 20 --min_support 1 --min_mapq 10 --min_read_len 1000 --report_readid {bam_sort} {ref} {vcf} {work_dir}".format(
            cuteSV=cuteSV, bam_sort=bam_merge_T_N, vcf=sv_vcf2_T_N, tr=tr_bed, ref=ref_fasta,work_dir = cuteSV_sample_dir) + '\n') # --minsvlen default 35bp
    out.write('''echo "$(date) 1. Finish to CuteSV: %s" ''' % sampleid + '\n')




## Run sh
stdout = cuteSV_script_ms.replace(".sh", ".o")
stderr = cuteSV_script_ms.replace(".sh", ".e")
for std in [stdout, stderr]:
    if os.path.exists(std):
        os.system("rm %s" % std)
# pools.apply_async(ms_bam, args=(script_ms,stdout,stderr))
ms_bam(cuteSV_script_ms,stdout,stderr)



#Filter somatic SV
sv_bed2_T_N = os.path.join(sniffles2_sample_dir, '%s_T_N_merge_minimap2_cuteSV_v2.bed' % sampleid)
cuteSV_vcftobed(sv_vcf2_T_N, sv_bed2_T_N, 3)








