
import os
import re
import json
import pandas as pd
from multiprocessing import Pool
from config import work_dir,chr2_fa,samtools,sniffles2,tr_bed,bgzip,tabix,truvari,work1_dir,sv_length_2,sortBed
from function import simulate_function,tumor_add_tag_function,ms_bam,bed_to_vcf



##  The general idea is to first simulate the geneline SV of blood and then simulate the SV of tumor, DEL INV DUP INS in different lengths (100/1000/2000/5000/10,000/50000/50000/10,000/10,000,000/100,000,000/) inside the same fq (excluding duplicate regions).
## blood file bulit

# if you do not want to do this Step,We give exclude_Repeat_masked_chr2.bed and exclude_Repeat_masked_1.bed in our 03.SV_lenth_emulate
Repeat_masked_bed = "mainChr_hg38_UCSC.RepeatMasked.Unselected.bed" # for UCSC
exclude_Repeat_masked = os.path.join(work_dir, "exclude_Repeat_masked_chr2.bed")
exclude_Repeat_masked_1 = os.path.join(work_dir, "exclude_Repeat_masked_1.bed")
cmd_1 = "grep  -w -E 'chr2'  %s  > %s" % (Repeat_masked_bed, exclude_Repeat_masked)
cmd_2 = "grep  Satelite|Retroposon|RC  %s  > %s" % (exclude_Repeat_masked, exclude_Repeat_masked_1)
os.system(cmd_1)
os.system(cmd_2)
print(cmd_2)



#filna filter 为 exclude_chr2.bed

# randomregion.r  chr2.dim.tsv  exclude_chr2.bed and geneline_random.bed We also give in our 03.SV_lenth_emulate

os.system(
    "Rscript randomregion.r -d chr2.dim.tsv -n 1000 -l 1000 -x exclude_chr2.bed -v 'deletion,inversion,insertion,tandem duplication,reciprocal translocation' -r '20:20:20:20:20' | {sortBed} > geneline_random.bed".format(sortBed=sortBed))  #Simulate geneline_random sv


TB = 'blood'
work_dir = os.path.join(work1_dir)
if not os.path.exists(work1_dir):
    os.makedirs(work1_dir)
HACk_bed = os.path.join(work1_dir, "geneline_random.bed")
hack_out = os.path.join(work1_dir, "%s_hack.out" % (TB))
chrom_dim_tsv = os.path.join(work1_dir, "%s_haplochroms.dim.tsv" % (TB))
maxdims_tsv = os.path.join(work1_dir, "%s_maxdims.tsv" % (TB))
shorts_laser_simple = os.path.join(work1_dir, "%s_shorts.laser.simple.bed" % (TB))
laser_1_out = os.path.join(work1_dir, "%s_laser.1.out" % (TB))
os.system("VISOR HACk -g %s -b %s -o %s" % (chr2_fa, HACk_bed, hack_out))
os.system("cut -f1,2 %s/*.fai %s.fai > %s" % (hack_out, chr2_fa, chrom_dim_tsv))
os.system("cat %s | sort  | awk '$2 > maxvals[$1] {lines[$1]=$0;"
          " maxvals[$1]=$2} END { for (tag in lines) print lines[tag] }' > %s" % (chrom_dim_tsv, maxdims_tsv))
os.system('''awk 'OFS=FS="\t"''{print $1, "1", $2, "80.0", "80.0"}' %s > %s ''' % (
    maxdims_tsv, shorts_laser_simple)) ###  80 60 50 "80.0", "100.0"   "80.0", "80.0"  "100.0", "100.0"
os.system("VISOR LASeR -g %s -s %s -b %s -o %s --threads 100 --coverage 30  --fastq --tag" % (
    chr2_fa, hack_out, shorts_laser_simple, laser_1_out))  # 80 60 50



pools = Pool(16)

for SV_length in sv_length_2:  #
    pools.apply_async(simulate_function, (SV_length,))
    # simulate_function(SV_length)

pools.close()
pools.join()




##  map fastq  minimap2 + sniffles

# blood_bam_dir = 'blood_laser.1.out'
# agg_tag_blood.sh is 03.SV_length_emulate
cmd_2 = '/bin/bash  agg_tag_blood.sh  {work_dir}/blood_laser.1.out  1 >  {work_dir}/blood_laser.1.out/blood_add_o.0 2 >  {work_dir}/blood_laser.1.out/blood_add_e.0'.format(work_dir=work1_dir)

os.system(cmd_2)

pools = Pool(16)

for SV_length in sv_length_2:
    # tumor_add_tag_function(SV_length)
    pools.apply_async(tumor_add_tag_function, (SV_length,))

pools.close()
pools.join()




for SV_length in sv_length_2:
    TB = 'SV_length_%s' % SV_length
    work_dir = os.path.join(work1_dir, 'tumor', 'SV_length_%s' % SV_length)
    laser_1_out = os.path.join(work_dir, "%s_laser.1.out" % (TB))
    bam_file = os.path.join(laser_1_out, "sim.srt.bam")
    #加tag
    bam_sort_tag_tumor = os.path.join(laser_1_out, "sim.srt.sort.tag.bam")
    bam_sort_tag_blood = os.path.join(
        "{work1_dir}/blood_laser.1.out/sim.srt.sort.tag.bam".format(work1_dir=work1_dir))
    sh_merge_sniffles = os.path.join(work_dir, "merge_sniffles.sh")
    sv_vcf2_repeat = os.path.join(work_dir, "merge_minimap2_sniffles_visor.vcf")
    bam_merge_repeat = os.path.join(work_dir, '{SV_length}_sim.srt.sort.tag.merge.bam'.format(SV_length=SV_length))
    if os.path.exists(bam_merge_repeat):
        os.system('rm %s'%bam_merge_repeat)
    # sniffles v2.0
    with open(sh_merge_sniffles, 'w') as out:
        out.write("#! /bin/bash" + '\n')
        out.write('''echo "$(date) 1. Start to merge bam" ''' + '\n')
        # if sampleid != "PAAD07":
        out.write("{samtools} merge -@ 40 -h {bam_N} {out_bam} {bam_N} {bam_T}".format(
            samtools=samtools, bam_N=bam_sort_tag_blood, bam_T=bam_sort_tag_tumor, out_bam=bam_merge_repeat) + '\n')
        out.write('''%s index -@ 25 %s \n''' % (samtools,bam_merge_repeat))
        out.write('''echo "$(date) 1. Finish to merge bam" ''' + '\n')
        out.write('''echo "$(date) 2. Start to sniffles2" ''' + '\n')
        # if sampleid != "PAAD07":
        out.write(
            "{sniffles} -i {bam_sort} -v {vcf} --tandem-repeats {tr} -t 50 --minsupport 1 --mapq 10 --min-alignment-length 1000 "
            "--output-rnames --allow-overwrite --long-ins-length 100000 --reference {ref}".format(
                sniffles=sniffles2, bam_sort=bam_merge_repeat, vcf=sv_vcf2_repeat, tr=tr_bed,
                ref=chr2_fa) + '\n')  # --minsvlen default 35bp

##  运行sh脚本

pools = Pool(16)

for SV_length in sv_length_2:
    TB = 'SV_length_%s' % SV_length
    work_dir = os.path.join(work1_dir, 'tumor', 'SV_length_%s' % SV_length)
    laser_1_out = os.path.join(work_dir, "%s_laser.1.out" % (TB))
    sh_merge_sniffles = os.path.join(work_dir, "merge_sniffles.sh")
    stdout = sh_merge_sniffles.replace(".sh", ".o")
    stderr = sh_merge_sniffles.replace(".sh", ".e")
    pools.apply_async(ms_bam, args=(sh_merge_sniffles, stdout, stderr))

pools.close()
pools.join()
del pools

##  Filter somatic SV
for SV_length in sv_length_2:
    TB = 'SV_length_%s' % SV_length
    work_dir = os.path.join(work1_dir, 'tumor', 'SV_length_%s' % SV_length)
    laser_1_out = os.path.join(work_dir, "%s_laser.1.out" % (TB))
    sv_vcf2_repeat = os.path.join(work_dir, "merge_minimap2_sniffles_visor.vcf")
    f = sv_vcf2_repeat
    print(f)
    with open(f, 'r') as fin:
        records = [x.strip().split("\t") for x in fin.readlines() if not re.search('##', x)]
    vcfDf = pd.DataFrame.from_records(records[1:])
    print(records[0])
    vcfDf.columns = records[0]
    vcfDf['RNAMES'] = vcfDf['INFO'].apply(lambda x: x.split("RNAMES=")[-1].split(";")[0])
    vcfDf['TR'] = vcfDf.apply(lambda x: [r for r in x['RNAMES'].split(",") if re.search('tumor', r)], axis=1)
    vcfDf['NR'] = vcfDf.apply(lambda x: [r for r in x['RNAMES'].split(",") if re.search('blood', r)], axis=1)
    vcfDf['Support'] = vcfDf.apply(lambda x: x['TR'] + x['NR'], axis=1)
    vcfDf['TS'] = vcfDf['TR'].apply(lambda x: len(x))
    vcfDf['NS'] = vcfDf['NR'].apply(lambda x: len(x))
    #Adjust the order of the columns in vcfDf to the order in the records according to the order in which the records are arranged
    vcfDf = vcfDf[(vcfDf['TS'] >= 3) & (vcfDf['NS'] == 0)]
    # records_2 = ['#CHROM','POS','POS2']+ records[0][2:]
    vcfDf['POS'] = vcfDf['POS'].astype(int)
    vcfDf = vcfDf.sort_values(by=['#CHROM', 'POS'])
    vcfDf_filter_TS_NS = vcfDf[records[0]]
    vcfDf_filter_TS_NS.to_csv(f.replace(".vcf", ".filter_TS_NS.no.head.vcf"), sep="\t", index=False)
    vcfDf_filter_TS_NS_no_head = f.replace(".vcf", ".filter_TS_NS.no.head.vcf")
    vcfDf_filter_TS_NS_head = f.replace(".vcf", ".filter_TS_NS.head.vcf")
    os.system(
        "cat vcf.head.txt %s > %s" % (
        vcfDf_filter_TS_NS_no_head, vcfDf_filter_TS_NS_head))


## The gold standard for making every length
for SV_length in sv_length_2:
    work_dir = os.path.join(work1_dir, 'tumor', 'SV_length_%s' % SV_length)
    tumor_random_bed = os.path.join(work_dir, "tumor_random.bed")
    df = pd.read_csv(tumor_random_bed, sep="\t", header=None)
    bed_to_vcf(tumor_random_bed, SV_length)
    output_file = tumor_random_bed.replace(".bed", "tumor.golde.no.head.vcf")
    output_file_head = tumor_random_bed.replace(".bed", "tumor.golde.head.vcf")
    os.system(
        "cat vcf.head.txt %s > %s" % (
        output_file, output_file_head))



for SV_length in sv_length_2:
    work_dir = os.path.join(work1_dir, 'tumor', 'SV_length_%s' % SV_length)
    tumor_random_bed = os.path.join(work_dir, "tumor_random.bed")
    base_path = tumor_random_bed.replace(".bed", "tumor.golde.head.vcf")
    sv_vcf2_repeat = os.path.join(work_dir, "merge_minimap2_sniffles_visor.vcf")
    comp_path = sv_vcf2_repeat.replace(".vcf", ".filter_TS_NS.head.vcf")
    output = os.path.join(work_dir, "truvari_result")
    if os.path.exists(output):
        os.system("rm -rf %s" % output)
    os.system('{bgzip} -c {comp_path}  > {comp_path}.gz'.format(comp_path=comp_path, bgzip=bgzip))
    os.system('{tabix} -p vcf  {comp_path}.gz'.format(tabix=tabix, comp_path=comp_path))
    os.system('{bgzip} -c {base_path}  > {base_path}.gz'.format(base_path=base_path, bgzip=bgzip))
    os.system('{tabix} -p vcf {base_path}.gz'.format(tabix=tabix, base_path=base_path))
    os.system('{truvari} bench -b {base_path}.gz -c {comp_path}.gz -o {output}  --pctseq 0'.format(truvari=truvari,
                                                                                                   base_path=base_path,
                                                                                                   comp_path=comp_path,
                                                                                                   output=output))
    #  --pctseq 0 i.e., without considering the sequence information



for SV_length in sv_length_2:
    work_dir = os.path.join(work1_dir, 'tumor', 'SV_length_%s' % SV_length)
    tumor_random_bed = os.path.join(work_dir, "tumor_random.bed")
    sv_vcf2_repeat = os.path.join(work_dir, "merge_minimap2_sniffles_visor.vcf")
    for sv_type in ['DEL', 'INS', 'DUP', 'INV']:
        no_head_base_path = tumor_random_bed.replace(".bed", "tumor.golde.no.head.vcf")
        no_head_comp_path = sv_vcf2_repeat.replace(".vcf", ".filter_TS_NS.no.head.vcf")
        base_path_sv_type = tumor_random_bed.replace(".bed", "tumor.golde.no.head.%s.vcf" % sv_type)
        comp_path_sv_type = sv_vcf2_repeat.replace(".vcf", ".filter_TS_NS.no.head.%s.vcf" % sv_type)
        no_head_base_path_df = pd.read_csv(no_head_base_path, sep="\t", header=0)
        no_head_base_path_df['Svtype'] = no_head_base_path_df['INFO'].apply(
            lambda x: x.split("SVTYPE=")[-1].split(";")[0])
        base_path_sv_type_df = no_head_base_path_df[no_head_base_path_df['Svtype'] == sv_type]
        base_path_sv_type_df = base_path_sv_type_df.drop(['Svtype'], axis=1)
        base_path_sv_type_df.to_csv(base_path_sv_type, sep="\t", header=True, index=False)
        no_head_comp_path_df = pd.read_csv(no_head_comp_path, sep="\t", header=0)
        no_head_comp_path_df['Svtype'] = no_head_comp_path_df['INFO'].apply(
            lambda x: x.split("SVTYPE=")[-1].split(";")[0])
        comp_path_sv_type_df = no_head_comp_path_df[no_head_comp_path_df['Svtype'] == sv_type]
        #del Svtype
        comp_path_sv_type_df = comp_path_sv_type_df.drop(['Svtype'], axis=1)
        comp_path_sv_type_df.to_csv(comp_path_sv_type, sep="\t", header=True, index=False)
        comp_path_sv_type_head = comp_path_sv_type.replace(".vcf", "head.%s.vcf" % sv_type)
        base_path_sv_type_head = base_path_sv_type.replace(".vcf", "head.%s.vcf" % sv_type)
        os.system(
            "cat vcf.head.txt %s > %s" % (
            comp_path_sv_type, comp_path_sv_type_head))
        os.system(
            "cat vcf.head.txt %s > %s" % (
            base_path_sv_type, base_path_sv_type_head))
        comp_path = comp_path_sv_type_head
        base_path = base_path_sv_type_head
        output = os.path.join(work_dir, "%s_truvari_result" % sv_type)
        if os.path.exists(output):
            os.system("rm -rf %s" % output)
        os.system('{bgzip} -c {comp_path}  > {comp_path}.gz'.format(comp_path=comp_path, bgzip=bgzip))
        os.system('{tabix} -p vcf  {comp_path}.gz'.format(tabix=tabix, comp_path=comp_path))
        os.system('{bgzip} -c {base_path}  > {base_path}.gz'.format(base_path=base_path, bgzip=bgzip))
        os.system('{tabix} -p vcf {base_path}.gz'.format(tabix=tabix, base_path=base_path))
        os.system('{truvari} bench -b {base_path}.gz -c {comp_path}.gz -o {output}  --pctseq 0'.format(truvari=truvari,
                                                                                                       base_path=base_path,
                                                                                                       comp_path=comp_path,
                                                                                                       output=output))
        #  --pctseq 0  i.e., without considering the sequence information




##  Statistical Recall and Precision
for SV_length in sv_length_2:
    work_dir = os.path.join(work1_dir, 'tumor', 'SV_length_%s' % SV_length)
    recall_pre = os.path.join(work1_dir, "all_all_recall_precision.txt")
    output = os.path.join(work_dir, "truvari_result")
    summary = os.path.join(output, "summary.json")
    with open(summary, 'r') as f1:
        data = json.load(f1)
    with open(recall_pre, 'a') as f2:
        f2.write('{}\t{}\t{}\t{}\n'.format(SV_length, data['recall'], data['precision'], data['f1']))



for SV_length in sv_length_2:
    work_dir = os.path.join(work1_dir, 'tumor', 'SV_length_%s' % SV_length)
    tumor_random_bed = os.path.join(work_dir, "tumor_random.bed")
    sv_vcf2_repeat = os.path.join(work_dir, "merge_minimap2_sniffles_visor.vcf")
    recall_pre = os.path.join(work1_dir, "all_sv_typeall_recall_precision.txt")
    for sv_type in ['DEL', 'INS', 'DUP', 'INV']:
        output = os.path.join(work_dir, "%s_truvari_result" % sv_type)
        summary = os.path.join(output, "summary.json")
        with open(summary, 'r') as f1:
            data = json.load(f1)
        with open(recall_pre, 'a') as f2:
            f2.write('{}\t{}\t{}\t{}\t{}\n'.format(SV_length,sv_type, data['recall'], data['precision'], data['f1']))




