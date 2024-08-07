

import os

import pandas as pd

from config import work1_dir,chr2_fa,sortBed


## Simulation of the tumor SV lenghth

def simulate_function(SV_length):
    # work1_dir = '/NAS/wg_fzt/benchmark/simu_somatic/INS_reads_random/new_result/01.LiuXuan_contral_SV_length'
    blood_bed = os.path.join(work1_dir, "geneline_random.bed")
    exclude_bed = os.path.join(work1_dir, "exclude_chr2.bed")
    work_dir = os.path.join(work1_dir, 'tumor', 'SV_length_%s' % SV_length)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    tumor_exclude_bed = os.path.join(work_dir, "tumor_exclude_chr2.bed")
    tumor_random_bed = os.path.join(work_dir, "tumor_random.bed")
    tumor_random_bed_2 = os.path.join(work_dir, "tumor_random_cat_blood.bed")
    os.system(
        "cat {blood_bed} {exclude_bed}  > {tumor_exclude_bed}".format(blood_bed=blood_bed, exclude_bed=exclude_bed,
                                                                      tumor_exclude_bed=tumor_exclude_bed))
    os.system(
        "Rscript randomregion.r -d chr2.dim.tsv -n 400 -l {SV_length} -x {tumor_exclude_bed} -v 'deletion,inversion,insertion,tandem duplication' -r '25:25:25:25' | {sortBed} > {tumor_random_bed}".format(
            SV_length=SV_length, tumor_exclude_bed=tumor_exclude_bed, tumor_random_bed=tumor_random_bed,sortBed=sortBed)) #randomregion.r chr2.dim.tsv are bulit for VISOR,we also giave in our 03.SV_lenth_emulate
    TB = 'SV_length_%s' % SV_length
    hack_out = os.path.join(work_dir, "%s_hack.out" % (TB))
    if os.path.exists(hack_out):
        os.system('rm -rf %s' % hack_out)
    chrom_dim_tsv = os.path.join(work_dir, "%s_haplochroms.dim.tsv" % (TB))
    maxdims_tsv = os.path.join(work_dir, "%s_maxdims.tsv" % (TB))
    shorts_laser_simple = os.path.join(work_dir, "%s_shorts.laser.simple.bed" % (TB))
    laser_1_out = os.path.join(work_dir, "%s_laser.1.out" % (TB))
    if os.path.exists(laser_1_out):
        os.system('rm -rf %s' % laser_1_out)
    os.system("cat %s %s > %s" % (tumor_random_bed, blood_bed, tumor_random_bed_2))
    os.system("VISOR HACk -g %s -b %s -o %s" % (chr2_fa, tumor_random_bed_2, hack_out))
    os.system("cut -f1,2 %s/*.fai %s.fai > %s" % (hack_out, chr2_fa, chrom_dim_tsv))
    os.system("cat %s | sort  | awk '$2 > maxvals[$1] {lines[$1]=$0;"
              " maxvals[$1]=$2} END { for (tag in lines) print lines[tag] }' > %s" % (chrom_dim_tsv, maxdims_tsv))
    os.system('''awk 'OFS=FS="\t"''{print $1, "1", $2, "80.0", "100.0"}' %s > %s ''' % (
        maxdims_tsv, shorts_laser_simple)) ###  80 60 50     "80.0", "100.0"     "80.0", "80.0"   "100.0", "100.0"
    os.system("VISOR LASeR -g %s -s %s -b %s -o %s --threads 30 --coverage 30  --fastq --tag" % (
    chr2_fa, hack_out, shorts_laser_simple, laser_1_out))



## tumor add tag
def tumor_add_tag_function(SV_length):
    TB = 'SV_length_%s' % SV_length
    work_dir = os.path.join(work1_dir, 'tumor', 'SV_length_%s' % SV_length)
    laser_1_out = os.path.join(work_dir, "%s_laser.1.out" % (TB))
    bam_file = os.path.join(laser_1_out, "sim.srt.bam")
    #add tag
    sh_tumor_add_tag = os.path.join(work1_dir, "agg_tag_tumor_1.sh") # in 03.SV_length_emulate
    sh_tumor_add_o = os.path.join(work_dir, "tumor_add_o.0")
    sh_tumor_add_e = os.path.join(work_dir, "tumor_add_e.e")
    cmd_1 = "/bin/bash {sh_tumor_add_tag} {work_dir}  1 > {sh_tumor_add_o} 2> {sh_tumor_add_e}".format(
        sh_tumor_add_tag=sh_tumor_add_tag, work_dir=laser_1_out, sh_tumor_add_o=sh_tumor_add_o,
        sh_tumor_add_e=sh_tumor_add_e)
    os.system(cmd_1)


def ms_bam(script_ms_x, stdout, stderr):
    os.system('/bin/bash {script} 1 > {out} 2 > {err}'.format(
        out=stdout, err=stderr, script=script_ms_x))



def bed_to_vcf(vcf_hg38_bed,SV_length ):
    df_bed = pd.read_csv(vcf_hg38_bed, sep="\t", header=None,names=['Chrom', 'start', 'end', 'svtype', 'seq', 'none'])
    df_bed['svtype'] = df_bed['svtype'].apply(lambda x:
                                      'DEL' if x =='deletion'
                                      else('INS' if x =='insertion'
                                           else ('INV' if x =='inversion'
                                                 else ('DUP' if x =='tandem duplication'
                                                       else 'BND'
                                                       )
                                                 )
                                           )
                                      )
    df = pd.DataFrame(columns=['Chrom1', 'Pos1', 'Chrom2', 'Pos2', 'SV_Size or breakpoints distance', 'SV_type', 'SV_id'])
    df['Chrom1'] = df_bed['Chrom']
    df['Pos1'] = df_bed['start']
    df['Chrom2'] = df_bed['Chrom']
    df['Pos2'] = df_bed['end']
    df['SV_Size or breakpoints distance'] = SV_length
    df['SV_type'] = df_bed['svtype']
    df['SV_id'] = df.apply(lambda x: 'SV' + str(x['Pos1']) + '_' + str(x['Pos2']), axis=1)
    # bulit vcf
    vcf_df = pd.DataFrame(columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'])
    vcf_df['#CHROM'] = df['Chrom1']
    vcf_df['POS'] = df['Pos1']
    vcf_df['ID'] = df['SV_id']
    vcf_df['REF'] = 'N'
    vcf_df['ALT'] = df.apply(
        lambda x: '<%s>' % x['SV_type'] if x['SV_type'] != 'BND' else 'N]{chr2}:{Pos2}]'.format(chr2=x['Chrom2'],
                                                                                                Pos2=x['Pos2']), axis=1)
    vcf_df['QUAL'] = '60'
    vcf_df['FILTER'] = 'GT'
    vcf_df['INFO'] = df.apply(lambda x: ';'.join(
        ['PRECISE', 'SVTYPE=%s' % x['SV_type'], 'SVLEN=%s' % x['SV_Size or breakpoints distance'],
         'END=%s' % (x['Pos2']), 'SUPPORT=10', 'RNAMES=1,2,3,4,5,6,7,8,9,10',
         'COVERAGE=15,28,28,30,32', 'STRAND=+-', 'AF=0.643', 'STDEV_LEN=11.508;STDEV_POS=79.536;SUPPORT_LONG=0']) if x[
                                                                                                                         'SV_type'] != 'BND' else ';'.join(
        ['PRECISE', 'SVTYPE=%s' % x['SV_type'], 'SUPPORT=4', 'RNAMES=1,2,3,4', 'COVERAGE=15,28,28,30,32', 'STRAND=+-',
         'AF=0.643', 'CHR2=%s' % str(x['Chrom2']), 'STDEV_POS=0.000']), axis=1)
    vcf_df['FORMAT'] = 'GT:GQ:DR:DV'
    vcf_df['SAMPLE'] = '0/0:60:27:1'
    # write vcf for Truvari
    #sort by '#CHROM', 'POS'
    vcf_df = vcf_df.sort_values(by=['#CHROM', 'POS'])
    output_file = vcf_hg38_bed.replace(".bed", "tumor.golde.no.head.vcf")
    vcf_df.to_csv(output_file, sep='\t', index=False)

