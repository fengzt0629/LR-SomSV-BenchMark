import os
import matplotlib.pyplot as plt
import pandas as pd
from multiprocessing import Pool

def bam_tag(script_tag,fq,new_ID,somatic_sample_dir,fq_dir,hg38_fa,stdout,stderr):
    cmd = '/bin/bash {script} {a} {b} {c} {d} {e} 1 > {out} 2 > {err}'.format(
        script=script_tag, a=fq, b=new_ID, c=somatic_sample_dir, d=fq_dir, e=hg38_fa,out=stdout, err=stderr)
    os.system(cmd)
    print(cmd)



def ms_bam(script_ms_x,stdout,stderr):
    os.system( '/bin/bash {script} 1 > {out} 2 > {err} '.format(
        out=stdout, err=stderr, script=script_ms_x))


def sniffles2_vcftobed(vcf, bed, support):
    """vcf to bed
    header: chrom start end (chr start end) svtype id length RE RNAMES IMPRECISE/PRECISE STD STRAND
    1-based to 0-based
    """
    bedout = open(bed, 'w')
    with open(vcf, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            else:
                content = line.split()
                chr1 = content[0]
                start = content[1]
                svid = content[2]
                info = content[7].split(';')
                fil = info[0]
                svtype = info[1].split('=')[1]
                if chr1 not in ['chrM', 'chrY']:
                    if svtype == 'BND':
                        # if ']' in content[4]:
                        #     content4 = content[4].split(']')
                        # else:
                        #     content4 = content[4].split('[')
                        # end = content4[1].split(':')[1]
                        end = content[4].replace('N', '').split(':')[1][:-1]
                        chr2 = info[-2].split('=')[1]
                        if chr2 not in ['chrM', 'chrY']:
                            svlen = 1
                            re = info[2].split('=')[1]
                            if int(re) >= support:
                                rnames = info[3].split('=')[1]
                                # std = info[-1].split('=')[1]
                                # suptype = info[-4].split('=')[1]
                                anotation = ['TRA', svid, str(abs(int(svlen))), re, rnames, fil]
                                newline1 = '\t'.join(
                                    [chr1, str(int(start) - 1), start, chr2, str(int(end) - 1), end] + anotation) + '\n'
                                newline2 = '\t'.join(
                                    [chr2, str(int(end) - 1), end, chr1, str(int(start) - 1), start] + anotation) + '\n'
                                newline = newline1 + newline2
                                bedout.write(newline)
                    else:
                        svlen = info[2].split('=')[1]
                        re = info[4].split('=')[1]
                        if int(re) >= support:
                            if abs(int(svlen)) >= 50:
                                rnames = info[5].split('=')[1]
                                # std = info[-1].split('=')[1]
                                # suptype = info[-4].split('=')[1]
                                end = info[3].split('=')[1]
                                anotation = ["NA","NA","NA",svtype, svid, str(abs(int(svlen))), re, rnames, fil]
                                newline = '\t'.join([chr1, str(int(start) - 1), end] + anotation) + '\n'
                                bedout.write(newline)
                    bedout.flush()
    bedout.close()


def cuteSV_vcftobed(vcf, bed, support):
    """vcf to bed
    header: chrom start end (chr start end) svtype id length RE RNAMES IMPRECISE/PRECISE STD STRAND
    1-based to 0-based
    """
    bedout = open(bed, 'w')
    with open(vcf, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            else:
                content = line.split()
                chr1 = content[0]
                start = content[1]
                svid = content[2]
                info = content[7].split(';')
                fil = info[0]
                svtype = info[1].split('=')[1]
                if chr1 not in ['chrM', 'chrY']:
                    if svtype == 'BND':
                        # if ']' in content[4]:
                        #     content4 = content[4].split(']')
                        # else:
                        #     content4 = content[4].split('[')
                        # end = content4[1].split(':')[1]
                        end = content[4].replace('N', '').split(':')[1][:-1]
                        # chr2 = info[-2].split('=')[1]
                        chr2 = content[4].replace('N','').split(':')[0][1:] # cuteSV_vcf
                        if chr2 not in ['chrM', 'chrY']:
                            svlen = 1
                            re = info[2].split('=')[1]
                            if int(re) >= support:
                                rnames = info[3].split('=')[1]
                                # std = info[-1].split('=')[1]
                                # suptype = info[-4].split('=')[1]
                                anotation = ['TRA', svid, str(abs(int(svlen))), re, rnames, fil]
                                newline1 = '\t'.join(
                                    [chr1, str(int(start) - 1), start, chr2, str(int(end) - 1), end] + anotation) + '\n'
                                newline2 = '\t'.join(
                                    [chr2, str(int(end) - 1), end,chr1, str(int(start) - 1), start] + anotation) + '\n'
                                newline = newline1 + newline2
                                bedout.write(newline)
                    else:
                        svlen = info[2].split('=')[1]
                        # re = info[4].split('=')[1] #sniffles2 vcf
                        re = info[6].split('=')[1] # cuteSV_vcf
                        if svtype  == 'DUP' or svtype == 'INV': # cuteSV_vcf
                            re = info[4].split('=')[1] # cuteSV_vcf
                        print(re)
                        # print(info)
                        if int(re) >= support:
                            if abs(int(svlen)) >= 50:
                                # rnames = info[5].split('=')[1] #sniffles2 vcf
                                rnames = info[-1].split('=')[1] # cuteSV_vcf
                                if svtype  == 'DUP' or svtype == 'INV': # cuteSV_vcf
                                    rnames = info[6].split('=')[1] # cuteSV_vcf
                                elif svtype == 'DEL':
                                    rnames = info[-2].split('=')[1] # cuteSV_vcf
                                end = info[3].split('=')[1]
                                anotation = ["NA","NA","NA",svtype, svid, str(abs(int(svlen))), re, rnames, fil]
                                newline = '\t'.join([chr1, str(int(start) - 1), end] + anotation) + '\n'
                                bedout.write(newline)
                    bedout.flush()
    bedout.close()

