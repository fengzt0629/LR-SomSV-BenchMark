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
