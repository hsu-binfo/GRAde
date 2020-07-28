#!/home/chialanghsu/miniconda3/bin/python
#################################################
## Basecalling by flappie, and demultiplex by qcat
#################################################
import sys
import os
import subprocess
import string

conda3 = os.environ['CONDA3']
indir = os.path.abspath(sys.argv[1])

root = os.path.dirname(indir)
fqdir = os.path.join(root, "fastq")

flappie = "/home/chialanghsu/usr/bin/flappie"
cmd = """%s ${fast5} > ${fastq}""" % flappie

if not os.path.exists(outdir):
    print("Error: %s does not exist" % outdir)
    sys.exit(1)

for sub in os.listdir(indir):
    for fast5 in os.listdir(os.path.join(indir, sub)):
        fastq = fast5.replace("fast5", "fastq")
        fast5 = os.path.join(indir, sub, fast5)
    fastq = os.path.join(fqdir, fastq)
    qcmd = string.Template(cmd).safe_substitute({
         "fast5": fast5,
         "fastq": fastq
         })
    #print qcmd
    subprocess.check_call(qcmd, shell = True)



