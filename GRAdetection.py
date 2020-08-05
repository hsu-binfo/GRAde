#from sqlalchemy.ext.declarative import declarative_base
#import vcf
import logging
import json
import argparse
#from lib.common import Common
import os
#from lib.common import MetaClass
import sys
import subprocess
import string
import gzip

#class mMinION(Common):
class GARdetection():
    #__metaclass__ = MetaClass
    #base = declarative_base()
	def __init__(self, args, path):
		self.outdir = args.outdir
		self.samplename = args.samplename
		self.fastq = args.fastq
		self.minReadLength = args.minReadLength
		self.maxReadLength = args.maxReadLength
		self.thread = args.thread
		self.errRate = args.correctedErrorRate
		self.path = path
		self.scripts = os.path.join(path, "scripts")
		self.db = os.path.join(path, "db")
		#self.sangerresults = sangerresults
		#self.confdir = "/data1/BINFO0039_GRA/scripts"

		with open(os.path.join(path, 'config.json')) as json_data_file:
			self.config = json.load(json_data_file)

		# set up logger
		self.logger = logging.getLogger(__name__)
		self.logger.setLevel(logging.INFO)
		# complete setup of logger
		handler = logging.FileHandler(self.outdir + "/pipeline.log", 'w')
		handler.setLevel(logging.INFO)
		formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
		handler.setFormatter(formatter)
		console = logging.StreamHandler()
		console.setLevel(logging.INFO)

		self.logger.addHandler(handler)
		self.logger.addHandler(console)

	def filter(self):
		#samtools view -h N4.trimmed.sorted.bam | awk 'length($10) > 1000 || $1 ~ /^@/' | samtools view -bS > N4.trimmed.sorted.filtered.bam
		tempFq = self.outdir + "/" + self.samplename + ".temp.fq"
		tempBam = self.outdir + "/" + self.samplename + ".uncorrected.bam"
		outfile = self.outdir + "/" + self.samplename + ".filtered.fastq"

		## mapping onto chr8 and excluding too short reads
		mapping = """${ngmlr} -r ${ref} -t ${thread} -q ${fq} -x ont | \
			${samtools} sort -@ ${thread} -O BAM > ${out}"""

		command = string.Template(mapping).safe_substitute({
			"ref": os.path.join(self.db, self.config["db"]["genome"]),
			"fq" : self.fastq,
			"len": self.minReadLength,
			"thread": self.thread,
			"out": tempBam,
			"ngmlr": self.config["software"]["ngmlr"],
			"samtools": self.config["software"]["samtools"]
			})
		print("[COMMEND]", command)
		subprocess.call(command, shell = True)

		## index temp bam file
		bamindex = "${samtools} index %s" % tempBam
		command = string.Template(bamindex).safe_substitute({
			"samtools": self.config["software"]["samtools"]
		})
		print("[COMMEND]", command)
		subprocess.call(command, shell = True)

		## retrive reads at the CYP11B1/CYP11B2 loci
		template = """${samtools} view -h ${bam} 'chr8:142872206-142918080' | \
			awk 'length($10) > ${minlen} || $1 ~ /^@/' | \
			awk 'length($10) < ${maxlen} || $1 ~ /^@/' | \
			${samtools} view -bS | \
			${bedtools} bamtofastq -i stdin -fq ${out}"""

		command = string.Template(template).safe_substitute({
			"bam": tempBam,
			"out": outfile,
			"minlen": self.minReadLength,
			"maxlen": self.maxReadLength,
			"samtools": self.config["software"]["samtools"],
			"bedtools": self.config["software"]["bedtools"]
			})

		subprocess.check_call(command, shell = True)

		#os.remove(tempBam)
		#os.remove(tempBam + ".bai")
		return outfile

	def read_correction(self, fastq):
		#correct reads
		template = """${canu} -correct -p ${samplename} -d ${outdir} \
			genomeSize=5k overlapper=mhap utgReAlign=true stopOnReadQuality=false \
			-nanopore-raw ${fq}"""

		command = string.Template(template).safe_substitute({
			"samplename": self.samplename,
			"outdir": self.outdir,
			"fq": fastq,
			"canu": self.config["software"]["canu"]
			})
		print("[COMMEND]",command)
		subprocess.check_call(command, shell = True)

		corrected = self.outdir + "/" + self.samplename + ".correctedReads.fasta.gz"
		return corrected

	def mapping(self, fastq):
		outfile = self.outdir + "/" + self.samplename + ".ngmlr.bam"
		mapping = """${ngmlr} -r ${ref} -t ${thread} -q ${fq} -x ont | \
			${samtools} sort -@ ${thread} -O BAM > ${out}"""

		command = string.Template(mapping).safe_substitute({
			"ref": os.path.join(self.db, self.config["db"]["genome"]),
			"fq" : fastq,
			"thread": self.thread,
			"out": outfile,
			"ngmlr": self.config["software"]["ngmlr"],
			"samtools": self.config["software"]["samtools"]
			})
		print(command)
		subprocess.check_call(command, shell = True)
		## index temp bam file
		subprocess.check_call("samtools index %s" % outfile, shell = True)

		return outfile

	def mapping_minimap2(self, fastq):
		outfile = self.outdir + "/" + self.samplename + ".minimap2.bam"
		mapping = """${minimap2} -ax map-ont -t ${thread} -a ${ref} ${fq} | \
			${samtools} sort -@ ${thread} -O BAM > ${out}"""

		command = string.Template(mapping).safe_substitute({
			"ref": os.path.join(self.db, self.config["db"]["mmi_index"]),
			"fq" : fastq,
			"thread": self.thread,
			"out": outfile,
			"minimap2": self.config["softwate"]["minimap2"],
			"samtools": self.config["software"]["samtools"]
			})
		print(command)
		subprocess.check_call(command, shell = True)
		## index temp bam file
		subprocess.check_call("samtools index %s" % outfile, shell = True)

		return outfile

	def read_correction_2nd(self, fasta):
		outfile = os.path.join(self.outdir, self.samplename + ".summary.txt")

		tempDir = os.path.join(self.outdir, "ssw_temp")
		aln1Dir = os.path.join(tempDir, "aln1")
		aln2Dir = os.path.join(tempDir, "aln2")
		originDir = os.path.join(tempDir, "read.origin")
		editDir = os.path.join(tempDir, "read.edit")
		for dn in [aln1Dir, aln2Dir, originDir, editDir]:
			if not os.path.exists(dn): os.makedirs(dn)

		tempID = 1
		with gzip.GzipFile(fasta, 'rb') as f:
			lines = []
			for line in f:
				lines.append(line.strip().decode())
				if len(lines) == 2:
					fa1 = os.path.join(originDir, "read."+str(tempID)+".fa")
					fa2 = os.path.join(editDir, "read."+str(tempID)+".fa")

					fh = open(fa1, "w")
					fh.write("%s\n%s\n" % (lines[0], lines[1]))
					fh.close()

					b1 = self.ssw(fa1, os.path.join(self.db, self.config["db"]["cyp11b1"]), os.path.join(aln1Dir, "read."+str(tempID)+".b1.aln"))
					b2 = self.ssw(fa1, os.path.join(self.db, self.config["db"]["cyp11b2"]), os.path.join(aln1Dir, "read."+str(tempID)+".b2.aln"))

					run_read_edit = """%s %s %s > %s""" % (os.path.join(self.scripts, self.config["scripts"]["edit_read"]), b1, b2, fa2)
					subprocess.check_call(run_read_edit, shell = True)

					b1 = self.ssw(fa2, os.path.join(self.db, self.config["db"]["cyp11b1"]), os.path.join(aln2Dir, "read."+str(tempID)+".b1.aln"))
					b2 = self.ssw(fa2, os.path.join(self.db, self.config["db"]["cyp11b2"]), os.path.join(aln2Dir, "read."+str(tempID)+".b2.aln"))

					lines = []
					tempID += 1

		run_multialigner = """%s %s > %s""" % (os.path.join(self.scripts, self.config["scripts"]["aligner"]), aln2Dir, outfile)
		subprocess.check_call(run_multialigner, shell = True)
		run_funsion_plot = """Rscript %s -f %s -o %s -s %s -b %s""" \
			% (os.path.join(self.scripts, self.config["scripts"]["fusionPlot"]), outfile, self.outdir, self.samplename, os.path.join(self.db, self.config["db"]["blackList"]))
		subprocess.check_call(run_funsion_plot, shell = True)
		os.system('rm -rf %s' % tempDir)

	def ssw(self, fa, ref, out):
		template = """${ssw_test} -c ${ref} ${fa} > ${out}"""
		cmd = string.Template(template).safe_substitute({
			"ref": ref,
			"fa": fa,
			"out": out,
			"ssw_test": self.config["software"]["ssw_test"]
			})
		# print("[COMMAND]", cmd)
		subprocess.call(cmd, shell = True)
		#print(cmd)
		return(out)

	def fqTofa(self, fq):
		fa = os.path.join(self.outdir, self.samplename + ".filtered.fasta")
		subprocess.check_call("sed -n '1~4s/^@/>/p;2~4p' %s > %s" % (fq, fa), shell = True)
		subprocess.check_call("gzip %s" % fa, shell = True)

		outfile = fa + ".gz"
		return(outfile)

def main():
	parser = argparse.ArgumentParser(description='mMinION - Run fastq file for a sample through mitochondraial analysis pipeline', usage='%(prog)s [options]', add_help=False)
	parser.add_argument('--samplename',help='[Required] The name of the sample being processed - files will take this name', required = True)
	parser.add_argument('--fastq',help='[Required] Full path to the Fastq file', required = True)
	parser.add_argument('--outdir',help='[Required] The output directory for all files generated', required = True)
	parser.add_argument('--correctedErrorRate', required = False, default = 0.8)
	parser.add_argument('--minReadLength', required = False, default = 3000)
	parser.add_argument('--maxReadLength', required = False, default = 5000)
	parser.add_argument('--thread', required = False, default = 8)
	parser.add_argument('--skipCorrection', required = False, action = 'store_true')

	toolDir = os.path.dirname(os.path.abspath(sys.argv[0]))

	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit(0)

	args = parser.parse_args()

	if not os.path.exists(args.outdir):
		os.makedirs(args.outdir)

	m = GARdetection(args, toolDir)
	fq = m.filter()

	if args.skipCorrection:
		fa = m.fqTofa(fq)
		m.read_correction_2nd(fa)
		bam = m.mapping(fq)
	else:
		corrected_fasta = m.read_correction(fq)
		m.read_correction_2nd(corrected_fasta)
		bam = m.mapping(corrected_fasta)

	#bam2 = m.mapping_minimap2(corrected_fastq)

if __name__ == '__main__':
	main()
