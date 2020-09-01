
import sys
import re
import os
import glob

path = os.path.dirname(os.path.abspath(sys.argv[0]))
source = os.path.join(path, "../db")

#fn = sys.argv[1]
b1_min = 84
b1_max = 3965

b2_min = 18
b2_max = 3906

readDir = sys.argv[1]
#index = int(sys.argv[2])

prog = re.compile(r'\w+\:\s+(\d+)\s+[\w\-]+')

def fasta_parsing(fn, trim_min, trim_max):
	seq = ''
	with open(fn) as f:
		next(f)
		for line in f:
			seq = seq + line.strip()
	return seq[trim_min-1:trim_max]

def aln_parsing(fn, trim_min, trim_max):
	start_pos = None
	stop_pos = None
	out_seq = [' '] * (trim_max - trim_min + 1)
	with open(fn) as f:
	    for _ in range(4):
	        next(f)
	    for target in f:
	        aln = f.next().strip()
	        query = f.next()
	        blank = f.next()
	        target_seq = re.search(r'\s+([ATCG-]+)\s+', target).group(1)
	        query_seq = re.search(r'\s+([ATCG-]+)\s+', query).group(1)
	        target_site = int(prog.search(target).group(1)) - 1
	        query_site = int(prog.search(query).group(1))

	        #print aln
	        for i in range(len(target_seq)):
	            #if target_seq[i] != "-": target_site += 1
	            if target_seq[i] != "-":
	            	target_site += 1
	            	if (target_site >= trim_min and target_site <= trim_max):
	            		out_seq[target_site - trim_min] = query_seq[i]

	return(out_seq)

def major_sites():
	cyp_aln = os.path.join(source, "cyp_4k.aln")
	b1_major = []
	b2_major = []
	with open(cyp_aln) as f:
	    for _ in range(4):
	        next(f)
	    for target in f:
	        aln = f.next().strip()
	        query = f.next()
	        blank = f.next()
	        b1_seq = re.search(r'\s+([ATCG-]+)\s+', target).group(1)
	        b2_seq = re.search(r'\s+([ATCG-]+)\s+', query).group(1)
	        b1_site = int(prog.search(target).group(1)) - 1
	        b2_site = int(prog.search(query).group(1))-1

	        #print aln
	        for i in range(len(b1_seq)):
	            if b1_seq[i] != "-":
	            	b1_site += 1
	            	if b1_seq[i] != b2_seq[i]:
	            		b1_major.append(b1_site)

	        for j in range(len(b2_seq)):
	        	if b2_seq[i] != "-":
	        		b2_site += 1
	        		if b2_seq[j] != b1_seq[j]:
	        			b2_major.append(b2_site)
	return(b1_major, b2_major)

b1_fa = fasta_parsing(os.path.join(source, 'CYP11B1_4k.fasta'), b1_min, b1_max)
b2_fa = fasta_parsing(os.path.join(source, 'CYP11B2_4k.fasta'), b2_min, b2_max)

b1_major, b2_major = major_sites()


print "Pos\tCount\tMis\tRate\tMajor\tGene"
### CYP11B1
b1_count = [0.0] * (b1_max - b1_min + 1)
b1_mismatch = [0.0] * (b1_max - b1_min + 1)
b1_rate = [0.0] * (b1_max - b1_min + 1)

for fn in glob.glob(os.path.join(readDir, '*.b1.aln')):
	s = aln_parsing(fn, b1_min, b1_max)
	for i in range(len(s)):
		if s[i] != ' ':
			b1_count[i] += 1
			if s[i] != b1_fa[i]:
				b1_mismatch[i] += 1

for j in range(len(b1_count)):
	pos = b1_min + j
	count = 1 if b1_count[j] == 0 else b1_count[j]
	mismatch = b1_mismatch[j]
	rate =  mismatch / count
	isMajor = 'Yes' if pos in b1_major else 'No'
	print "%i\t%i\t%i\t%.3f\t%s\t%s" %  (pos, count, mismatch, rate, isMajor, "CYP11B1")

b2_count = [0.0] * (b2_max - b2_min + 1)
b2_mismatch = [0.0] * (b2_max - b2_min + 1)
b2_rate = [0.0] * (b2_max - b2_min + 1)

#### CYP11B2
for fn in glob.glob(os.path.join(readDir, '*.b2.aln')):
	s = aln_parsing(fn, b2_min, b2_max)
	for i in range(len(s)):
		if s[i] != ' ':
			b2_count[i] += 1
			if s[i] != b2_fa[i]:
				b2_mismatch[i] += 1

for j in range(len(b2_count)):
	pos = b2_min + j
	count = 1 if b2_count[j] == 0 else b2_count[j]
	mismatch = b2_mismatch[j]
	rate =  mismatch / count
	isMajor = 'Yes' if pos in b2_major else 'No'
	print "%i\t%i\t%i\t%.3f\t%s\t%s" %  (pos, count, mismatch, rate, isMajor, "CYP11B2")
