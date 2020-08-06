
import sys
import re

b1 = sys.argv[1]
b2 = sys.argv[2]

prog = re.compile(r'\w+\:\s+(\d+)\s+[\w\-]+')

B1_edit_table = {}
with open(b1) as f:
    for _ in range(4):
        next(f)
    for target in f:
        aln = f.next().strip()
        query = f.next()
        blank = f.next()
        target_seq = re.search(r'\s+([ATCG-]+)\s+', target).group(1)
        query_seq = re.search(r'\s+([ATCG-]+)\s+', query).group(1)
        target_site = int(prog.search(target).group(1)) - 1
        query_site = int(prog.search(query).group(1)) - 1
        for i in range(len(aln)):
            #if target_seq[i] != "-": target_site += 1
            if query_seq[i] != "-": 
            	query_site = int(query_site)
            	query_site += 1
            else:
            	query_site += 0.1

            origin = query_seq[i]
            modified = target_seq[i]
            pattern = origin + modified
            if origin == "-": origin = ''
            if modified == "-": modified = ''
            if aln[i] == "*" or aln[i] == " ":
            	#print query_site, origin, modified, pattern
            	B1_edit_table[query_site] = {"origin": origin, "modified": modified, "pattern": pattern}

out_seq = ''
with open(b2) as f:
    for _ in range(4):
        next(f)
    for target in f:
        aln = f.next().strip()
        query = f.next()
        blank = f.next()
        target_seq = re.search(r'\s+([ATCG-]+)\s+', target).group(1)
        query_seq = re.search(r'\s+([ATCG-]+)\s+', query).group(1)
        target_site = int(prog.search(target).group(1)) - 1
        query_site = int(prog.search(query).group(1)) - 1
        for i in range(len(aln)):
            #if target_seq[i] != "-": target_site += 1
            if query_seq[i] != "-": 
            	query_site = int(query_site)
            	query_site += 1
            else:
            	query_site += 0.1

            origin = query_seq[i]
            modified = target_seq[i]
            pattern = origin + modified

            if origin == "-": origin = ''
            if modified == "-": modified = ''

            if query_site in B1_edit_table and pattern == B1_edit_table[query_site]["pattern"]:
            	out_seq += modified
            else:
            	out_seq += origin

print ">temp_fasta"
print out_seq

#diff_match(b2, "B1", db, reg)
#diff_match(b1, "B2", db, reg)

