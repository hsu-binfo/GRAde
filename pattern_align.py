#!/usr/bin/python

import sys
import re

b1 = sys.argv[1]
b2 = sys.argv[2]
b12 = "/home/chialanghsu/nanopore/db/cyp.aln"

prog = re.compile(r'\w+\:\s+(\d+)\s+[\w\-]+')

align = {}
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
            	query_site += 1
            if target_seq[i] != "-":
                target_site += 1

            query_base = query_seq[i]
            if query_base != "-":
                align[query_site] = [query_base, target_seq[i], 142873254 + target_site - 1, '', '']
            #origin = query_seq[i]
            #modified = target_seq[i]
            #pattern = origin + modified
            #if origin == "-": origin = ''
            #if modified == "-": modified = ''
            #if aln[i] == "*" or aln[i] == " ":
            	#print query_site, origin, modified, pattern
            	#B1_edit_table[query_site] = {"origin": origin, "modified": modified, "pattern": pattern}

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
                query_site += 1
            if target_seq[i] != "-":
                target_site += 1

            query_base = query_seq[i]
            if query_base != "-" and query_site in align:
                align[query_site][3] = target_seq[i]
                align[query_site][4] = 142910559 + target_site - 1

current_seq = "N"
b1_count = 3
b2_count = 3
min_offset = 5

b1_break = ''
b2_break = ''
for k, v in align.iteritems():
    if v[0] == v[1] and v[0] != v[3]:
        #b2_break = v[4]
        current_seq = "B1"
        if b1_count < min_offset:
            b1_count += 1
        if b2_count > 0:
            b2_count -= 1
    elif v[0] == v[3] and v[0] != v[1]:
        #b1_break = v[2]
        current_seq = "B2"
        if b2_count < min_offset:
            b2_count += 1
        if b1_count > 0:
            b1_count -= 1
    #elif v[0] != v[1] and v[0] != v[3]:
    #    current_seq = "N"

    print "%i\t%s\t%s\t%s\t%s\t%s\t%s\t%i\t%i" % (k, v[0], v[1], str(v[2]), v[3], str(v[4]), current_seq, b1_count, b2_count)

    if b1_count != b2_count:
        current_seq = "B1" if b1_count > b2_count else "B2"

#diff_match(b2, "B1", db, reg)
#diff_match(b1, "B2", db, reg)

