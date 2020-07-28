#!/usr/bin/python

import sys
import re

infile = "cyp.aln"

prog = re.compile(r'\w+\:\s+(\d+)\s+[\w\-]+')

B1_genome = 142873253
B2_genome = 142910558

print "B1\tB1N\tB2\tB2N\tB1_genome\tB2_genome"
with open(infile) as f:
    for _ in range(4):
        next(f)
    for target in f:
        aln = f.next().strip()
        query = f.next()
        blank = f.next()

        target_seq = re.search(r'\s+([ATCG-]+)\s+', target).group(1)
        query_seq = re.search(r'\s+([ATCG-]+)\s+', query).group(1)

        #print target_seq

        target_site = int(prog.search(target).group(1)) - 1
        query_site = int(prog.search(query).group(1)) - 1

        for i in range(len(aln)):
            if target_seq[i] != "-": target_site += 1
            if query_seq[i] != "-": query_site += 1

            if aln[i] == "*" or aln[i] == " ":
                print "%i\t%s\t%i\t%s\t%i\t%i" % (target_site, target_seq[i], query_site, query_seq[i], B1_genome + target_site, B2_genome + query_site)



