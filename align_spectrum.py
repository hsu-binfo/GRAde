#!/usr/bin/python

import sys
import re

b1 = sys.argv[1]
b2 = sys.argv[2]

def createDB():
    dbFile = "/home/chialanghsu/nanopore/db/cyp.diff.db"

    target = {"B1":{}, "B2": {}}
    region = {"B1":[], "B2": []}
    with open(dbFile) as f:
        next(f)
        for line in f:
            data = line.strip().split("\t")
            target["B1"][data[0]] = {'target': data[1], 'query': data[3], 'genome': data[5]}
            target["B2"][data[2]] = {'target': data[3], 'query': data[1], 'genome': data[4]}
            region["B1"].append(data[0])
            region["B2"].append(data[2])

    return([target, region])


window = 20
offset = 5

db, reg = createDB()
prog = re.compile(r'\w+\:\s+(\d+)\s+[\w\-]+')

def diff_match(fn, tag, db, reg):
    q = "B1" if tag == "B2" else "B2"
    db = db[q]
    reg = reg[q]
    pre_index = 0
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
            query_site = int(prog.search(query).group(1)) - 1
            for i in range(len(aln)):
                if target_seq[i] != "-": target_site += 1
                if query_seq[i] != "-": query_site += 1

                if aln[i] == "*" or aln[i] == " ":
                    if str(target_site) in db and \
                        db[str(target_site)]['target'] == target_seq[i] and \
                        db[str(target_site)]['query'] == query_seq[i]:

                        genome_pos = db[str(target_site)]['genome']
                        current_index = reg.index(str(target_site))
                        dist = current_index - pre_index
                        pre_index = current_index
                        print "%i\t%s\t%s\t%s\t%s\t%i\t%s" % (query_site, tag, target_site, target_seq[i], query_seq[i], dist, genome_pos)

diff_match(b2, "B1", db, reg)
diff_match(b1, "B2", db, reg)

