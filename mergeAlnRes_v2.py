#!/usr/bin/python

import sys
import collections
import numpy

end_index = int(sys.argv[1])
prefix = sys.argv[2]

record = {}
b1_dist = {}
b2_dist = {}

diff_file = "/home/chialanghsu/nanopore/db/cyp.diff.db"
diff_db = []
with open(diff_file) as f:
    next(f)
    for line in f:
        data = line.strip().split()
        data.append(0)
        data.append(0)
        diff_db.append(data)

for i in range(1, end_index):
    fn = "read." + str(i) + ".aln"

    ## ignoring no-fusion reads.
    counter = {'B1': 0, 'B2': 0, 'N': 0, 'M': 0}
    with open(fn) as f:
        for line in f:
            data = line.strip().split("\t")
            counter[data[6]] += 1
    total = counter["B1"] + counter["B2"] + 0.0
    ratio = counter["B1"] / total

    if (ratio > 0.2 and ratio < 0.8):
        with open(fn) as f:
            for line in f:
                data = line.strip().split("\t")
                if data[6] == "B1":
                    if data[3] not in b1_dist: b1_dist[data[3]] = 0
                    b1_dist[data[3]] += 1
                elif data[6] == "B2":
                    if data[5] not in b2_dist: b2_dist[data[5]] = 0
                    b2_dist[data[5]] += 1

hd = open(prefix + ".align_dist.txt", "w")
hd.write("B1_gpos\tB1_nt\tB2_gpos\tB2_nt\tB1_pos\tB2_pos\tB1_freq\tB2_freq\n")
with open(diff_file) as f:
    next(f)
    for line in f:
        data = line.strip().split()
        b1_count = 0
        b2_count = 0
        if data[4] in b1_dist: b1_count = b1_dist[data[4]]
        if data[5] in b2_dist: b2_count = b2_dist[data[5]]
        hd.write("%s\t%s\t%s\t%s\t%s\t%s\t%i\t%i\n" % (data[0], data[1], data[2], data[3], data[4], data[5], b1_count, b2_count))

