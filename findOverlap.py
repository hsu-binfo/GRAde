#!/usr/bin/python

import sys

resFile = sys.argv[1]
diffFile = sys.argv[2]

res = {}
with open(resFile) as f:
    for line in f:
        data = line.strip().split("\t")
        uid = data.pop(0)
        res[uid] = data

index = 1
with open(diffFile) as f:
    for line in f:
        data = line.strip().split("\t")
        if data[4] in res and data[5] in res:
            print "%i\t%s\t%s\t%s" % (index, "\t".join(data), "\t".join(res[data[4]]), "\t".join(res[data[5]]))
        elif data[4] in res:
            print "%i\t%s\t%s\tNA\tNA\tNA" % (index, "\t".join(data), "\t".join(res[data[4]]))
        elif data[5] in res:
            print "%i\t%s\tNA\tNA\tNA\t%s" % (index, "\t".join(data), "\t".join(res[data[5]]))
        index += 1

