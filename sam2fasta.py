#!/usr/bin/python

import sys

infile = sys.argv[1]

with open(infile) as f:
    for line in f:
        if not line.startswith("@"):
            data = line.strip().split("\t")
            print ">%s\n%s" % (data[0], data[9])
