# #!/usr/bin/python

# import sys
# import collections
# import numpy

# end_index = int(sys.argv[1])
# prefix = sys.argv[2]

# record = {}
# b1_dist = {}
# b2_dist = {}
# for i in range(1, end_index):
#     fn = "read." + str(i) + ".aln"

#     ## ignoring no-fusion reads.
#     counter = {'B1': 0, 'B2': 0, 'N': 0}
#     with open(fn) as f:
#         for line in f:
#             data = line.strip().split("\t")
#             counter[data[6]] += 1

#     total = counter["N"] + counter["B1"] + counter["B2"] + 1.0
#     if (counter["B1"]/total > 0.2 and counter["B1"]/total < 0.8):
#         with open(fn) as f:
#             for line in f:
#                 data = line.strip().split("\t")
#                 if data[3] != "":
#                     if data[3] not in b1_dist:
#                         b1_dist[data[3]] = {"B1": 0, "B2": 0, "N": 0}
#                     b1_dist[data[3]][data[6]] += 1
#                 if data[5] != "":
#                     if data[5] not in b2_dist:
#                         b2_dist[data[5]] = {"B1": 0, "B2": 0, "N": 0}
#                     b2_dist[data[5]][data[6]] += 1



# hd = open(prefix + ".b1_dist.txt", "w")
# hd.write("Position\tB1\tB2\tUnknown\n")
# keylist = b1_dist.keys()
# keylist.sort()
# for k in keylist:
#     v = b1_dist[k]
#     hd.write("%s\t%i\t%i\t%i\n" % (k, v["B1"], v["B2"], v["N"]))
# hd.close()

# hd = open(prefix + ".b2_dist.txt", "w")
# hd.write("Position\tB1\tB2\tUnknown\n")
# keylist = b2_dist.keys()
# keylist.sort()
# for k in keylist:
#     v = b2_dist[k]
#     hd.write("%s\t%i\t%i\t%i\n" % (k, v["B1"], v["B2"], v["N"]))
# hd.close()
