#!/usr/bin/env python

#script to sum columns and output error rates from usearch9 mapping to reference and custom
#user field of -userfield ql+ids+mism+diffs

import sys

total = 0
count = 0
identical = 0
mismatches = 0
indels = 0

with open(sys.argv[1], 'rU') as input:
    for line in input:
        count += 1
        line = line.replace('\n', '')
        cols = line.split('\t')
        total += int(cols[0])
        identical += int(cols[1])
        mismatches += int(cols[2])
        ind = int(cols[3]) - int(cols[2])
        indels += ind

subs_errors = mismatches / float(total) * 100
indel_errors = indels / float(total) * 100

print "------------------------------"
print "Summarizing %s" % sys.argv[1]
print "------------------------------"
print "Aligned reads: "+ "{0:,}".format(count) + " ("+ "{0:,}".format(total)+" bp)"
print "Substitution errors: "+ "{0:.4f}%".format(subs_errors)
print "InDel errors: "+ "{0:.4f}%".format(indel_errors)
     
        