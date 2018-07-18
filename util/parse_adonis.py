#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *
import sys
import os
from natsort import natsorted
import pandas as pd

#try to get all results from folder of phyloseq_R results

#Adonis looks like this currently,  first is adonis test, then betadisper

'''

Call:
adonis(formula = rp ~ metadata[[treatments[y]]], permutations = 9999) 

Permutation: free
Number of permutations: 9999

Terms added sequentially (first to last)

                           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
metadata[[treatments[y]]]   1     0.868 0.86802  4.2181 0.02987 0.0192 *
Residuals                 137    28.193 0.20579         0.97013         
Total                     138    29.061                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 9999

Response: Distances
           Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
Groups      1 0.00137 0.0013715 0.0648   9999 0.8026
Residuals 137 2.89871 0.0211585                     
'''

#the name of the test is in the filename, and variable is also in the filename
'''
mhc.location.adonis-raupcrick.txt
mhc.DNA_quant_mean.adonis-braycurtis.txt

'''
    
#pass name of folder to script
input = sys.argv[1]
if not os.path.isdir(input):
    print("Error: script is expecting a folder as input")
    sys.exit(1)

Results = {}
headers = []
samples = []
for file in natsorted(os.listdir(input)):
    adonisp = '1.00'
    betap = '1.00'
    if 'adonis-' in file:
        basename = file.split('.')[0]
        variable = file.split('.')[1]
        distance = file.split('.')[2].replace('adonis-', '')
        if not distance in headers:
            headers.append(distance)
        if not basename in samples:
            samples.append(basename)
        file = os.path.join(input, file)
        with open(file, 'rU') as filein:
            for line in filein:
                line = line.replace('\n', '')
                if line.startswith('metadata'):
                    cols = line.split(' ')
                    if '*' in cols[-1] or cols[-1] == '.':
                        adonisp = cols[-2]
                    else:
                        adonisp = cols[-1]
                if line.startswith('Groups'):
                    cols = line.split(' ')
                    if '*' in cols[-1] or cols[-1] == '.':
                        betap = cols[-2]
                    else:
                        betap = cols[-1]

        if not variable in Results:
            Results[variable] = [float(adonisp), float(betap)]
        else:
            Results[variable].append(float(adonisp))
            Results[variable].append(float(betap))



df = pd.DataFrame.from_dict(Results, orient="index")
cols = pd.MultiIndex.from_product([samples, headers, ['Adonis', 'Betadisper']])
df.columns = cols
transposed = df.transpose()

if len(sys.argv) > 2:
    if sys.argv[2]:
        transposed.to_csv(sys.argv[2])
else:
    print(transposed.to_string())

