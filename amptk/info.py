#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import sys
import os
import pandas as pd
from amptk import amptklib
from amptk.__version__ import __version__


def getVersion():
    git_version = amptklib.git_version()
    if git_version:
        version = __version__+'-'+git_version
    else:
        version = __version__
    return version


def main():
    parentdir = os.path.join(os.path.dirname(amptklib.__file__))
    db_list = []
    okay_list = []
    search_path = os.path.join(parentdir, 'DB')
    for file in os.listdir(search_path):
        file_data = []
        if file.endswith(".udb"):
            if file.startswith('.'):
                continue
            okay_list.append(file)
            info_file = file + '.txt'
            file_data.append(file.rstrip('.udb'))
            with open(os.path.join(search_path, info_file), 'r') as info:
                for line in info:
                    line = line.strip()
                    if '\t' in line:
                        file_data += line.split('\t')
                    else:
                        file_data += line.split(' ')
            db_list.append(file_data)

    if len(db_list) < 1:
        db_print = "No DB configured, run 'amptk install' or 'amptk database' command."
    else:
        df = pd.DataFrame(db_list)
        df.columns = ['DB_name', 'DB_type', 'FASTA', 'Fwd Primer',
                      'Rev Primer', 'Records', 'Source', 'Version', 'Date']
        dfsort = df.sort_values(by='DB_name')
        db_print = dfsort.to_string(index=False, justify='center')
    print('------------------------------')
    print('Running AMPtk v {:}'.format(getVersion()))
    print('------------------------------')
    print('Taxonomy Databases Installed: {:}'.format(os.path.join(
        parentdir, 'DB')))
    print('------------------------------')
    print(db_print)
    print('------------------------------')


if __name__ == "__main__":
    main()
