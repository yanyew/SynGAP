#!/usr/bin/env python
# -*- coding:utf-8 -*-
import datetime
import os
import sys
import tarfile


def initdb(args):
    script_dir, filename = os.path.split(os.path.abspath(sys.argv[0]))
    if str(args.sp) != 'plant' and str(args.sp) != 'animal':
        print('The species type of masterdb should be plant or animal\n')
        sys.exit()
    masterdb_dir = str(script_dir) + '/bin/masterdb/' + str(args.sp)
    if not os.path.exists(masterdb_dir):
        os.makedirs(masterdb_dir, exist_ok=True)
    tar = tarfile.open(str(args.file), "r:gz")
    gz_file_list = tar.getnames()
    masterdb_dir_file_list = \
        [f for f in os.listdir(masterdb_dir) if os.path.isfile(os.path.join(masterdb_dir, f))]
    splist = masterdb_dir + '/species.list'
    global sp_list
    if os.path.exists(splist):
        if sorted(masterdb_dir_file_list) == sorted(gz_file_list):
            with open(splist, 'r') as file:
                lines = file.readlines()
            sp_list = [line.split('\t')[0] for line in lines]
            for sp in sp_list:
                if os.path.exists(masterdb_dir + '/' + sp + '.fa') and os.path.exists(masterdb_dir + '/' + sp + '.gff3'):
                    continue
                else:
                    tar.extractall(masterdb_dir)
                    tar.close()
                    print('Import masterbd done! You can check the species list in \033[0;35m' + splist + '\033[0m')
                    sys.exit()
            print("The masterdb has been imported.\n"
                  "Skip `\033[0;35msyngap initdb\033[0m`.")
        else:
            tar.extractall(masterdb_dir)
            tar.close()
            print('Import masterbd done! You can check the species list in \033[0;35m' + splist + '\033[0m')
            sys.exit()
    else:
        tar.extractall(masterdb_dir)
        tar.close()
        print('Import masterbd done! You can check the species list in \033[0;35m' + splist + '\033[0m')
        sys.exit()