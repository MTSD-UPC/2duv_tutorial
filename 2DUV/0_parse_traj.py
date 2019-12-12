
#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2018 Hao Ren <renh.cn@gmail.com>
# altered by Qian Zhang 2019 <zhangqian.allen@gmail.com>
#
# Distributed under terms of the LGPLv3 license.

"""
Parse the MD trajectory in pdb format.
Each snapshot will be splitted into protein and solvent.
The protein part will be written into two files in xyz and pdb formats, respectively.
The solvent part will only be written into an xyz file 'job.sol'

We use mini-batch for swarm of snapshots along a long trajectory.
"""
import os
import argparse

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--i',help='input file name',type=str)
parser.add_argument('--s',help='batch size/default=10',type=int,default=10)
parser.add_argument('--b',help='number of snapshot',type=int,default=None)
parser.add_argument('--protein',type=int)
parser.add_argument('--water',type=int)
parser.add_argument('--ions',type=int)
args = parser.parse_args()


batch_size = args.s 
output_dir = 'snapshots'
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

n_batches = args.b // batch_size

def write_xyz(fname, pdblines):
    with open(fname, 'w') as fh:
        fh.write('{}\n\n'.format(len(pdblines)))
        newlines = []
        for l in pdblines:
            label = l[13]
            xyz = l[28:56]
            newlines.append(label + xyz + '\n')
        fh.writelines(newlines)
    return


with open(args.i, 'r') as fh:
    for i_batch in range(n_batches):
        batch_path = '{}/{:05d}'.format(output_dir, (i_batch+1)*10)
        if not os.path.exists(batch_path):
            os.mkdir(batch_path)
        for i_sample in range(batch_size):
            # read the i_snap-th snapshot, and write into a subdirectory
            i_snap = i_batch * batch_size + i_sample
            snap_path = '{}/{:05d}'.format(batch_path, i_snap)
            if not os.path.exists(snap_path):
                os.mkdir(snap_path)
            protein = []
            solvent = []

            if not os.path.exists(snap_path):
                os.mkdir(snap_path)
            for i in range(5):
                fh.readline() # 5-line header
            for i in range(args.protein):
                # 4672 lines for protein, no HEM
                protein.append(fh.readline())
            #for i in range(188):
                # 188 lines for HEM
            #    fh.readline()
            for i in range(args.water):
                # lines from 4678 to 129612 for solvent
                solvent.append(fh.readline())
            for i in range(args.ions):
                fh.readline() # Cl ions
            for i in range(2):
                fh.readline() # 2 lines tail

            fh_jobpdb = open('{}/job.pdb'.format(snap_path), 'w')
            fh_jobpdb.writelines(protein)
            write_xyz('{}/job.xyz'.format(snap_path), protein)
            write_xyz('{}/job.sol'.format(snap_path), solvent)
            print("Processed snapshot #{}".format(i_snap))
print('Done successfully!')


