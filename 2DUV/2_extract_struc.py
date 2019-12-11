#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2018 Hao Ren <renh.cn@gmail.com>
# Modified by Qian Zhang  2019 <zhangqian.allen@gmail.com>
# Distributed under terms of the LGPLv3 license.

import re
import os 
import argparse
import numpy as np
import json
import tarfile 
import sys

parser = argparse.ArgumentParser(description='manul to this scipt')
parser.add_argument('-f','--file',help='input pdb file name',type=str)
parser.add_argument('-c','--chains',help='number of chains in pdb file; default=4',type=int,default=4)
parser.add_argument('-b','--snapshots',help='number of snapshot',type=int,default=None)
parser.add_argument('-p','--pattern',help='secondary structure pattern(str)',type=str)
parser.add_argument('-s','--ssfile',help='do_dssp file,default=ss.dat',type=str,default='ss.dat')
args = parser.parse_args()


def extract(line, pattern):
    result = []
    ah = re.findall(pattern, l) # find all alpha helices with length >= 10
    #print([len(x) for x in ah])
    start = 0
    for ih in ah:
        start = l.find(ih, start)
        end = start + len(ih)
        #print(start, len(ih), end)
        
        result.append([start+1, end, len(ih)]) # add 1 to count from 1
        # now 'end' is the index just the end of the segment
        start = end+1
    return result


def findIndex(chains, start, end):
    # chain_A, chain_B, chain_C, chain_D = chains
    seq = ''.join(chains)
    # seq = ''
    # for chain in chains:
    #     seq+=chain
    # find aromatic residue numbers
    aromatic = []
    i = 0
    for res in seq:
        if res in "FYW":
            aromatic.append((i, res))
        i += 1
    #print(aromatic)
    #print(len(aromatic))
    len_chains = [len(chain) for chain in chains]
    range_chains = np.cumsum(len_chains)
    #print(range_chains)
    for i, l in enumerate(list(range_chains)):
        #print(i, l)
        if end < l:
            break

    in_chain = i
    start_in_chain = start
    end_in_chain = end
    for i in range(in_chain):
        start_in_chain -= len_chains[i]
        end_in_chain -= len_chains[i]
    #print('your sequence resides in chain {}, from {} to {}'.format(
    #    'ABCD'[in_chain], start_in_chain, end_in_chain))

    # Find peptide exciton numbers
    peptide_excitons = [2*(x-1) for x in len_chains]
    base = 0
    for i in range(in_chain):
        base += peptide_excitons[i]
    #print(base)
    #print(start_in_chain, end_in_chain)
    pep_start_index = base + start_in_chain * 2 - 1
    pep_end_index = base + (end_in_chain-1) * 2
    #print('peptide indices from {} to {}'.format(pep_start_index, pep_end_index))

    # find aromatic exciton numbers
    # first count number of aromatic residues before this sequence
    n_before = 0
    for res in aromatic:
        if res[0] < start:
            n_before += 1
        else:
            break
    #print('There are {} aromatic residues before this sequence'.format(n_before))
    base = sum(peptide_excitons) + n_before * 4

    n_aromatic = 0
    for i in range(start, end+1):
        if seq[i] in "FWY": # if this residue is aromatic
            n_aromatic += 1
    #print('There are {} aromatic residues.'.format(n_aromatic))
    aromatic_start_index = base + 1
    aromatic_end_index = aromatic_start_index + n_aromatic * 4

    indicies = np.concatenate(
            (np.arange(pep_start_index, pep_end_index+1),
                np.arange(aromatic_start_index, aromatic_end_index))
            )

    return((in_chain, start_in_chain, end_in_chain),indicies)


if __name__ == "__main__":

    # step1. get sequence of protein
    pdb_file = args.file
    n_sanps = args.snapshots
    chains = []
    chain_name = [chr(i) for i in range(65,65+args.chains)] 
    aa_codes = {
        'ALA':'A','CYS':'C','ASP':'D','GLU':'E',
        'PHE':'F','GLY':'G','HIS':'H','LYS':'K',
        'ILE':'I','LEU':'L','MET':'M','ASN':'N',
        'PRO':'P','GLN':'Q','ARG':'R','SER':'S',
        'THR':'T','VAL':'V','TYR':'Y','TRP':'W'}
    with open(pdb_file,'r') as f1:
        lines = f1.readlines()
    for c_name in chain_name:
        seqi = ''
        for line in lines:
            colums = line.split()
            if colums[0] == 'SEQRES' and colums[2]==c_name:
                for resname in colums[4:]:
                    seqi += aa_codes[resname]
        chains.append(seqi)
    np.save('chains.npy',chains)

    #step2. extract secondary structural fragments
    fh = open(args.ssfile, 'r')
    pattern = args.pattern
    # pattern = 'HHHHHHHHHH*'  #represent that find all alpha helices with length >= 10

    n_residues = int(fh.readline())
    for i_snap in range(n_sanps):
        if i_snap % 10 == 0:
            print('working on snapshot #{}'.format(i_snap))
        batch_dir = 'snapshots/{:05d}'.format((i_snap // 10)*10 + 10)
        snap_dir = '{}/{:05d}'.format(batch_dir, i_snap)
        helix_path = '{}/helices'.format(snap_dir)
        if not os.path.exists(helix_path):
            os.mkdir(helix_path)
        out_fh = open('{}/helix.dat'.format(helix_path), 'w')
        out_fh.write('#start end nmodes chain start_c end_c mode-indicies ->\n')

        l = fh.readline().strip()
        l = ''.join(l.split('='))
        ah = extract(l, pattern)
        for ih in ah:
            indicies = findIndex(chains, ih[0], ih[1])
            inchain = indicies[0]
            out_fh.write('{:5d}{:5d}{:7d} {:>3s}{:8d}{:7d} '.format(
                ih[0], ih[1], len(indicies[1]), 'ABCD'[inchain[0]], inchain[1],
                inchain[2]
                ))
            fmt = '{:5d}' * len(indicies[1]) + '\n'
            out_fh.write(fmt.format(*(indicies[1])))

    fh.close()
    print('Extracting secondary structural fragment have been done successfully!')

