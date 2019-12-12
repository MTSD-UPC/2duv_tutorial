#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2018 Hao Ren <renh.cn@gmail.com>
# Modified by Qian Zhang 2019 <zhangqian.allen@gmail.com>
# 
# Distributed under terms of the LGPLv3 license.

"""
Read the alpha helix information in snapshot_dir/helices/helix.dat
create a subdirectory for each helix segment in the helices directory,
and write the following info:
    1. segment sequence and positions
    2. Edipl.dat and Mdipl.dat
    3. coord.dat
    4. Hamil.dat read from Hamil_ext.dat of the whole system

"""
import os
import numpy as np
import json
import tarfile
import argparse
import sys

parser = argparse.ArgumentParser(description='manul to this script')
parser.add_argument('-s','--start',help='start number of snapshot',type=int)
parser.add_argument('-e','--end',help='end number of snapshot',type=int)
parser.add_argument('-d','--dir-name',help='dir name of secondary structure',type=str)
args = parser.parse_args()

start = args.start
end = args.end
stru_dir = args.dir_name

chains = np.load('chains.npy').tolist()
seq = ''.join(chains)
#print([len(x) for x in chains])
#print(len(seq))

def parse_seg(seg, i_seg, snap_dir):
    seg = seg.split()
    seg_info = {}
    start, end = seg[:2]
    seg_info['start'] = int(start)
    seg_info['end'] = int(end)
    
    nmodes = int(seg[2])
    seg_info['num_modes'] = nmodes

    seg_info['chain'] = seg[3]
    start_in_chain, end_in_chain = [int(x) for x in seg[4:6]]
    if seg_info.get('chain', 'X') in 'BD':
        start_in_chain += 1
        end_in_chain += 1
    seg_info['start_in_chain'] = start_in_chain
    seg_info['end_in_chain'] = end_in_chain
    modes = [int(x) for x in seg[6:]]
    assert(len(modes) == nmodes)
    seg_info['modes'] = modes

    seg_info['i_seg'] = i_seg
    seg_info['snap_dir'] = snap_dir
    seg_dir = '{}/{}/seg{:02d}'.format(snap_dir,stru_dir,i_seg)
    seg_info['seg_dir'] = seg_dir
    return seg_info


def create_coord(seg_info, protein):
    # snap_dir = seg_info.get('snap_dir')
    # i_seg = seg_info.get('i_seg')
    

    chain = seg_info.get('chain')
    start_in_chain = seg_info.get('start_in_chain')
    end_in_chain = seg_info.get('end_in_chain')
    seq_range = range(start_in_chain, end_in_chain + 1)

    labels = []
    coords = []
    for l in protein:
        i_res = int(l[23:26]) 
        i_chain = l[21]
        if i_chain != chain:
            continue
        if i_res in seq_range:
            labels.append(l[13])
            coords.append([float(x) for x in l[29:54].split()])
    with open('{}/coord.dat'.format(seg_info.get('seg_dir')), 'w') as fh:
        fh.write('#start\n')
        for c in coords:
            fh.write('{:16.8f} {:16.8f} {:16.8f}\n'.format(*c))
        fh.write('#end')
    return

def create_dipole(seg_info, dipoles, label):
    modes = np.array(seg_info.get('modes')) - 1
    dip = dipoles[modes]
    with open('{}/{}dipl.dat'.format(
        seg_info.get('seg_dir'), label
        ), 'w') as fh:
        fh.write('#start\n')
        for d in dip:
            fh.write('{:16.8f} {:16.8f} {:16.8f}\n'.format(*d))
        fh.write('#end')
    return

def create_Hamil(seg_info, Hamil):
    modes = np.array(seg_info.get('modes')) - 1
    N = len(modes)
    H_sub = Hamil[np.ix_(modes, modes)]
    with open('{}/Hamil.dat'.format(seg_info.get('seg_dir')), 'w') as fh:
        fh.write('#start\n')
        for i in range(N):
            fmt = '{:16.8f}' + ' {:17.8f}' * i + '\n'
            fh.write(fmt.format(*(H_sub[i,:(i+1)])))
        fh.write('#end')
    return


if __name__ == '__main__':
    for i_snap in range(int(start),int(end)+1):
        print('working on snapshot #', i_snap)
        batch_dir = 'snapshots/{:05d}'.format(i_snap // 10 * 10 + 10)
        snap_dir = '{}/{:05d}'.format(batch_dir, i_snap)
        helices_dir = '{}/{}'.format(snap_dir,stru_dir)
        helices_dat = '{}/helix.dat'.format(helices_dir)

        with open('{}/job.pdb'.format(snap_dir), 'r') as fh:
            protein = fh.readlines()
        edip = np.loadtxt('{}/Edipl.dat'.format(snap_dir))
        mdip = np.loadtxt('{}/Mdipl.dat'.format(snap_dir))
        nmode = int(np.loadtxt('{}/numberofmode.txt'.format(snap_dir),dtype=object)[1])
        Htar = tarfile.open('{}/H.tgz'.format(snap_dir))
        Htar.extract('Hamil_ext.dat')
        Hamil = np.zeros([nmode, nmode])
        with open('Hamil_ext.dat', 'r') as fh:
            lines = fh.readlines()[1:-1]
        os.remove('Hamil_ext.dat')
        for i in range(nmode):
            Hamil[i,:(i+1)] = [float(x) for x in lines[i].split()]
        

        with open(helices_dat, 'r') as fh:
            helices = fh.readlines()[1:]
        for i_seg, seg in enumerate(helices):
            seg_info = parse_seg(seg, i_seg, snap_dir)
            if not os.path.exists(seg_info['seg_dir']):
                os.mkdir(seg_info['seg_dir'])
            with open('{}/seg.json'.format(seg_info.get('seg_dir')), 'w') as fh:
                json.dump(seg_info, fh)
            #print(seg_info)
            with open('{}/numberofmodes.dat'.format(seg_info['seg_dir']), 'w') as fh:
                fh.write('{}\n'.format(seg_info.get('num_modes')))
            create_coord(seg_info, protein)
            create_dipole(seg_info, edip, 'E')
            create_dipole(seg_info, mdip, 'M')
            create_Hamil(seg_info, Hamil)
            
