
#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2018 Hao Ren <renh.cn@gmail.com>
# Modified by Qian Zhang 2019 <zhangian.allen@gmail.com>
# Distributed under terms of the LGPLv3 license.

"""
Please confirm that input_dir and Spectron dir
"""
import os
import numpy as np

def get_lines(fname):
    with open(fname, 'r') as fh:
        lines = fh.readlines()
    return lines

input_dir = '/home2/qzhang/test/inputs'
cdla_template = get_lines('{}/input_CDLA.inp'.format(input_dir))
xxxx_template = get_lines('{}/input_KI_FUV_xxxx.inp'.format(input_dir))
xxxy_template = get_lines('{}/input_KI_FUV_xxxy.inp'.format(input_dir))
#zzzz_template = get_lines('{}/input_KI_FUV_zzzz.inp'.format(input_dir))
KI_template = {'xxxx':xxxx_template, 'xxxy':xxxy_template}
        #'zzzz':zzzz_template}

segs = [x for x in os.listdir('.') if x.startswith('seg')]
for seg in segs:
    os.chdir(seg)
    #print(seg)
    with open('numberofmodes.dat') as fh:
        nmodes = int(fh.readline())

    cdla_lines = [l for l in cdla_template]
    cdla_lines[17] = cdla_lines[17].replace('NMODS', '{}'.format(nmodes))
    with open('input_CDLA.inp', 'w') as fh:
        fh.writelines(cdla_lines)

    for conf in ['xxxx', 'xxxy']:
        lines = [l for l in KI_template.get(conf)]
        lines[19] = lines[19].replace('NMODS', '{}'.format(nmodes))
        with open('input_KI_FUV_{}.inp'.format(conf), 'w') as fh:
            fh.writelines(lines)
    os.system('/opt/spectron/2.7/spectron2 -i input_CDLA.inp > cdla.log')
    os.system('/opt/spectron/2.7/spectron2 -i input_KI_FUV_xxxx.inp > xxxx.log')
    os.system('/opt/spectron/2.7/spectron2 -i input_KI_FUV_xxxy.inp > xxxy.log')
    #os.system('/opt/spectron/2.7/spectron2 -i input_KI_FUV_zzzz.inp > zzzz.log')
    cd = np.loadtxt('sig-CD.dat')[:,1]
    la = np.loadtxt('sig-LA.dat')[:,1]
    np.savez_compressed('LACD.npz', la=la, cd=cd)
    os.remove('sig-CD.dat')
    os.remove('sig-LA.dat')

    xxxx = np.loadtxt('FUV-KI-xxxx.dat')[:,4:]
    xxxy = np.loadtxt('FUV-KI-xxxy.dat')[:,4:]
    np.savez_compressed('2dfuv.npz', xxxx=xxxx, xxxy=xxxy)
    os.remove('FUV-KI-xxxx.dat')
    os.remove('FUV-KI-xxxy.dat')





    os.remove('spec_cal.sm')
    os.chdir('..')




