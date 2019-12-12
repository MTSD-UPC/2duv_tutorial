'''
@Description: This script can cal H
@Author: Qian Zhang 
@Date: 2019-09-09 11:25:57
@LastEditTime: 2019-12-11 20:39:04
@LastEditors: Please set LastEditors
'''

'''
本脚本用于计算片段的2DUV光谱，使用服务器的bench队列，
两个参数分别是需要计算的始末snapshots序号。
队列可以通过修改 4_sub_calspectra.pbs 文件进行修改
'''
import os
import sys
import numpy as np 
import shutil


start = sys.argv[1]
end = sys.argv[2]

try:
    stru_path = np.load('stru_path.npy')
except: 
    #step1. get stru_path
    path_now = os.getcwd()
    dir_path = path_now+'/snapshots/'
    stru_path = []
    list_dir = os.listdir(dir_path)
    for dir1 in list_dir:
        dir_path1 = dir_path + dir1
        for dir2 in os.listdir(dir_path1):
            dir_path2 = dir_path1 + '/' + dir2
            stru_path.append(dir_path2+'/helices')
    np.save('stru_path.npy',stru_path)

path_now = os.getcwd()
select_path = np.array(stru_path)[int(start):int(end)+1]
for path in select_path:
    os.chdir(path_now)
    shutil.copyfile('4_sub_calspectra.pbs',path+'/sub_calspectra.pbs')
    shutil.copyfile('4_cal_spectra.py',path+'/cal_spectra.py')
    os.chdir(path)
    os.system('qsub sub_calspectra.pbs')

print('All works have been submitted!')
        
