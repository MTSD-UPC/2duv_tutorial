'''
@Author: Qian Zhang
@Date: 2019-12-10 22:07:07
@LastEditTime: 2019-12-11 19:50:17
@LastEditors: Please set LastEditors
@Description: Cal. Hamilton Matrix
'''
'''
此文件放在蛋白目录下，要求当前目录中有文件夹snapshots/
'''

import numpy as np 
import os 
import sys 
import shutil

start = sys.argv[1]  #传入第一个参数，开始snapshot序号
end = sys.argv[2]    #传入的第二个参数，结束的snapshot序号

try:
    path_list = np.load('path_list.npy')
except:
    # 生成path——list：生成所有的snapshot的路径文档
    path_now = os.getcwd()
    dir_path = path_now + '/'+'snapshots/'
    path_list = []
    list_dir = os.listdir(dir_path)
    for dir1 in list_dir:
        dir_path1 = dir_path+dir1
        for dir2 in os.listdir(dir_path1):
            dir_path2 = dir_path1+'/'+dir2
            path_list.append(dir_path2)

    path_list.sort()
    np.save('path_list.npy',path_list)

select_path = np.array(path_list)[int(start):int(end)+1]
genH_path = os.getcwd()+'/1_bench_genH.sh'
for path in select_path:
    shutil.copyfile(genH_path,path+'/bench_genH.sh')
    os.chdir(path)
    os.system('qsub bench_genH.sh')
    print('submit #{}'.format(path[-5:]))

print('All works have been submitted!')
