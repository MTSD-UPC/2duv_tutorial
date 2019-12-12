## 2DUV计算简明教程

本文档不讲原理，只记录计算流程，快速计算出大量二维紫外光谱。

这里假设我们需要计算hemoglobin蛋白的10个snapshot，然后提取a-helix片段计算2DUV光谱。

**注：本文档所使用的绝大多数脚本是基于任老师的脚本修改**


### 确认文件

`2DUV/`目录下为计算2DUV所需脚本文件，`GromacsFile/`目录下为跑MD模拟所需要的文件，将所有脚本下载：

```bash
git clone https://github.com/MTSD-UPC/2duv_tutorial.git
```

### 确认环境

- Gromacs

- EHEF
- SPECTRON 
- Python 3.X

在你的home ~/目录下尽量添加以下内容到.bashrc文件（因为避免后面哪个脚本忘记写声明了...）：

``` bash
export PATH=/opt/openmpi/1.6.5-intel2011/bin:$PATH 
source /opt/gromacs/4.6.7/bin/GMXRC
export EHEF='/home/ren/EHEF'  #改成自己的目录下的也可以
```

确认SPECTRON软件的安装位置：`/opt/spectron/2.7/spectron2`，不是的话需要改脚本``4_cal_spectra.py``

执行 `source .bashrc`

### MD跑轨迹

这里以Gromacs为例，我的蛋白的PDB ID为[1hda](https://www.rcsb.org/structure/1hda)。

首先打开1hda.pdb 删除里面的所有以HEATEM开头的行数，然后保存文件。

将pdb文件复制到`GromacsFile/`目录下。

```bash

pdb2gmx -f 1hda.pdb -o 1hda.gro -water spc 
8 #选择charmm27

editconf -f 1hda.gro -o new_box.gro -c -d 1.0 -bt cubic

genbox -cp new_box.gro -cs spc216.gro -o solv.gro -p topol.top

grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr

genion -s ions.tpr -o solv_ions.gro -pname NA -nname CL -neutral

```

到目前为止，你会看到输出提示告诉你`Will try to add 0 NA ions and 2 CL ions  `,意思是需要添加2个CL，0个NA，在输入`12`,选择Water那个选项。

这时文件夹里有很多文件，找到`topo.top`文件，打开拖到最后进行修改，将SOL 的个数减2，增加CL的个数为2，NA因为是0就不需要添加了。在本例中，我的top文件修改如下：

> ```ba
> [ molecules ]
> ; Compound        #mols
> Protein_chain_A     1
> Protein_chain_B     1
> Protein_chain_C     1
> Protein_chain_D     1
> SOL                21717          #之前是21719
> CL                  2           #这一行是新加的，之前没有的。
> ```

接着输命令：

```bash
 grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr
 
 mdrun -v -deffnm em   #这一步有点慢
 
 grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
 
```

跑nvt 撰写脚本：

```shell
#!/bin/bash
#PBS -N n-ra
##PBS -M renh@upc.edu.cn
##PBS -m abe
##PBS -o /path/to/stdout
##PBS -e /path/to/stderr

#PBS -q bench
#PBS -l nodes=1:ppn=20
#PBS -l walltime=9:00:00

echo -n 'We work on:'
cat $PBS_NODEFILE
cd $PBS_O_WORKDIR

source /opt/gromacs/4.6.7/bin/GMXRC

mdrun -nt 20 -pin on -deffnm nvt
```

提交任务跑nvt。

修改`npt.mdp`文件的相关选项如下：

> ```bash
> nsteps		= 5000		; 2 * 5000 = 10 ps
> dt		    = 0.002		; 2 fs
> ; Output control
> nstxout		= 500		; save coordinates every 1.0 ps
> nstvout		= 500		; save velocities every 1.0 ps
> nstenergy	= 500		; save energies every 1.0 ps
> nstlog		= 500		; update log file every 1.0 ps
> ```

解释： `nsteps`总步长，共500步，每一步2fs，共10ps，每隔`nstxout=500`步取一个snapshot，即每隔1ps取一个snapshot，而我们想要10个snapshot的话，就需要10ps。**注**：一般默认每隔1000fs就认为两个结构完全不同，所以一般都会设置每隔1000fs取一个snapshot。

跑完nvt之后，接着输命令合成npt文件：

```bash
grompp -f npt.mdp -c nvt.gro -p topol.top -o npt.tpr
```

接着写脚本跑npt:(和nvt相比就把nvt替换掉npt就行)

```shell
#!/bin/bash
#PBS -N n-ra
##PBS -M renh@upc.edu.cn
##PBS -m abe
##PBS -o /path/to/stdout
##PBS -e /path/to/stderr

#PBS -q bench
#PBS -l nodes=1:ppn=20
#PBS -l walltime=9:00:00

echo -n 'We work on:'
cat $PBS_NODEFILE
cd $PBS_O_WORKDIR

source /opt/gromacs/4.6.7/bin/GMXRC

mdrun -nt 20 -pin on -deffnm npt
```

等着跑完接着输命令：

```bash

trjconv -s npt.tpr -f npt.trr -pbc whole -o trj_pre.pdb
```

**注意** 这里会提示输出的类型，同时需要记录几个数据，很重要。

1. Protein 总共有**8788**个原子        
2. Water 有 **65151**个原子
3. 离子 有 **2** 个 （正负离子相加）

```bash
0  #System 都输出来               
```

至此，trj_pre.pdb 就是你输出来的包含10个snapshot的轨迹文件了！（**注意**： 其实我们实际输出了11个snapshot，因为默认输出了初始构型（0s 的时候的构型），这个无所谓）

### 利用do_dssp生成二级结构文件

```bash

do_dssp -f npt.trr -s npt.tpr -o ss.xpm -ssdump ss.dat
```

会生成两个文件`ss.xpm`与`ss.dat`，打开`ss.xpm`会显示每个二级结构的标识符含义，比如：

> ```bash
> "~  c #FFFFFF " /* "Coil" */,
> "S  c #008000 " /* "Bend" */,
> "T  c #FFFF00 " /* "Turn" */,
> "H  c #0000FF " /* "A-Helix" */,
> "I  c #800080 " /* "5-Helix" */,
> "G  c #808080 " /* "3-Helix" */,
> ```

`ss.dat`是我们需要用于接下来计算的内容。

### 计算2DUV光谱

本过程将会用到大量的脚本，将之前的pdb数据库源文件`1hda.pdb`、生成的轨迹文件`trj_pre.pdb`、二级结构文件`ss.dat`拷贝到`2DUV/`目录下。

#### 拆分轨迹文件

脚本名称：`0_parse_traj.py`

脚本使用：可以使用`python 0_parse_traj.py -h` 查看使用说明，以本例为例：

```bash

python 0_parse_traj.py --i trj_pre.pdb --b 10 --protein 8788 --water 65151 --ions 2
```

这样会在本目录下生成`snapshots/`文件夹，并且把每个snapshot都拆分成`job.pdb`、`job.sol`、`job.xyz`三个文件，等计算完毕，进入下一步。

#### 计算Hamilton矩阵

脚本名称：`1_cal_Hamil.py`、`1_bench_genH.sh`

使用方法：`python 1_cal_Hamil.py [start] [end]` 【start】、【end】分别代表想要计算的开始与最后一个的snapshot序号。

**注意：以上脚本必须放在同一文件夹下，不要改名字**

*可以根据需要更改`1_bench_genH.sh`文件中的队列信息，但不要更改使用CPU的核数，每个任务单核即可*

```bash

python 1_cal_Hamil.py 0 9  #这行命令代表会提交序号为0~9的snapshot的计算任务共10个任务给节点计算
```

按目前chimy的bench节点的算力来看，计算一个snapshot需要20分钟左右，在空闲的情况下，可以同时计算20个snapshot。

#### 提取二级结构片段

脚本名称：`2_extract_struc.py`

脚本使用：`python 2_extract_srtuc.py -h` 可查看参数帮助，以本次为例：

```bash

python 2_extract_struc.py -f 1hda.pdb -c 4 -b 10 -p HHHHHHHHHH* -s ss.dat

```

> ```bash
> #-f pdb文件名
> #-c 蛋白的链的个数
> #-b 总snapshot数目
> #-p 二级结构的pattern,HHHHHHHHHH* 代表查找10个残基以上的helix片段。
> #-s do_dssp 生成的二级结构文件
> ```

#### 提取片段的Hamilton

**注意本小节需要计算完Hamilton矩阵才可以继续，需确认待提取的snapshot已经计算完毕Hamilton矩阵（确认snapshot文件夹内生成`H.tgz`等系列文件）**。

脚本名称：`3_extract_Hamil.py`

脚本使用：`python 3_extract_Hamil.py [start] [end]` ，[start]、[end] 参数代表想要提取的snpshot的起始与终止序号，以本例为例：

```bash

python 3_extract_Hamil.py 0 9
#为了节约时间，本脚本可以指定计算snapshots的范围，
#在需要计算很多snapshot时，往往Hamilton矩阵计算较慢，需要排队。
#因此可以先对已经计算完Hamilton矩阵的snapshot应用本脚本。
```

#### 计算2DUV光谱

脚本名称：`4_run_calspectra.py`、`4_sub_calspectra.pbs`、`4_cal_spectra.py`

脚本使用：`4_sub_calspectra.pbs`中可设置队列名称，默认的是**bench**队列（short队列经常报错不知道为啥）

​					`4_cal_spectra.py` 中设置inputs文件夹的**绝对路径**，与spectron软件的默认**绝对路径**，默认如下：

> #4_cal_spectra.py
>
> ```python
> ...
> 
> input_dir = '/home2/qzhang/test/inputs'    # about 20 lines
> 
> ...
> # about 45 lines
> ​	os.system('/opt/spectron/2.7/spectron2 -i input_CDLA.inp > cdla.log')
> 
> ​	os.system('/opt/spectron/2.7/spectron2 -i input_KI_FUV_xxxx.inp > xxxx.log')
> 
> ​	os.system('/opt/spectron/2.7/spectron2 -i input_KI_FUV_xxxy.inp > xxxy.log')
> ```

​				`python 4_run_calspectron.py [start] [end]` 【start】【end】参数指定需要计算snapshot的开始与结束序号，以本次为例，计算0~9 snapshot中的所有片段：

```bash

python 4_run_calspectra.py 0 9
```

提交完所有snapshot的作业即可，这个过程计算较慢，一般如果一个snapshot有20个片段需要计算一个小时。

