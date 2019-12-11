#!/bin/bash
#PBS -N 2DUV-EHEF-Cal
##PBS -M renh@upc.edu.cn
##PBS -m abe
##PBS -o /path/to/stdout
##PBS -e /path/to/stderr

#PBS -q bench 
#PBS -l nodes=1:ppn=1
#PBS -l walltime=9:00:00

echo -n 'We work on:'
cat $PBS_NODEFILE
cd $PBS_O_WORKDIR




# genH.sh
# Copyright (C) 2018 Hao Ren <renh.cn@gmail.com>
#
# Distributed under terms of the LGPLv3 license.
#


EHEF=/home/ren/EHEF

# copy chromophores.dat
cp $EHEF/dichrocalc/params/chromophores.dat .

# generate dcinput
$EHEF/dcinput/dcinput -sc job.pdb &> s3.log

# Extract system info
echo job.pdb | $EHEF/fluct/extract_peptide_RCSBpdb_July2009 &> s4.log

# convert
cp job.sol o0axyz.dat
$EHEF/fluct/readbaseinfo-new &> s5.log
mv ttaxyz.dat ttaxyz.sol

cp job.xyz o0axyz.dat
$EHEF/fluct/readbaseinfo-new &> s52.log

# map atoms with database
echo job.inp | $EHEF/fluct/FindAmino_FromDichInp_mapsequence_2010Jan-ren &> s6.log

# copy
cp $EHEF/amino_acid/pep/NMA*.TEP .
cp $EHEF/amino_acid/aromatic/*.TEP .

# calc coupling
echo job.inp | $EHEF/fluct/xyz2Amino_FromDichInp_Pep_aromatic_coup_2010Jan-ren &> s8.log

# gen
$EHEF/dichrocalc/dichrocalc -i job.inp -p $EHEF/dichrocalc/params/ -v -d 5 &> s9.log

# rename
mv pep_input_PtTrEsp_Update_Hamil.dat job-pep_input_PtTrEsp_Update_Hamil.dat
mv ami_input_PtTrEsp_Update_Hamil.dat job-ami_input_PtTrEsp_Update_Hamil.dat
echo job.inp | $EHEF/fluct/Read_DichOutput_TrEspcoupAll_2_SpectronInput_Aromatic_2010Jan &> s10.log

# compress
tar zcf H.tgz Hamil*


rm -f *.TEP job-* chromophores.dat o0a* tta* Hamil* job.dbj *.log

filename=$PBS_O_WORKDIR
echo $PBS_O_WORKDIR>/home2/qzhang/2duv_max/hemoglobin/Donefile/${filename:0-5}.txt

