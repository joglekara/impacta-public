#!/bin/sh
#PBS -S /bin/sh
#PBS -N gnt
#PBS -l procs=16,walltime=100:00:00
#PBS -l pmem=3800mb
#PBS -A agrt1_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -M archisj@umich.edu
#PBS -m abe
#PBS -j oe
#PBS -V
echo "I ran on:"
cat $PBS_NODEFILE

#############################################################################################
# EDIT BELOW HERE 
#############################################################################################

export ROOT_DIR=/scratch/agrt_flux/archisj/devilray/impray/
export DATA_DIR=/scratch/agrt_flux/archisj/th_out/
export INPUTFILE=imstdin

## OPTIONAL - CAN BE BLANK, no spaces
export RUNTITLE=wtgen

#############################################################################################

export DIR_NAME=impacta_${RUNTITLE}_${PBS_JOBID%.nyx.engin.umich.edu}

mkdir ${DATA_DIR}/${DIR_NAME}
cp -f ${ROOT_DIR}/bin/impacta ${DATA_DIR}/${DIR_NAME}/.
cp -f ${ROOT_DIR}/input/${INPUTFILE} ${DATA_DIR}/${DIR_NAME}/imstdin
#cp -fr ${ROOT_DIR}/bin/impacta_data_in ${DATA_DIR}/${DIR_NAME}/.

# Write down the version
cp ${ROOT_DIR}/version.txt ${DATA_DIR}/${DIR_NAME}/version.txt

# Create root directory
cd ${DATA_DIR}/${DIR_NAME}

# Run code
mpirun -np 16 ./impacta -read imstdin
#-renorm_temp 0.01
