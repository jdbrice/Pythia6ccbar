#!/bin/bash

cd $SHARED_SCRATCH/jdb12/pythia_histo/

echo "JOBID $1"
source /projects/geurts/jb31/ROOT/bin/thisroot.sh

# export ROOTSYS=/projects/geurts/jb31/ROOT
# export PATH=$ROOTSYS/bin:$PATH
# export LD_LIBRARY_PATH=$ROOTSYS/lib/root
export LD_LIBRARY_PATH=/projects/geurts/jb31/ROOTbuild/pythia6:$LD_LIBRARY_PATH


trandom="$(od -vAn -N4 -tu4 < /dev/urandom | tr -d '[:space:]')"
echo "seed = $trandom"
/home/jdb12/work/Pythia6HistoMaker/GENERATOR 4 100000 $trandom >& log_$trandom.log
