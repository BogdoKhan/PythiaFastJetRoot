#!/bin/sh
source /home/krizek/setenv.sh
cd /nethome/krizek/PYTHIA8/pythia8183/PythiaFastJetIntro
./compile_pythia.sh pysec

#run analysis
./pysec.exe 5 0.4 
