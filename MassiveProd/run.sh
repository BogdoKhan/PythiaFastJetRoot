#!/bin/bash

pythiaTune=14
nasada=0

pthatMin=("4"  "11" "21" "36" "56" "84"  "117" "156" "200" "249")  #hard bin
pthatMax=("11" "21" "36" "56" "84" "117" "156" "200" "249" "1000")   #hard bin
TTLow=("6" "12" "20")      #TT bin
TTHigh=("7" "20" "30")     #TT bin 
jetR=("0.4")     #jet R


#pthatMin=("0")  #hard bin
#pthatMax=("0")   #hard bin
#TTLow=("20" )      #TT bin
#TTHigh=("30")     #TT bin 
#jetR=("0.4")     #jet R


############################################
nhb=${#pthatMin[@]}
nTT=${#TTLow[@]}
njetR=${#jetR[@]}

ipp=0
#for hbin in `seq 1 40`; do
for hbin in `seq 1 150`; do
   for (( ihb=0; ihb<${nhb}; ihb++ )); do
      for (( itt=0; itt<${nTT}; itt++ )); do
         for (( ijr=0; ijr<${njetR}; ijr++ )); do
             rndNum=$((nasada+ipp))
             echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
             echo "./pygen.exe $pythiaTune $rndNum ${TTLow[$itt]} ${TTHigh[$itt]} ${pthatMin[$ihb]} ${pthatMax[$ihb]} ${jetR[$ijr]}"
             echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
             ./pygen.exe $pythiaTune $rndNum ${TTLow[$itt]} ${TTHigh[$itt]} ${pthatMin[$ihb]} ${pthatMax[$ihb]} ${jetR[$ijr]}
             ipp=$((ipp+1)) 
          done
       done
   done
done
