#!/bin/bash

#### Important note! this code only works with the  (module load anaconda/2_4.3.0 )!! Lukas scripts depent on it!
#Calculate all the rates for the ions. Order parameter is already calculated and stores in folder /home/tb/secruz/DATA/PhD/P1_all_atom/step2/nadine-data/TREMD_RNA-TIP3P/analysis/python/averaged-2/${ion}
#mkdir averaged2
cd averaged2/
#ion in LIO NIO KIO CIO 
ion=$1
#ion OP O6G N7G
atom=$2
echo '####################################################'
echo '######## Calculating'${ion} '#######################'
echo '####################################################'
  mkdir ${ion}/${atom}
  rm ${ion}/${atom}/*
  cp ../analyze-kin-General-mono.py  ${ion}/${atom}/ 
  cp ../eval-start.py ${ion}/${atom}/
  cp ../eval-lifetime.py ${ion}/${atom}/

  cd ${ion}/${atom}/
  mkdir RESULTS
  python analyze-kin-General.py  Ba2 OP (O6G or N7A) '/home/tb/secruz/DATA/PhD/P1_all_atom/step2/rates-micro-macro/Ba2/rates-all/order-parameter/'  '/home/tb/secruz/DATA/PhD/P1_all_atom/step2/rates-micro-macro/bash-files/index-temp/replica_temp_divalent.dat'
  python eval-start1.py ${ion}
  python eval-lifetime1.py ${ion}

  cd ../../

done
cd ../
echo 'FINISH!'	
