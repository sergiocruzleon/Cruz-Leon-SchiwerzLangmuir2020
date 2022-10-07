DIR='/home/tb/secruz/DATA/PhD/P1_all_atom/step2/rates-micro-macro'
REP_T_DIV='/home/tb/secruz/DATA/PhD/P1_all_atom/step2/rates-micro-macro/bash-files/index-temp/replica_temp_divalent.dat'
REP_T_MONO='/home/tb/secruz/DATA/PhD/P1_all_atom/step2/rates-micro-macro/bash-files/index-temp/replica_temp_monovalent.dat'
RESULTS_ALL='/home/tb/secruz/DATA/PhD/P1_all_atom/step2/rates-micro-macro/RESULTS_ALL'

index=0
#r_dagger=[LIO, NIO, KIO CIO]
MN7_r_dagger=(0.28 0.35 0.35 0.4)
MN7_r_dagger2=(0.35 0.45 0.5 0.5)

MOP_r_dagger=(0.25 0.30 0.40 0.40)
MOP_r_dagger2=(0.40 0.40 0.50 0.50)

MO6_r_dagger=(0.30 0.30 0.4 0.40) #MO6_r_dagger=(0.29 0.30 0.4 0.45)
MO6_r_dagger2=(0.40 0.5 0.5 0.55) #MO6_r_dagger2=(0.35 0.42 0.5 0.55)

max_num=(25000 25000 25000  25000)

for ion in 'LIO' 'NIO' 'KIO' 'CIO';
do
	cd ${DIR}
	cd ${ion}/rates-all/
	pwd
	#ORDER PARAMETER
	mkdir order-parameter
	cp ${DIR}/bash-files/mono_order_parameter_N7A.sh . 
	cp ${DIR}/bash-files/mono_order_parameter_O6.sh .
	cp ${DIR}/bash-files/mono_order_parameter_OP.sh .
	cp ${DIR}/bash-files/averaged-order-parameter.py .
	cp ${DIR}/bash-files/averaged-order-parameter2.py .
	# Calculate order parameters	
	echo '####################################################'
	echo '######## Order parameter for '${ion} 'OP ###############'
	echo '####################################################'

	echo '#### Min dist: ' ${MOP_r_dagger[${index}]}
	echo '#### Max dist: '  ${MOP_r_dagger2[${index}]}

	sh mono_order_parameter_OP.sh ${ion} ${ion}/dist order-parameter ${MOP_r_dagger[${index}]} ${MOP_r_dagger2[${index}]}   #dist_bound, #dist_unbound

	echo '####################################################'
	echo '######## Order parameter for '${ion} ' N7A ###############'
	echo '####################################################'

	echo '#### Min dist: ' ${MN7_r_dagger[${index}]}
	echo '#### Max dist: '  ${MN7_r_dagger2[${index}]}

	sh mono_order_parameter_N7A.sh ${ion} ${ion}/dist order-parameter ${MN7_r_dagger[${index}]} ${MN7_r_dagger2[${index}]}  #dist_bound, #dist_unbound

	echo '####################################################'
	echo '######## Order parameter for '${ion} ' O6G ###############'
	echo '####################################################'


	echo '#### Min dist: ' ${MO6_r_dagger[${index}]}
	echo '#### Max dist: '  ${MO6_r_dagger2[${index}]}

	sh mono_order_parameter_O6.sh ${ion} ${ion}/dist order-parameter ${MO6_r_dagger[${index}]} ${MO6_r_dagger2[${index}]}   #dist_bound, #dist_unbound
	
	
	#KINETICS

	mkdir kinetics
	for atom in 'OP' 'N7A' 'O6G'; do	
	cd ${DIR}/${ion}/rates-all/kinetics/
		mkdir ${atom}/
		cp ${DIR}/bash-files/analyze-kin-General-mono.py  ${atom}/ 
	  	cp ${DIR}/bash-files/eval-start.py ${atom}/
	  	cp ${DIR}/bash-files/eval-lifetime.py ${atom}/
		cd ${atom}/
		mkdir RESULTS
		echo '####################################################'
		echo '######## Kinetics of '${ion} ${atom}'#################'
		echo '####################################################'
		pwd
      	        python analyze-kin-General-mono.py  ${ion} ${atom} ${DIR}/${ion}/rates-all/order-parameter  ${REP_T_MONO} ${max_num[${index}]}
   	        python eval-start.py ${ion}
                python eval-lifetime.py ${ion}

		# Write results 
		mkdir ${RESULTS_ALL}/${ion}
		mkdir ${RESULTS_ALL}/${ion}/${atom}
		cp RESULTS/* ${RESULTS_ALL}/${ion}/${atom}/
	done
	let "index= index+1"
	 
done


#############################################################################################################
####### DIVALENT IONS
#############################################################################################################

#r_dagger=[Ca, Sr, Ba]
#N7_r_dagger=(0.28248405 0.30048465 0.326974)
#N7_r_dagger2=(0.3996819 0.4118152 0.4395652)
N7_r_dagger=(0.370 0.350 0.35)
N7_r_dagger2=(0.42 0.46 0.50)

#OP_r_dagger=(0.297 0.298 0.328)
#OP_r_dagger2=(0.36 0.38 0.38)
OP_r_dagger=(0.25 0.30 0.30)
OP_r_dagger2=(0.40 0.55 0.5)

#O6_r_dagger=(0.2375002 0.27829285 0.2832782 0.30178365)
#O6_r_dagger2=(0.344744 0.39695535 0.3922331 0.41143885)
O6_r_dagger=(0.30 0.38 0.37)   # O6_r_dagger=(0.30 0.38 0.37)
O6_r_dagger2=(0.45 0.55 0.50) # O6_r_dagger2=(0.45 0.55 0.50)

let "index=0"


max_num=(250000 250000 25000)
for ion in 'Ca2' 'Sr2' 'Ba2'; 
#for ion in 'Ca2'; 
do
	cd ${DIR}
	cd ${ion}/rates-all/
	pwd
	
	#ORDER PARAMETER
	mkdir order-parameter
	cp ${DIR}/bash-files/div_order_parameter_N7A.sh . 
	cp ${DIR}/bash-files/div_order_parameter_O6.sh .
	cp ${DIR}/bash-files/div_order_parameter_OP.sh .
	cp ${DIR}/bash-files/averaged-order-parameter.py .
	cp ${DIR}/bash-files/averaged-order-parameter2.py .
	# Calculate order parameters	
	echo '####################################################'
	echo '######## Order parameter for '${ion} 'OP ###############'
	echo '####################################################'
	echo '#### Min dist: ' ${OP_r_dagger[${index}]}
	echo '#### Max dist: '  ${OP_r_dagger2[${index}]}
	#sh div_order_parameter_OP.sh ${ion} ${ion}/dist order-parameter 0.35 1.0  #dist_bound, #dist_unbound
	sh div_order_parameter_OP.sh ${ion} ${ion}/dist order-parameter ${OP_r_dagger[${index}]} ${OP_r_dagger2[${index}]}  #dist_bound, #dist_unbound 
	echo '####################################################'
	echo '######## Order parameter for '${ion} ' N7A ###############'
	echo '####################################################'
	echo '#### Min dist: ' ${N7_r_dagger[${index}]}
	echo '#### Max dist: '  ${N7_r_dagger2[${index}]}
	sh div_order_parameter_N7A.sh ${ion} ${ion}/dist order-parameter  ${N7_r_dagger[${index}]} ${N7_r_dagger2[${index}]}   #dist_bound, #dist_unbound
	echo '####################################################'
	echo '######## Order parameter for '${ion} ' O6G ###############'
	echo '####################################################'
	echo '#### Min dist: ' ${O6_r_dagger[${index}]}
	echo '#### Max dist: '  ${O6_r_dagger2[${index}]}
	sh div_order_parameter_O6.sh ${ion} ${ion}/dist order-parameter  ${O6_r_dagger[${index}]} ${O6_r_dagger2[${index}]}   #dist_bound, #dist_unbound
		
	#KINETICS

	mkdir kinetics
	for atom in 'OP' 'N7A' 'O6G'; do	
	cd ${DIR}/${ion}/rates-all/kinetics/
		mkdir ${atom}/
		cp ${DIR}/bash-files/analyze-kin-General-div.py  ${atom}/ 
	  	cp ${DIR}/bash-files/eval-start.py ${atom}/
	  	cp ${DIR}/bash-files/eval-lifetime.py ${atom}/
		cd ${atom}/
		mkdir RESULTS
		echo '####################################################'
		echo '######## Kinetics of '${ion} ${atom}'#################'
		echo '####################################################'
		pwd
      	        python analyze-kin-General-div.py  ${ion} ${atom} ${DIR}/${ion}/rates-all/order-parameter  ${REP_T_DIV} ${max_num[${index}]}
   	        python eval-start.py ${ion}
                python eval-lifetime.py ${ion}

		# Write results 
		mkdir ${RESULTS_ALL}/${ion}
		mkdir ${RESULTS_ALL}/${ion}/${atom}
		cp RESULTS/* ${RESULTS_ALL}/${ion}/${atom}/
	done
	let "index= index+1"
done

