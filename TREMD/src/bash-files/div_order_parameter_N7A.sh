#Run all cases

# Use:
#sh run_all_OP.sh <ion-name> <folder-distances>  <folder-output> <cut-on> <cut_off>
#example: sh  sh order_parameter_N7A LIO LIO/dist order-parameter 0.25 1.0

#mkdir $3/$1
#path_input='/home/tb/secruz/PhD/P1_all_atom/step2/nadine-data/TREMD_RNA-TIP3P/analysis/LIO/dist-demux/'
#path_output='/home/tb/secruz/PhD/P1_all_atom/step2/nadine-data/TREMD_RNA-TIP3P/analysis/python'
#mkdir LI-OP

path_input=$2
path_output=$3

atom='N7A'
mkdir ${path_output}/$1/${atom}

for((i=0; i<20; i++)); do
        echo 'Calculating order parameter' $i
	python averaged-order-parameter.py ${path_input}/dist_MD_${i}_$1-N7A.xvg   ${path_input}/dist_MD_${i}_$1-N7A.xvg $4 $5 ${path_output}/$1/${atom}/N7A_${i}.xvg
        echo '...'
done
