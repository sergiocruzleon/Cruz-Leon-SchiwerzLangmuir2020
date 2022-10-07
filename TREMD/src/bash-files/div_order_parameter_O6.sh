#Run all cases

# Use:
#sh run_all_OP.sh <ion-name> <folder-distances>  <folder-output> <cut-on> <cut_off>
#example: sh run_all_averaged_OP.sh LIO ../LIO/dist-demux averaged 0.25 1.0

#mkdir ./$3/$1
#path_input='/home/tb/secruz/PhD/P1_all_atom/step2/nadine-data/TREMD_RNA-TIP3P/analysis/LIO/dist-demux/'
#path_output='/home/tb/secruz/PhD/P1_all_atom/step2/nadine-data/TREMD_RNA-TIP3P/analysis/python'
#mkdir LI-OP

path_input=$2
path_output=$3
atom='O6G'
mkdir ${path_output}/$1/${atom}

for((i=0; i<20; i++)); do
        echo 'Calculating order parameter' $i
	python averaged-order-parameter2.py ${path_input}/dist_MD_${i}_$1-O6G.xvg   ${path_input}/dist_MD_${i}_$1-O6G.xvg $4 $5 ${path_output}/$1/${atom}/O6G_${i}.xvg
done
