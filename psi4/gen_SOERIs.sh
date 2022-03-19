for molecule in $(echo "BH CH2 H2O HF NH")
do
	for basis in $(echo "STO-3G 3-21G 6-31G 6-31+G 6-31++G 6-31G* 6-31G** 6-311G")
	do
		shortid=$(echo $basis | sed "s/-//g" | sed "s/+/plus/g" | sed "s/*/star/g")
		python ./calc_SOERI_${molecule}.py ${basis} > ${molecule}_${shortid}.dat
	done
done
