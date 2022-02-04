for i in "3-21G" "6-311G" "6-31G" "6-31G*" "6-31G**" "6-31+G" "6-31++G" "STO-3G" ;
do
	python ./calc_energies_BH.py $i > energies_BH_$(echo $i | sed "s/-//g" | sed "s/*/star/g" | sed "s/+/plus/g").dat
done

for i in "3-21G" "6-311G" "6-31G" "6-31G*" "6-31G**" "6-31+G" "6-31++G" ;
do
	python ./calc_energies_CH2.py $i > energies_CH2_$(echo $i | sed "s/-//g" | sed "s/*/star/g" | sed "s/+/plus/g").dat
done

for i in "3-21G" "6-311G" "6-31G" "6-31G*" "6-31G**" "6-31+G" "6-31++G" ;
do
	python ./calc_energies_H2O.py $i > energies_H2O_$(echo $i | sed "s/-//g" | sed "s/*/star/g" | sed "s/+/plus/g").dat
done

for i in "3-21G" "6-311G" "6-31G" "6-31G*" "6-31G**" "6-31+G" "6-31++G" ;
do
	python ./calc_energies_HF.py $i > energies_HF_$(echo $i | sed "s/-//g" | sed "s/*/star/g" | sed "s/+/plus/g").dat
done

for i in "3-21G" "6-311G" "6-31G" "6-31G*" "6-31G**" "6-31+G" "6-31++G" ;
do
	python ./calc_energies_NH.py $i > energies_NH_$(echo $i | sed "s/-//g" | sed "s/*/star/g" | sed "s/+/plus/g").dat
done
