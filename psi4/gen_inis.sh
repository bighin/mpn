for i in *.dat ;
do
	for j in 1 2 3 ;
	do
		jnext=$(($j+1))
		cat template.ini | sed "s/DATAFILE/$i/g" | sed "s/MINORDER/$j/g" | sed "s/MAXORDER/$jnext/g" > run10_$(echo $i | sed "s/psi4\///g" | sed "s/\.dat//g")_${j}${jnext}.ini ;
	done
done
