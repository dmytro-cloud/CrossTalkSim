for i in {1..5}
do
	for j in 1 25 50 75
	do
		if [ $j -eq 1 ]
		then
			python3 electron_holes_prob.py ../output/HPK_full_0.${j}_d${i}.root HPK_0.000${j}00_d0.00${i}000.txt
		else
			python3 electron_holes_prob.py ../output/HPK_full_0.${j}_d${i}.root HPK_0.000${j}0_d0.00${i}000.txt
		fi
	done
	python3 electron_holes_prob.py ../output/HPK_full_1.00_d${i}.root HPK_0.001000_d0.00${i}000.txt
done

