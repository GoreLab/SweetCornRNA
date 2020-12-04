#!/bin/bash

for i in `seq 0 10 22580`; do
	j01=$(($i + 1))
	nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args $j01 v1.1 PH207 > out.txt &

	j02=$(($i + 2))
	nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args $j02 v1.1 PH207 > out.txt &

	j03=$(($i + 3))
	nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args $j03 v1.1 PH207 > out.txt &

	j04=$(($i + 4))
	nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args $j04 v1.1 PH207 > out.txt &

	j05=$(($i + 5))
	nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args $j05 v1.1 PH207 > out.txt &

	j06=$(($i + 6))
	nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args $j06 v1.1 PH207 > out.txt &

	j07=$(($i + 7))
	nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args $j07 v1.1 PH207 > out.txt &

	j08=$(($i + 8))
	nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args $j08 v1.1 PH207 > out.txt &

	j09=$(($i + 9))
	nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args $j09 v1.1 PH207 > out.txt &

	j10=$(($i + 10))
	nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args $j10 v1.1 PH207 > out.txt &

	j11=$(($i + 11))
	nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args $j11 v1.1 PH207 > out.txt &

	j12=$(($i + 12))
	nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args $j12 v1.1 PH207 > out.txt &

	j13=$(($i + 13))
	nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args $j13 v1.1 PH207 > out.txt &

	j14=$(($i + 14))
	nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args $j14 v1.1 PH207 > out.txt &

	j15=$(($i + 15))
	nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args $j15 v1.1 PH207 > out.txt &

	j16=$(($i + 16))
	nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args $j16 v1.1 PH207 > out.txt &

	j17=$(($i + 17))
	nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args $j17 v1.1 PH207 > out.txt &

	j18=$(($i + 18))
	nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args $j18 v1.1 PH207 > out.txt &

	j19=$(($i + 19))
	nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args $j19 v1.1 PH207 > out.txt &

	j20=$(($i + 20))
	nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args $j20 v1.1 PH207 > out.txt

done


nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args 22601 v1.1 PH207 > out.txt &
nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args 22602 v1.1 PH207 > out.txt &
nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args 22603 v1.1 PH207 > out.txt &
nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args 22604 v1.1 PH207 > out.txt &
nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args 22605 v1.1 PH207 > out.txt &
nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args 22606 v1.1 PH207 > out.txt &
nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args 22607 v1.1 PH207 > out.txt &
nohup R --vanilla --slave < 1.3.1-BLUE_each.R --args 22608 v1.1 PH207 > out.txt
