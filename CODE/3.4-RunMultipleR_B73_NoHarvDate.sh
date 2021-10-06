#!/bin/bash

for i in `seq 0 30 18840`; do
	j01=$(($i + 1))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j01 B73 > out.txt &

	j02=$(($i + 2))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j02 B73 > out.txt &

	j03=$(($i + 3))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j03 B73 > out.txt &

	j04=$(($i + 4))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j04 B73 > out.txt &

	j05=$(($i + 5))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j05 B73 > out.txt &

	j06=$(($i + 6))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j06 B73 > out.txt &

	j07=$(($i + 7))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j07 B73 > out.txt &

	j08=$(($i + 8))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j08 B73 > out.txt &

	j09=$(($i + 9))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j09 B73 > out.txt &

	j10=$(($i + 10))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j10 B73 > out.txt &

	j11=$(($i + 11))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j11 B73 > out.txt &

	j12=$(($i + 12))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j12 B73 > out.txt &

	j13=$(($i + 13))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j13 B73 > out.txt &

	j14=$(($i + 14))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j14 B73 > out.txt &

	j15=$(($i + 15))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j15 B73 > out.txt &

	j16=$(($i + 16))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j16 B73 > out.txt &

	j17=$(($i + 17))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j17 B73 > out.txt &

	j18=$(($i + 18))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j18 B73 > out.txt &

	j19=$(($i + 19))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j19 B73 > out.txt &

	j20=$(($i + 20))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j20 B73 > out.txt &

	j21=$(($i + 21))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j21 B73 > out.txt &

	j22=$(($i + 22))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j22 B73 > out.txt &

	j23=$(($i + 23))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j23 B73 > out.txt &

	j24=$(($i + 24))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j24 B73 > out.txt &

	j25=$(($i + 25))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j25 B73 > out.txt &

	j26=$(($i + 26))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j26 B73 > out.txt &

	j27=$(($i + 27))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j27 B73 > out.txt &

	j28=$(($i + 28))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j28 B73 > out.txt &

	j29=$(($i + 29))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j29 B73 > out.txt &

	j30=$(($i + 30))
	nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args $j30 B73 > out.txt

done

nohup R --vanilla --slave < 3.4-BLUE_each_NoHarvDate.R --args 18871 B73 > out.txt &
