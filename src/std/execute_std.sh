#!/bin/bash
#52 periodos, 108 instancias

form=std

for id in {1..108} #$(seq 1)
do
	python3 clsr_math_std_mip.py c52_${id}.txt >> report/out_${form}_c52_${id}.txt
done
