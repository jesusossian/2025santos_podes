from pathlib import Path
import os
import leitura as ler
import optimization as opt
import relax_fix as rf
import fix_optimize as fo
import gera_particoes as gera
import numpy as np
import pandas as pd
import sys
from datetime import datetime, date
import time

# parameters
file_name = sys.argv[1]
tam_part = int(sys.argv[2]) # tamanho das particoes
num_fix = int(sys.argv[3]) # qtd variaveis fixadas 

USE_FOP = True

result_path   = Path('result/')
instance_path = Path('../../data/c1sifa')

def timer(start_time=None):
	if not start_time:
		start_time = datetime.now()
		return start_time
	elif start_time:
		temp_sec = (datetime.now() - start_time).total_seconds()
		return temp_sec

def main():
    # read instances
	N, PP, PR, FP, FR, HR, HP, D, R, C = ler.leitura_instance(os.path.join(instance_path, file_name))

	yp_sol = np.zeros(N)
	yr_sol = np.zeros(N)

	SD = (np.zeros((N,N))).tolist()
	SR = (np.zeros((N,N))).tolist()
	for  i in range(N):
		SD[i][i] = D[i]
		SR[i][i] = R[i]
		for j in range(i+1, N):
			SD[i][j] = SD[i][j-1] + D[j]
			SR[i][j] = SR[i][j-1] + R[j]

	subset = gera.gera_particoes(N, tam_part, num_fix)
		
	start_rf = timer()
	
	for conj in subset:
		rf_obj, yp_sol, yr_sol = rf.relax_fix(conj, yp_sol, yr_sol, N, PP, PR, FP, FR, HR, HP, D, R, SD, SR, C)

	temp_rf = timer(start_rf)
	
	if USE_FOP == True:
		
		start_fo = timer()
		
		for conj in subset:
			fo_obj, yp_sol, yr_sol = fo.fix_optimize(conj, yp_sol, yr_sol, N, PP, PR, FP, FR, HR, HP, D, R, SD, SR, C)
			
		temp_fo = timer(start_fo)

	temp_total = timer(start_rf)
	
	if USE_FOP == True:
		arquivo = open(os.path.join(result_path,'clsr_std_math_result.txt'),'a')
		arquivo.write(file_name+';'
			+str(round(rf_obj,2))+';'
			+str(round(temp_rf,2))+';'
			+str(round(fo_obj,2))+';'
			+str(round(temp_fo,2))+';'
			+str(round(temp_total,2))
			+'\n')
		arquivo.close()
	else :
		arquivo = open(os.path.join(result_path,'clsr_std_rf_result.txt'),'a')
		arquivo.write(file_name+';'
			+str(round(rf_obj,2))+';'
			+str(round(temp_rf,2))+';'
			+str(round(temp_fo,2))+';'
			+str(round(temp_total,2))
			+'\n')
		arquivo.close()

if __name__== "__main__" :
	main()
