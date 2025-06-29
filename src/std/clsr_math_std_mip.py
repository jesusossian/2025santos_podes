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
import itertools

# parameters
file_name = sys.argv[1]
tam_part = int(sys.argv[2]) # tamanho das particoes
num_fix = int(sys.argv[3]) # qtd variaveis fixadas 

# paths
result_path = Path('result/')
instance_path = Path('../../data/c1sifa')

def main():
    
    # read instances
	N, PP, PR, FP, FR, HP, HR, D, R, C = ler.leitura_instance(os.path.join(instance_path,file_name))
	
	yp_val = np.zeros(N)
	yr_val = np.zeros(N)

	SD = (np.zeros((N,N))).tolist()
	SR = (np.zeros((N,N))).tolist()
	for  i in range(N):
		SD[i][i] = D[i]
		SR[i][i] = R[i]
		for j in range(i+1, N):
			SD[i][j] = SD[i][j-1] + D[j]
			SR[i][j] = SR[i][j-1] + R[j]
		
	subset = gera.gera_particoes(N, tam_part, num_fix)
	
	start_time = time.time()
	for conj in subset:
		rf_objval, yp_val, yr_val = rf.relax_fix(conj, yp_val, yr_val, N, PP, PR, FP, FR, HP, HR, D, R, SD, SR, C)

	rf_rtime = time.time() - start_time

	start_time = time.time()
	for conj in subset:
		fo_objval, yp_val, yr_val = fo.fix_optimize(conj, yp_val, yr_val, N, PP, PR, FP, FR, HP, HR, D, R, SD, SR, C)

	fo_rtime = time.time() - start_time

	objval, objbound, mgap, rtime, ncount, tmp = opt.clsr_std(N, PP, PR, FP, FR, HP, HR, D, R, SD, SR, C, yp_val, yr_val)
	
	time_total = rf_rtime + fo_rtime + rtime

	arquivo = open(os.path.join(result_path,'clsr_std_math_mip_result.txt'),'a')
	arquivo.write(file_name+';'
	 	+str(round(objval,2))+';'
	 	+str(round(objbound,2))+';'
	 	+str(round(mgap,2))+';'
	 	+str(round(time_total,2))+';'
	 	+str(round(ncount,2))+';'
	 	+str(round(tmp,2))
	 	+'\n')
	arquivo.close()

if __name__== "__main__" :
	main()
