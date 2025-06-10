from pathlib import Path
import os
import leitura as ler
import optimization as opt
import relax_fix as rf
import fix_optimize as fop
import gera_particoes as gera
import numpy as np
import pandas as pd
import sys
from datetime import datetime, date
import time

#parameters
file_name = sys.argv[1]
tam_particao = int(sys.argv[2]) 
num_fix = int(sys.argv[3])

USE_FOP = True

result_path   = Path('result/')
instance_path = Path('../../data/csifa')

#global variables
objval = 0
objbound = 0
nodecount = 0
runtime = 0
mipgap = 0

def timer(start_time=None):
	if not start_time:
		start_time = datetime.now()
		return start_time
	elif start_time:
		temp_sec = (datetime.now() - start_time).total_seconds()
		return temp_sec

def main():
    # read instances
	N, PP, PR, FP, FR, HR, HP, D, R, C = ler.leitura_instance(os.path.join(instance_path,file_name))

	xp_sol = [0]*N
	xr_sol = [0]*N
	sp_sol = [0]*N
	sr_sol = [0]*N
	yp_sol = [0]*N
	yr_sol = [0]*N

	rf_xp_sol = [0]*N
	rf_xr_sol = [0]*N
	rf_sp_sol = [0]*N
	rf_sr_sol = [0]*N
	rf_yp_sol = [0]*N
	rf_yr_sol = [0]*N

	fop_xp_sol = [0]*N
	fop_xr_sol = [0]*N
	fop_sp_sol = [0]*N
	fop_sr_sol = [0]*N
	fop_yp_sol = [0]*N
	fop_yr_sol = [0]*N

	SD = (np.zeros((N,N))).tolist()
	SR = (np.zeros((N,N))).tolist()
	for  i in range(N):
		SD[i][i] = D[i]
		SR[i][i] = R[i]
		for j in range(i+1, N):
			SD[i][j] = SD[i][j-1] + D[j]
			SR[i][j] = SR[i][j-1] + R[j]
			
	subset = gera.gera_particoes(N, tamanho_particao=tam_particao, num_par_fix=num_fix)
		
	start_rf = timer()
	
	for conj in subset:
		rf_obj, rf_xp_sol, rf_xr_sol, rf_sp_sol, rf_sr_sol, rf_yp_sol, rf_yr_sol, \
			rf_bestbound, rf_numnode, rf_gap, rf_elapsed = rf.relax_fix(conj, rf_yp_sol, \
				rf_yr_sol, N, PP, PR, FP, FR, HR, HP, D, R, SD, SR, C)

	temp_rf = timer(start_rf)

	fop_obj, fop_xp_sol, fop_xr_sol, fop_sp_sol, fop_sr_sol, fop_yp_sol, fop_yr_sol = \
	    rf_obj, rf_xp_sol, rf_xr_sol, rf_sp_sol, rf_sr_sol, rf_yp_sol, rf_yr_sol 
	
	temp_opt = 0.0

	if USE_FOP == True:
		
		start_opt = timer()
		
		for conj in subset:
			fop_obj, fop_xp_sol, fop_xr_sol, fop_sp_sol, fop_sr_sol, fop_yp_sol, fop_yr_sol, \
				fop_bestbound, fop_numnode, fop_gap, fop_elapsed = \
			fop.fix_and_optimize(conj, fop_yp_sol, fop_yr_sol, N, PP, PR, FP, FR, HR, HP, \
				D, R, SD, SR, C, fop_xp_sol, fop_xr_sol, fop_sp_sol, fop_sr_sol)
			
		temp_opt = timer(start_opt)

	temp_total = timer(start_rf)
	
	if USE_FOP == True:
		arquivo = open(os.path.join(result_path,'clsr_std_relax_and_opt_table.txt'),'a')
		arquivo.write(file_name+';'
			+str(round(rf_obj,2))+';'
			+str(round(temp_rf,2))+';'
			+str(round(fop_obj,2))+';'
			+str(round(temp_opt,2))+';'
			+str(round(temp_total,2))
			+'\n')
		arquivo.close()
	else :
		arquivo = open(os.path.join(result_path,'clsr_std_relax_fix_table.txt'),'a')
		arquivo.write(file_name+';'
			+str(round(rf_obj,2))+';'
			+str(round(temp_rf,2))+';'
			+str(round(temp_opt,2))+';'
			+str(round(temp_total,2))
			+'\n')
		arquivo.close()

if __name__== "__main__" :
	main()
