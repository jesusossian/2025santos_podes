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

file_name = sys.argv[1]
tam_particao = int(sys.argv[2])
num_fix = int(sys.argv[3])

USE_FOP = True# Se usa o fix and optimize


RESULT_PATH   = Path('result/')
INSTANCE_PATH = Path('../../data/c1sifa')

from datetime import *
def timer(start_time=None):
	if not start_time:
		start_time = datetime.now()
		return start_time
	elif start_time:
		temp_sec = (datetime.now() - start_time).total_seconds()
		return temp_sec

def main():

	N, PP, PR, FP, FR, HR, HP, D, R ,C= ler.leitura_instance(os.path.join(INSTANCE_PATH,file_name))

	rf_zsp_sol = (np.zeros((N,N))).tolist()
	rf_zsr_sol = (np.zeros((N,N))).tolist()
	rf_zr_sol = (np.zeros((N,N))).tolist()
	rf_l_sol = [0]*N
	rf_yp_sol = [0]*N
	rf_yr_sol = [0]*N

	fo_zsp_sol = (np.zeros((N,N))).tolist()
	fo_zsr_sol = (np.zeros((N,N))).tolist()
	fo_zr_sol = (np.zeros((N,N))).tolist()
	fo_l_sol = [0]*N
	fo_yp_sol = [0]*N
	fo_yr_sol = [0]*N

	zsp_sol = (np.zeros((N,N))).tolist()
	zsr_sol = (np.zeros((N,N))).tolist()
	zr_sol = (np.zeros((N,N))).tolist()
	l_sol = [0]*N
	yp_sol = [0]*N
	yr_sol = [0]*N

	SD = (np.zeros((N,N))).tolist()
	SR = (np.zeros((N,N))).tolist()
	for  i in range(N):
		SD[i][i] = D[i]
		SR[i][i] = R[i]
		for j in range(i+1, N):
			SD[i][j] = SD[i][j-1] + D[j]
			SR[i][j] = SR[i][j-1] + R[j]
	
	subset = gera.gera_particoes(N,tamanho_particao = tam_particao,num_par_fix= num_fix)

	start_rf = timer()
	for conj in subset:
		rf_obj,rf_zsp_sol,rf_zsr_sol,rf_zr_sol,rf_l_sol,rf_yp_sol,rf_yr_sol,rf_bestbound,rf_numnode,rf_gap,rf_elapsed = rf.relax_fix(conj,rf_yp_sol,rf_yr_sol,N,PP,PR,FP,FR,HR,HP,D,R,SD,SR,C)

	temp_rf = timer(start_rf)

	fo_zsp_sol,fo_zsr_sol,fo_zr_sol,fo_l_sol,fo_yp_sol,fo_yr_sol = rf_zsp_sol,rf_zsr_sol,rf_zr_sol,rf_l_sol,rf_yp_sol,rf_yr_sol

	temp_opt = 0.0

	if USE_FOP == True:

		fo_melhor_obj,fo_zsp_melhor_sol,fo_zsr_melhor_sol,fo_zr_melhor_sol,fo_l_melhor_sol,fo_yp_melhor_sol,fo_yr_melhor_sol, fo_melhor_bestbound, fo_melhor_numnode,fo_melhor_gap,fo_melhor_elapsed = rf_obj,rf_zsp_sol,rf_zsr_sol,rf_zr_sol,rf_l_sol,rf_yp_sol,rf_yr_sol, rf_bestbound, rf_numnode,rf_gap,rf_elapsed
		start_opt = timer()
		for conj in subset:
			fo_obj,fo_zsp_sol,fo_zsr_sol,fo_zr_sol,fo_l_sol,fo_yp_sol,fo_yr_sol,fo_bestbound,fo_numnode,fo_gap,fo_elapsed = fop.fix_and_optimize(conj,fo_yp_sol,fo_yr_sol,N, PP, PR, FP, FR, HR, HP, D, R, SD,SR,C,fo_zsp_sol,fo_zsr_sol,fo_zr_sol,fo_l_sol)

			if fo_obj <= fo_melhor_obj:
				fo_melhor_obj,fo_zsp_melhor_sol,fo_zsr_melhor_sol,fo_zr_melhor_sol,fo_l_melhor_sol,fo_yp_melhor_sol,fo_yr_melhor_sol, fo_melhor_bestbound, fo_melhor_numnode,fo_melhor_gap,fo_melhor_elapsed = fo_obj,fo_zsp_sol,fo_zsr_sol,fo_zr_sol,fo_l_sol,fo_yp_sol,fo_yr_sol, fo_bestbound, fo_numnode,fo_gap,fo_elapsed

		temp_opt = timer(start_opt)
	
	temp_total = timer(start_rf)
			
	if USE_FOP == True:

		arquivo = open(os.path.join(RESULT_PATH,'clsr_sp_relax_and_opt_table.txt'),'a')
		arquivo.write(file_name+';'
			+str(round(rf_obj,2))+';'
			+str(round(temp_rf,2))+';'
			+str(round(fo_obj,2))+';'
			+str(round(temp_opt,2))+';'
			+str(round(temp_total,2))
			+'\n'
		)
		arquivo.close()
	else :
		arquivo = open(os.path.join(RESULT_PATH,'clsr_sp_relax_fix_table.txt'),'a')
		arquivo.write(file_name+';'
			+str(round(rf_obj,2))+';'
			+str(round(temp_rf,2))+';'
			+str(round(rf_numnode,2))+';'
			+str(round(temp_total,2))
			+'\n'
		)
		arquivo.close()

if __name__== "__main__" :
	main()
