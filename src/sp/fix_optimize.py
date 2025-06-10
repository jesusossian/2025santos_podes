import gurobipy as gp
from gurobipy import GRB
import numpy as np

MAX_CPU_TIME = 3600.0
EPSILON = 0.000001

def fix_and_optimize(particoes,yp_sol,yr_sol,N,PP,PR,FP,FR,HR,HP,D,R,SD,SR,C,zsp_sol,zsr_sol,zr_sol,l_sol):

	zsp_sol1 = (np.zeros((N,N))).tolist()
	zsr_sol1 = (np.zeros((N,N))).tolist()
	zr_sol1 = (np.zeros((N,N))).tolist()
	l_sol1 = [0]*N
	yp_sol1 = [0]*N
	yr_sol1 = [0]*N

	objval = 0.0
	bestbound = 0.0
	numnode = 0.0
	elapsed = 0.0
	gap = 0.0

	CSP = (np.zeros((N,N))).tolist()
	CSR = (np.zeros((N,N))).tolist()
	CR = (np.zeros((N,N))).tolist()
	CL = [0]*N

	for i in range(N):
		for j in range(i,N):
			CR[i][j]  = sum(HR[t]*SR[i][t] for t in range(i,j))
			CSP[i][j] = PP[i] * SD[i][j] + sum(HP[t]*SD[t+1][j] for t in range(i,j))
			CSR[i][j] = PR[i] * SD[i][j] + sum(HP[t]*SD[t+1][j] for t in range(i,j))

	for i in range(N):
		CL[i] = sum(HR[j]*SR[i][j] for j in range(i,N))
	
	try:

		# Create a new model
		model = gp.Model("CLSR")

		# Create variables
	
		zsp = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="z_sp")
		zsr = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="z_sr")
		zr = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="z_r")
		l = model.addVars(list(range(N)), lb = 0.0, ub = 1.0, vtype=GRB.CONTINUOUS, name="l")
		yp = model.addVars(list(range(N)), vtype=GRB.BINARY, name="yp")
		yr = model.addVars(list(range(N)), vtype=GRB.BINARY, name="yr")
  
		for i in range(N):
			if i not in particoes:
				yp[i].lb = yp_sol[i]
				yp[i].ub = yp_sol[i]
				yr[i].lb = yr_sol[i]
				yr[i].ub = yr_sol[i]
			else:
				yp[i].start = yp_sol[i]
				yr[i].start = yr_sol[i]
		for i in range(N):
			l[i].start = l_sol[i]
			for j in range(i,N):
				zsp[i,j].start = zsp_sol[i][j]
				zsr[i,j].start = zsr_sol[i][j]
				zr[i,j].start = zr_sol[i][j]

		model.update()

		# set objective

		obj = None
		for i in range(N):
			obj += yp[i]*FP[i] + yr[i]*FR[i] + l[i]*CL[i] 
			for j in range(i,N):
				obj += zsp[i,j]*CSP[i][j] + zsr[i,j]*CSR[i][j]+zr[i,j]*CR[i][j] 

		model.setObjective(obj, sense = GRB.MINIMIZE)

		# add constraints
		model.addConstr(
			gp.quicksum(zsp[0,j] + zsr[0,j] for j in range(N)) ==1
		)
		
		model.addConstrs(
			gp.quicksum(zsp[i,t-1] + zsr[i,t-1] for i in range(t)) - 
			gp.quicksum(zsp[t,j] + zsr[t,j] for j in range(t, N)) == 0  for t in range(1,N) 
		)
				
		model.addConstrs(
			gp.quicksum(zsp[t,j] for j in range(t,N)) <= yp[t] for t in range(N)
		)
			
		model.addConstrs(
			gp.quicksum(zsr[t,j] for j in range(t,N)) <= yr[t] for t in range(N)
		)
			
		model.addConstr(
			gp.quicksum(zr[0,j] for j in range(N)) + l[0]==1
		)
					
		model.addConstrs(
			gp.quicksum(zr[i,t-1] for i in range(0,t)) == 
			gp.quicksum(zr[t,j]  for j in  range(t,N)) + l[t] for t in range(1,N)
		)       
				
		model.addConstrs(
			gp.quicksum(zr[i,t] for i in range(0,t+1)) <= yr[t] for t in range(N)
		)    
			
		model.addConstrs(
			gp.quicksum(SR[i][t]*zr[i,t] for i in range(t+1) ) ==
			gp.quicksum(SD[t][j]*zsr[t,j] for j in range(t,N)) for t in range(N)
		)
		
		model.addConstrs(
			gp.quicksum(SD[t][k]*zsp[t,k] for k in range(t,N)) + 
			gp.quicksum(SD[t][k]*zsr[t,k] for k in range(t,N)) <= C for t in range(N)
		)
	 
		# parameters 
		model.setParam(GRB.Param.TimeLimit,MAX_CPU_TIME)
		model.setParam(GRB.Param.MIPGap,EPSILON)
		model.setParam(GRB.Param.Threads,1)
		#model.setParam(GRB.Param.Cuts,-1)
		#model.setParam(GRB.Param.Presolve,-1)

		# optimize model
		model.optimize()

		objval = model.ObjVal
		objbound = model.ObjBound
		nodecount = model.NodeCount
		runtime = model.Runtime
		mipgap = model.MIPGap

		for i in range(N):
			for j in range(i,N):
				zsp_sol1[i][j] = zsp[i,j].X
				zsr_sol1[i][j] = zsr[i,j].X
				zr_sol1[i][j] = zr[i,j].X

		l_sol1  = [l[i].X for i in range(N)]
		yp_sol1 = [yp[i].X for i in range(N)]
		yr_sol1 = [yr[i].X for i in range(N)]

	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))

	return objval, zsp_sol1, zsr_sol1, zr_sol1, l_sol1, yp_sol1, yr_sol1, objbound, nodecount, mipgap, runtime
