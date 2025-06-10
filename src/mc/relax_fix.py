import gurobipy as gp
from gurobipy import GRB
import numpy as np
MAX_CPU_TIME = 3600.0
EPSILON = 0.000001

def relax_fix(particoes,yp_sol ,yr_sol,N, PP, PR, FP, FR, HR, HP, D, R, SD,SR,C):

	xp_sol1 = [0]*N
	xr_sol1 = [0]*N
	wp_sol1 = (np.zeros((N,N))).tolist()
	wr_sol1 = (np.zeros((N,N))).tolist()
	vor_sol1 = (np.zeros((N,N))).tolist()
	yp_sol1 = [0]*N
	yr_sol1 = [0]*N

	objval = 0.0
	bestbound = 0.0
	numnode = 0.0
	elapsed = 0.0
	gap = 0.0

	CP = [0]*N
	CR = [0]*N
	KP = 0
	KR = 0 

	for i in range(N):
		CP[i] = PP[i]

		for j in range(i,N):
			CP[i] = CP[i] + HP[j]

	for i in range(N):
		CR[i] = PR[i]
		for j in range(i,N):
			CR[i] += HP[j]
		for j in range(i,N):
			CR[i] -= HR[j]

	for i in range(N):
		AUX = 0 
		for j in range(i+1):
			AUX += D[j]

		KP += HP[i]*AUX

	for i in range(N):
		AUX = 0

		for j in range(i+1):
			AUX += R[j]

		KR += HR[i]*AUX

	K = KR - KP

	try:

		# create model
		model = gp.Model("CLSR")

		# create variables
		wp = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="wp")
		wr = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="wr")
		vor = model.addVars([(i,j) for i in range(N) for j in range(i,N)], vtype=GRB.CONTINUOUS, name="vor")
		xp = model.addVars(list(range(N)), vtype=GRB.CONTINUOUS, name="xp")
		xr = model.addVars(list(range(N)), vtype=GRB.CONTINUOUS, name="xr")
		yp = model.addVars(list(range(N)), vtype=GRB.BINARY, name="yp")
		yr = model.addVars(list(range(N)), vtype=GRB.BINARY, name="yr")
		
		for i in range(N):
			if i > max(particoes) :
				yp[i].VType = gp.GRB.CONTINUOUS
				yp[i].lb    = 0
				yp[i].ub    = 1
				yr[i].VType = gp.GRB.CONTINUOUS
				yr[i].lb    = 0
				yr[i].ub    = 1
			elif i in particoes:
				yp[i].VType = gp.GRB.BINARY
				yp[i].lb = 0
				yp[i].ub = 1
				yr[i].VType = gp.GRB.BINARY
				yr[i].lb = 0
				yr[i].ub = 1
			else:
				yp[i].lb = yp_sol[i]
				yp[i].ub = yp_sol[i]
				yr[i].lb = yr_sol[i]
				yr[i].ub = yr_sol[i]
	
		model.update()

		# set objective
		obj = 0
		for i in range(N):
			obj += xp[i]*CP[i]
			obj += yp[i]*FP[i]
			obj += xr[i]*CR[i]
			obj += yr[i]*FR[i]
		obj += K

		# add constraints
		model.setObjective(obj, sense = GRB.MINIMIZE)
		
		# add constraints
		for i in range(N):
			ctr = 0.0
			for j in range(i+1):
				ctr += wp[j,i]
				ctr += wr[j,i]
			model.addConstr(ctr >= D[i])

		for i in range(N):
			ctr = 0.0
			for j in range(i+1):
				ctr+= vor[j,i]
			for j in range(i,N):
				ctr+=(-wr[i,j])
			model.addConstr(ctr == 0)

		for i in range(N):
			ctr =0.0
			for j in range(i,N):
				ctr += vor[i,j]
			model.addConstr(ctr <= R[i])

		for i in range(N):
			for j in range(i,N):
				model.addConstr(wp[i,j] + yp[i]*(-D[j]) <=0)

		for i in range(N):
			for j in range(i,N):
				model.addConstr(wr[i,j] + yr[i]*(-min(SR[0][i],D[j])) <=0)

		for i in range(N):
			for j in range(i,N):
				model.addConstr(vor[i,j] + yr[j]*(-R[i]) <= 0)

		for i in range(N):
			ctr = 0.0 
			ctr += xp[i]
			for j in range(i,N):
				ctr += (-wp[i,j])
			model.addConstr(ctr == 0)

		for i in range(N):
			ctr = 0.0
			ctr += xr[i]
			for j in range(i,N):
				ctr += (-wr[i,j])
			model.addConstr(ctr == 0)

		model.addConstrs(
			xp[i] + xr[i] <= C for i in range(N)
		)
		
		# Parameters 
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
	
		xp_sol1  = [xp[i].X for i in range(N)]
		xr_sol1  = [xr[i].X for i in range(N)]
		for i in range(N):
			for j in range(i,N):
				wp_sol1[i][j] = wp[i,j].X
				wr_sol1[i][j] = wr[i,j].X
				vor_sol1[i][j] = vor[i,j].X
		yp_sol1  = [yp[i].X for i in range(N)]
		yr_sol1  = [yr[i].X for i in range(N)]

	except gp.GurobiError as e:
		print('Error code ' + str(e.errno) + ': ' + str(e))

	return objval, xp_sol1, xr_sol1, wp_sol1, wr_sol1, vor_sol1, yp_sol1, yr_sol1, objbound, nodecount, mipgap, runtime
