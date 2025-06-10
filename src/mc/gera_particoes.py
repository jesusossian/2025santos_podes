
def gera_particoes(N,tamanho_particao=5,num_par_fix=2,semente = 5,indice_geracao = 1):
	tam_jane = tamanho_particao - num_par_fix
	subset = []
	 
	if indice_geracao == 1:
		for i in range(0,N,tam_jane):
			if i + tamanho_particao > N:
				subset.append([k for k in range(i,N)])
			else:
				subset.append([k for k in range(i,i+tamanho_particao)])
	else:
		from random import sample
		import random
		random.seed(semente)
		lista_periodos = [k  for k in range(N)]
		while True:
			 
			sub_conj = []
			 
			if len(lista_periodos) >= tamanho_particao:
				sub_conj = sample(lista_periodos,tamanho_particao)
			else:
				sub_conj = lista_periodos[:]
			subset.append(sub_conj)
			for i in sub_conj:
				lista_periodos.remove(i)
			if len(lista_periodos)==0:
				break
			
	return subset