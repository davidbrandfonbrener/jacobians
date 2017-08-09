import numpy as np
from time import time
from multiprocessing import Pool
from sage.all import *


# ith column, ith row removed
def minor(arr, i):
    return arr[np.array(range(i)+range(i+1, arr.shape[0]))[:,np.newaxis],
               np.array(range(i)+range(i+1,arr.shape[1]))]

# function to create a random graph G_n,p with cyclic Jacobian and check 
# 	(1) if a fixed delta_xy generates Jac(G)
#	(2) if there exists a delta_xy that generates Jac(G)
# uses the theorem from our paper and just checks gcd(Jac(G), Jac(G_1))
def test_if_dxy_generator(n):
    set_random_seed()
    p = .5
    G = graphs.RandomGNP(n,p)
    L = np.array(G.laplacian_matrix())
    redL = matrix(minor(L, 0))
    snf, P, _ = redL.smith_form()
    
    # generate connected graph with cyclic Jacobian 
    while not G.is_connected() or snf[-2, -2] != 1:
        G = graphs.RandomGNP(n,.5)
        L = np.array(G.laplacian_matrix())
        redL = matrix(minor(L, 0))
        snf, P, _ = redL.smith_form()
    
    # indicator variables of generators
    fixedxy = 0
    existsxy = 0

    # set a = |Jac(G)|
    a = redL.determinant()

    # loop over possible xy to check for generators
    for i in range(n):
        for j in range(n):
           if i < j : 
           		# calculate b = |Jac(G_1)|
                m = matrix(minor(L, i))
                if L[i,j] == 0:
                	m[j-1, j-1] += 1
               	else:
               		m[j-1, j-1] -= 1
                b = m.determinant()
                
                if gcd(a, b) == 1:
                    existsxy = 1
                    # check the fixed generator
                    if i == 0 and j == 1:
                    	fixedxy = 1
                   	break
        if existsxy == 1:
        	break

    return [existsxy, fixedxy]
    

# Test Conjecture that if G is biconnected then exponent(Jac(G)) >= n
def test_dxy_order_n(n):
    set_random_seed()
    p = .5
    G = graphs.RandomGNP(n,p)
    L = np.array(G.laplacian_matrix())
    redL = matrix(minor(L, 0))
    snf, P, _ = redL.smith_form()
    
    # generate biconnected graph 
    while not G.is_biconnected():
        G = graphs.RandomGNP(n,.5)
        L = np.array(G.laplacian_matrix())
        redL = matrix(minor(L, 0))
        snf, P, _ = redL.smith_form()
    
    # indicator variables of order >= n
    dxy_order_n = 0

    # find inverse of reduced Laplacian
    m = redL**(-1)

    # loop over possible xy to check orders
    # x if fixed at 0
    for j in range(1,n):
        v = [0] * (n-1)
        v[j-1] = 1
        v = vector(v) 

        # calculate order delta_xy
        order = lcm([denominator(k) for k in v*m])

        if order >= 1:
        	dxy_order_n = 1
        	break

    return dxy_order_n

    
if __name__ == "__main__":
 
	p = Pool(6)
	
	# values of n to loop over
	n_vals = [5, 10, 20, 40]

	# number of trials to conduct for each n
	trials = 100000

	with open('jacobians_results.txt', 'w') as f:
		for n in n_vals:
			t1 = time()
			# conduct trials
			out = np.array(p.map(test_if_dxy_generator, [n]*trials))
			out2 = np.array(p.map(test_dxy_order_n, [n]*trials))

			# write output 
			f.write('\n n:' + str(n) + '\n exists generator prob:' + str((1.0 * sum(out[:,0])) / len(out[:,0])) + '\n fixed generator prob:' + str((1.0 * sum(out[:,1])) / len(out[:,1])) + '\n')
			f.write(' order n dxy prob:' + str((1.0 * sum(out2)) / len(out2)) +'\n')

			t2 = time()
			print t2 - t1


