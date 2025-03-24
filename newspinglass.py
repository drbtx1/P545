# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 21:16:07 2025

@author: dbtx
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

I = np.identity(2)
sz = np.array([[1,0],[0,-1]])
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])

transverse_field = sz
longitudinal_field = sx
N = 7
#h = 1
periodic = False

if periodic == False:
    coupling = np.random.uniform(low = 0, high = 1, size = (N-1,))
if periodic == True:    
    coupling = np.random.uniform(low = 0, high = 1, size = (N,))

#coupling = np.ones(N-1)

#takes number of sites N and boolean periodic, returns sum of tensor products 
#describing nearest neighbor interactions
def build_J(N,periodic, longitudinal_field, coupling):
    J_operator = longitudinal_field
    tensor_sum_J = 0
    #cross-site interactions
    for dim in range(N):
        #create list of operators
        operators = []
        for index in range(0,N): 
            operators.append(I)  #set all to identity and adjust appropriates sits below
     #for each iteration of loop, except last one, set that site and the next one to J_operator         
        if dim != N-1:           
            operators[dim] = J_operator*coupling[dim]
            operators[dim+1] = J_operator    #don't double count J
    #only set last site and the next (first) one to  if the periodic bc are needed        
        if dim == N-1 and periodic == True:
            operators[dim] = J_operator*coupling[dim]
            operators[0] = J_operator
            
        #create tensor product of operators     
        product = operators[0]
        for index in range(1,N):
            product = np.kron(product, operators[index])
    #check not to add extraneous identity in the case that loop is in last iteration
    #and periodic bc are not needed        
        if dim != N-1 or periodic == True:     
            tensor_sum_J = tensor_sum_J + product
    return tensor_sum_J        
            
#takes number of sites N and boolean periodic, returns sum of tensor products 
#describing localized term in  transverse direction at each site
def build_h(N, periodic, transverse_field):
    #intra-site   
    h_operator = transverse_field
    tensor_sum_h = 0
    for dim in range(N):
        #create list of operators
        operators = []
        for index in range(0,N): 
            operators.append(I)
        operators[dim] = h_operator  
        product = operators[0]
        for index in range(1,N):
            product = np.kron(product, operators[index])
        tensor_sum_h = tensor_sum_h + product   
    return tensor_sum_h 
  

                
#takes number of sites N, energy scaling J and h, and boolean periodic (boundary conditions)
#return lists of eigenvalues and eigenfunctions for corresponding Hamiltonian 
def findHamiltonian(N,h,coupling, periodic = False,longitudinal_field = sx, transverse_field = sz):   
    tensor_sum_J = build_J(N, periodic,longitudinal_field, coupling)
    tensor_sum_h = build_h(N, periodic, transverse_field)
    H = -tensor_sum_J - h*tensor_sum_h  
    eigenvalues, eigenvectors = np.linalg.eig(H)
    #eigenvectors = np.transpose(eigenvectors)
    return eigenvalues, eigenvectors.astype(float)




def findThermalDistribution(eigenvalues, beta):
    Z = 0
    prob = np.zeros(len(eigenvalues))
    for val in range(len(eigenvalues)):
        Z += np.exp(-beta*eigenvalues[val]) 
        #print(Z)    
    for val in range(len(eigenvalues)):
        prob[val] = np.exp(-beta*eigenvalues[val])/Z
    return prob   
    

#takes N total sites, and nth specific site, and a Pauli operator
#returns tensor product of pauli on n and identity operator on all other sites
#this is used to expectation value of spin in various directions on site n 
def expectationOperator(n,N, pauli = sz):      #N sites, nth location
    operator = []
    for i in range(N):
        operator.append(I)
    operator[n] = pauli
    product = operator[0]
    for index in range(1,N):
        product = np.kron(product, operator[index])
    return product   

pauliOperatorProducts = []
for n in range(N):
    pauliOperatorProducts.append(expectationOperator(n,N))
#print(pauliOperatorProducts)


#expects a tensor product, and an eigenvector  
#returns expectation value for that operator for the given eigenvector
def findExpectationValue(operator, eigenvector):
    #return np.matmul(np.matmul(eigenvector.conj().T,pauliOperatorProduct), eigenvector)
    return np.matmul(eigenvector.conj().T,np.matmul(operator, eigenvector))   

#ffor each eigenvector find average of <M> for each site, then average
#weight each eigenvector by thermal distribution
#plot that vs B

#returns M for single sites given a specific eigenvector
def localM(vec):
    locM = np.empty(N)    
    for site in range(N):
        locM[site] = findExpectationValue(pauliOperatorProducts[site], vec)
    print(locM)    
    return locM


#vals, vecs = findHamiltonian(N,0.5, coupling)
#thermalweights = findThermalDistribution(vals, 10) 
#weighted_vecs = np.empty((len(vals), len(vals)))
#for idx in range(len(vals)):
#    weighted_vecs[:,idx] = vecs[:,idx]*thermalweights[idx]


def globalM(vec):
    return sum(localM(vec)/N)    


#print(globalM(vecs[:,0]))

#print(vecs[:,3])    
#print(vals)
def overallM(N,h,coupling, B):
    sumM = 0
    vals, vecs = findHamiltonian(N,h, coupling)
    thermalweights = findThermalDistribution(vals, B) 
    for idx in range(len(vals)):
        sumM += globalM(vecs[:,idx])*thermalweights[idx]
    return sumM/N
B = 10


H = []
M = []
for h in np.arange(-4,4,0.1):
    H.append(h)
    M.append(overallM(N, h, coupling, B))
plt.plot(H,M)    
plt.show
    
Blist = []
M = []
h = 1
for B in np.arange(0.1,10, 0.1):
    Blist.append(B)
    M.append(overallM(N, h, coupling, B))
plt.plot(Blist,M)    
plt.xlabel("Beta")      
plt.ylabel("<M>")       
        
plt.show()    


                         