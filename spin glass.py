# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 18:36:32 2025

@author: dbtx
"""

#code uses standard Hamiltonian for Ising model with and addition spin that can 
#move and couple with various pairs of spins


import numpy as np
#import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
#from sortev import *
#from sortcolumns import *
#from sortvecs import *
#from scipy.linalg import eig
from sympy import Matrix
#import tensorflow as tf

#Following lines make plots look a little more Latex-like
#mpl.rcParams['mathtext.fontstyle'] = 'cm' 
mpl.rcParams['font.family'] = 'serif' # or 'sans-serif' or 'monospace'
mpl.rcParams['font.serif'] = 'cmr10'
mpl.rcParams['font.sans-serif'] = 'cmss10'
mpl.rcParams['font.monospace'] = 'cmtt10'
mpl.rcParams["axes.formatter.use_mathtext"] = True # to fix the minus signs

I = np.identity(2)
sz = np.array([[1,0],[0,-1]])
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])
splus = (sx + 1j*sy)/2
sminus = (sx - 1j*sy)/2
transverse_field = sz
longitudinal_field = sx
N = 3
h = 1
periodic = False
if periodic == False:
    coupling = np.random.uniform(low = 0, high = 1, size = (N-1,))
if periodic == True:    
    coupling = np.random.uniform(low = 0, high = 1, size = (N,))
print(coupling)
coupling = np.ones(2**N)
print(coupling)
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
            operators[dim] = J_operator*couping[dim]
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

vals, vecs = findHamiltonian(N, 0, coupling)
#print(vals)

def findThermalDistribution(eigenvalues, beta):
    Z = 0
    prob = np.zeros(len(eigenvalues))
    for val in range(len(eigenvalues)):
        Z += np.exp(-beta*eigenvalues[val]) 
        #print(Z)    
    for val in range(len(eigenvalues)):
        prob[val] = np.exp(-beta*eigenvalues[val])/Z
    return prob   
    
#print(findThermalDistribution(vals, 1))    





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


#expects a list of tensor products (give expectation values at each site), and 
#a list of eigenvectors  
#returns a list of expectation values (corresponding to each eigenvector) for that site
def findExpectationValue(operator, eigenvector):
    #return np.matmul(np.matmul(eigenvector.conj().T,pauliOperatorProduct), eigenvector)
    return np.matmul(eigenvector.conj().T,np.matmul(operator, eigenvector))     



def averageExpectationValue(pauliOperatorProducts, vecs, weight):
    #print(pauliOperatorProducts)
    count = len(vecs[:, 0])
    sum = 0
    for idx in range(count):
        for site in range(N):
            sum += weight[idx]*findExpectationValue(pauliOperatorProducts[site], vecs[:,idx] )
    return sum/count     
H = []
M = []
for h in np.arange(0,2, 0.1):
    H.append(h)
    vals, vecs = findHamiltonian(N, h, coupling)
    print(vecs[0])
    weight = coupling# np.zeros(2**N)
    #weight[0] = 1
    M.append(averageExpectationValue(pauliOperatorProducts, vecs, weight))
    
plt.plot(H,M)    












    
def findAndSaveMagnetization(N,periodic,eV,J,maxh,starth,steps,\
                             longitudinal_field, transverse_field, direction_of_M):
    h_per_J = []
    operator_list = []
    avgExpValue = []
    #realPart = []
    #imagPart = []

#build list of operators to find expectation value of nth site
#nth element in each list holds tensor product of sigma operator at nth site
#and all other sites holding identity operator
    for site in range(0,N):
        operator_list.append(expectationOperator(site,N,direction_of_M))
     
    last = []   
    for step in range(steps): 
        #hold expectation values at nth site 
        expectation_at_site = [] 
        h = starth + (maxh-starth)*step/steps
        #print(h)
        eigenvalues, eigenvectors = findHamiltonian(N,J,h,periodic,longitudinal_field, transverse_field)
        #assure sorted in decreasing order
        idx = eigenvalues.argsort()[::-1]   
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:,idx]
        #startstep = 5
        #if step >= startstep:
        #    eigenvectors = sortEigenvectorsasarray(eigenvectors,last)
        #print(eigenvectors)    
        last = eigenvectors 
        #print(eigenvectors[:,eV])
        #print(eigenvalues[:,eV])
        for site in range(N):
            expectation_at_site.append(findExpectationValue(operator_list[site], eigenvectors[:,eV]))
            #expectation_at_site.append(findExpectationValue(operator_list[site], eigenvectors[eV]))
            #print(findExpectationValue(operator_list[site], eigenvectors[:,eV]))
        avg = sum(expectation_at_site)/len(expectation_at_site)   
        #print(avg)
        h_per_J.append(h/J)
        avgExpValue.append(avg)  #average across N sites, all with same eigenvector
        #realPart.append(avg.real)
        #print(avg.real)
        #imagPart.append(avg.imag)
        #print(avg.imag)
        
        
        
    #data = avgExpValue    
    title = "N=" +str(N)+",eigenvector = " + str(eV) + ", periodic =" + str(periodic)
    '''filetype = '.csv'
    dest = 'C:/Users/dabuch/Ising Data/Magnetization/PeriodicTrue/JxhzMx/N' + str(N) + '/' + title + filetype
    
    data = {"h/J": h_per_J, "Expectation Value averaged over all sites": avgExpValue}
    df = pd.DataFrame(data)
    df.to_csv(dest)'''
    #mpl.rcParams['text.usetex'] = True
    plt.xlabel("h/J")
    plt.ylabel("M")#" (averaged over all sites)")
    #plt.title(title)   
    plt.plot(h_per_J, avgExpValue)
    #plt.legend(["Real","Imaginary"])
    #
    #plt.plot(h_per_J, realPart)
    #plt.plot(h_per_J, imagPart)
    #image_filetype = '.png'
    #dest = 'C:/Users/dbtx/.spyder-py3/' + title +  image_filetype
    #plt.savefig(dest)
    #plt.show()
    
    
    
            
                    
        