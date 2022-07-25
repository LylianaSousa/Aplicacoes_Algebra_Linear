#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 00:25:19 2022

@author: lyliana
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 23:04:45 2022

@author: lyliana Myllena Santos de Sousa
NUSP: 11223740

Aplicações da Algebra linear  - Estudo da Solução Numérica Exata de uma Equação Diferencial Ordinária de segunda ordem

"""
import numpy as np

#Definir uma função que me da os n-valores de x que usaremos durante a aproximação
def vetor_x(a,b,n,h):
    x = np.zeros(n-1)
    for j in range(1,n):
        x[j-1] = a + (j*h)
    
    return x

#Definir uma função que me da os n-valores de f(x) que usaremos durante a aproximação
def vetor_b(h,n,x):
    b = np.zeros(n-1)
    for j in range(1,n):
        X = x[j-1] 
        b[j-1] = (h**2)*((np.cos(X) - (np.sin(X))**2)*np.exp(np.cos(X)))
              
    return b

#Definir uma função que cria a nossa matriz tridiagonal que representa a discretização da EDO
def matriz_discretizada(n,b,x):
    A = np.zeros((n-1,n-1))
    for linha in range(1,n):
        if linha == 1: 
            A[linha- 1][linha - 1] = 2
            A[linha- 1][linha] = - 1
            
        if linha == n-1:
            A[linha- 1][linha - 2] = - 1
            A[linha- 1][linha-1] = 2
            
        if linha>1 and linha<n-1: 
            A[linha- 1][linha - 2] = - 1
            A[linha- 1][linha - 1] = 2
            A[linha- 1][linha] = - 1
    return A

#Definir uma função que altera a nossa matriz Ã, de acordo com a condição de contorno escolhida
def condiçao_contorno(n,b):
    
    u_a = np.exp(1)#float(input("Insira f(X), com a variável x maiuscula: "))
    u_b = np.exp(1)#float(input("Insira f(X), com a variável x maiuscula: "))
    b[0]= b[0] + u_a 
    b[n-2] = b[n-2] + u_b 
    
    return b

def soluçao_exata(n,x):
    u_exata = np.zeros(n-1)
    for j in range(1,n):
        X = x[j-1] 
        u_exata[j-1] = np.exp(np.cos(X))
              
    return u_exata
 
if __name__ == '__main__':
    
    a = 0#input("a? ")
    b = np.pi #input("b? ")
    n_vet = [64,128,256,512,1024,2048,4096,8192,16384,32768]
    
    for i in range(0,len(n_vet) -1):
        n = n_vet[i]
        h = (b-a)/n
        x = vetor_x(a,b,n,h)
        B = vetor_b(h,n,x)
        A = matriz_discretizada(n,B,x)
        B = condiçao_contorno(n,B)
        u = np.linalg.solve(A, B)
        u_exata = soluçao_exata(n,x)
        norma_euclidiana = np.linalg.norm(u_exata - u, 2)
        norma_infinita = np.linalg.norm(u_exata - u, np.inf)
        
        i += 1
        n2 = n_vet[i]
        h2 = (b-a)/n2
        x = vetor_x(a,b,n2,h2)
        B = vetor_b(h2,n2,x)
        A = matriz_discretizada(n2,B,x)
        B = condiçao_contorno(n2,B)
        u = np.linalg.solve(A, B)
        u_exata = soluçao_exata(n2,x)
        norma_infinita_meio = np.linalg.norm(u_exata - u, np.inf)
        
        p = np.absolute(np.log2(norma_infinita/norma_infinita_meio))
        
        print(f"{n} & {h:.3e} & {norma_euclidiana:.3e} & {norma_infinita:.3e} & {p:.3e} \\\\ \n")
        i += 1