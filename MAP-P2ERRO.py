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
def vetor_fx(n,x):
    equation = "X+1"#input("Insira f(X), com a variável x maiuscula: ")
    fx = np.zeros(n-1)
    for j in range(1,n):
        X = x[j-1] 
        fx[j-1] = eval(equation)
              
    return fx

#Definir uma função que cria a nossa matriz tridiagonal que representa a discretização da EDO
def matriz_discretizada(n,fx,x):
    A = np.zeros((n-1,n))
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
        A[linha- 1][n-1] = fx[linha-1]
    return A

#Definir uma função que altera a nossa matriz Ã, de acordo com a condição de contorno escolhida
def condiçao_contorno(n,A):
    u_a = 10#float(input("Insira f(X), com a variável x maiuscula: "))
    u_b = 10#float(input("Insira f(X), com a variável x maiuscula: "))
    A[0][n-1] = A[0][n-1] + u_a 
    A[n-2][n-1] = A[n-2][n-1] + u_b 
    
    return A

#Definir uma função que reduz a matriz para se encaixar no algoritmo de thomas
def matriz_thomas_algoritmo(n,A):
    A_esc = np.zeros((n-1,4))
    for linha in range(1,n):
        if linha == 1: 
            A_esc[linha- 1][0] = 2
            A_esc[linha- 1][1] = - 1 
            
        if linha == n-1:
            A_esc[linha- 1][1] = - 1
            A_esc[linha- 1][2] = 2
            
        if linha>1 and linha<n-1: 
            A_esc[linha- 1][0] = - 1
            A_esc[linha- 1][1] = 2
            A_esc[linha- 1][2] = - 1
        A_esc[linha-1][3] = A[linha- 1][n-1]
    return A_esc

#Definir uma função que soluciona a nossa matriz pelo Algoritmo de Thomas
def thomas_algoritmo(n,A_esc):
    
    A_esc[0, 2] = A_esc[0, 2]/A_esc[0, 1]
    A_esc[0, 3] = A_esc[0, 3]/A_esc[0, 1]
    
    for i in range(1, n-1):
        ptemp = A_esc[i, 1] - (A_esc[i, 0] * A_esc[i-1, 2])
        A_esc[i, 2] /= ptemp
        A_esc[i, 3] = (A_esc[i, 3] - (A_esc[i, 0] * A_esc[i-1, 3]))/ptemp
    
    x = [0 for i in range(1,n)]
    x[n-2] = A_esc[n-2, 3]
    
    for i in range(1,n-2):
        x[n-2-i] = A_esc[n-2-i, 3] - (A_esc[n-2-i, 2]*x[n-2 -(i+1)])
    return x
 
if __name__ == '__main__':
    
    n = 4#input("n? ")
    a = 0#input("a? ")
    b = 6#input("b? ")
    h = (b-a)/n
    x = vetor_x(a,b,n,h)
    fx = vetor_fx(n,x)
    A = matriz_discretizada(n,fx,x)
    print(A)
    A = condiçao_contorno(n,A)
    A_esc = matriz_thomas_algoritmo(n,A)
    from sklearn.preprocessing import scale
    A = scale( A, axis=0, with_mean=True, with_std=True, copy=True )
    print(A)