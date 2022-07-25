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
def condiçao_contorno(n,fx):
    u_a = 10#float(input("Insira f(X), com a variável x maiuscula: "))
    u_b = 10#float(input("Insira f(X), com a variável x maiuscula: "))
    fx[0]= fx[0] + u_a 
    fx[n-2] = fx[n-2] + u_b 
    
    return fx
 
if __name__ == '__main__':
    
    n = 4#input("n? ")
    a = 0#input("a? ")
    b = 6#input("b? ")
    h = (b-a)/n
    x = vetor_x(a,b,n,h)
    fx = vetor_fx(n,x)
    A = matriz_discretizada(n,fx,x)
    print("\n Matriz A ")
    print(A)
    fx = condiçao_contorno(n,fx)
    print("\n Vetor b")
    print(fx)
    u = np.linalg.solve(A, fx)
    print("\n Vetor u")
    print(u)