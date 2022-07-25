#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  3 15:08:15 2022

@author: lyliana Myllena Santos de Sousa
NUSP: 11223740

Aplicações da Algebra linear  - Estudo de Algoritmos de Solução Numérica de Sistemas Lineares de Equações.

"""

#bibliotecas que usarei ao longo da análise
import numpy as np
import time
from scipy import linalg

#Função que cria a Matrix Quadrado Magico
def quadradomagico(n) :
    if (n > 0 and n % 2 == 0) :
        #  Use to collect result
        matrix = [[0] * (n) for _ in range(n)]
        #  Some auxiliary variable
        k = int(n / 2)
        row = 0
        col = int(k / 2)
        if ((int(n / 2)) % 2 != 0) :
            #  When n is Even and its divisible by 2 reminder is not zero
            #  6,10,14,18,22,26..(4*n +2)
            #  Exmple (6/2 = 3) and its reminder not zero
            #  Execute loop (n/2)*2 times.
            value = 1
            while (value <= k * k) :
                #  Add element in top left grid
                matrix[row][col] = value
                #  Add element in bottom right grid
                matrix[k + row][k + col] = value + (k * k)
                #  Add element in top right grid
                matrix[row][k + col] = value + (k * k * 2)
                #  Add element in bottom left grid
                matrix[k + row][col] = value + (k * k * 3)
                col += 1
                row -= 1
                if (value % k == 0) :
                    col -= 1
                    row += 2
                else :
                    if (col == k) :
                        col = 0
                    elif (row < 0) :
                        row += k
                value += 1
        else :
            #  Doubly Even order magic square of n
            #  (4*n)
            row = 0
            while (row < n) :
                col = 0
                while (col < n) :
                    matrix[row][col] = (n * row) + col + 1
                    col += 1
                row += 1
            row = 0
            while (row < int(n / 4)) :
                col = 0
                while (col < int(n / 4)) :
                    #  Top Left corner
                    matrix[row][col] = (n * n + 1) - matrix[row][col]
                    col += 1
                row += 1
            row = 0
            while (row < int(n / 4)) :
                col = 3 * (int(n / 4))
                while (col < n) :
                    #  Top right corner
                    matrix[row][col] = (n * n + 1) - matrix[row][col]
                    col += 1
                row += 1
            row = int(3 * n / 4)
            while (row < n) :
                col = 0
                while (col < int(n / 4)) :
                    #  Bottom Left corner
                    matrix[row][col] = (n * n + 1) - matrix[row][col]
                    col += 1
                row += 1
            row = int(3 * n / 4)
            while (row < n) :
                col = int(3 * n / 4)
                while (col < n) :
                    #  Bottom right corner
                    matrix[row][col] = (n * n + 1) - matrix[row][col]
                    col += 1
                row += 1
            row = int(n / 4)
            while (row < int(3 * n / 4)) :
                col = int(n / 4)
                while (col < int(3 * n / 4)) :
                    #  Centre elements
                    matrix[row][col] = (n * n + 1) - matrix[row][col]
                    col += 1
                row += 1
    elif (n > 0 and n % 2 != 0) :   
        matrix = np.zeros((n,n), dtype=int)

        passo = 1
        i, j = 0, n//2
        
        while passo <= n**2:
            matrix[i, j] = passo
            passo += 1
            newi, newj = (i-1) % n, (j+1)% n
            if matrix[newi, newj]:
                i += 1
            else:
                i, j = newi, newj
    
    #Exibindo resultados
    row = 0
    
    while (row < n) :
            col = 0
            while (col < n) :
                print(" ",matrix[row][col], end ="")
                col += 1
            print()
            row += 1
    k = (int(n * (n * n + 1) / 2))
    
    return np.array(matrix)

#Fatoração QR/ Classical Gram-Schmidt
def fatoração_QR(A):
    """Classical Gram-Schmidt (CGS) algorithm"""
    m, n = A.shape
    R = np.zeros((n, n))
    Q = np.zeros((m, n))
    v = np.zeros((m, n))
    for j in range(0,n):
        v[:,j] = A[:,j]
        for i in range(0,j-1):
            R[i, j] = np.dot(Q[i,:].T, A[:,j])
            v[:,j]  = v[:,j]  - np.dot(R[i, j], Q[i,:])
        R[j, j] = linalg.norm(v[:, j])
        Q[:, j] = v[:, j] / R[j,j]
    return Q, R

#Fatoração Gram-Schmidt Ortogonalizado
def fatoração_GS(A):
    """Modified Gram-Schmidt (CGS) algorithm"""
    m, n = A.shape
    R = np.zeros((n, n))
    Q = np.zeros((m, n))
    v = np.zeros((m, n))
    for i in range(0,n):
        v[:,i] = A[:,i]
    for i in range(0,n):
            R[i, i] = linalg.norm(v[:, i])
            Q[:,i]  = v[:,i]/R[i, i]
            for j in range(i+1,n-1):
                R[i, j] = np.dot(Q[i,:].T, v[:, j])
                v[:, j] = v[:, j] - np.dot(R[i, j], Q[i,:])
    return Q, R

#Função que fornece a Matriz de HouseHolder
def householder(a):
    #find prependicular vector to mirror
    u = a / (a[0] + np.copysign(np.linalg.norm(a), a[0]))
    u[0] = 1
    H = np.eye(a.shape[0])
    #finding Householder projection
    H -= (2 / np.dot(u, u)) * np.matmul(u[:, None], u[None, :])
    return H

#Fatoração por Reflexão de Householder
def fatoração_RH(A):
    """Householder reflection algorithm"""
    m, n = A.shape
    Q = np.eye(m)
    for i in range(n - (m == n)):
        H = np.eye(m)
        H[i:, i:] = householder(A[i:, i])
        Q = np.matmul(Q,H)
        A = np.matmul(H,A)
    return Q, A
    

if __name__ == '__main__':
    
    n = 5#int(input("n? "))
    A = quadradomagico(n)
    

    # Q, R = linalg.qr(A)
    # print('-----------------')
    # print('Resultado')
    # print('-----------------\n Q = ')
    # print(Q)
    # print('-----------------\n R = ')
    # print(R)
    # print('-----------------\n A = ')
    # print(np.matmul(Q,R))
    
    Q_fatoração_QR, R_fatoração_QR = fatoração_QR(A)
    # print('-----------------')
    # print('Fatoração QR')
    # print('-----------------\n Q = ')
    # print(Q_fatoração_QR)
    # print('-----------------\n R = ')
    # print(R_fatoração_QR)
    # print('-----------------\n A = ')
    A_fatoração_QR = np.matmul(Q_fatoração_QR,R_fatoração_QR)
    norma_fatoração_QR = linalg.norm(A - A_fatoração_QR)
    norma_qualidade_QR = linalg.norm(np.identity(n) - np.matmul(Q_fatoração_QR.T,Q_fatoração_QR))
    print("norma_fatoração_QR = ", norma_fatoração_QR)
    print("norma_qualidade_QR = ", norma_qualidade_QR)
    
    Q_fatoração_GS, R_fatoração_GS = fatoração_GS(A)
    # print('-----------------')
    # print('Fatoração GS')
    # print('-----------------\n Q = ')
    # print(Q_fatoração_GS)
    # print('-----------------\n R = ')
    # print(R_fatoração_GS)
    # print('-----------------\n A = ')
    A_fatoração_GS = np.matmul(Q_fatoração_GS,R_fatoração_GS)
    norma_fatoração_GS = linalg.norm(A - A_fatoração_GS)
    norma_qualidade_GS = linalg.norm(np.identity(n) - np.matmul(Q_fatoração_GS.T,Q_fatoração_GS))
    print("norma_fatoração_GS = ", norma_fatoração_GS)
    print("norma_qualidade_GS = ", norma_qualidade_GS)
    
    Q_fatoração_RH, R_fatoração_RH = fatoração_RH(A)
    # print('-----------------')
    # print('Fatoração RH')
    # print('-----------------\n Q = ')
    # print(Q_fatoração_RH)
    # print('-----------------\n R = ')
    # print(R_fatoração_RH)
    # print('-----------------\n A = ')
    A_fatoração_RH = np.matmul(Q_fatoração_RH,R_fatoração_RH)
    norma_fatoração_RH = linalg.norm(A - A_fatoração_RH)
    norma_qualidade_RH = linalg.norm(np.identity(n) - np.matmul(Q_fatoração_RH.T,Q_fatoração_RH))
    print("norma_fatoração_RH = ", norma_fatoração_RH)
    print("norma_qualidade_RH = ", norma_qualidade_RH)