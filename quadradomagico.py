#  Python 3 program for
#  Magic Squares of Even Order
import numpy as np
#  Find magic square of Even Order
def magicSquare(n) :
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
    #  Display given size
    print("\n  Magic square of size (",n, "X",n,")")
    #  Display result
    row = 0
    
    while (row < n) :
            col = 0
            while (col < n) :
                print(" ",matrix[row][col], end ="")
                col += 1
            print()
            row += 1
    k = (int(n * (n * n + 1) / 2))
    print("  Sum of each rows and columns is",k)
    
    return matrix

if __name__=="__main__":   
    #  Test
    print(magicSquare(7))
    print(magicSquare(4))