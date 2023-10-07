import numpy as np
import sys
import matplotlib.pyplot as plt


def finite_difference(h, c):
    row_num = int((1 / h) + 1)
    clmn_num = int((1 / h) + 1)
    zero_matr = np.zeros((row_num, clmn_num), dtype='float64')

    print(zero_matr)

    num_of_equation = (int(1 / h) - 1) * (int(1 / h) - 1)
    solution_matr = np.zeros((num_of_equation, num_of_equation), dtype='float64')
    print(solution_matr)
    current = 0
    for j in range(1, clmn_num - 1):
        for i in range(1, row_num - 1):
            zero_matr[i, j] = current
            solution_matr[current, current] = 4 + c * h * h
            if i != row_num - 2:
                solution_matr[current, current + 1] = -1
            if i != 1:
                solution_matr[current, current - 1] = -1
            if j != clmn_num - 2:
                solution_matr[current, current + (row_num - 2)] = -1
            if j != 1:
                solution_matr[current, current - (row_num - 2)] = -1
            current += 1
    print(solution_matr)
    plt.spy(solution_matr, markersize=7, marker='.', color='blue')
    plt.grid()
    plt.show()


    with open('solution_matr.txt', 'w') as f:
        print(solution_matr, file=f)

finite_difference(0.25, 0.1)
