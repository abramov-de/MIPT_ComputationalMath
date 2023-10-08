import numpy as np
import matplotlib.pyplot as plt
import time
from scipy import linalg


def finite_difference(h, c):
    row_num = int((1 / h) + 1)
    clmn_num = int((1 / h) + 1)
    zero_matr = np.zeros((row_num, clmn_num), dtype='float64')

    num_of_equation = (int(1 / h) - 1) * (int(1 / h) - 1)
    coef_matr = np.zeros((num_of_equation, num_of_equation), dtype='float64')

    flag = 0
    for j in range(1, clmn_num - 1):
        for i in range(1, row_num - 1):
            zero_matr[i, j] = flag
            coef_matr[flag, flag] = 4 + c * (h ** 2)
            if i != int(1 / h) - 1:
                coef_matr[flag, flag + 1] = -1
            if i != 1:
                coef_matr[flag, flag - 1] = -1
            if j != int(1 / h) - 1:
                coef_matr[flag, flag + (int(1 / h) - 1)] = -1
            if j != 1:
                coef_matr[flag, flag - (int(1 / h) - 1)] = -1
            flag += 1

    # with open('coeff_matr.txt', 'w') as f:
    #     print(coef_matr, file=f)

    right_part = np.empty(num_of_equation)
    right_part.fill(h ** 2)

    return coef_matr, right_part, clmn_num, row_num, zero_matr


coef_matr, right_part, clmn_num, row_num, zero_matr = finite_difference(1/150, 0.1)
plt.spy(coef_matr, markersize=7, marker='.', color='blue')
plt.grid()
plt.title('Visualization of coefficient matrix.')
plt.show()

print('Matrix of coefficients: \n', coef_matr, '\n')
print('Right part of matrix: \n', right_part, '\n')


def lu_solver(coef_matr, right_part):
    start = time.time()

    lu, piv = linalg.lu_factor(coef_matr)
    x = linalg.lu_solve((lu, piv), right_part)

    finish = time.time()

    print('Solutions:\n', x, '\n')

    flag = 0
    for j in range(1, clmn_num - 1):
        for i in range(1, row_num - 1):
            zero_matr[i, j] = x[flag]
            flag += 1
    solution_matr = zero_matr

    print('Solutions matrix:\n', solution_matr, '\n')
    print('Solution time: ', finish - start)

    plt.imshow(solution_matr, cmap='rainbow')
    plt.show()

    flag = 0
    grid_size = []
    times = []
    h_denominator = 4
    while flag == 0:
        start = time.time()
        finite_difference(1 / h_denominator, 0.1)
        finish = time.time()
        grid_size.append(h_denominator + 1)
        times.append(finish - start)
        if h_denominator < 150:
            h_denominator += 10
        else:
            break
        print(h_denominator + 1)
        print(finish - start)

    plt.plot(grid_size, times)
    plt.xlabel('Size of grid')
    plt.ylabel('Solution time')
    plt.grid()
    plt.show()

lu_solver(coef_matr, right_part)
