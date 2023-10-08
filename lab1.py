import numpy as np
import matplotlib.pyplot as plt


def finite_difference(h, c):
    row_num = int((1 / h) + 1)
    clmn_num = int((1 / h) + 1)
    zero_matr = np.zeros((row_num, clmn_num), dtype='float64')

    # print(zero_matr)

    num_of_equation = (int(1 / h) - 1) * (int(1 / h) - 1)
    coef_matr = np.zeros((num_of_equation, num_of_equation), dtype='float64')
    # print(coef_matr)
    current = 0
    for j in range(1, clmn_num - 1):
        for i in range(1, row_num - 1):
            zero_matr[i, j] = current
            coef_matr[current, current] = 4 + c * (h ** 2)
            if i != row_num - 2:
                coef_matr[current, current + 1] = -1
            if i != 1:
                coef_matr[current, current - 1] = -1
            if j != clmn_num - 2:
                coef_matr[current, current + (row_num - 2)] = -1
            if j != 1:
                coef_matr[current, current - (row_num - 2)] = -1
            current += 1
    # print(coef_matr)
    plt.spy(coef_matr, markersize=7, marker='.', color='blue')
    plt.grid()
    plt.show()

    with open('solution_matr.txt', 'w') as f:
        print(coef_matr, file=f)

    right_part = np.empty(num_of_equation)
    right_part.fill(h ** 2)

    # print(right_part)
    # print(type(right_part))


finite_difference(0.25, 0.1)
