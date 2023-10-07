import numpy as np


def finite_difference(h):
    row_num = int((1 / h) + 1)
    clmn_num = int((1 / h) + 1)
    zero_matr = np.zeros((row_num, clmn_num), dtype='int32')

    print(zero_matr)

    num_of_equation = (int(1 / h) - 1) * (int(1 / h) - 1)
    solution_matr = np.zeros((num_of_equation, num_of_equation), dtype='int32')
    print(solution_matr)


finite_difference(0.25)
