import numpy as np
import matplotlib.pyplot as plt


def iter_solver(alpha, beta, tau, callback):
    residual = np.zeros((10000, 1000))
    u_solution = np.zeros(1000)
    rk_r0_arr = []

    f = np.zeros(1000)
    for i in range(len(f)):
        if (i >= 494) & (i <= 504):
            f[i] = 1
        else:
            f[i] = 0

    flag = 0
    while (True):
        if flag >= 10000:
            break
        residual[flag] = f - callback(alpha, beta, u_solution)
        u_solution = u_solution + tau * residual[flag]

        rk_r0_arr = np.append(rk_r0_arr, np.linalg.norm(residual[flag]) / np.linalg.norm(residual[0]))
        if np.linalg.norm(residual[flag]) / np.linalg.norm(residual[0]) <= 0.0001:
            break

        flag += 1

    return u_solution, rk_r0_arr, flag


def mult_matr(alpha, beta, u):
    res_mult = np.zeros(1000)
    for i in range(1000):
        if i == 0:
            res_mult[0] = (2 + beta) * u[0] + (-1) * u[1]
        elif i == 999:
            res_mult[999] = (2 + alpha) * u[999] + (-1) * u[998]
        else:
            res_mult[i] = (2 + alpha) * u[i] + (-1) * u[i - 1] + (-1) * u[i + 1]
    return res_mult


precond = np.zeros(1000)
precond[0] = 12
for i in range(1, 1000):
    precond[i] = 2.01


def iter_with_precond(alpha, beta, precond):
    f = np.zeros(1000)
    for i in range(len(f)):
        if (i >= 494) & (i <= 504):
            f[i] = 1
        else:
            f[i] = 0
    rk_r0_arr = []
    r0 = np.linalg.norm(f)
    rk_r0 = r0
    x_k = np.zeros(1000)
    precondition = np.reciprocal(precond)
    flag = 0
    while rk_r0 > 0.0001:
        Ax_k = mult_matr(alpha, beta, x_k)
        x_k_next = x_k - np.multiply(precondition, mult_matr(alpha, beta, x_k) - f)
        rk_r0 = np.linalg.norm(Ax_k - f) / r0
        rk_r0_arr += [rk_r0]
        x_k = x_k_next
        flag += 1
    return x_k_next, rk_r0_arr, flag


def conjgrad(A, b, x):
    r = b - np.dot(A, x)
    p = r
    rsold = np.dot(np.transpose(r), r)
    r_arr = []
    flag = 0
    for i in range(len(b)):
        Ap = np.dot(A, p)
        alpha = rsold / np.dot(np.transpose(p), Ap)
        x = x + np.dot(alpha, p)
        r = r - np.dot(alpha, Ap)
        rsnew = np.dot(np.transpose(r), r)
        if np.linalg.norm(rsnew) < 0.0001:
            break
        p = r + (rsnew / rsold) * p
        rsold = rsnew
        r_arr += [rsold]
        flag += 1
    return x, r_arr, flag


A = np.zeros((1000, 1000))
A[0][0] = 12
for i in range(0, 1000):
    A[i][i] = 2.01
    if i <= 998:
        A[i + 1][i] = -1
        A[i][i + 1] = -1

f = np.zeros(1000)
for i in range(len(f)):
    if (i >= 494) & (i <= 504):
        f[i] = 1
    else:
        f[i] = 0

x = np.zeros(1000)

solution1, residuals1, iter_number1 = iter_solver(0.01, -3, 0.165, mult_matr)
# print(solution1)
# print(len(residuals1))
# print(iter_number1)
iter_number1_arr = range(iter_number1)
# print(len(iter_number1_arr))
solution2, residuals2, iter_number2 = iter_with_precond(0.01, 10, precond)
# print(solution2)
# print(residual2)
# print(iter_number2)
iter_number2_arr = range(iter_number2)

solution3, residuals3, iter_number3 = conjgrad(A, f, x)
# print(solution3)
# print(residual3)
# print(iter_number3)
iter_number3_arr = range(iter_number3)

x = range(1000)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(x, solution1, s=47, c='g', marker='^', label='SI')
ax1.scatter(x, solution2, s=10, c='m', marker='s', label='SI + diag. precond.')
ax1.scatter(x, solution3, s=6, c='k', marker='.', label='CG')
plt.legend(loc='upper right')
plt.grid()
plt.show()

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.scatter(iter_number1_arr, residuals1, s=1, c='g', marker='^', label='SI')
ax2.scatter(iter_number2_arr, residuals2, s=1, c='m', marker='s', label='SI + diag. precond.')
ax2.scatter(iter_number3_arr, residuals3, s=7, c='k', marker='.', label='CG')
plt.legend(loc='upper right')
plt.grid()
plt.yscale("log")
plt.ylim([0.000001, 100])
plt.show()
