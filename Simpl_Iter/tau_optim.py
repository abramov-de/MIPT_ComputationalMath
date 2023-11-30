import numpy as np

A = np.zeros((1000, 1000))

A[0][0] = 2 + (-3)
for i in range(1, 1000):
    A[i][i] = 2 + 0.01

for i in range(0, 999):
    A[i+1][i] = -1
    A[i][i+1] = -1

# w, v = np.linalg.eig(A)

norm_A = np.linalg.norm(A)
rev_A = np.linalg.inv(A)
norm_rev_A = np.linalg.norm(rev_A)


print(norm_A * norm_rev_A)
