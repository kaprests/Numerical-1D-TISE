import numpy as np
from matplotlib import pyplot as plt

L = 1
n = 100
m = 1
h_bar = 1
delta_x = L/(n+1)
x_vec = np.linspace(0, L, n)

H = np.zeros([n, n])

V = lambda x: 0

for i in range(n):
	# Fills main diagonal
	H[i][i] = ((h_bar**2)/(m*(delta_x**2))) + V(x[i])
for i in range(n-1):
	# Fills second diagonals below and above main
	H[i+1][i] = H[i][i +1] = (h_bar**2)/(2*m*(delta_x**2))
	

