import numpy as np
from scipy.linalg import eigh_tridiagonal
from matplotlib import pyplot as plt

L = 1
n = 100
m = 1
h_bar = 1
delta_x = L/(n+1)
x_vec = np.linspace(0, L, n)

#H = np.zeros([n, n])
main_diag = np.zeros(n)
off_diag = np.zeros(n-1)

# Box potential
V = lambda x: 0

for i in range(n):
	# Fills main diagonal
	#H[i][i] = ((h_bar**2)/(m*(delta_x**2))) + V(x_vec[i])
	main_diag[i] = ((h_bar**2)/(m*(delta_x**2))) + V(x_vec[i])

for i in range(n-1):
	# Fills second diagonals below and above main
	#H[i+1][i] = H[i][i +1] = (h_bar**2)/(2*m*(delta_x**2))
	off_diag[i] = -(h_bar**2)/(2*m*(delta_x**2))
	

energies, wave_funcs = eigh_tridiagonal(main_diag, off_diag)
wave_funcs = wave_funcs.T


plt.plot(x_vec, wave_funcs[0])
plt.show()
