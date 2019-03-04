import numpy as np
from scipy.linalg import eigh_tridiagonal
from matplotlib import pyplot as plt

L = 2.1
n = 1000
m = 1
h_bar = 1
delta_x = L/(n+1)
x_vec = np.linspace(0, L, n)

#H = np.zeros([n, n])
main_diag = np.zeros(n)
off_diag = np.zeros(n-1)

# Box potential
V_box = lambda x: 0

# Single well potential
def single_well(x, V0=-3000, w=0.1, x_0=L/2 - 0.1/2):
	if x >= x_0 and x <= x_0 + w:
		return V0
	else:
		return 0*x

# Double well potential

# Many wells potential


def analyze(V):
	V_vals = np.zeros(n)
	for i in range(n):
		V_vals[i] = V(x_vec[i])


	for i in range(n):
		# Fills main diagonal
		#H[i][i] = ((h_bar**2)/(m*(delta_x**2))) + V(x_vec[i])
		main_diag[i] = ((h_bar**2)/(m*(delta_x**2))) + V_vals[i]

	for i in range(n-1):
		# Fills second diagonals below and above main
		#H[i+1][i] = H[i][i +1] = (h_bar**2)/(2*m*(delta_x**2))
		off_diag[i] = -(h_bar**2)/(2*m*(delta_x**2))

	energies, wave_funcs = eigh_tridiagonal(main_diag, off_diag)
	wave_funcs = wave_funcs.T

	for i in range(5):
		plt.plot([0, L], [energies[i]]*2)
		plt.plot(x_vec, V_vals)
		plt.plot(x_vec, energies[i] + wave_funcs[i]*1000)
	plt.show()

analyze(V_box)
analyze(single_well)

