import numpy as np
from scipy.linalg import eigh_tridiagonal
from matplotlib import pyplot as plt


# Box potential
V_box = lambda x: 0*x

# n well potential
def well(x, w, n_w, b,V0=-3000):
	n = len(x)
	V_mid_ex = np.array([[V0]*int(w*n) + [0]*int(b*n)]*n_w).flatten()
	V_mid = V_mid_ex[:-int(b*n)]
	
	m = len(V_mid)
	diff = n - m
	len_b = diff//2
	len_f = diff - len_b
	V_front = np.array([0]*len_f)
	V_back = np.array([0]*len_b)

	V_temp = np.append(V_front, V_mid)
	V_vec = np.append(V_temp, V_back)

	return V_vec


def analyze(V, n_w):
	w = 0.1
	b = 0.1
	L = (20 + n_w)*w + (n_w -1)*b
	n = 2500
	m = 1
	h_bar = 1
	delta_x = L/(n+1)
	x_vec = np.linspace(0, L, n)


	#H = np.zeros([n, n])
	main_diag = np.zeros(n)
	off_diag = np.zeros(n-1)
		
	if n_w == 0:
		V_vals = V(x_vec)

	else:
		V_vals = V(x_vec, w, n_w, b)
		print(V_vals)

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

analyze(V_box, 0)
analyze(well, 5)

