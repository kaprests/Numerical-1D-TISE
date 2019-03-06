import numpy as np
from scipy.linalg import eigh_tridiagonal
from matplotlib import pyplot as plt

w = 0.1
b = 0.05
n = 1000
m = 1
h_bar = 1
num_well = 4

# Box potential
V_box = lambda x: 0*x

# n well potential
def well(x, w, n_w, b, L, V0=-3000):
	n = len(x)
	V_mid_ex = np.array([[V0]*int(w*n/L) + [0]*int(b*n/L)]*n_w).flatten()
	V_mid = V_mid_ex[:-int(b*n/L)]
	
	m = len(V_mid)
	diff = n - m
	len_b = diff//2
	len_f = diff - len_b
	V_front = np.array([0]*len_f)
	V_back = np.array([0]*len_b)

	V_temp = np.append(V_front, V_mid)
	V_vec = np.append(V_temp, V_back)

	return V_vec


def analyze(n_w, iter_lim = 5):
	L = (20 + n_w)*w + (n_w -1)*b
	n = 1000
	x_vec = np.linspace(0, L, n)
	delta_x = L / (n + 1)

	#H = np.zeros([n, n])
	main_diag = np.ones(n)
	off_diag = np.ones(n-1)
	if n_w == 0:
		V_vals = V_box(x_vec)

	else:
		V_vals = well(x_vec, w, n_w, b, L)





	#H = np.zeros([n, n])
	main_diag = np.zeros(n)
	off_diag = np.zeros(n-1)




	main_diag *= ((h_bar**2)/(m*(delta_x**2)))
	main_diag += V_vals
	off_diag *= -(h_bar**2)/(2*m*(delta_x**2))
	
	energies, wave_funcs = eigh_tridiagonal(main_diag, off_diag)
	wave_funcs = wave_funcs.T
	return energies, wave_funcs, iter_lim, L, V_vals, n_w


def band_widths(lower, upper):
	x_vec = np.linspace(lower, upper, upper-lower)
	bw = np.zeros([3, upper-lower])
	for i in range(lower, upper):
		energies, _,_,_,_ = analyze(i)
		for j in range(3):
			bw[j][i-lower] = np.absolute(energies[i*j] - energies[(j+1)*i -1])


	for i in range(3):
		plt.plot(x_vec, bw[i])
	plt.savefig("bw.png")
	plt.show()


def plot_wave_funcs(energies, wave_funcs, iter_lim, L, V_vals, n_w):
	x_vec = np.linspace(0, L, len(V_vals))
	for i in range(iter_lim): 
		plt.plot([0, L], [energies[i]]*2)
		plt.plot(x_vec, V_vals)
		plt.plot(x_vec, energies[i] + wave_funcs[i]*3000)
	plt.savefig("figW" + str(n_w) + "E" + str(iter_lim) + ".png")
	plt.show()


#energies, wave_funcs, iter_lim, L, V_vals, n_w = analyze(8, 24)
#plot_wave_funcs(energies, wave_funcs, iter_lim, L, V_vals, n_w)
#energies, wave_funcs, iter_lim, L, V_vals, n_w = analyze(9, 27)
#plot_wave_funcs(energies, wave_funcs, iter_lim, L, V_vals, n_w)
energies, wave_funcs, iter_lim, L, V_vals, n_w = analyze(1, 3)
plot_wave_funcs(energies, wave_funcs, iter_lim, L, V_vals, n_w)

band_widths(2, 10)
