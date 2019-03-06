import numpy as np
from scipy.linalg import eigh_tridiagonal
from matplotlib import pyplot as plt

w = 0.1
b = 0.05
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

	return V_vec, n


def analyze(n_w, iter_lim = 5):
	n = n_w*200 + 500
	L = (20 + n_w)*w + (n_w -1)*b
	n = 5000
	x_vec = np.linspace(0, L, n)
	delta_x = L / (n + 1)

	#H = np.zeros([n, n])
	main_diag = np.ones(n)
	off_diag = np.ones(n-1)
		
	if n_w == 0:
		V_vals = V_box(x_vec)

	else:
		V_vals, _ = well(x_vec, w, n_w, b, L)

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
		print("heyhey: ", energies, "\n")
		for j in range(3):
			bw[j][i-lower] = np.absolute(energies[i*j] - energies[(j+1)*i -1])
			#print("yoo: ", energies[(j+1)*i - 1])
			#print(bw[j][i - lower])
			#print(i, " øvre E: ", (j+1)*i - 1, "nedre E:", i*j, "diff: ", bw[j][i-lower])
			print(i, "lower, upper:",energies[i*j], " ", energies[(j+1)*i - 1])


	for i in range(3):
		plt.plot(x_vec, bw[i])
	plt.savefig("bw.png")
	plt.show()

def band_widths_noob(lower, upper):
	x_vec = np.linspace(2, 10, 9)
	bw_1 = np.zeros(9)
	bw_2 = np.zeros(9)
	bw_3 = np.zeros(9)

	for i in range(9):
		n_w = i + 2
		energies, a, z, c, d, e = analyze(n_w, n_w * 3)
		bw_1[i] = energies[n_w - 1] - energies[0]
		#print(n_w, "diff", bw_1[i], " energies[n_w-1] ", energies[n_w - 1], " energies[0] ", energies[0])
		#print(energies[:31])
		bw_2[i] = energies[2*n_w - 1] - energies[n_w]
		bw_3[i] = energies[3*n_w - 1] - energies[2*n_w]

	plt.plot(x_vec, bw_1)
	plt.plot(x_vec, bw_2)
	plt.plot(x_vec, bw_3)

	plt.savefig("bwnoob.png")
	plt.show()

def plot_wave_funcs(energies, wave_funcs, iter_lim, L, V_vals, n_w):
	x_vec = np.linspace(0, L, len(V_vals))
	for i in range(iter_lim): 
		plt.plot([0, L], [energies[i]]*2)
		plt.plot(x_vec, V_vals)
		plt.plot(x_vec, energies[i] + wave_funcs[i]*3000)
	plt.savefig("figW" + str(n_w) + "E" + str(iter_lim) + ".png")
	plt.show()


E_b, w_b, iter_lim_b, L_b, V_vals_b = analyze(0)
plot_wave_funcs(E_b, w_b, iter_lim_b, L_b, V_vals_b)
energies, wave_funcs, iter_lim, L, V_vals = analyze(1)
plot_wave_funcs(energies, wave_funcs, iter_lim, L, V_vals)

band_widths(1, 50)
