import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve
from scipy.constants import hbar, m_e, eV
from scipy.linalg import eigh_tridiagonal

hbar *=1e9 #Plancks constant, scaled for nm as recomended by mr. MAN ;)
L_well = 1 #nm - length/width of each potential well
V0 = 3 #eV - depth of the wells
n_well = 10 # number of datapoints per well
n_bar = 5 # number of datapoints per barrer
fact = (hbar**2)/(m_e*eV) # precalculated factor to lessen float operations


# Constructs a potential with num_w number of wells
def potential(num_w, V0=-V0):
	V1 = np.array([[V0]*n_well + [0]*n_bar]).flatten() # one well with barrier to the right
	V_mid = np.tile(V1, num_w) 
	V_front = np.zeros(n_well*10)
	V_back = np.zeros(n_well*10 - n_bar)
	V = np.append(V_front, V_mid)
	V = np.append(V, V_back)
	well_voids = max(num_w -1, 0)
	L = num_w*L_well + well_voids*n_bar*(L_well/n_well) + 20*L_well
	return V, L # returns the final potential and length of the system
	

# Calculates eigenenergies and eigenvector(wave functions)
def analyze(num_w, add_E_lvls = 0):
	V, L = potential(num_w)
	n = len(V)
	delta_x = L/(n+1)
	x_vec = np.linspace(0, L, n)
	if num_w == 0:
		num_E_lvls = 3 + add_E_lvls
	else:
		num_E_lvls = 3*num_w + add_E_lvls
	
	main_diag = np.ones(n) * fact/(delta_x**2) + V
	off_diag = np.ones(n-1) * -fact/(2*delta_x**2)

	energies, wave_funcs = eigh_tridiagonal(main_diag, off_diag)
	wave_funcs = wave_funcs.T

	return energies, wave_funcs, x_vec, num_E_lvls, L, V


def plot_well_wave_funcs(num_w, add_E_lvls=0):
	energies, wave_funcs, x_vec, iter_lim, L, V = analyze(num_w, add_E_lvls)
	scale = V0
	if num_w == 0:
		scale = 1
	for i in range(iter_lim):
		plt.plot([0,L], [energies[i]]*2)
		plt.plot(x_vec, V)
		plt.plot(x_vec, energies[i] + wave_funcs[i]*scale)
	plt.savefig("figW" + str(num_w) + "E" + str(iter_lim) + ".pdf")
	plt.show()


# Computes and plots energy band widths for potentials with lower to upper number of wells
def band_widths(lower, upper):
	x_vec = np.linspace(lower, upper, upper-lower)
	band_widths = np.zeros([3, upper - lower])
	for i in range(lower, upper):
		energies,_,_,_,_,_ = analyze(i)
		for j in range(3):
			band_widths[j][i - lower] = np.absolute(energies[i*j] - energies[(j+1)*i -1])

	for i in range(3):
		plt.plot(x_vec, band_widths[i])
	#plt.savefig("bw.pdf")
	plt.show()


# Numerically computes the energies from the analytic solutions of TISE and prints them
def trancendent_sol():
	a = L_well/2
	
	z0 = (a/hbar)*np.sqrt(2*m_e*V0*eV)
	
	RS = lambda z : np.sqrt((z0/z)**2 -1)
	RS_asym = lambda z : -1/(np.sqrt((z0/z)**2 -1))
	LS = lambda z : np.tan(z)

	z_tan = np.linspace(0.1,4.5, 1000)
	z = np.linspace(0.1,4.05, 1000)
	for i in range(len(z_tan)):
		if np.absolute(LS(z_tan[i])) > 9:
			z_tan[i] = np.nan

	f = lambda z : RS(z) - LS(z)
	f_asym = lambda z : RS_asym(z) - LS(z)
	
	z_zeros = np.array([1.2, 3.5])
	z_zeros = fsolve(f, z_zeros)

	z_zeros_asym = 2
	z_zeros_asym = fsolve(f_asym, z_zeros_asym)
	
	z_zeros = np.array([z_zeros[0], z_zeros_asym[0], z_zeros[1]])

	E = lambda z : ((z*hbar)**2)/((a**2)*2*m_e*eV) - V0 # Okay, still not perf. div by a**2 arbitrarily to get the energies to match, but someone should figure this out later.
	for i in range(len(z_zeros)):
		print("Energi #" + str(i+1) + ": ", E(z_zeros[i]))

	energies,_,_,_,_,_ = analyze(1)
	for i in range(3):
		print("Num energies: ", energies[i])

	plt.ylim(-5, 5)
	plt.plot(z_tan, LS(z_tan))
	plt.plot(z, RS(z))
	plt.plot(z, RS_asym(z))
	#plt.savefig("trancend_plot.pdf")
	plt.show()
	
	#return z_zeros


plot_well_wave_funcs(0, 10)
plot_well_wave_funcs(1)
plot_well_wave_funcs(10)
plot_well_wave_funcs(50)
trancendent_sol()
band_widths(2, 50)
