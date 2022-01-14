###############################################
def setupMatplotlib(height=8.,width=6.):
	set1=brew.get_map('Set1','Qualitative',8).mpl_colors
	dark=brew.get_map('Dark2','Qualitative',8).mpl_colors
	paired=brew.get_map('Paired','Qualitative',12).mpl_colors
	reds=brew.get_map('Reds','Sequential',8).mpl_colors
	blues=brew.get_map('Blues','Sequential',9,reverse='True').mpl_colors
	spec=brew.get_map('Spectral','Diverging',8).mpl_colors
	plt.rcParams['xtick.direction']='out'
	plt.rcParams['ytick.direction']='out'
	plt.rcParams['lines.linewidth']= 2.0
	plt.rcParams['lines.color']= 'black'
	plt.rcParams['legend.frameon']=True
	plt.rcParams['font.family'] = 'serif'
	plt.rcParams['legend.fontsize']=10
	plt.rcParams['font.size'] = 12
	plt.rcParams['axes.color_cycle']=set1
	# For ipython notebook display set default values.
	#plt.rcParams['lines.markersize'] = 12
	plt.rcParams['figure.figsize'] = (height,width)
	plt.rcParams['grid.linewidth'] = 1

	# General settings used by display and print contexts.
	plt.rcParams['axes.axisbelow'] = True
	grid_line_color = '0.5'
	plt.rcParams['grid.color'] = grid_line_color
	plt.rcParams['grid.linestyle'] = '-'

###############################################
def commonFormat(ax_el,centerx=None,centery=None,location = 'lower right'):
	#ax_el.set_xlim(0,0.08)
	#ax_el.grid(True)
	#nur einfache Achsen, kein Rahmen
	ax_el.spines['top'].set_visible(0)
	ax_el.spines['right'].set_visible(0)
	ax_el.spines['bottom'].set_linewidth(0.5)
	ax_el.spines['left'].set_linewidth(0.5)
	ax_el.xaxis.set_ticks_position('bottom')
	ax_el.yaxis.set_ticks_position('left')
	if ((centerx is not None) and (centery is not None)):
		ax_el.spines['left'].set_position(('data', centerx))
		ax_el.spines['bottom'].set_position(('data', centery))
		ax_el.spines['right'].set_position(('data', centerx - 1))
		ax_el.spines['top'].set_position(('data', centery - 1))
	# Shink current axis's height by 10% on the bottom
	ax_el.legend(loc = location)#='lower right')
	#box = ax_el.get_position()
	#ax_el.set_position([box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.9])
	# Put a legend below current axis
	#ax_el.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=2)

def plot_all(sol,name,index=0):
	#Plot results
	setupMatplotlib(10.,4.)
	fig, ax = plt.subplots(nrows=1,ncols=2)
	plt.subplots_adjust(wspace=0.5)
	sol_elastic_strain = -0.5*(2.*sol.get_strain().T[0]+1.0-sol.get_inelastic_tens().T[0,2])
	#print(0.5*(sol.get_inelastic_tens().T[0,2]-1.0)-sol.get_strain().T[0])
	
	#stress-strain
	#ax[0].set_xlim(0.0, 0.20)
	#ax[0].set_ylim(0.10, 40)
	#print(E_11_ana-sol.get_strain().T[0])
	ax[0].plot(-exp_strain11,exp_stress,ls='-',label='Experimental',color='r')#,marker = 'o')
	ax[0].plot(-sol.get_strain().T[0],-sol.get_stress().T[0],label='$C_1$ = 10.0 ',ls='--',color='b')
	#ax[0].plot(-sol2.get_strain().T[0],-sol2.get_stress().T[0],label='$C_1$ = 20.0',ls='--',color='g')
	#ax[0].plot(-sol3.get_strain().T[0],-sol3.get_stress().T[0],label='$C_1$ = 30.0',ls='--',color='b')
	commonFormat(ax[0],0.,0.,'upper left')
	ax[0].set_xlabel('$E_{11}$')
	ax[0].set_ylabel('$S_{11}$ (kPa)')
	#print(len(sol.get_strain()))
	# time-stress
	#ax[1].set_xlim(0.0, 16000.)
	#ax[1].set_ylim(0.10, 40)
	ax[1].plot(exp_time,exp_stress,ls='-',label='Experimental',color='r')#,marker = 'o')
	ax[1].plot(sol.get_time(),-sol.get_stress().T[0],ls='--',label='$C_1$ = 10.0',color='b')
	#ax[1].plot(sol2.get_time(),-sol2.get_stress().T[0],ls='--',label='$C_1$ = 20.0',color='g')
	#ax[1].plot(sol3.get_time(),-sol3.get_stress().T[0],ls='--',label='$C_1$ = 30.0',color='b')
	commonFormat(ax[1],0.,0.,'upper right')
	ax[1].set_xlabel('$Time$ (min)')
	ax[1].set_ylabel('$S_{11}$ (kPa)')
	#
	#time-strain
	#ax[2].set_xlim(0.0, 1000.)
	#ax[2].set_ylim(0.0, -0.05)
	#ax[2].plot(exp_time,-exp_strain,ls='-',label='Experimental',color='b')#,marker = 'o')
	#ax[2].plot(sol.get_time(),-sol.get_strain().T[0],ls='-',label='C_1 = 0.5',color='r')
	#ax[2].plot(sol2.get_time(),-sol2.get_strain().T[0],ls='--',label='C_1 = 3.0',color='g')
	#ax[2].plot(sol3.get_time(),-sol3.get_strain().T[0],ls='-.',label='C_1 = 5.0',color='b')
	#commonFormat(ax[2],0.,0.,'upper right')
	#ax[2].set_xlabel('$Time (min)$')
	#ax[2].set_ylabel('$E_{zz}$')
	
	# elastic strain - plastic strain
	#ax[2].set_xlim(0.0, 1000.)
	#ax[2].set_ylim(0.0, -0.05)
	#ax[2].plot(exp_time,-exp_strain,ls='-',label='Experimental',color='b')#,marker = 'o')
	#ax[2].plot(-sol.get_strain().T[0],sol_elastic_strain,ls='-',label='C_1 = 0.5',color='r')
	#ax[2].plot(sol2.get_time(),-sol2.get_strain().T[0],ls='--',label='C_1 = 3.0',color='g')
	#ax[2].plot(sol3.get_time(),-sol3.get_strain().T[0],ls='-.',label='C_1 = 5.0',color='b')
	#commonFormat(ax[2],0.,0.,'upper right')
	#ax[2].set_xlabel('$E_{11}$')
	#ax[2].set_ylabel('$E_{11e}$')
	
	
	#print(Sig_11_ana-sol.get_stress().T[0])
	#print(np.min(Sig_11_ana - sol.get_stress().T[0]))
	#plastic strain-time
	#ax[3].plot(sol.get_time(),sol.get_inelastic_tens()[:,index][:,0] + sol.get_inelastic_scal().T[0]/3.,label='xx pl')\

	#ax[3].plot(sol.get_time(),sol.get_inelastic_tens()[:,index][:,1] + sol.get_inelastic_scal().T[0]/3.,label='yy pl')
	#ax[3].plot(sol.get_time(),sol.get_inelastic_tens()[:,index][:,2] + sol.get_inelastic_scal().T[0]/3.,label='zz pl')
	#ax[3].plot(sol.get_time(),sol.get_inelastic_tens()[:,index][:,3],label='xy pl')
	#ax[3].plot(sol.get_time(),sol.get_inelastic_tens()[:,index][:,4],label='yz pl')
	#ax[3].plot(sol.get_time(),sol.get_inelastic_tens()[:,index][:,5],label='xz pl')
	#ax[3].plot(sol.get_time(),sol.get_inelastic_scal().T[0],label='e_vol')
	#commonFormat(ax[3],0.,0.)
	#ax[3].set_xlabel('time')
	#ax[3].set_ylabel('plastic strain')
	fig.tight_layout()
	fig.savefig(name)
	return None


#########################################
#array to remember which positions are stress (or strain) controlled and (do not) require global iteration
def get_control_positions(loading,string):
	res = np.empty(0)
	for i in range(len(loading)):
		if (loading[i] == string):
			res = np.append(res,i)
	return res

def calc_tangent(dt,eps_t,eps_tdt,sig_t,deps_dt_t,hist_tens_t,hist_scal_t,material):
	pert = 1.e-8
	Jac_num = np.zeros((6,6))
	for i in range(6):
		for j in range(6):
			upper = 1.*eps_tdt
			lower = 1.*eps_tdt
			upper[j] += pert
			lower[j] -= pert
			res_up = material.calc_new_stress(dt,eps_t,upper,sig_t,deps_dt_t,hist_tens_t,hist_scal_t)
			res_low = material.calc_new_stress(dt,eps_t,lower,sig_t,deps_dt_t,hist_tens_t,hist_scal_t)
			up = res_up[i]
			low = res_low[i]
			Jac_num[i,j] = (up-low)/(2.*pert)
	return Jac_num

def convergence_check(residual,increment,global_tolerance):
	error_residual = np.linalg.norm(residual)
	error_increment = np.linalg.norm(increment)
	if (np.maximum(error_increment,error_residual) > global_tolerance):
		return False
	else:
		return True

def global_Newton(load,material,name="solution"):
	#for boundary conditions/global NR algorithm:
	loading_scheme = load[0][3]
	#set up solution tracker
	sol = solution_3D(material.get_number_of_internal_tensors(),material.get_number_of_internal_scalars(),name)
	pos_e = get_control_positions(loading_scheme,'e')
	pos_s = get_control_positions(loading_scheme,'s')
	#Kelvin mapping operator
	Kelv = np.array([1.,1.,1.,np.sqrt(2.),np.sqrt(2.),np.sqrt(2.)])
	total_time = 0
	#Loop over all load steps
	for i in range(1,len(load)):
		dt = load[i][1]
		local_time = 0.
		#Convert Load to Kelvin mapped coordinates
		load[i][2]*=Kelv
		#determine load increment vector
		incr = np.array((load[i][2] - load[i-1][2])/(load[i][0]) * dt)
		#Loop over all time increments in a load case
		while ((load[i][0] - local_time) >= DBL_EPS):
			#increment time measures
			local_time += dt
			total_time += dt
			#internal stress from previous increment
			sig_int = np.array(sol.get_stress()[len(sol.get_stress())-1]) #np.array necessary, otherwise address (overwrites values)
			sig_t = np.array(sol.get_stress()[len(sol.get_stress())-1])
			sig_ext = np.array(sig_int)#Attention: critical if previous increment not converged
			#get previous strain
			eps_t = np.array(sol.get_strain()[len(sol.get_strain())-1])
			eps_tdt = np.array(sol.get_strain()[len(sol.get_strain())-1])
			deps_dt_t = np.array(sol.get_strain_rate()[len(sol.get_strain_rate())-1])
			hist_tens_t = np.array(sol.get_inelastic_tens()[len(sol.get_inelastic_tens())-1])
			hist_scal_t = np.array(sol.get_inelastic_scal()[len(sol.get_inelastic_scal())-1])
			for j in range(6):
				if (load[i][3][j] == 'e'):
					eps_tdt[j] = eps_t[j] + incr[j] #put in bc and determine total strain
				else:
					sig_ext[j] += incr[j] ##external stress increased by driving force
			#perform trial stress integration
			sig_int = material.calc_new_stress(dt,eps_t,eps_tdt,sig_t,deps_dt_t,hist_tens_t,hist_scal_t)
			#determine initial residual
			residual = sig_ext-sig_int
			residual = np.delete(residual,pos_e,0)
			eps_inc = 1.*residual
			iter = 0
			#perform global (equilibrium) iteration
			while (not (convergence_check(residual,eps_inc,1.e-12)) and iter < ITER_MAX):
				iter += 1
				stiffn = material.get_dsig_dE()
				#stiffn_num = calc_tangent(dt,eps_t,eps_tdt,sig_t,deps_dt_t,hist_tens_t,hist_scal_t,material)
				#if (np.abs(stiffn_num - stiffn) > 1.e-4).any():
					#print("Deviation larger than tolerance detected.\n Maximum deviation between GLOBAL Jac_num and Jac_analyt is: \n")
					#print(np.max(np.abs(stiffn_num - stiffn)))
					#print(stiffn_num - stiffn)
				#print("Global Jacobian:")
				#print("Maximum absolute deviation ", np.max(np.abs(stiffn_num - stiffn)))
				#print("Maximum relative deviation ", np.nanmax(np.abs(np.divide(stiffn_num - stiffn,stiffn))))
				#print(stiffn)
				#build global Jacobian out of stiffness matrix
				#For Newton it seems we don't need to modify the RHS since the corresponding entry in
				#the solution vector is zero (no increment, since we have the solution already)
				K = np.delete(np.delete(stiffn,pos_e,1),pos_e,0)
				#get increment
				eps_inc = np.linalg.solve(K,residual)
				#increment strain only where no bc applied

				#scaling = np.minimum(0.0001/np.max(eps_inc),1.)
				#print(scaling)
				for j in range(len(eps_inc)):
						eps_tdt[pos_s[j]] += eps_inc[j]#*scaling
				#perform stress integration
				sig_int = material.calc_new_stress(dt,eps_t,eps_tdt,sig_t,deps_dt_t,hist_tens_t,hist_scal_t)
				#print("numerical global\n",calc_tangent(dt,eps_t,eps_tdt,sig_t,deps_dt_t,hist_tens_t,hist_scal_t,material))
				#sig_int = material.calc_new_stress(dt,eps_t,eps_tdt,sig_t,deps_dt_t,hist_tens_t,hist_scal_t)
				#print("extracted global\n",material.get_dsig_dE())
				#sys.exit()
				#compute residual
				residual = sig_ext-sig_int
				residual = np.delete(residual,pos_e,0)
			#append solution upon convergence
			#print(iter, " global iterations\n")#, sig_ext, " external stress\n", sig_int, " internal stress\n")
			sol.append_state(total_time, np.array([eps_tdt]), np.array([(eps_tdt-eps_t)/dt]), np.array([sig_int]),
					material.get_new_inelastic_tensors(),material.get_new_inelastic_scalars())

	#reverse Kelvin mapping and go to "normal" tensor coordinates.
	sol.Kelvin_to_normal(material.get_number_of_internal_tensors())
	return sol

def global_Newton_delta(load,material,name="solution_delta"):
	#for boundary conditions/global NR algorithm:
	loading_scheme = load[0][3]
	#set up solution tracker
	sol = solution_3D(material.get_number_of_internal_tensors(),material.get_number_of_internal_scalars(),name,True) # true value for internal tensor, added in viscoelastic model
	pos_e = get_control_positions(loading_scheme,'e')
	pos_s = get_control_positions(loading_scheme,'s')
	#Kelvin mapping operator
	Kelv = np.array([1.,1.,1.,np.sqrt(2.),np.sqrt(2.),np.sqrt(2.)])
	total_time = 0
	#Loop over all load steps
	for i in range(1,len(load)):
		dt = load[i][1]
		local_time = 0.
		#Convert Load to Kelvin mapped coordinates
		load[i][2]*=Kelv
		#Get new loading scheme
		pos_e = get_control_positions(load[i][3],'e')
		pos_s = get_control_positions(load[i][3],'s')
		#determine load increment vector
		incr = np.array(load[i][2]/load[i][0] * dt)
		#Loop over all time increments in a load case
		while ((load[i][0] - local_time) >= DBL_EPS):
			#increment time measures
			local_time += dt
			total_time += dt
			#internal stress from previous increment
			sig_int = np.array(sol.get_stress()[len(sol.get_stress())-1]) #np.array necessary, otherwise address (overwrites values)
			sig_t = np.array(sol.get_stress()[len(sol.get_stress())-1])
			sig_ext = np.array(sig_int)#Attention: critical if previous increment not converged
			#get previous strain
			eps_t = np.array(sol.get_strain()[len(sol.get_strain())-1])
			eps_tdt = np.array(sol.get_strain()[len(sol.get_strain())-1])
			deps_dt_t = np.array(sol.get_strain_rate()[len(sol.get_strain_rate())-1])
			hist_tens_t = np.array(sol.get_inelastic_tens()[len(sol.get_inelastic_tens())-1])
			hist_scal_t = np.array(sol.get_inelastic_scal()[len(sol.get_inelastic_scal())-1])
			for j in range(6):
				if (load[i][3][j] == 'e'):
					eps_tdt[j] = eps_t[j] + incr[j] #put in bc and determine total strain
				else:
					sig_ext[j] += incr[j] #external stress increased by driving force
			#perform trial stress integration
			sig_int = material.calc_new_stress(dt,eps_t,eps_tdt,sig_t,deps_dt_t,hist_tens_t,hist_scal_t)
			#determine initial residual
			residual = sig_ext-sig_int
			residual = np.delete(residual,pos_e,0)
			eps_inc = 1.*residual
			iter = 0
			#perform global (equilibrium) iteration
			while (not (convergence_check(residual,eps_inc,1.e-12)) and iter < ITER_MAX):
				iter += 1
				stiffn = material.get_dsig_dE()
				#build global Jacobian out of stiffness matrix
				#For Newton it seems we don't need to modify the RHS since the corresponding entry in
				#the solution vector is zero (no increment, since we have the solution already)
				K = np.delete(np.delete(stiffn,pos_e,1),pos_e,0)
				#get increment
				eps_inc = np.linalg.solve(K,residual)
				#increment strain only where no bc applied
				for j in range(len(eps_inc)):
						eps_tdt[pos_s[j]] += eps_inc[j]
				#perform stress integration
				sig_int = material.calc_new_stress(dt,eps_t,eps_tdt,sig_t,deps_dt_t,hist_tens_t,hist_scal_t)
				#compute residual
				residual = sig_ext-sig_int
				residual = np.delete(residual,pos_e,0)
			#append solution upon convergence
			#print(iter, " global iterations\n")#, sig_ext, " external stress\n", sig_int, " internal stress\n")
			sol.append_state(total_time, np.array([eps_tdt]), np.array([(eps_tdt-eps_t)/dt]), np.array([sig_int]),
					material.get_new_inelastic_tensors(),material.get_new_inelastic_scalars())
		load[i][2]/=Kelv #This brings the GL strain kelvin mapping back to normal

	#print('self.__epsv',float(material.get_new_inelastic_tensors().T[0]))
	#reverse Kelvin mapping and go to "normal" tensor coordinates.
	sol.Kelvin_to_normal(material.get_number_of_internal_tensors())
	return sol


def write_solution(sol,name='material_results.txt'):
	file = open('material_results.txt','w')
	file.write("Time,eps_xx,eps_yy,eps_zz,eps_xy,eps_yz,eps_xz,sig_xx,sig_yy,sig_zz,sig_xy,sig_yz,sig_xz\n")
	for i in range(len(sol.get_time())):
		string = str(sol.get_time()[i]) + ' ' + ' '.join(map(str,sol.get_strain()[i])) + ' ' + ' '.join(map(str,sol.get_stress()[i])) + "\n"
		file.write(string.replace(' ',','))
	file.close()
	return None

def set_material(C1,D2,alpha,C1v,D2v,alphav,eta,C1v2,D2v2,alphav2,eta2,C1p,D2p,alphap,cp):
	#choose material
	material = viscoelastoplastic(C1,D2,alpha,C1v,D2v,alphav,eta,C1v2,D2v2,alphav2,eta2,C1p,D2p,alphap,cp)#Parameters in MPa
	return material
	
def set_material_ana(stretch,C1,vStretch,C1v,D2v,eta,vStretch2,C1v2,D2v2,eta2):
	material_ana = ana_viscoelastic(stretch,C1,vStretch,C1v,D2v,eta,vStretch2,C1v2,D2v2,eta2)
	return material_ana

def set_MC(E):
	nu = 0.17 #guessed
	G = E/(2.*(1.+nu))
	K = E/(3.*(1.-2.*nu))
	
	return material
	
	
	
###########################################################################################################

# objective function for material model fit

def objective(params, load, exp):
	#ld = lc.exp_ld()
	ld = lc.experiment_ld()
	mat = viscoelastoplastic(params['C1'].value,10.,params['alpha'].value,0.001,0.001,0.,0.001,0.001,0.001,0.,0.001,params['C1p'].value,60.,params['alphap'].value,params['cp'].value)
	mod = global_Newton_delta(ld,mat)
	plot_all(mod,"fit_model.pdf")
	f = mod.get_stress().T[0] - exp
	return f
	
def test_fit():
	#load = lc.exp_ld()
	load = lc.experiment_ld()
	stress_exp = exp.get_stress()
	print(stress_exp.shape)
	#plot_all(exp,"test_fit_experiment.pdf")
	params = Parameters()
	#	(Name, Value, Vary, Min, Max, Expr)
	params.add_many(('C1', 3.0, True, None, None, None),
			('alpha', 0.0, True, None, None, None),
			('C1p', 0.03, True, None, None, None),
			('alphap', 0., True, None, None, None),
			('cp',0.04, True, None, None, None))
			#('D2', 20., True, None, None, None),
			#('alpha', 0.2, True, None, None, None),
			#('C1v', 1.0, True, None, None, None),
			#('D2v', 1., True, None, None, None),
			#('alphav', 0.2, True, None, None, None),
			#('eta1',0.001, True, None, None, None),
			#('C1v2', 1.0, True, None, None, None),
			#('D2v2', 1., True, None, None, None),
			#('alphav2', 0.2, True, None, None, None),
			#('eta2',0.001, True, None, None, None),
			#('C1p', 10.0, True, None, None, None),
			#('D2p', 20., True, None, None, None),
			#('alphap', 0.2, True, None, None, None),
			#('cp',50.0, True, None, None, None))
	# do fit with leastsq
	result = minimize(objective, params, args=(load,stress_exp),method='leastsq')
	# calculate final result
	final = result.residual
	report_fit(params)
	return result
			

###########################################################################################################
#main
import numpy as np
import matplotlib.pyplot as plt
import brewer2mpl as brew
import sys
#from lmfit import minimize, Parameters, Parameter, report_fit

#If classes don't work, create an empty __init__.py in the subdirectory
from Classes.solution_3D import solution_3D #import solution class (time, strain and stress series)
from Classes.Minkley_CS import viscoelastoplastic #Minkley with corner smoothing, scaled residuals
#from Classes.analy import ana_viscoelastic
from Classes.exp import experiment as exp


exp_stress = exp.get_stress()
exp_strain11 = exp.get_strain_11()
exp_time = exp.get_time()


setupMatplotlib()

#Close all figures
plt.close('all')

#Global variables
DBL_EPS = np.finfo(np.float).eps * 1.e4 #double tolerance
ITER_MAX = 50

#choose load case
import load_cases as lc
#load1 = lc.relax()
#load3 = lc.ramp_relax() #
#load1 = lc.relax(1.) #ramp_unload_relax, the first load case
#load2 = lc.relax(2.)
#load3 = lc.relax(3.)
ld = lc.experiment_uld()


#load = lc.loading_unloading(1.,0.85,-4.45614e-5,4.45614e-5)
#load2 = lc.loading_unloading(1.,0.85,-4.45614e-4,4.45614e-4)
#load1 = lc.relax_steps(1.,0.80,10., 0.1,-0.06) #ramp_unload_relax, the first load case
#load2 = lc.relax_steps(1.,0.90,10., 0.1,-0.06)
#load3 = lc.relax_steps(1.,0.95,10., 0.1,-0.06)
#
#material = set_material(0.04,0.,0.,0.05,0.,0.,0.08,0.05,0.,0.,0.09)
#material1 = set_material(1.0,0.,0.,0.5,0.5,0.,25.,0.5,0.,0.,5.0)
#material2 = set_material(1.0,0.,0.,0.5,1.0,0.,25.,0.5,0.,0.,5.0)
#material3 = set_material(1.0,0.,0.,0.5,2.0,0.,25.,0.5,0.,0.,5.0)
#material3 = set_material(0.04,0.,0.,2.0,0.,0.,0.08,0.06,0.,0.,0.09)
# material test for viscoelastoplastic
#material1 = set_material(1.0,2.0,0.5,1.0,0.5,0.5,2.0,1.,0.5,0.5,4.0,1.,0.5,0.0,10.)
material1 = set_material(9.0,500.0,0.0,  0.0001,0.0001,0.0,0.001,0.0001,0.0001,0.0,0.001, 70., 500., 00., 0.1)
#material2 = set_material(15.,4000.0,0.,  0.0001,0.0001,0.0,0.001,0.0001,0.0001,0.0,0.001, 80., 1000.,  0., 0.08)
#material3 = set_material(20.,4000.0,0.,  0.0001,0.0001,0.0,0.001,0.0001,0.0001,0.0,0.001, 80., 1000.,  0., 0.08)



#v_num = []
#v_num += [material.get_new_inelastic_tensors()]
#print(v_num)

################################################################################################
# analytical solution #

#t_ana = np.linspace(0,10.,1010)
#stretch_ana = 0.95
#material_ana = set_material_ana(stretch_ana,0.04,1.0,0.05,0.,0.08,1.,0.06,0.,0.09)
#E_11_ana = material_ana.GL_11(t_ana,stretch_ana)
#print('analytical', len(E_11_ana)) # debugging
#Sig_11_ana = material_ana.pk_11(t_ana,stretch_ana,0.05,0.08,0.06,0.09)



###############################################################

#sol = global_Newton_delta(load4, material)
#
sol1 = global_Newton_delta(ld,material1)
#sol2 = global_Newton_delta(ld,material2)
#sol3 = global_Newton_delta(ld,material3)
#
#sol4 = global_Newton_delta(load4,material)
# viscoelastoplastic
#sol = global_Newton_delta(load,material1)
#sol2 = global_Newton_delta(load2,material1)

plot_all(sol1,"ch6_para_ep_exp_uld.pdf",1)
#plot_all(sol1,sol2,"Para_two maxwells.pdf",1) # nramp = 4, nunld = 4; 4 holdings for ramping and 4 holdings for unloading
#plot_all(sol,sol2,"3rdtrial.pdf",1)
#plot_all(sol,sol2,"test.pdf",1)

##################################################################################################




#fit = test_fit()
















