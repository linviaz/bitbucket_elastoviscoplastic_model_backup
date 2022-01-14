#

import numpy as np
DBL_EPS = np.finfo(np.float).eps #* 1.0e4 #double tolerance
import sys

class viscoelastoplastic: 
	def __init__(self,Shear=1., Bulk=1., Shape= 0., vShear = 1., vBulk = 1., vShape =0.5, viscosity =1.,vShear2 = 1., vBulk2 = 1., vShape2 =0.5, viscosity2 =1., pShear = 1., pBulk = 1., pShape =0.5, pMultiplier =1.):
		# modified neo-hooke
		self.__C1 = Shear
		self.__D2 = Bulk
		self.__alpha = Shape
		# 1st Maxwell
		self.__C1v = vShear # C1v parameter for the viscosity 
		self.__D2v = vBulk #D2v parameter for the voscisity
		self.__eta = viscosity #viscosity parameter eta_v
		self.__alphav = vShape #shape change parameter for the viscous free Helmholtz energy function tensor
		# 2nd Maxwell
		self.__C1v2 = vShear2 #number 2 indicates the second Maxwell element
		self.__D2v2 = vBulk2
		self.__eta2 = viscosity2 
		self.__alphav2 = vShape2
		# elasto-plastic
		self.__C1p = pShear
		self.__D2p = pBulk
		self.__alphap = pShape
		self.__cp = pMultiplier
		#
		#
		self.__PD = self.__build_PD() #projection for deviatoric part
		self.__PS = self.__build_PS() #projection for spherical part
		self.__ivec = np.array([1.0,1.0,1.0,0.0,0.0,0.0]) #6D vector of unit tensor
		self.__ident = np.identity(6) #6x6 identity matrix
		#
		self.__dsig_dE = self.__build_D(self.__C1,self.__D2)
		self.__n_inel_tens = 3 #number of 6D state vectors
		self.__epsv = np.array([[[1.,1.,1.,0.,0.,0.]]*self.__n_inel_tens]) #viscoelastic strain tensor C_v and plastic strain tensor C_p
		self.__n_inel_scal = 0 #number of state scalars
		self.__int_tol = 1.0e-14
		self.__max_it = 20
		np.set_printoptions(precision=2,linewidth=250,suppress=False) #for debugging
	
	def get_number_of_internal_tensors(self):
		return self.__n_inel_tens
	
	def get_number_of_internal_scalars(self):
		return self.__n_inel_scal
	
	def calc_new_stress(self,dt,eps_t0,eps_t1,sig_t0,deps_dt_t0,int_tens,int_scal):
		#print(eps_t1)
		CGv_i = 1.0*int_tens[0]
		#Cv_11 = np.array(CGv_i[0]) #for debug
		CGv_i2 = 1.0*int_tens[1]
		##print(Cv_11)
		#print(int_tens[0])
		CGp_i = 1.0*int_tens[2]
		#
		sig_i = self.__viscoelastoplastic(eps_t1,CGv_i,CGv_i2,CGp_i)
		#sig_i = 1.* sig_t0
		#print(sig_t0, sig_i)
		#sys.exit()
		#first evaluation
		#print(CGv_i[0])
		###########visco##################		
		#
		#print(self.residual(dt,np.array([0.]*6),np.array([0.]*6),np.array([1.,1.,1.,0.,0.,0.]),np.array([1.,1.,1.,0.,0.,0.])))
		res = self.residual(dt,eps_t0,eps_t1,sig_i,int_tens,CGv_i,CGv_i2,CGp_i)
		#print(res)
		K = self.jacobian(dt,eps_t0,eps_t1,int_tens,CGv_i,CGv_i2,CGp_i)
		counter = 0
		delta = res #only for convergence check
		#local loop
		while ((np.linalg.norm(res) > self.__int_tol) and (np.linalg.norm(delta) > self.__int_tol) and (counter < self.__max_it)):
			counter += 1
			#compute Jacobian
			K = self.jacobian(dt,eps_t0,eps_t1,int_tens,CGv_i,CGv_i2,CGp_i)
			#K = self.num_jacobian(dt,eps_t1,sig_i,int_tens,CGv_i)
			#self.test_jacobian(dt,eps_t0,eps_t1,sig_i,int_tens,CGv_i,CGv_i2,CGp_i)
			#compare
			#solve for increment
			delta = np.linalg.solve(K,-res)
			#update variables
			#print(K-self.num_jacobian(dt,eps_t1,sig_i,int_tens,CGv_i))
			sig_i += np.array(delta[0:6])
			CGv_i += np.array(delta[6:12])
			CGv_i2 += np.array(delta[12:18])
			CGp_i += np.array(delta[18:24])
			#print(Cv_11, counter)
			#
			#compute new residual
			res = self.residual(dt,eps_t0,eps_t1,sig_i,int_tens,CGv_i,CGv_i2,CGp_i)
		#
		#find consistent algorithmic tangent
		CG = self.__CG(eps_t1)
		CGt = self.__CG(eps_t0)
		CGi = self.__inverse(CG)
		#
		Cvi = self.__inverse(CGv_i)
		C_ev = self.__C_ev(CG,CGv_i)
		#
		Cvi2 = self.__inverse(CGv_i2)
		C_ev2 = self.__C_ev(CG,CGv_i2)
		#print('CGv_i numerical:',Cv_11,CGv_i[0])
		#
		Cpi = self.__inverse(CGp_i)
		C_ep = self.__C_ev(CG,CGp_i)
		#
		I3dPsi_dI3 = 2.*self.__D2*np.log(self.__I3(CG))-self.__C1*self.__expA(self.__alpha, CG)
		ddPsi_ddI1 = self.__C1*self.__expA(self.__alpha, CG)*self.__alpha
		I3ddPsi_dI1dI3 = - self.__C1*self.__expA(self.__alpha,CG)*self.__alpha
		I3I3ddPsi_ddI3 = self.__C1*(1.+self.__alpha)*self.__expA(self.__alpha, CG) + 2.*self.__D2*(1.- np.log(self.__I3(CG)))
		#
		vdPsi_dI1 = self.__C1v*self.__expA(self.__alphav,C_ev)
		vI3dPsi_dI3 = 2.*self.__D2v*np.log(self.__I3(C_ev))-self.__C1v*self.__expA(self.__alphav, C_ev)
		vddPsi_ddI1 = self.__C1v*self.__expA(self.__alphav, C_ev)*self.__alphav
		vI3ddPsi_dI1dI3 = - self.__C1v*self.__expA(self.__alphav,C_ev)*self.__alphav
		vI3I3ddPsi_ddI3 = self.__C1v*(1.+self.__alphav)*self.__expA(self.__alphav, C_ev) + 2.*self.__D2v*(1.- np.log(self.__I3(C_ev)))
		#
		v2dPsi_dI1 = self.__C1v2*self.__expA(self.__alphav2,C_ev2)
		v2I3dPsi_dI3 = 2.*self.__D2v2*np.log(self.__I3(C_ev2))-self.__C1v2*self.__expA(self.__alphav2, C_ev2)
		v2ddPsi_ddI1 = self.__C1v2*self.__expA(self.__alphav2, C_ev2)*self.__alphav2
		v2I3ddPsi_dI1dI3 = - self.__C1v2*self.__expA(self.__alphav2,C_ev2)*self.__alphav2
		v2I3I3ddPsi_ddI3 = self.__C1v2*(1.+self.__alphav2)*self.__expA(self.__alphav2, C_ev2) + 2.*self.__D2v2*(1.- np.log(self.__I3(C_ev2)))	
		#
		pdPsi_dI1 = self.__C1p*self.__expA(self.__alphap,C_ep)
		pI3dPsi_dI3 = 2.*self.__D2p*np.log(self.__I3(C_ep))-self.__C1p*self.__expA(self.__alphap, C_ep)
		pddPsi_ddI1 = self.__C1p*self.__expA(self.__alphap, C_ep)*self.__alphap
		pI3ddPsi_dI1dI3 = - self.__C1p*self.__expA(self.__alphap,C_ep)*self.__alphap
		pI3I3ddPsi_ddI3 = self.__C1p*(1.+self.__alphap)*self.__expA(self.__alphap, C_ep) + 2.*self.__D2p*(1.- np.log(self.__I3(C_ep)))	
		#
		dG_de_1 = -4./self.__C1*(ddPsi_ddI1*np.outer(self.__ivec,self.__ivec)+I3ddPsi_dI1dI3*(np.outer(CGi,self.__ivec)+np.outer(self.__ivec,CGi))+(I3I3ddPsi_ddI3+I3dPsi_dI3)*np.outer(CGi,CGi)- I3dPsi_dI3*self.__s_odot_s(CGi) + vddPsi_ddI1*np.outer(Cvi,Cvi) + vI3ddPsi_dI1dI3*(np.outer(CGi,Cvi)+np.outer(Cvi,CGi)) + (vI3I3ddPsi_ddI3+vI3dPsi_dI3)*np.outer(CGi,CGi) - 		vI3dPsi_dI3*self.__s_odot_s(CGi) + v2ddPsi_ddI1*np.outer(Cvi2,Cvi2) + v2I3ddPsi_dI1dI3*(np.outer(CGi,Cvi2)+np.outer(Cvi2,CGi)) + (v2I3I3ddPsi_ddI3+v2I3dPsi_dI3)*np.outer(CGi,CGi) - v2I3dPsi_dI3*self.__s_odot_s(CGi) + pddPsi_ddI1*np.outer(Cpi,Cpi) + pI3ddPsi_dI1dI3*(np.outer(CGi,Cpi)+np.outer(Cpi,CGi)) + (pI3I3ddPsi_ddI3+pI3dPsi_dI3)*np.outer(CGi,CGi) - 		pI3dPsi_dI3*self.__s_odot_s(CGi))
		dG_de_2 = -8./self.__eta*(vddPsi_ddI1*np.outer(CG,Cvi)+vI3ddPsi_dI1dI3*(np.outer(CG,CGi)+np.outer(CGv_i,Cvi))+vdPsi_dI1*self.__ident+(vI3I3ddPsi_ddI3+vI3dPsi_dI3)*np.outer(CGv_i,CGi))
		dG_de_3 = -8./self.__eta2*(v2ddPsi_ddI1*np.outer(CG,Cvi2)+v2I3ddPsi_dI1dI3*(np.outer(CG,CGi)+np.outer(CGv_i2,Cvi2))+v2dPsi_dI1*self.__ident+(v2I3I3ddPsi_ddI3+v2I3dPsi_dI3)*np.outer(CGv_i2,CGi))
		#print("Analytical\n",dG_de_1)
		#print("Numerical\n",self.check_jacobian(eps_t1))
		#print("the norm value",np.linalg.norm(CG-CGt))
		dG_de_4 = np.where(np.linalg.norm(CG-CGt)==0.,0.,-4.*self.__cp/(np.linalg.norm(CG-CGt)*dt)*np.outer((pdPsi_dI1*CG + pI3dPsi_dI3*CGp_i), (CG-CGt)) - 4.*self.__cp*np.linalg.norm(CG-CGt)/dt *(pddPsi_ddI1*np.outer(CG,Cpi)+pI3ddPsi_dI1dI3*(np.outer(CG,CGi)+np.outer(CGp_i,Cpi))+pdPsi_dI1*self.__ident+(pI3I3ddPsi_ddI3+pI3dPsi_dI3)*np.outer(CGp_i,CGi)))
		#
		dG_de = np.append(dG_de_1, np.append(dG_de_2,np.append(dG_de_3,dG_de_4,0),0),0)
		dz_de = np.linalg.solve(K,-1.0*dG_de)
		#print("Analytical and numerical tangent dr/de difference\n",dG_de - self.check_tangent(dt,eps_t0,eps_t1,sig_i,int_tens,CGv_i,CGv_i2,CGp_i))
		#print("Analytical and numerical tangent dr/de maximum difference\n",np.max(dG_de - self.check_tangent(dt,eps_t0,eps_t1,sig_i,int_tens,CGv_i,CGv_i2,CGp_i)))
		#print("Numerical\n",self.check_tangent(dt,eps_t1,sig_i,int_tens,CGv_i))
		#sys.exit()
		#extract tangent moduli and add hydrostatic extension
		self.__dsig_dE = self.__C1*dz_de[0:6]
		#print(self.__dsig_dE, self.__dsig_dE.shape)
		#print(counter, " local iterations")
		#update solution variables
		sig_i *= self.__C1 #redimensionalise
		#sys.exit()
		#print('sig_i in cal_new _stress',sig_i)
		self.__epsv = np.append([[CGv_i]],np.append([[CGv_i2]],[[CGp_i]],1),1)
		return sig_i

	#function to test analytical tangent
	def test_jacobian(self,dt,eps_t0,eps_t1,sig_i,int_tens,CGv_i,CGv_i2,CGp_i):
		num_tol = 1.e-8
		Jac_analyt = self.jacobian(dt,eps_t0,eps_t1,int_tens,CGv_i,CGv_i2,CGp_i)
		Jac_num = self.num_jacobian(dt,eps_t0,eps_t1,sig_i,int_tens,CGv_i,CGv_i2,CGp_i)
		if (np.abs(Jac_num - Jac_analyt) > num_tol).any():
			print("Deviation larger than tolerance detected.\n Maximum deviation between Jac_num and Jac_analyt is: \n")
			print(np.max(np.abs(Jac_num - Jac_analyt)))
			print(Jac_num - Jac_analyt)
			print('jacobian_num',Jac_num)
			print('analytical jacobian',Jac_analyt)
			sys.exit()
		else:
			print("Analytical and numerical viscoelastic Jacobians: Maximum deviation",np.max(np.abs(Jac_num - Jac_analyt)))
		return None
				
		

	#derivative check
	def num_jacobian(self,dt,eps_t0,strain,sig_i,int_tens,CGv_i,CGv_i2,CGp_i):
		res = np.zeros((24,24))
		eps = 1.e-8
		for i in range(24):
			for j in range(24):
				if (j<6):
					top = 1.*sig_i
					top[j] += eps
					bot = 1.*sig_i
					bot[j] -= eps	
					res[i,j] = (self.residual(dt,eps_t0,strain,top,int_tens,CGv_i,CGv_i2,CGp_i)[i] - self.residual(dt,eps_t0,strain,bot,int_tens,CGv_i,CGv_i2,CGp_i)[i])/(2.*eps)
				elif (j<12):
					top = 1.*CGv_i
					top[j-6] += eps
					bot = 1.*CGv_i
					bot[j-6] -= eps
					res[i,j] = (self.residual(dt,eps_t0,strain,sig_i,int_tens,top,CGv_i2,CGp_i)[i] - self.residual(dt,eps_t0,strain,sig_i,int_tens,bot,CGv_i2,CGp_i)[i])/(2.*eps)
				elif (j<18):
					top = 1.*CGv_i2
					top[j-12] += eps
					bot = 1.*CGv_i2
					bot[j-12] -= eps
					res[i,j] = (self.residual(dt,eps_t0,strain,sig_i,int_tens,CGv_i,top,CGp_i)[i] - self.residual(dt,eps_t0,strain,sig_i,int_tens,CGv_i,bot,CGp_i)[i])/(2.*eps)
				elif (j<24):
					top = 1.*CGp_i
					top[j-18] += eps
					bot = 1.*CGp_i
					bot[j-18] -= eps
					res[i,j] = (self.residual(dt,eps_t0,strain,sig_i,int_tens,CGv_i,CGv_i2,top)[i] - self.residual(dt,eps_t0,strain,sig_i,int_tens,CGv_i,CGv_i2,bot)[i])/(2.*eps)	
		return res

	# r(S, CGv_i, CGv_i2)
	def residual(self,dt,eps_t0,eps_t1,sig_i,int_tens,CGv_i, CGv_i2,CGp_i):
		cauchy_green = self.__CG(eps_t1)
		cauchy_green_t = self.__CG(eps_t0)
		#print(eps_t1)
		C_ev = self.__C_ev(cauchy_green,CGv_i) # for calculating C_ev invariants only
		C_ev2 = self.__C_ev(cauchy_green,CGv_i2) #second maxwell element
		#print(self.__I3(C_ev))
		C_ep = self.__C_ev(cauchy_green,CGp_i)
		#
		vdPsi_dI1 = self.__C1v*self.__expA(self.__alphav, C_ev)
		vI3dPsi_dI3 = self.__D2v*2.*np.log(self.__I3(C_ev))-self.__C1v*self.__expA(self.__alphav, C_ev)
		#
		v2dPsi_dI1 = self.__C1v2*self.__expA(self.__alphav2, C_ev2)
		v2I3dPsi_dI3 = self.__D2v2*2.*np.log(self.__I3(C_ev2))-self.__C1v2*self.__expA(self.__alphav2, C_ev2)
		#
		pdPsi_dI1 = self.__C1p*self.__expA(self.__alphap, C_ep)
		pI3dPsi_dI3 = self.__D2p*2.*np.log(self.__I3(C_ep))-self.__C1p*self.__expA(self.__alphap, C_ep)
		#
		G_sig = sig_i - self.__viscoelastoplastic(eps_t1,CGv_i, CGv_i2,CGp_i) 
		G_Cv = (CGv_i - int_tens[0])/dt - 4./self.__eta *(vdPsi_dI1*cauchy_green+vI3dPsi_dI3*CGv_i) # here is the reason of discrepency btwn analytical and numerical solutions
		G_Cv2 = (CGv_i2 - int_tens[1])/dt - 4./self.__eta2 *(v2dPsi_dI1*cauchy_green+v2I3dPsi_dI3*CGv_i2) #analytical use former step, numerical use initial value
		G_Cp = (CGp_i - int_tens[2])/dt - 2.*self.__cp*np.linalg.norm(cauchy_green - cauchy_green_t)/dt*(pdPsi_dI1*cauchy_green + pI3dPsi_dI3*CGp_i) # plastic RCG tensor rate residual
		result = np.append(G_sig, np.append(G_Cv, np.append(G_Cv2,G_Cp)))
		return result
	
	def check_tangent(self,dt,eps_t0,eps_t1,sig_i,int_tens,CGv_i,CGv_i2,CGp_i):
		res = np.zeros((24,6)) # of the 24, 6 are stress, 6 are cv_i, 6 are cv_i2, and 6 are cp
		eps = 1.e-8
		for i in range(24):
			for j in range(6):
				top = 1.*eps_t1
				top[j] += eps
				bot = 1.*eps_t1
				bot[j] -= eps
				res[i,j] = (self.residual(dt,eps_t0,top,sig_i,int_tens,CGv_i,CGv_i2,CGp_i)[i] - self.residual(dt,eps_t0,bot,sig_i,int_tens,CGv_i,CGv_i2,CGp_i)[i])/(2.*eps)
		return res
				
				

	def __viscoelastoplastic(self,eps_t1,CGv_i, CGv_i2,CGp_i): # output the dimensionless 2nd pk stress 
		cauchy_green = self.__CG(eps_t1)
		CGi = self.__inverse(cauchy_green)
		#
		C_ev = self.__C_ev(cauchy_green,CGv_i) 
		Cvi = self.__inverse(CGv_i)
		#
		C_ev2 = self.__C_ev(cauchy_green,CGv_i2) 
		Cvi2 = self.__inverse(CGv_i2)
		#print(CGv_i)
		#print(cauchy_green,"\n",self.__I3(cauchy_green))
		#
		C_ep = self.__C_ev(cauchy_green,CGp_i)
		Cpi = self.__inverse(CGp_i)
		#
		dPsi_dI1 = self.__C1*self.__expA(self.__alpha, cauchy_green)
		I3dPsi_dI3 = self.__D2*2.*np.log(self.__I3(cauchy_green))-self.__C1*self.__expA(self.__alpha, cauchy_green)
		#
		vdPsi_dI1 = self.__C1v*self.__expA(self.__alphav, C_ev)
		vI3dPsi_dI3 = self.__D2v*2.*np.log(self.__I3(C_ev))-self.__C1v*self.__expA(self.__alphav, C_ev)
		#
		v2dPsi_dI1 = self.__C1v2*self.__expA(self.__alphav2, C_ev2)
		v2I3dPsi_dI3 = self.__D2v2*2.*np.log(self.__I3(C_ev2))-self.__C1v2*self.__expA(self.__alphav2, C_ev2)
		#
		pdPsi_dI1 = self.__C1p*self.__expA(self.__alphap, C_ep)
		pI3dPsi_dI3 = self.__D2p*2.*np.log(self.__I3(C_ep))-self.__C1p*self.__expA(self.__alphap, C_ep)
		#
		return 2./self.__C1*(dPsi_dI1*self.__ivec + I3dPsi_dI3*CGi + vdPsi_dI1*Cvi + vI3dPsi_dI3*CGi + v2dPsi_dI1*Cvi2 + v2I3dPsi_dI3*CGi + pdPsi_dI1*Cpi + pI3dPsi_dI3*CGi) 
		

	def __CG(self,GL):
		return 2*GL + self.__ivec

	#viscoelastic Jacobian ####################################################
	def jacobian(self,dt,eps_t0,eps_t1,int_tens,CGv_i,CGv_i2,CGp_i):
		#
		CG = self.__CG(eps_t1)
		CGt = self.__CG(eps_t0)
		CGi = self.__inverse(CG)
		#
		Cvi = self.__inverse(CGv_i)
		C_ev = self.__C_ev(CG,CGv_i)
		#
		Cvi2 = self.__inverse(CGv_i2)
		C_ev2 = self.__C_ev(CG,CGv_i2)
		#
		Cpi = self.__inverse(CGp_i)
		C_ep = self.__C_ev(CG,CGp_i)
		#
		dPsi_dI1 = self.__C1*self.__expA(self.__alpha, CG)
		I3dPsi_dI3 = self.__D2*2.*np.log(self.__I3(CG))-self.__C1*self.__expA(self.__alpha, CG)
		#
		vdPsi_dI1 = self.__C1v*self.__expA(self.__alphav, C_ev)
		vI3dPsi_dI3 = self.__D2v*2.*np.log(self.__I3(C_ev))-self.__C1v*self.__expA(self.__alphav, C_ev)
		vddPsi_ddI1 = self.__C1v*self.__alphav*self.__expA(self.__alphav, C_ev)
		vddPsi_dI1dI3 = - self.__C1v*self.__alphav*self.__expA(self.__alphav,C_ev)/self.__I3(C_ev)
		vddPsi_ddI3 = self.__C1v*(1.+self.__alphav)*self.__expA(self.__alphav,C_ev)/((self.__I3(C_ev))**2.) + 2.*self.__D2v*(1.-np.log(self.__I3(C_ev)))/((self.__I3(C_ev))**2) 
		#
		v2dPsi_dI1 = self.__C1v2*self.__expA(self.__alphav2, C_ev2)
		v2I3dPsi_dI3 = self.__D2v2*2.*np.log(self.__I3(C_ev2))-self.__C1v2*self.__expA(self.__alphav2, C_ev2)
		v2ddPsi_ddI1 = self.__C1v2*self.__alphav2*self.__expA(self.__alphav2, C_ev2)
		v2ddPsi_dI1dI3 = - self.__C1v2*self.__alphav2*self.__expA(self.__alphav2,C_ev2)/self.__I3(C_ev2)
		v2ddPsi_ddI3 = self.__C1v2*(1.+self.__alphav2)*self.__expA(self.__alphav2,C_ev2)/((self.__I3(C_ev2))**2.) + 2.*self.__D2v2*(1.-np.log(self.__I3(C_ev2)))/((self.__I3(C_ev2))**2) 
		#
		pdPsi_dI1 = self.__C1p*self.__expA(self.__alphap, C_ep)
		pI3dPsi_dI3 = self.__D2p*2.*np.log(self.__I3(C_ep))-self.__C1p*self.__expA(self.__alphap, C_ep)
		pddPsi_ddI1 = self.__C1p*self.__alphap*self.__expA(self.__alphap, C_ep)
		pddPsi_dI1dI3 = - self.__C1p*self.__alphap*self.__expA(self.__alphap,C_ep)/self.__I3(C_ep)
		pddPsi_ddI3 = self.__C1p*(1.+self.__alphap)*self.__expA(self.__alphap,C_ep)/((self.__I3(C_ep))**2.) + 2.*self.__D2p*(1.-np.log(self.__I3(C_ep)))/((self.__I3(C_ep))**2) 
		#
		K_11 = self.__ident
		K_12 = 2./self.__C1*np.dot((vddPsi_ddI1*np.outer(Cvi,CG)+ vddPsi_dI1dI3*self.__I3(C_ev)*(np.outer(Cvi,CGv_i)+np.outer(CGi,CG))+(self.__I3(C_ev)*vddPsi_ddI3*self.__I3(C_ev)+ vI3dPsi_dI3)*np.outer(CGi,CGv_i)+vdPsi_dI1*self.__ident), self.__s_odot_s(Cvi)) #+vdPsi_dI1*self.__s_odot_s(Cvi)
		K_13 = 2./self.__C1*np.dot((v2ddPsi_ddI1*np.outer(Cvi2,CG)+ v2ddPsi_dI1dI3*self.__I3(C_ev2)*(np.outer(Cvi2,CGv_i2)+np.outer(CGi,CG))+(self.__I3(C_ev2)*v2ddPsi_ddI3*self.__I3(C_ev2)+ v2I3dPsi_dI3)*np.outer(CGi,CGv_i2)+v2dPsi_dI1*self.__ident), self.__s_odot_s(Cvi2))
		K_14 = 2./self.__C1*np.dot((pddPsi_ddI1*np.outer(Cpi,CG)+ pddPsi_dI1dI3*self.__I3(C_ep)*(np.outer(Cpi,CGp_i)+np.outer(CGi,CG))+(self.__I3(C_ep)*pddPsi_ddI3*self.__I3(C_ep)+ pI3dPsi_dI3)*np.outer(CGi,CGp_i)+pdPsi_dI1*self.__ident), self.__s_odot_s(Cpi))
		K_21 = np.zeros((6,6))
		K_22 = (1./dt - 4./self.__eta*vI3dPsi_dI3)*self.__ident + 4./self.__eta*np.dot((vddPsi_ddI1*np.outer(CG,CG)+self.__I3(C_ev)*vddPsi_dI1dI3*(np.outer(CG,CGv_i)+np.outer(CGv_i,CG))+ (vddPsi_ddI3*self.__I3(C_ev)*self.__I3(C_ev)+vI3dPsi_dI3)*np.outer(CGv_i,CGv_i)), self.__s_odot_s(Cvi))
		K_23 = np.zeros((6,6))
		K_24 = np.zeros((6,6))
		K_31 = np.zeros((6,6))
		K_32 = np.zeros((6,6))
		K_33 = (1./dt - 4./self.__eta2*v2I3dPsi_dI3)*self.__ident + 4./self.__eta2*np.dot((v2ddPsi_ddI1*np.outer(CG,CG)+self.__I3(C_ev2)*v2ddPsi_dI1dI3*(np.outer(CG,CGv_i2)+np.outer(CGv_i2,CG))+ (v2ddPsi_ddI3*self.__I3(C_ev2)*self.__I3(C_ev2)+v2I3dPsi_dI3)*np.outer(CGv_i2,CGv_i2)), self.__s_odot_s(Cvi2))
		K_34 = np.zeros((6,6))
		K_41 = np.zeros((6,6))
		K_42 = np.zeros((6,6))
		K_43 = np.zeros((6,6))
		K_44 = self.__ident*(1.-2.*self.__cp*np.linalg.norm(CG - CGt)*pI3dPsi_dI3)/dt + 2.*self.__cp*np.linalg.norm(CG-CGt)/dt*np.dot((pddPsi_ddI1*np.outer(CG,CG)+self.__I3(C_ep)*pddPsi_dI1dI3*(np.outer(CG,CGp_i)+np.outer(CGp_i,CG))+(pddPsi_ddI3*self.__I3(C_ep)*self.__I3(C_ep)+pI3dPsi_dI3)*np.outer(CGp_i,CGp_i) ), self.__s_odot_s(Cpi))
		#
		#print(K_11.shape,K_12.shape,K_21.shape,K_22.shape)
		#sys.exit()
		row1 = np.append(K_11,np.append(K_12,np.append(K_13,K_14,1),1),1)
		row2 = np.append(K_21,np.append(K_22,np.append(K_23,K_24,1),1),1)
		row3 = np.append(K_31,np.append(K_32,np.append(K_33,K_34,1),1),1)
		row4 = np.append(K_41,np.append(K_42,np.append(K_43,K_44,1),1),1)
		K_full = np.append(row1,np.append(row2,np.append(row3,row4,0),0),0)
		return K_full
	
	def __build_D(self, G, K):
		#Hilfskonstanten
		lam = (3.*K-2.*G)/3.
		#linear elastic moduli
		D = np.array([[lam+2.*G,lam,lam,0.,0.,0.],
			[lam,lam+2.*G,lam,0.,0.,0.],
			[lam,lam,lam+2.*G,0.,0.,0.],
			[0.,0.,0.,2.*G,0.,0.],
			[0.,0.,0.,0.,2.*G,0.],
			[0.,0.,0.,0.,0.,2.*G]])
		return D
	
	#projects 6D vector into vector of the deviatoric tensor coordinates
	def __build_PD(self):
		I = np.eye(3)
		N = np.zeros((3,3))
		fullI = np.ones((3,3))
		PD = np.append(np.append(I - 1./3.*fullI,N,1),np.append(N,I,1),0)
		return PD
	
	#projects 6D vector into vector of the spherical tensor coordinates
	def __build_PS(self):
		N = np.zeros((3,3))
		fullI = np.ones((3,3))
		PS = np.append(np.append(1./3.*fullI,N,1),np.append(N,N,1),0)
		return PS
	
	#dtermines trace of vector of tensor coordinates
	def __I1(self,vec):
		return vec[0]+vec[1]+vec[2]
	
	#determine effective stress
	def __sig_eff(self,vec):
		return np.sqrt(3.*self.__J2(vec))
	

	#return first invariant
	#def __I1(self,sig):
		#return self.__trace(sig)
	
	#return second deviatoric invariant
	def __J2(self,sig):
		sig_dev = np.dot(self.__PD,sig)
		return 0.5*np.dot(sig_dev,sig_dev)

	#return third deviatoric invariant
	def __J3(self,sig):
		sig_dev = np.dot(self.__PD,sig)
		return np.linalg.det(self.__vec_to_tens(sig_dev))

	def __I3(self,A):
		return np.linalg.det(self.__vec_to_tens(A))
########################################################################

	def __C_ev(self, CG, CGv_i): # this is only valid for calculating the invariants, not actual C_ev
		Cvi = self.__inverse(CGv_i)
		Cvi_tens = self.__vec_to_tens(Cvi)
		CG_tens = self.__vec_to_tens(CG)
		CCvi_tens = np.dot(CG_tens,Cvi_tens)
		C_ev = np.array([CCvi_tens[0,0],CCvi_tens[1,1],CCvi_tens[2,2],np.sqrt(2.)*CCvi_tens[0,1],\
				np.sqrt(2.)*CCvi_tens[1,2],np.sqrt(2.)*CCvi_tens[0,2]])
		return C_ev
	


	def __expA(self,alphav,vec):
		I1_ev = self.__I1(vec)
		I3_ev = self.__I3(vec)
		expA = np.exp(alphav*(I1_ev - np.log(I3_ev)-3.))
		return expA


########################################################################


	#invert a tensor (takes and returns the Kelvin vector)
	def __inverse(self,vec):
		tens = self.__vec_to_tens(vec)
		if (np.abs(np.linalg.det(tens)) < DBL_EPS):
			inverse = np.zeros((3,3))
		else:
			inverse = np.linalg.inv(tens)
		return np.array([inverse[0,0],inverse[1,1],inverse[2,2], np.sqrt(2.) * inverse[0,1], \
			np.sqrt(2.) * inverse[1,2], np.sqrt(2.) * inverse[0,2]])

	#finds principal values and return array sorted in ascending order
	def __eigenvalues(self,vec):
		tens = self.__vec_to_tens(vec)
		eig=np.linalg.eigvals(tens) #linalg.eigvals(tens): compute the eigenvalues of tens
		eig.sort()#sort ascending
		eig = eig[::-1] #reverse --> descending
		return eig

	
	#map vector back to tensor coordinate matrix (also removes Kelvin mapping factors)
	def __vec_to_tens(self,vec):
		return np.array([[vec[0],vec[3]/np.sqrt(2),vec[5]/np.sqrt(2)], \
			[vec[3]/np.sqrt(2),vec[1],vec[4]/np.sqrt(2)], \
			[vec[5]/np.sqrt(2),vec[4]/np.sqrt(2),vec[2]]])

	#construct \sigma^{-1} \odot \sigma^{-1} matrix
	def __s_odot_s(self,sigma_inv):
		#
		odot = np.zeros((6,6))
		odot[0,0] = sigma_inv[0]*sigma_inv[0]
		odot[0,1] = odot[1,0] = sigma_inv[3]*sigma_inv[3]/2.
		odot[0,2] = odot[2,0] = sigma_inv[5]*sigma_inv[5]/2.
		odot[0,3] = odot[3,0] = sigma_inv[0]*sigma_inv[3]
		odot[0,4] = odot[4,0] = sigma_inv[3]*sigma_inv[5]/np.sqrt(2.)
		odot[0,5] = odot[5,0] = sigma_inv[0]*sigma_inv[5]
		#
		odot[1,1] = sigma_inv[1]*sigma_inv[1]
		odot[1,2] = odot[2,1] = sigma_inv[4]*sigma_inv[4]/2.
		odot[1,3] = odot[3,1] = sigma_inv[3]*sigma_inv[1]
		odot[1,4] = odot[4,1] = sigma_inv[1]*sigma_inv[4]
		odot[1,5] = odot[5,1] = sigma_inv[3]*sigma_inv[4]/np.sqrt(2.)
		#
		odot[2,2] = sigma_inv[2]*sigma_inv[2]
		odot[2,3] = odot[3,2] = sigma_inv[5]*sigma_inv[4]/np.sqrt(2.)
		odot[2,4] = odot[4,2] = sigma_inv[4]*sigma_inv[2]
		odot[2,5] = odot[5,2] = sigma_inv[5]*sigma_inv[2]
		#
		odot[3,3] = sigma_inv[0]*sigma_inv[1] + sigma_inv[3]*sigma_inv[3]/2.
		odot[3,4] = odot[4,3] = sigma_inv[3]*sigma_inv[4]/2. + sigma_inv[5]*sigma_inv[1]/np.sqrt(2.)
		odot[3,5] = odot[5,3] = sigma_inv[0]*sigma_inv[4]/np.sqrt(2.) + sigma_inv[3]*sigma_inv[5]/2.
		#
		odot[4,4] = sigma_inv[1]*sigma_inv[2] + sigma_inv[4]*sigma_inv[4]/2.
		odot[4,5] = odot[5,4] = sigma_inv[3]*sigma_inv[2]/np.sqrt(2.) + sigma_inv[5]*sigma_inv[4]/2.
		#
		odot[5,5] = sigma_inv[0]*sigma_inv[2] + sigma_inv[5]*sigma_inv[5]/2.
		return odot

	
	#return tangent moduli
	def get_dsig_dE(self):
		return self.__dsig_dE
	
	def get_new_inelastic_tensors(self):
		return self.__epsv
	
	def get_new_inelastic_scalars(self):
		return None
	
	
	
	
	
