
import numpy as np
from scipy.ndimage import measurements
import os
import time

class IPfilm:

	def __init__(self,Nrows,Ncol,theta):

		#Parameters
		self.Nrows = int(Nrows) #row
		self.Ncol = int(Ncol) #col
		self.theta = (theta*np.pi)/180. # convertion from angles to radians is needed.
		self.g = 9.81
		self.dimensional_constant = 1#0.140/self.Nrows
		self.rho_nw = 1.225
		self.rho_w = 1.205
		self.delta_rho = self.rho_nw-self.rho_w
		self.L = 1
		self.path='/home/monem/Documents/python/IPfilm/'
		self.pathsimulation = self.path+"Nrows"+str(self.Nrows)+"_Ncol"+str(self.Ncol)+"_theta"+str(self.theta)+"/"
		
		#Element specification.
		self.pore = 1
		self.bead = 2
		self.bond = 3
		self.invaded_pore = 0
		self.invaded_bond = 5
		self.film = 4

		if not os.path.exists(self.pathsimulation):
			os.makedirs(self.pathsimulation)

	def structure(self):

		L,g,Nrows,Ncol,delta_rho,theta = self.L,self.g,self.Nrows,self.Ncol,self.delta_rho, self.theta
		BL = np.zeros([Nrows+4,(Ncol+4)],int)
		PL = np.zeros([Nrows+4,Ncol+4],float)
		Lattice_only_bonds = np.zeros([Nrows,Ncol/2],int)

		#Fyller ytterpunktene med verdien 9.
		BL[:,-2:] = BL[-2:,:] = BL[:2,:] = BL[:,:2] = 9;
		
		#Bondary condition.
		for i in range(2,Nrows+2,1):
			if i%2 == 1:
				BL[i,2:-2:2] = self.bead
				BL[i,3:-2:2] = self.bond

			else:
				BL[i,2:-2:2] = self.bond
				BL[i,3:-2:2] = self.pore

		if BL[2,3] == 1:
			BL[2,3:-2:2] = self.invaded_pore

		#Random seed number
		np.random.seed(1)

		#Filling the presssure lattice.
		for i in xrange(Nrows+4):
			for j in xrange(Ncol+4):
				PL[i,j]= np.random.uniform() - self.dimensional_constant*(Nrows-i)*g*delta_rho*np.sin(theta)


		self.BL = BL
		self.PL = PL
		self.termination_criteria = sum(BL[-5,:])
		self.Lattice_only_bonds = Lattice_only_bonds

		return BL,PL

	def invasion(self):

		Nrows, Ncol = self.Nrows, self.Ncol
		#setting up the lattices and pressure list.
		BL,PL = self.structure()
		SL = self.bond2site(BL)
		pressure = []
		count = 1
		start_time = time.time()

		tmp = 1
		while tmp <= 30:

			Pmax=1; itmp=0; jtmp=0; ibond=0; jbond=0
			#Looping trough all the bonds in the bond matrix, checking the vertical situation of the inner points.
			#lw, num = measurements.label(SL)

			for i in range(2,(Nrows+2),1):
				for j in range(2,(Ncol+2),1):

					#Downward invasion
					if (BL[i,j] == self.bond and BL[i-1,j] == self.invaded_pore and BL[i+1,j] == self.pore and PL[i,j]<Pmax):
						Pmax = PL[i,j]
						itmp =i+1; jtmp=j; ibond= i; jbond=j

					#Invasion towards the right
					if (BL[i,j] == self.bond and BL[i,j-1] == self.invaded_pore and BL[i,j+1] == self.pore and PL[i,j]<Pmax):
						Pmax = PL[i,j]
						itmp = i; jtmp=j+1; ibond = i; jbond=j

					#Invasion towards the left
					if (BL[i,j] == self.bond and BL[i,j+1] == self.invaded_pore and BL[i,j-1] == self.pore and PL[i,j]<Pmax):
						Pmax = PL[i,j]
						itmp = i; jtmp=j-1; ibond = i; jbond=j

					#print Pmax
					
			BL[itmp,jtmp] = self.invaded_pore
			BL[ibond,jbond] = self.invaded_bond

			#Logging pressure
			pressure.append(Pmax)
			
			#Filmflow
			self.evaluatefilm(itmp,jtmp,ibond,jbond)
			
			#Cutoff
			if sum(BL[-3,:]) != self.termination_criteria:
				break

			count += 1
			tmp += 1

		self.BL = BL
		film = self.filmthresh()
		SL = self.bond2site(BL)
		Lattice_only_bonds = self.onlybonds(BL)

		#write to file
		np.savetxt(self.pathsimulation + "pressure"+str(Nrows)+"_"+str(Ncol)+".txt",[pressure])
		np.savetxt(self.pathsimulation + "SL"+str(Nrows)+"_"+str(Ncol)+".txt",SL)
		np.savetxt(self.pathsimulation + "BL"+str(Nrows)+"_"+str(Ncol)+".txt",BL)
		np.savetxt(self.pathsimulation + "film"+str(Nrows)+"_"+str(Ncol)+".txt",film)

		return BL, SL,film, Lattice_only_bonds

	def bond2site(self,BL):

		Nrows, Ncol = self.Nrows, self.Ncol
		SL = np.zeros([int(Nrows/2),int(Ncol/2)],int)

		count = 0
		for i in range(2,Nrows+1,2):			
			if i%2 != 1:
				SL[count,:] = BL[i,3:-2:2]
				count += 1
			pass

		return SL

	def site2bond(self,SL):

		Nrows,Ncol = self.Nrows,self.Ncol
		BL = self.BL

		count = 0
		for i in range(1,Nrows+1,1):
			if i%2 != 1:
				BL[i,3:-2:2] = SL[count,:]
				count+=1

		return BL

	def evaluatefilm(self,itmp,jtmp,ibond,jbond):

		BL = self.BL

		if (jtmp-jbond) < 0:
			#invasion to the left.
			if (BL[itmp,jtmp] == BL[itmp+2,jtmp]):
				self.BL[itmp+1,jtmp]  = self.film
			elif (BL[itmp,jtmp] == BL[itmp-2,jtmp]):
				self.BL[itmp-1,jtmp]  = self.film
			elif (BL[itmp,jtmp] == BL[itmp,jtmp-2]):
				self.BL[itmp,jtmp-1]  = self.film

		elif (itmp-ibond) < 0:
			#invasion from the bottom.
			if (BL[itmp,jtmp] == BL[itmp,jtmp+2]):
				self.BL[itmp,jtmp+1] = self.film
			elif (BL[itmp,jtmp] == BL[itmp,jtmp-2]):
				self.BL[itmp,jtmp-1] = self.film
			elif (BL[itmp,jtmp] == BL[itmp-2,jtmp]):
				self.BL[itmp-1,jtmp] = self.film

		elif (jtmp-jbond) > 0:
			#invasion to the right.
			if (BL[itmp,jtmp] == BL[itmp+2,jtmp]):
				self.BL[itmp+1,jtmp]  = self.film
			elif (BL[itmp,jtmp] == BL[itmp,jtmp+2]):
				self.BL[itmp-1,jtmp]  = self.film
			elif (BL[itmp,jtmp] == BL[itmp-2,jtmp]):
				self.BL[itmp,jtmp+1]  = self.film

	def onlybonds(self,BL):

		Nrows,Ncol = self.Nrows,self.Ncol
		Lattice_only_bonds = self.Lattice_only_bonds
		BL = self.BL
		counti = 0; countj = 0

		for i in range(2,Nrows+2,1):

			if BL[i,2] == 3 or BL[i,2] == 5:
				for j in range(2,Ncol+2,2):
					Lattice_only_bonds[counti,countj] = BL[i,j]
					countj = countj+1

			if BL[i,3] == 3 or BL[i,3] == 5:
				for j in range(3,Ncol+2,2):
					Lattice_only_bonds[counti,countj] = BL[i,j]
					countj = countj+1

			countj =  0
			counti += 1

		return Lattice_only_bonds

	def plotting(self,BL, SL,film,Lattice_only_bonds):

		import matplotlib.pylab as plt

		plt.figure(1)
		#plt.imshow(BL[2:-2,2:-2], interpolation='nearest',cmap = 'jet')
		plt.imshow(BL,interpolation='nearest',cmap = 'tab10')
		plt.title('pore = 1, bead = 2, bond = 3, film = 4, invaded_pore = 0')
		plt.colorbar()
		
		plt.figure(2)
		plt.imshow(SL,cmap='Greys', interpolation='nearest')
		plt.title("SL")
		plt.colorbar()

		plt.figure(5)
		plt.imshow(Lattice_only_bonds,interpolation='nearest',cmap = 'jet')
		plt.title('Lattice_only_bonds')
		plt.colorbar()

		plt.figure(3)
		plt.imshow(film[2:-2,2:-2], interpolation='nearest')
		plt.title("film")
		plt.colorbar()

		plt.figure(4)
		plt.imshow(self.PL,cmap='Greys', interpolation='nearest')
		plt.title("SL")
		plt.colorbar()
		plt.show()

	def filmthresh(self):

		Nrows,Ncol = self.Nrows,self.Ncol
		BL = self.BL

		film = np.zeros([Nrows+4,(Ncol+4)],int)

		for i in range(2,(Nrows+2),1):
			for j in range(2,(Ncol+2),1):

				if BL[i,j] >= 4 or BL[i,j] == 2:
					film[i,j] = 1;
				else:
					film[i,j] = 0;

		return film


def execute():

	#Checking the execution time.
	import time
	start_time = time.time()
	c=IPfilm(Nrows=6*3,Ncol=6*2,theta=30) # Nrow and Ncol must be even numbers.
	c.structure()
	BL, SL,film, Lattice_only_bonds = c.invasion()
	print Lattice_only_bonds
	c.plotting(BL,SL,film,Lattice_only_bonds)

	print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
	execute()