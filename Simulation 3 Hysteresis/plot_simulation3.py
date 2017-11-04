import numpy
from numpy import linspace
import matplotlib.pyplot as plt
import h5py
import matplotlib.cm as cm


def plot(N, n, T, J, trial):
	plt.figure()
	plt.ylabel('Magnetization per site')
	plt.xlabel('B')
	plt.axhline(0, color='black')
	plt.axvline(0, color='black')
	plt.grid(True, which='both')

	h5f = h5py.File('data/Periodic-boundaries/N='+str(N)+'/Trial '+str(trial)+'/hysteresis,N='+str(N)+',n='+str(n)+',t='+str(T)+',j='+str(J)+'.h5', 'r')
	h5f = h5py.File('data/Periodic-boundaries/N=75/Trial 1 (short equilibration)/hysteresis,N=75,n=50,t=2.3,j=1.h5', 'r')
	Ms = h5f['Ms'][:]
	Bs = h5f['Bs'][:]

	Ms = numpy.insert(Ms, 0, 0)
	Bs = numpy.insert(Bs, 0, 0)

	

	sign = 1.0

	h5f.close()
	plot1 = plt.plot(Bs, Ms)
	plt.scatter(Bs, Ms)	
	plt.axis([-0.3, 0.3, -1.1, 1.1])
	# plt.savefig('data/Periodic-boundaries/N='+str(N)+'/Trial '+str(trial)+'/hysteresis,N='+str(N)+',n='+str(n)+',t='+str(T)+',j='+str(J)+'_connected.png')
	plt.show()
	# plt.savefig('data/Periodic-boundaries/N=75/Trial 1 (short equilibration)/hysteresis,N=75,n=50,t=2.3,j=1_connected.png')

N=75
trial=3
n=50
J=1
ts=[2.1, 1.7, 1.8, 2.0] 	

# for trial in range(3):
# for t in ts:
# 	plot(N, n, t, J, trial+1) 

# plot(N, 75, 2.1, J, trial+1) 
# plot(N, 100, 1.7, J, trial+1) 
# plot(N, 100, 1.8, J, trial+1) 
plot(N, n, 2.3, J, trial+1) 