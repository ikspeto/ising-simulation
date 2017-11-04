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
	Ms = h5f['Ms'][:]
	Bs = h5f['Bs'][:]

	Ms = numpy.insert(Ms, 0, 0)
	Bs = numpy.insert(Bs, 0, 0)

	sign = 1.0

	h5f.close()
	plot1 = plt.plot(Bs, Ms)
	plt.scatter(Bs, Ms)	
	plt.axis([-0.7, 0.7, -1.5, 1.5])
	# plt.savefig('data/Periodic-boundaries/N='+str(N)+'/Trial '+str(trial)+'/hysteresis,N='+str(N)+',n='+str(n)+',t='+str(T)+',j='+str(J)+'_connected.png')
	plt.show()

N=25
trial=1
n=100
J=0.75
ts=[1.2, 1.3, 1.4, 1.5] 	

# for trial in range(3):
# 	for t in ts:
plot(N, n, ts[1], J, trial) 