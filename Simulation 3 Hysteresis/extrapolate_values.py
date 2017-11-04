import numpy
from numpy import linspace
import matplotlib.pyplot as plt
import h5py
import matplotlib.cm as cm



def extrapolate(N, n, T, J, trial):
	h5f = h5py.File('data/Periodic-boundaries/N='+str(N)+'/Trial '+str(trial)+'/hysteresis,N='+str(N)+',n='+str(n)+',t='+str(t)+',j='+str(J)+'.h5', 'r')
	Ms = h5f['Ms'][:]
	Bs = h5f['Bs'][:]
	Ms = numpy.insert(Ms, 0, 0)
	Bs = numpy.insert(Bs, 0, 0)
	h5f.close()


	i = 0	
	while(Ms[i] < 0.85):
		i += 1
	while(Ms[i] > 0.85):
		i += 1
	H2 = Bs[i]

	while(Ms[i] > -0.85):
		i += 1
	H1 = Bs[i]

	Hc = -(H2 - numpy.float32(H2-H1)/2.0)

	data_file = open('data/Periodic-boundaries/N='+str(N)+'/Trial '+str(trial)+'/hysteresis,N='+str(N)+',n='+str(n)+',t='+str(t)+',j='+str(J)+'.txt', "w")

	data_file.write("H1 = " + str(H1)+'\n')
	data_file.write("H2 = " + str(H2)+'\n')
	data_file.write("Hc = " + str(numpy.around(Hc, decimals=2))+'\n')

	data_file.close()

N=10
trial=1
n=100
J=0.75
ts=[1.2, 1.3, 1.4, 1.5]

for trial in range(3):
	for t in ts:
		extrapolate(N, n, t, J, trial+1) 