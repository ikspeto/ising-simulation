import numpy
from numpy import linspace
import matplotlib.pyplot as plt
import h5py
import matplotlib.cm as cm
# import matplotlib.colors as colors



def plot_simulation(N, B, n, J):

	# h5f = h5py.File('data/b='+str(B)+',j=0,n=50,t=0.0-5.0.h5', 'r')
	h5f = h5py.File('j=0,n=100,T=1,B=-5.h5', 'r')
	Ms = h5f['Ms'][:]
	Bs = h5f['Bs'][:]
	h5f.close()

	jet = plt.get_cmap('plasma') 

	colors = [ cm.jet(x) for x in range(0500) ]

	#ax = plt.subplots()

	ax = plt.figure()
	plt.ylabel('Magnetization per site')
	plt.xlabel('H')
	
	plt.scatter(Bs, Ms)
	plt.axhline(0, color='black')
	plt.axvline(0, color='black')
	plt.grid(True, which='both')
	# ax.set_aspect('equal')
	# ax.grid(True, which='both')
	# ax.axhline(y=0, color='k')
	# ax.axvline(x=0, color='k')
	LMs = Ms[48:52]
	LBs = Bs[48:52]
	m,b = numpy.polyfit(LBs, LMs, 1)

	residuals = LMs - (LBs*m+b)
	print "RMSE",(numpy.sum(residuals**2)/(residuals.size-2))**0.5
	print 'm = '+str(m)
	line_x = [-10, 10]
	line_y = [-10*m+b, 10*m+b]	

	plt.plot(line_x, line_y)
	plt.axis([-6, 6, -1.5, 1.5])
	plt.show()

	# plt.savefig('data/MvsT,N='+str(N)+',n='+str(n)+',j='+str(J)+'.png')
	plt.savefig('j=0,n=100,T=1,B=-5.png')

	# plt.figure()
	# plt.ylabel('1/Magnetization per site')
	# plt.xlabel('T')

	# for i in range(Ms.size):
	# 	Ms[i] = 1/Ms[i]

	# plt.scatter(Ts, Ms, color=colors[200])
	# plt.axhline(0, color='black')
	# plt.axvline(0, color='black')
	# plt.grid(True, which='both')
	# plt.savefig('data/1-MvsT,N='+str(N)+',n='+str(n)+',j='+str(J)+'.png')

	# plt.figure()
	# plt.ylabel('Magnetization per site')
	# plt.xlabel('1/T')

	# for i in range(Ms.size):
	# 	Ms[i] = 1/Ms[i]
	# 	Ts[i] = 1/Ts[i]

	# plt.scatter(Ts, Ms, color=colors[200])
	# plt.axhline(0, color='black')
	# plt.axvline(0, color='black')
	# plt.grid(True, which='both')
	# plt.savefig('data/1-MvsT,N='+str(N)+',n='+str(n)+',j='+str(J)+'.png')

	# plt.show()

N = 50
B = 0.1
n = 100
J = 0
plot_simulation(N, B, n, J)