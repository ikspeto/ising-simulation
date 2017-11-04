#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy
from numpy import linspace
import matplotlib.pyplot as plt
import h5py
import matplotlib.cm as cm
# import matplotlib.colors as colors



def plot_simulation(N, B, n, J):

	h5f = h5py.File('data/N=50,b=0.3,j=0,n=50,t=0.1-2.0.h5', 'r')
	Ms = h5f['Ms'][:]
	Ms1 = h5f['Ms'][:]
	Ts = h5f['Ts'][:]
	h5f.close()

	h5f = h5py.File('data/N=50,b=0.3,j=0,n=50,t=2.0-4.0.h5', 'r')
	Ms2 = h5f['Ms'][:]
	Ts2 = h5f['Ts'][:]
	h5f.close()

	jet = plt.get_cmap('plasma') 

	colors = [ cm.jet(x) for x in range(0500) ]

	# #ax = plt.subplots()

	ax = plt.figure()
	plt.ylabel('Magnetization per site')
	plt.xlabel('T')
	
	plt.scatter(Ts, Ms, color=colors[300])
	plt.scatter(Ts2, Ms2, color=colors[300])
	plt.axhline(0, color='black')
	plt.axvline(0, color='black')
	plt.grid(True, which='both')
	# # ax.set_aspect('equal')
	# # ax.grid(True, which='both')
	# # ax.axhline(y=0, color='k')
	# # ax.axvline(x=0, color='k')
	
	# # plt.savefig('data/MvsT,N='+str(N)+',n='+str(n)+',j='+str(J)+'.png')

	# plt.figure()
	# plt.ylabel('Magnetization per site')
	# plt.xlabel('1/T')

	# LTs = Ts[1:48:1]
	# LMs = Ms[1:48:1]
	# for i in range(LMs.size):
	# 	LTs[i] = 1/LTs[i]

	# plt.scatter(LTs, LMs, color=colors[200])
	# plt.axhline(0, color='black')
	# plt.axvline(0, color='black')
	# plt.grid(True, which='both')

	# print "M vs 1/T"
	# m,b = numpy.polyfit(LTs, LMs, 1)

	# residuals = Ms - (Ts*m+b)
	# sssum = 0
	# NN = 0
	# i = 0
	# for a in residuals:
	# 	if a > 1000:
	# 		residuals[i] = 0
	# 		NN+=1
	# 	sssum += residuals[i]**2
	# sssum /= numpy.float32(residuals.size-2-NN)
	# sssum = sssum**0.5
	# print sssum
	# print "RMSE",(numpy.sum(residuals**2)/(residuals.size-2))**0.5

	# print 'm = '+str(m)
	# line_x = [-10, 10]
	# line_y = [-10*m+b, 10*m+b]	

	# plt.plot(line_x, line_y)

	# # plt.savefig('data/1-MvsT,N='+str(N)+',n='+str(n)+',j='+str(J)+'.png')




	print "-------------------------------------------------------"
	print "x vs T"
	plt.figure()
	plt.ylabel(r'1/$\chi$') 
	plt.xlabel('T')

	for i in range(Ms.size):
		Ms[i] = 1/(Ms[i]/0.3)
		Ts[i] = Ts[i]
		Ms2[i] = 1/(Ms2[i]/0.3)	

	LMs = Ms[1:20:1]
	LTs = Ts[1:20:1]

	m,b = numpy.polyfit(LTs, LMs, 1)

	residuals = LMs - (LTs*m+b)
	print "RMSE",(numpy.sum(residuals**2)/(residuals.size-2))**0.5
	print 'm = '+str(m)
	line_x = [-10, 10]
	line_y = [-10*m+b, 10*m+b]	

	# plt.plot(line_x, line_y)

	plt.scatter(LTs, LMs, color=colors[200])
	plt.scatter(Ts2, Ms2, color=colors[200])
	plt.axhline(0, color='black')
	plt.axvline(0, color='black')
	plt.grid(True, which='both')

	#plt.savefig('data/MvsinvT,N='+str(N)+',n='+str(n)+',j='+str(J)+'.png')

	plt.show()

N = 15
B = 0.2
n = 20
J = 0
plot_simulation(N, B, n, J)

# sssum = 0
# 	NN = 0
# 	i = 0
# 	for a in residuals:
# 		sssum += a**2
# 	sssum /= numpy.float32(residuals.size-2-NN)
# 	sssum = sssum**0.5
# 	print sssum