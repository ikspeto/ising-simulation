import numpy
import matplotlib.pyplot as plt
import h5py

def ising(N, T, B):
	# Initialization
	n = 2100000
	xs = [1, 0, -1, 0]
	ys = [0, 1, 0, -1]

	M = 0 # Magnetization
	E = 0 # Energy

	J = 1
	k = 1
	
	beta = 1.0/numpy.float32(k*T)

	grid = numpy.random.random_integers(0, 1, (N, N))
	
	for i in range(N):
		for j in range(N):
			if(grid[i, j] == 0):
				grid[i, j] -= 1

	# Initial Energy
	for i in range(N):
		for j in range(N):
			ssum = 0
			for x in range(4):
				i2 = i + xs[x]
				j2 = j + ys[x]
				
				if(i2 < 0):
					i2 = N-1
				elif(i2 >= N):
					i2 = 0
				if(j2 < 0):
					j2 = N-1
				elif(j2 >= N):
					j2 = 0
				ssum += grid[i2, j2]
			E += -J * grid[i, j] * (ssum)/2 - B * grid[i, j]
	# Initial Magnetization
	for i in range(N):
		for j in range(N):
			M += grid[i, j]

	# Simulation
	for z in range(n):
		#T = E
		i = numpy.random.random_integers(0, N-1, 1)[0]
		j = numpy.random.random_integers(0, N-1, 1)[0]

		ssum = 0
		
		for x in range(4):
			i2 = i + xs[x]
			j2 = j + ys[x]
			if(i2 < 0):
				i2 = N-1
			elif(i2 >= N):
				i2 = 0
			if(j2 < 0):
				j2 = N-1
			elif(j2 >= N):
				j2 = 0
			ssum += grid[i2, j2]

		delta_e = 2*J*grid[i, j]*(ssum) + 2 * B * grid[i, j]
		
		if(delta_e <= 0):
			grid[i, j] = -grid[i, j]
			E += delta_e
			M += 2*grid[i, j]
	
		elif(numpy.random.random(1)[0] < numpy.exp(-beta*delta_e)):#exps[delta_e]): #numpy.exp(-k*T*delta_e)):
			grid[i, j] = -grid[i, j]
			E += delta_e
			M += 2*grid[i, j]
	return (M, E)
	
def metropolis(B, n, T1, T2):
	N = 50
	T = T1
	step = numpy.float32((T2-T1)/n)
	print step
	
	print('b='+str(B)+',j=1,n=50,t='+str(T1)+'-'+str(T2)+'.h5')
	h5f = h5py.File('b='+str(B)+',j=1,n=50,t='+str(T1)+'-'+str(T2)+'.h5', 'w')

	Ms = numpy.empty([n], dtype='f')
	Es = numpy.empty([n], dtype='f')
	Ts = numpy.empty([n], dtype='f')
	
	for i in range(n):
		#T = numpy.random.random(1)[0]*5
		res = ising(N, T, B)
		Ms[i] = res[0]/numpy.float32(N*N)
		Es[i] = res[1]/numpy.float32(N*N)
		Ts[i] = T
		T += step

	h5f.create_dataset('Ms', data=Ms)
	h5f.create_dataset('Es', data=Es)
	h5f.create_dataset('Ts', data=Ts)
	h5f.close()


	# h5f = h5py.File('data.h5','r')
	# b = h5f['dataset_1'][:]
	# h5f.close()

	# plt.figure()
	# plt.ylabel('Magnetization per site')
	# plt.xlabel('Temperature')
	# plt.scatter(Ts, Ms)
	# plt.savefig('b='+str(5.0)+',j=1,n=50,t=2.0.png')
	# plt.figure()
	# plt.ylabel('Energy per site')
	# plt.xlabel('Temperature')
	# plt.scatter(Ts, Es)
	# plt.figure()
	# plt.ylabel('Magnetization per site')
	# plt.xlabel('Energy per site')
	# plt.scatter(Es, Ms)
	#plt.axis([0, 5, -3.5, 3.5])	
	#plt.show()



metropolis(5.0, 200, 17.0, 27.0)
metropolis(0.5, 200, 10.0, 27.0)
metropolis(1.0, 200, 10.0, 27.0)