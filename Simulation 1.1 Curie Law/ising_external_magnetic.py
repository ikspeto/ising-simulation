import numpy
import matplotlib.pyplot as plt
import h5py

def ising(N, T, B):
	# Initialization
	n = 10000000
	xs = [1, 0, -1, 0]
	ys = [0, 1, 0, -1]
	exps = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

	M = 0 # Magnetization
	E = 0 # Energy

	print B

	J = 0
	k = 1
	
	beta = 1.0/numpy.float32(T)

	exps[0] = numpy.exp(-beta*0)
	exps[2] = numpy.exp(-beta*2)
	exps[6] = numpy.exp(-beta*6)
	exps[8] = numpy.exp(-beta*8)
	#print(exps[2])

	grid = numpy.random.random_integers(0, 1, (N, N))
	#grid = numpy.ones((N, N), dtype=numpy.int)
	
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
	#plt.ion()
	#im = plt.imshow(grid, cmap='gray', interpolation="nearest")	

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
		#print(delta_e)
		
		if(delta_e <= 0):
			grid[i, j] = -grid[i, j]
			E += delta_e
			M += 2*grid[i, j]
		elif(numpy.random.random(1)[0] < numpy.exp(-beta*delta_e)): #exps[delta_e]): #numpy.exp(-k*T*delta_e)):
			grid[i, j] = -grid[i, j]
			E += delta_e
			M += 2*grid[i, j]

		#M = sum(sum(grid))
		#E = -sum(sum(delta_e)) / 2
		#Ms[z] = numpy.float32(M)/numpy.float32(N*N)
		#Es[z] = numpy.float32(E)/numpy.float32(N*N)
		
		
		#if (z%2000 == 0):
			#im.set_data(grid)
			#plt.pause(0.0001)

		#plt.subplot(211)
		#plt.imshow(grid)
		#plt.subplot(212)
		#plt.imshow(grid, cmap='Greys', interpolation='nearest')
		# plt.savefig('blkwht.png')
	M1 = 0
	for i in range(N):
		for j in range(N):
			M1 += grid[i, j]

	M2 = 0
	for i in range(N):
		for j in range(N):
			M2 += grid[i, j]*grid[i, j]

	#return (M2/numpy.float32(N*N)-(M1/numpy.float32(N*N))*(M1/numpy.float32(N*N)), E)
	return (M1, E)
	
	#print(grid)
	#im = plt.imshow(grid, cmap='gray', interpolation="nearest")	
	#plt.show()
	#print(grid)

	
	
	
def metropolis():
	N = 25
	n = 100
	T = 1
	B = -5
	step = (-numpy.float32(B)*2.0)/numpy.float32(n)

	Ms = numpy.empty([n], dtype='f')
	Es = numpy.empty([n], dtype='f')
	Ts = numpy.empty([n], dtype='f')
	Bs = numpy.empty([n], dtype='f')

	h5f = h5py.File('j=0,n='+str(n)+',N='+str(N)+'T='+str(T)+',B='+str(B)+'.h5', 'w')

	for i in range(n):
		#T = numpy.random.random(1)[0]*5
		#B = numpy.random.random(1)[0]
		res = ising(N, T, B)
		Ms[i] = res[0]/numpy.float32(N*N)
		#Ms[i] = res[0]/numpy.float32(T)
		Es[i] = res[1]/numpy.float32(N*N)
		Ts[i] = T
		Bs[i] = B
		B += step

	h5f.create_dataset('Ms', data=Ms)
	h5f.create_dataset('Es', data=Es)
	h5f.create_dataset('Ts', data=Ts)
	h5f.create_dataset('Bs', data=Bs)
	h5f.close()

	m,b = numpy.polyfit(Bs, Ms, 1)

	line_x = [-10, 10]
	line_y = [-10*m+b, 10*m+b]

	plt.figure()
	plt.ylabel('Magnetization per site')
	plt.xlabel('B')
	plt.scatter(Bs, Ms)
	# plt.plot(line_x, line_y)
	# plt.figure()
	# plt.ylabel('Energy per site')
	# plt.xlabel('Temperature')
	# plt.scatter(Ts, Es)
	# plt.figure()
	# plt.ylabel('Magnetization per site')
	# plt.xlabel('Energy per site')
	# plt.scatter(Es, Ms)
	#plt.axis([0, 5, -3.5, 3.5])	
	plt.show()

	

	# print m

	# C = m*T
	# Tc = T-C/m

	# f = open('Ms_vs_Bs_data.txt', 'w')
	# for i in range(n):
	# 	f.write(str(Bs[i]) + " " + str(Ms[i]) + "\n")
	# f.close()



metropolis()