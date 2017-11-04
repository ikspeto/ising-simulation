import numpy
import matplotlib.pyplot as plt
import h5py

# constants
xs = [1, 0, -1, 0]
ys = [0, 1, 0, -1]
J = 0.75 # Coupling strength
k = 1 # k Boltzman

# global variables
 # Initialize grid with random 1s and 0s
	
def initialize_grid(N):
	grid = numpy.random.random_integers(0, 1, (N, N)) # Initialize grid with random 1s and 0s
	#grid = numpy.ones((N, N), dtype=numpy.int) # Initialize grid with 1s

def calculate_energy(N, B):
	E = 0
	for i in range(N):
		for j in range(N):
			ssum = 0
			for x in range(4):
				i2 = i + xs[x]
				j2 = j + ys[x]
				
				if(i2 < 0 or i2 >= N or j2 < 0 or j2 >= N):
					continue
				
				ssum += grid[i2, j2] 
			E += -J * grid[i, j] * (ssum)/2 - B * grid[i, j]
	return E

def calculate_magnetization(N):
	M = 0
	for i in range(N):
		for j in range(N):
			M += grid[i, j]
	return M

def ising(N, T, B):
	# Initialization
	n_iterations = 1000000 # How many times the step is done

	M = 0 # Magnetization
	E = 0 # Energy

	beta = 1.0/numpy.float32(T) # 1/T
	

	# Initial Energyx
	E = calculate_energy(N, B)
	# Initial Magnetization
	M = calculate_magnetization(N)
	# Simulation
	for z in range(n_iterations):
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
		
		# if (z%2000 == 0):
		# 	im.set_data(grid)
		# 	plt.pause(0.0001)

		#plt.subplot(211)
		#plt.imshow(grid)
		#plt.subplot(212)
		#plt.imshow(grid, cmap='Greys', interpolation='nearest')
		# plt.savefig('blkwht.png')
	
	return (M, E)
	
def metropolis(N, T, n, B1, trial):
	n_iterations = n
	B = 0 # From -20 to +20
	step = numpy.float32(B1*5)/(numpy.float32(n))
	
	Ms = numpy.empty([n_iterations], dtype='f')
	Es = numpy.empty([n_iterations], dtype='f')
	Ts = numpy.empty([n_iterations], dtype='f')
	Bs = numpy.empty([n_iterations], dtype='f')
	#initialize_grid(N)

	for i in range(N): # Make 0s to -1s
		for j in range(N):
			if(grid[i, j] == 0):
				grid[i, j] -= 1

	plt.figure()
	plt.ylabel('Magnetization per site')
	plt.xlabel('B')
	plt.ion()
	plt.show()

	print ('data/N='+str(N)+'/Trial '+str(trial)+'/hysteresis,N='+str(N)+',n='+str(n)+',t='+str(T)+',j='+str(J)+'.png')
	print ('data/N='+str(N)+'/Trial '+str(trial)+'/hysteresis,N='+str(N)+',n='+str(n)+',t='+str(T)+',j='+str(J)+'.h5')

	initial_magnetization = calculate_magnetization(N)
	initial_energy = calculate_energy(N, B)
	plt.scatter(B, initial_magnetization/numpy.float32(N*N))

	for i in range(int(n/5)):
		res = ising(N, T, B)
		Ms[i] = res[0]/numpy.float32(N*N)
		Es[i] = res[1]/numpy.float32(N*N)
		Ts[i] = T
		Bs[i] = B
		plt.scatter(B, Ms[i])
		plt.pause(0.0001)
		B += step
	for i in range(int(n/5), int((3*n)/5)):
		res = ising(N, T, B)
		Ms[i] = res[0]/numpy.float32(N*N)
		Es[i] = res[1]/numpy.float32(N*N)
		Ts[i] = T
		Bs[i] = B
		plt.scatter(B, Ms[i])
		plt.pause(0.0001)
		B -= step
	for i in range(int((3*n)/5), n):
		res = ising(N, T, B)
		Ms[i] = res[0]/numpy.float32(N*N)
		Es[i] = res[1]/numpy.float32(N*N)
		Ts[i] = T
		Bs[i] = B
		plt.scatter(B, Ms[i])
		plt.pause(0.0001)
		B += step

	plt.savefig('data/N='+str(N)+'/Trial '+str(trial)+'/hysteresis,N='+str(N)+',n='+str(n)+',t='+str(T)+',j='+str(J)+'.png')

	Ms = numpy.insert(Ms, 0, initial_magnetization/numpy.float32(N*N))
	Bs = numpy.insert(Bs, 0, 0)
	Es = numpy.insert(Es, 0, initial_energy/numpy.float32(N*N))

	h5f = h5py.File('data/N='+str(N)+'/Trial '+str(trial)+'/hysteresis,N='+str(N)+',n='+str(n)+',t='+str(T)+',j='+str(J)+'.h5', 'w')
	h5f.create_dataset('Ms', data=Ms)
	h5f.create_dataset('Es', data=Es)
	h5f.create_dataset('Bs', data=Bs)
	h5f.close()

	# m,b = numpy.polyfit(Bs, Ms, 1)

	# line_x = [-10, 10]
	# line_y = [-10*m+b, 10*m+b]

	# plt.figure()
	# plt.ylabel('Magnetization per site')
	# plt.xlabel('B')
	# plt.scatter(Bs, Ms)
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
	# plt.show()

	

	# print m

	# C = m*T
	# Tc = T-C/m

	# f = open('Ms_vs_Bs_data.txt', 'w')
	# for i in range(n):
	# 	f.write(str(Bs[i]) + " " + str(Ms[i]) + "\n")
	# f.close()

grid = numpy.random.random_integers(0, 1, (50, 50))
# im = plt.imshow(grid, cmap='gray', interpolation="nearest")	
for i in range(2, 4):
		  # N,   t,   n,   J, trial
	metropolis(5, 1.2, 100, 0.6, i+1)
	metropolis(5, 1.3, 100, 0.6, i+1)
	metropolis(5, 1.4, 100, 0.6, i+1)
	metropolis(5, 1.5, 100, 0.6, i+1)

