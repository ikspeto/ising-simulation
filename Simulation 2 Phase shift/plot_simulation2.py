import numpy
from numpy import linspace
import matplotlib.pyplot as plt
import h5py
import matplotlib.cm as cm
# import matplotlib.colors as colors




# x = numpy.arange(10)
# ys = [i+x+(i*x)**2 for i in range(10)]
# #colors = iter(cm.rainbow(numpy.linspace(0, 1, len(ys)))) # plt.scatter(x, y, color=next(colors))

jet = plt.get_cmap('plasma') 

colors = [ cm.jet(x) for x in range(0500) ]

# cNorm  = colors.Normalize(vmin=0, vmax=4)
# scalarMap = cm.ScalarMappable(norm=cNorm	, cmap=jet)


# h5f = h5py.File('data/b=5.0,j=1,n=50,t=4-7.h5', 'r')
# Ms = h5f['Ms'][:]
# Ts = h5f['Ts'][:]
# h5f.close()

plt.figure()
plt.ylabel('Magnetization per site')
plt.xlabel('Temperature')
# plt.scatter(Ts, Ms, color=colors[0])

# h5f = h5py.File('data/b=5.0,j=1,n=50,t=1.0-4.0.h5', 'r')
# Ms = h5f['Ms'][:]
# Ts = h5f['Ts'][:]
# h5f.close()
# plt.scatter(Ts, Ms, color=colors[0])

# h5f = h5py.File('data/b=5.0,j=1,n=50,t=7.0-12.0.h5', 'r')
# Ms = h5f['Ms'][:]
# Ts = h5f['Ts'][:]
# h5f.close()
# plt.scatter(Ts, Ms, color=colors[0])

# h5f = h5py.File('data/b=5.0,j=1,n=50,t=12.0-17.0.h5', 'r')
# Ms = h5f['Ms'][:]
# Ts = h5f['Ts'][:]
# h5f.close()
# plt.scatter(Ts, Ms, color=colors[0])

# h5f = h5py.File('data/b=5.0,j=1,n=50,t=17.0-27.0.h5', 'r')
# Ms = h5f['Ms'][:]
# Ts = h5f['Ts'][:]
# h5f.close()
# plt.scatter(Ts, Ms, color=colors[0])





# h5f = h5py.File('data/b=1.0,j=1,n=50,t=1-4.h5', 'r')
# Ms = h5f['Ms'][:]
# Ts = h5f['Ts'][:]
# h5f.close()
# plt.scatter(Ts, Ms, color=colors[100])

# h5f = h5py.File('data/b=1.0,j=1,n=50,t=4-10.0.h5', 'r')
# Ms = h5f['Ms'][:]
# Ts = h5f['Ts'][:]
# h5f.close()
# plt.scatter(Ts, Ms, color=colors[100])

# h5f = h5py.File('data/b=1.0,j=1,n=50,t=10.0-27.0.h5', 'r')
# Ms = h5f['Ms'][:]
# Ts = h5f['Ts'][:]
# h5f.close()
# plt.scatter(Ts, Ms, color=colors[100])




# h5f = h5py.File('data/b=0.5,j=1,n=50,t=1-4.h5', 'r')
# Ms = h5f['Ms'][:]
# Ts = h5f['Ts'][:]
# h5f.close()
# plt.scatter(Ts, Ms, color=colors[200])

# h5f = h5py.File('data/b=0.5,j=1,n=50,t=4-10.0.h5', 'r')
# Ms = h5f['Ms'][:]
# Ts = h5f['Ts'][:]
# h5f.close()
# plt.scatter(Ts, Ms, color=colors[200])

# h5f = h5py.File('data/b=0.5,j=1,n=50,t=10.0-27.0.h5', 'r')
# Ms = h5f['Ms'][:]
# Ts = h5f['Ts'][:]
# h5f.close()
# plt.scatter(Ts, Ms, color=colors[200])


h5f = h5py.File('data/b=0.0,j=1,n=50,t=1-4.h5', 'r')
Ms = h5f['Ms'][:]
Ts = h5f['Ts'][:]
# for i in range(Ms.size):
# 	Ms[i] = abs(Ms[i])
h5f.close()
print Ms.size
plt.scatter(Ts, Ms, color=colors[300])

plt.show()

# h5f = h5py.File('data.h5','r')
# b = h5f['dataset_1'][:]
# h5f.close()

# plt.figure()
# plt.ylabel('Magnetization per site')
# plt.xlabel('Temperature')
# plt.scatter(Ts, Ms)
# plt.savefig('b=5.0,j=1,n=50,t=2.0.png')
# plt.figure()
# plt.ylabel('Energy per site')
# plt.xlabel('Temperature')
# plt.scatter(Ts, Es)
# plt.figure()
# plt.ylabel('Magnetization per site')
# plt.xlabel('Energy per site')
# plt.scatter(Es, Ms)
#plt.axis([0, 5, -3.5, 3.5])	