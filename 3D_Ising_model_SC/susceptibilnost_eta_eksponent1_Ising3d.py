import pyalps
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyalps.plot
import numpy as np
import pyalps.fit_wrapper as fw

#prepare the input parameters
parms = []
for l in [2,4,6,8,10,12]:
    for t in [4.45,4.46,4.47,4.48,4.49,4.50,4.51,4.52,4.53,4.54,4.55,4.56,4.57,4.58]:
        parms.append(
            { 
              'LATTICE'        : "square lattice", 
              'T'              : t,
              'J'              : 1 ,
              'THERMALIZATION' : 20000,
              'SWEEPS'         : 1000000,
              'UPDATE'         : "cluster",
              'MODEL'          : "Ising",
              'L'              : l
            }
    )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm7a',parms)
pyalps.runApplication('spinmc',input_file,Tmin=5)
# use the following instead if you have MPI
#pyalps.runApplication('spinmc',input_file,Tmin=5,MPI=4)

pyalps.evaluateSpinMC(pyalps.getResultFiles(prefix='parm7a'))

data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm7a'),'Connected Susceptibility')

connected_susc = pyalps.collectXY(data,x='T',y='Connected Susceptibility',foreach=['L'])

Tc=4.511
a=1.6

#make a data collapse of the connected susceptibility as a function of (T-Tc)/Tc:
for d in connected_susc:
    d.x -= Tc
    d.x = d.x/Tc
    l = d.props['L']
    d.x = d.x * pow(float(l),a)

two_minus_eta=1.96 #your estimate
for d in connected_susc:
    l = d.props['L']
    d.y = d.y/pow(float(l),two_minus_eta)

plt.figure()
pyalps.plot.plot(connected_susc)
plt.xlabel('$L^a(T-T_c)/T_c, Tc=4.511, a=1.6$')
plt.ylabel(r'$L^{\gamma/\nu}\chi_c,\gamma/\nu=$ %.4s' % two_minus_eta)
plt.title('3D Ising model')
plt.savefig("figure37.eps",dpi=300)
