import pyalps
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyalps.plot
import numpy as np
import pyalps.fit_wrapper as fw

#prepare the input parameters
parms = []
for l in [4,6,8,12,16]:
    for t in [9.70,9.71,9.72,9.73,9.74,9.75,9.76,9.77,9.78,9.79,9.80,9.81,9.82,9.83,9.84,9.85,9.86,9.87,9.88]:
        parms.append(
            { 
              'LATTICE'        : "bcc",
              'LATTICE_LIBRARY' : "bcc.xml", 
              'T'              : t,
              'J'              : 1 ,
              'THERMALIZATION' : 10000,
              'SWEEPS'         : 100000,
              'UPDATE'         : "cluster",
              'MODEL'          : "Ising",
              'L'              : l
            }
    )

#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm7m',parms)
pyalps.runApplication('spinmc',input_file,Tmin=5)
# use the following instead if you have MPI
#pyalps.runApplication('spinmc',input_file,Tmin=5,MPI=4)

pyalps.evaluateSpinMC(pyalps.getResultFiles(prefix='parm7m'))

data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm7m'),'Connected Susceptibility')

connected_susc = pyalps.collectXY(data,x='T',y='Connected Susceptibility',foreach=['L'])

Tc=9.79
a=2.2

#make a data collapse of the connected susceptibility as a function of (T-Tc)/Tc:
for d in connected_susc:
    d.x -= Tc
    d.x = d.x/Tc
    l = d.props['L']
    d.x = d.x * pow(float(l),a)

two_minus_eta=2.7 #your estimate
for d in connected_susc:
    l = d.props['L']
    d.y = d.y/pow(float(l),two_minus_eta)

plt.figure()
pyalps.plot.plot(connected_susc)
plt.xlabel('$L^a(T-T_c)/T_c, Tc=4.511, a=1.6$')
plt.ylabel(r'$L^{\gamma/\nu}\chi_c,\gamma/\nu=$ %.4s' % two_minus_eta)
plt.title('3D Ising model, FCC')
plt.savefig("figure83.eps",dpi=300)
