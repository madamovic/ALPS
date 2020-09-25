import pyalps
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyalps.plot
import numpy as np
import pyalps.fit_wrapper as fw
from math import sqrt

#prepare the input parameters
parms = []
for j2 in [0.,1.]:
    for t in np.linspace(0.01,4.0,40):
        parms.append(
            { 
              'LATTICE'        : "coupled ladders", 
              'local_S'        : 1.0,
              'ALGORITHM'      : 'loop',
              'SEED'           : 0,
              'T'              : t,
              'J0'             : 1 ,
              'J1'             : 1,
              'J2'             : j2,
              'THERMALIZATION' : 5000,
              'SWEEPS'         : 50000, 
              'MODEL'          : "spin",
              'L'              : 8,
              'W'              : 8
            }
    )
    
#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm8b',parms)
pyalps.runApplication('loop',input_file)



data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm8b'),['Staggered Susceptibility','Susceptibility','Staggered Magnetization'])
susc1=pyalps.collectXY(data,x='T',y='Susceptibility', foreach=['J2'])
susc2=pyalps.collectXY(data,x='T',y='Staggered Susceptibility', foreach=['J2'])
mag=pyalps.collectXY(data,x='T',y='Staggered Magnetization', foreach=['J2'])


plt.figure()
pyalps.plot.plot(mag)
plt.xlabel(r'$T$')
plt.ylabel(r'$\chi$')
plt.title('Staggered Magnetization')
plt.legend(loc='best')
plt.savefig('2D_heisenberg_ladder_4.eps',dpi=400)

plt.figure()
pyalps.plot.plot(susc1)
plt.xlabel(r'$T$')
plt.ylabel(r'$\chi$')
plt.title('Susceptibility')
plt.legend(loc='best')
plt.savefig('2D_heisenberg_ladder_5.eps',dpi=400)

plt.figure()
pyalps.plot.plot(susc2)
plt.xlabel(r'$T$')
plt.ylabel(r'$\chi_{s}$')
plt.title('Staggered Susceptibility')
plt.legend(loc='best')
plt.savefig('2D_heisenberg_ladder_6.eps',dpi=400)