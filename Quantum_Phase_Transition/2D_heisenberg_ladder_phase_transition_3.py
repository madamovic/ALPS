# LxL

import pyalps
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyalps.plot
import numpy as np

#prepare the input parameters
parms = []
for l in [8,10,12,16]:
    for j2 in np.linspace(0.3,0.4,100):
        parms.append(
            { 
              'LATTICE'        : "coupled ladders", 
              'local_S'        : 0.5,
              'ALGORITHM'      : 'loop',
              'SEED'           : 0,
              'BETA'           : 2*l,
              'J0'             : 1 ,
              'J1'             : 1,
              'J2'             : j2,
              'THERMALIZATION' : 5000,
              'SWEEPS'         : 50000, 
              'MODEL'          : "spin",
              'L'              : l,
              'W'              : l
            }
    )
    
#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm8d',parms)
pyalps.runApplication('loop',input_file)

data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm8d'),['Binder Ratio of Staggered Magnetization','Stiffness'])

binder=pyalps.collectXY(data,x='J2',y='Binder Ratio of Staggered Magnetization', foreach=['L'])
stiffness =pyalps.collectXY(data,x='J2',y='Stiffness', foreach=['L'])

for q in stiffness:
    q.y = q.y*q.props['L']

#make plot    
plt.figure()
pyalps.plot.plot(stiffness)
plt.xlabel(r'$J2$')
plt.ylabel(r'Stiffness $\rho_s L$')
plt.legend(loc='best')
plt.savefig('2D_heisenberg_ladder_12.eps',dpi=400)

plt.figure()
pyalps.plot.plot(binder)
plt.xlabel(r'$J_2$')
plt.ylabel(r'$g(m_s)$')
plt.title('coupled ladders')
plt.legend(loc='best')
plt.savefig('2D_heisenberg_ladder_13.eps',dpi=400)