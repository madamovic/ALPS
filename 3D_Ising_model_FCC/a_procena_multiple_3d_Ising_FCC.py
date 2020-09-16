
# coding: utf-8

# In[ ]:

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
    for t in np.linspace(0.0,12.0,120):
        parms.append(
            { 
              'LATTICE'        : "fcc",
              'LATTICE_LIBRARY' : "fcc.xml", 
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
input_file = pyalps.writeInputFiles('parm7e',parms)
pyalps.runApplication('spinmc',input_file,Tmin=5)


pyalps.evaluateSpinMC(pyalps.getResultFiles(prefix='parm7e'))

#load the susceptibility and collect it as function of temperature T
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm7e'),'Binder Cumulant')

binder_u4 = pyalps.collectXY(data,x='T',y='Binder Cumulant',foreach=['L'])

Tc=9.79 

for a in np.linspace(1.5,2.8,13):
  s=binder_u4
  for d in s:
      d.x -= Tc
      d.x = d.x/Tc
      l = d.props['L']
      d.x = d.x * pow(float(l),a)
      plt.figure()
      pyalps.plot.plot(s)
      plt.xlabel('$L^a(T-T_c)/T_c, a=%.2f$' % a)
      plt.ylabel('Binderov kumulant U4 $g$')
      plt.title('3D Izingov model, BCC')
      plt.legend(loc='best')
      plt.savefig("figura_a_nu_FCC%d.eps"%(a*10-14),dpi=300)
#     


