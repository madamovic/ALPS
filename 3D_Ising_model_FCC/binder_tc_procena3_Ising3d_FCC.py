
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
input_file = pyalps.writeInputFiles('parm7d',parms)
pyalps.runApplication('spinmc',input_file,Tmin=5)


pyalps.evaluateSpinMC(pyalps.getResultFiles(prefix='parm7d'))

#load the susceptibility and collect it as function of temperature T
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm7d'),'Binder Cumulant')

binder_u4 = pyalps.collectXY(data,x='T',y='Binder Cumulant',foreach=['L'])
# 
# #make a plot of the Binder cumulant:

# 
#perform a data collapse of the Binder cumulant: 
Tc=8.79 #your estimate
#a=1  #your estimate
# 
for d in binder_u4:
    d.x -= Tc
    d.x = d.x/Tc
#    l = d.props['L']
#    d.x = d.x * pow(float(l),a)
#     
plt.figure()
pyalps.plot.plot(binder_u4)
plt.xlabel('$t=(T-T_c)/T_c, T_c=8.79$')
plt.ylabel('Binderov kumulant U4 $g$')
plt.title('3D Izingov model, FCC')
plt.legend(loc='best')
plt.savefig('figure77.eps',dpi=300)


