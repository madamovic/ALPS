
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
for l in [2,4,6,8,10,12]:
    for t in np.linspace(0.0,6.0,60):
        parms.append(
            { 
              'LATTICE'        : "simple cubic lattice", 
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
input_file = pyalps.writeInputFiles('parm7i',parms)
pyalps.runApplication('spinmc',input_file,Tmin=5)
# use the following instead if you have MPI
#pyalps.runApplication('spinmc',input_file,Tmin=5,MPI=4)

pyalps.evaluateSpinMC(pyalps.getResultFiles(prefix='parm7i'))

#load the susceptibility and collect it as function of temperature T
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm7i'),'Binder Cumulant')

binder_u4 = pyalps.collectXY(data,x='T',y='Binder Cumulant',foreach=['L'])

# 
for d in binder_u4:
    d.x = np.around(d.x,1)


f = open('binderdata_rounded_t_3d_Ising.txt','w')
f.write(pyalps.plot.convertToText(binder_u4))
f.close()
