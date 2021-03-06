
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
    for t in [9.70,9.71,9.72,9.73,9.74,9.75,9.76,9.77,9.78,9.79,9.80,9.81,9.82,9.83,9.84,9.85,9.86,9.87,9.88]:
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
input_file = pyalps.writeInputFiles('parm7p',parms)
pyalps.runApplication('spinmc',input_file,Tmin=5)
# use the following instead if you have MPI
#pyalps.runApplication('spinmc',input_file,Tmin=5,MPI=4)

pyalps.evaluateSpinMC(pyalps.getResultFiles(prefix='parm7p'))

#load the susceptibility and collect it as function of temperature T
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm7p'),'|Magnetization|')
magnetization_abs = pyalps.collectXY(data,x='T',y='|Magnetization|',foreach=['L'])

Tc=9.79
a=2.2

#make a data collapse of the |magnetization| as a function of (T-Tc)/Tc
for d in magnetization_abs:
    d.x -= Tc
    d.x = d.x/Tc
    l = d.props['L']
    d.x = d.x * pow(float(l),a)


beta_over_nu=0.715 #your estimate    

for d in magnetization_abs:
    l = d.props['L']
    d.y = d.y / pow(float(l),-beta_over_nu)
#

plt.figure()
pyalps.plot.plot(magnetization_abs)
plt.xlabel('$L^a T-Tc/Tc')
plt.ylabel(r'Magnetizacija $|m|L^\beta/\nu, \beta/\nu=$ %.4s' % beta_over_nu)
plt.title('3D Izingov model, FCC')
plt.savefig("figure86.eps",dpi=300)

