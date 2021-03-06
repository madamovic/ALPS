
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
    for t in [6.28,6.29,6.30,6.31,6.32,6.33,6.34,6.35,6.36,6.37,6.38,6.39,6.40,6.41]:
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
input_file = pyalps.writeInputFiles('parm7l',parms)
pyalps.runApplication('spinmc',input_file,Tmin=5)
# use the following instead if you have MPI
#pyalps.runApplication('spinmc',input_file,Tmin=5,MPI=4)

pyalps.evaluateSpinMC(pyalps.getResultFiles(prefix='parm7l'))

#load the susceptibility and collect it as function of temperature T
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm7l'), 'Specific Heat')
spec_heat = pyalps.collectXY(data,x='T',y='Specific Heat',foreach=['L'])


# 
#make a fit of the specific heat as a function of L:
sh_mean=[]
for q in spec_heat:
    sh_mean.append(np.array([d.mean for d in q.y]))
#     
peak_sh = pyalps.DataSet()
peak_sh.props = pyalps.dict_intersect([q.props for q in spec_heat])
peak_sh.y = np.array([np.max(q) for q in sh_mean])
peak_sh.x = np.array([q.props['L'] for q in spec_heat])
 
sel = np.argsort(peak_sh.x)
peak_sh.y = peak_sh.y[sel]
peak_sh.x = peak_sh.x[sel]
# 
pars = [fw.Parameter(1), fw.Parameter(1)]
f = lambda self, x, pars: pars[0]()*np.power(x,pars[1]())
fw.fit(None, f, pars, peak_sh.y, peak_sh.x)
prefactor = pars[0].get()
alpha_nu = pars[1].get()
# 
plt.figure()
plt.plot(peak_sh.x, f(None, peak_sh.x, pars))
pyalps.plot.plot(peak_sh)
plt.xlabel('$L$')
plt.ylabel('Specificna toplota $c_v(T_c)$')
plt.title(r'3D Izingov model, BCC, $\alpha=$ %.4s' % alpha_nu)
plt.savefig("figure59.eps",dpi=300)


