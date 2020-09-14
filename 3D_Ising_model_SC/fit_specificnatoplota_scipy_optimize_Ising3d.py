
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
    for t in [4.45,4.46,4.47,4.48,4.49,4.50,4.51,4.52,4.53,4.54,4.55,4.56,4.57,4.58]:
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


from scipy import optimize

def test_funk(x,a,b):
    return a*x**b

params,params_covariance=optimize.curve_fit(test_funk,peak_sh.x,peak_sh.y)
# 
plt.figure()
plt.scatter(peak_sh.x,peak_sh.y,label='Podaci',color='b')
plt.plot(peak_sh.x,f(None, peak_sh.x, pars),label='ALPS',color='g')
plt.plot(peak_sh.x,test_funk(peak_sh.x,params[0],params[1]),label='Scipy',color='r')
plt.xlabel('$L$')
plt.ylabel('Specificna toplota $c_v(T_c)$')
plt.title(r'$\alpha_{ALPS}=$ %.13s, $\alpha_{SCIPY}=$ %.13s' % (alpha_nu,params[1]))
plt.legend(loc='upper left')
plt.savefig("figure44.eps",dpi=300)




