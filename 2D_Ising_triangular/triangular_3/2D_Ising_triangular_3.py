import pyalps
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyalps.plot
#%matplotlib inline
import numpy as np
from scipy import optimize
from scipy import interpolate
import pyalps.fit_wrapper as fw


parms = []
for l in [8,12,24,48,64,128]: 
    for t in [3.56,3.57,3.58,3.59,3.60,3.61,3.62,3.63,3.64,3.65,3.66]:
        parms.append(
            { 
              'LATTICE'        : "triangular lattice", 
              'T'              : t,
              'J'              : 1 ,
              'THERMALIZATION' : 20000,
              'SWEEPS'         : 100000,
              'UPDATE'         : "cluster",
              'MODEL'          : "Ising",
              'L'              : l
            }
    )
#write the input file and run the simulation
input_file = pyalps.writeInputFiles('parm7b',parms)
pyalps.runApplication('spinmc',input_file,Tmin=5)
# use the following instead if you have MPI
#pyalps.runApplication('spinmc',input_file,Tmin=5,MPI=2)

pyalps.evaluateSpinMC(pyalps.getResultFiles(prefix='parm7b'))

#load the susceptibility and collect it as function of temperature T
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm7b'),['|Magnetization|', 'Connected Susceptibility', 'Specific Heat', 'Binder Cumulant', 'Binder Cumulant U2'])
magnetization_abs = pyalps.collectXY(data,x='T',y='|Magnetization|',foreach=['L'])
connected_susc = pyalps.collectXY(data,x='T',y='Connected Susceptibility',foreach=['L'])
spec_heat = pyalps.collectXY(data,x='T',y='Specific Heat',foreach=['L'])
binder_u4 = pyalps.collectXY(data,x='T',y='Binder Cumulant',foreach=['L'])
binder_u2 = pyalps.collectXY(data,x='T',y='Binder Cumulant U2',foreach=['L'])


numeratorfigs=1

#
#make a fit of the connected susceptibility as a function of L:
cs_mean=[]
for q in connected_susc:
    cs_mean.append(np.array([d.mean for d in q.y]))


peak_cs = pyalps.DataSet()
peak_cs.props = pyalps.dict_intersect([q.props for q in connected_susc])
peak_cs.y = np.array([np.max(q) for q in cs_mean])
peak_cs.x = np.array([q.props['L'] for q in connected_susc])
# 



sel = np.argsort(peak_cs.x)


peak_cs.y = peak_cs.y[sel]
peak_cs.x = peak_cs.x[sel]

listafit=[]

for i in range(0,len(peak_cs.y)):
    listafit.append(np.array([peak_cs.x[i],peak_cs.y[i]]))

np.savetxt('susceptibilnost_podaci_fitovanje_2D_Ising_triangular.txt',listafit,delimiter=' ')

pars = [fw.Parameter(1), fw.Parameter(1)]


f = lambda self, x, pars: pars[0]()*np.power(x,pars[1]())
fw.fit(None, f, pars, peak_cs.y, peak_cs.x)
prefactor = pars[0].get()
gamma_nu = pars[1].get()


# 
plt.figure()
plt.plot(peak_cs.x, f(None, peak_cs.x, pars))
pyalps.plot.plot(peak_cs)
plt.xlabel('$L$')
plt.ylabel('Susceptibilnost $\chi_c(T_c)$')
plt.title(r'2D Izingov model, trougaona resetka, $\gamma/\nu=$ %.4s' % gamma_nu)
plt.savefig("figure_triangular_fit%d.eps"%(numeratorfigs),dpi=300)

numeratorfigs+=1

# OK OK OK OK OK OK 
#############################################
## SCIPY FIT


def test_funk(x,a,b):
    return a*x**b

params,params_covariance=optimize.curve_fit(test_funk,peak_cs.x,peak_cs.y)


# 
plt.figure()
plt.scatter(peak_cs.x,peak_cs.y,label='Podaci',color='b')
plt.plot(peak_cs.x,f(None, peak_cs.x, pars),label='ALPS',color='g')
plt.plot(peak_cs.x,test_funk(peak_cs.x,params[0],params[1]),label='Scipy',color='r')
plt.xlabel('$L$')
plt.ylabel('Susceptibilnost $\chi_c(T_c)$')
plt.title(r'$\gamma/\nu_{ALPS}=$ %.13s, $\gamma/\nu_{SCIPY}=$ %.13s' % (gamma_nu,params[1]))
plt.legend(loc='upper left')
plt.savefig("figure_triangular_fit%d.eps"%(numeratorfigs),dpi=300)

numeratorfigs+=1

# OK OK OK OK OK OK OK
#################################################################################3
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
plt.title(r'2D Izingov model, trougaona resetka, $\alpha/\nu=$ %.4s' % alpha_nu)
plt.savefig("figure_triangular_fit%d.eps"%(numeratorfigs),dpi=300)

numeratorfigs+=1




# OK OK OK OK OK OK OK OK OK OK
############################################
# SCIPY FIT



params,params_covariance=optimize.curve_fit(test_funk,peak_sh.x,peak_sh.y)
# 
plt.figure()
plt.scatter(peak_sh.x,peak_sh.y,label='Podaci',color='b')
plt.plot(peak_sh.x,f(None, peak_sh.x, pars),label='ALPS',color='g')
plt.plot(peak_sh.x,test_funk(peak_sh.x,params[0],params[1]),label='Scipy',color='r')
plt.xlabel('$L$')
plt.ylabel('Specificna toplota $c_v(T_c)$')
plt.title(r'$\alpha/\nu_{ALPS}=$ %.13s, $\alpha/\nu_{SCIPY}=$ %.13s' % (alpha_nu,params[1]))
plt.legend(loc='upper left')
plt.savefig("figure_triangular_fit%d.eps"%(numeratorfigs),dpi=300)

numeratorfigs+=1
##
# OK OK OK OK OK OK OK OK OK 
#################################################################################################
