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
input_file = pyalps.writeInputFiles('parm7k',parms)
pyalps.runApplication('spinmc',input_file,Tmin=5)
# use the following instead if you have MPI
#pyalps.runApplication('spinmc',input_file,Tmin=5,MPI=4)

pyalps.evaluateSpinMC(pyalps.getResultFiles(prefix='parm7k'))

#load the susceptibility and collect it as function of temperature T
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm7k'),'Connected Susceptibility')

connected_susc = pyalps.collectXY(data,x='T',y='Connected Susceptibility',foreach=['L'])

 
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

np.savetxt('susceptibilnost_podaci_fitovanje_3d_Ising_BCC.txt',listafit,delimiter=' ')

pars = [fw.Parameter(1), fw.Parameter(1)]


f = lambda self, x, pars: pars[0]()*np.power(x,pars[1]())
fw.fit(None, f, pars, peak_cs.y, peak_cs.x)
prefactor = pars[0].get()
gamma_nu = pars[1].get()



from scipy import optimize

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
plt.title('$\gamma_{ALPS}=$ %.13s, $\gamma_{SCIPY}=$ %.13s' % (gamma_nu,params[1]))
plt.legend(loc='upper left')
plt.savefig("figure66.eps",dpi=300)




