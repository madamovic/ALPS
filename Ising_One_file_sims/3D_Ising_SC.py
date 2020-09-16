#
#3D Ising model on simple cubic lattice
#
#
import pyalps
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyalps.plot
#%matplotlib inline
import numpy as np



numeratorfigs=24

#prepare the input parameters
parms = []
for l in [2,4,6,8,10,12]: 
    for t in np.linspace(0.01,6.0,60):
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
input_file = pyalps.writeInputFiles('parm7a',parms)
pyalps.runApplication('spinmc',input_file,Tmin=5)
# use the following instead if you have MPI
#pyalps.runApplication('spinmc',input_file,Tmin=5,MPI=2)

pyalps.evaluateSpinMC(pyalps.getResultFiles(prefix='parm7a'))

#load the susceptibility and collect it as function of temperature T
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm7a'),['|Magnetization|', 'Connected Susceptibility', 'Specific Heat', 'Binder Cumulant', 'Binder Cumulant U2'])
magnetization_abs = pyalps.collectXY(data,x='T',y='|Magnetization|',foreach=['L'])
connected_susc = pyalps.collectXY(data,x='T',y='Connected Susceptibility',foreach=['L'])
spec_heat = pyalps.collectXY(data,x='T',y='Specific Heat',foreach=['L'])
binder_u4 = pyalps.collectXY(data,x='T',y='Binder Cumulant',foreach=['L'])
binder_u2 = pyalps.collectXY(data,x='T',y='Binder Cumulant U2',foreach=['L'])



#make plots
plt.figure()
pyalps.plot.plot(magnetization_abs)
plt.xlabel('Temperatura $T$')
plt.ylabel('Magnetizacija $|m|$')
plt.title('3D Izingov model')
plt.legend(loc='best')
plt.savefig("figure%d.eps"%(numeratorfigs),dpi=300)

numeratorfigs+=1

plt.figure()
pyalps.plot.plot(connected_susc)
plt.xlabel('Temperatura $T$')
plt.ylabel('Susceptibilnost $\chi$')
plt.title('3D Izingov model')
plt.legend(loc='best')
plt.savefig("figure%d.eps"%(numeratorfigs),dpi=300)

numeratorfigs+=1

plt.figure()
pyalps.plot.plot(spec_heat)
plt.xlabel('Temperatura $T$')
plt.ylabel('Specificna toplota $c_v$')
plt.title('3D Izingov model')
plt.legend(loc='best')
plt.savefig("figure%d.eps"%(numeratorfigs),dpi=300)


numeratorfigs+=1

plt.figure()
pyalps.plot.plot(binder_u4)
plt.xlabel('Temperatura $T$')
plt.ylabel('Binderov kumulant U4 $g$')
plt.title('3D Izingov model')
plt.legend(loc='best')
plt.savefig("figure%d.eps"%(numeratorfigs),dpi=300)

numeratorfigs+=1

plt.figure()
pyalps.plot.plot(binder_u2)
plt.xlabel('Temperatura $T$')
plt.ylabel('Binderov kumulant U2 $g$')
plt.title('3D Izingov model')
plt.legend(loc='best')
plt.savefig("figure%d.eps"%(numeratorfigs),dpi=300)

numeratorfigs+=1

# OK OK OK OK OK OK 
###################################################################


f = open('binderdata_3d_Ising.txt','w')
f.write(pyalps.plot.convertToText(binder_u4))
f.close()


# OK OK OK OK OK OK
#############################################################
# ROUND
r=binder_u4

for d in r:
    d.x = np.around(d.x,1)


f = open('binderdata_rounded_t_3d_Ising.txt','w')
f.write(pyalps.plot.convertToText(r))
f.close()

# OK OK OK OK OK OK 
####################################################################
# REDOSLED
red=binder_u4

for d in red:
    d.x=np.around(d.x,1)


f=open('binderdata_rounded_t_3d_Ising.txt','w')
f.write(pyalps.plot.convertToText(red))
f.close()


lvrednost=np.array([q.props['L'] for q in red])
print(lvrednost)


sel=np.argsort(lvrednost)
print(sel)


red=np.array(red)



red=red[sel]






s=open('binderdata_rounded_t_redosled_3d_Ising.txt','w')
s.write(pyalps.plot.convertToText(red))
s.close()




###################################################################################################################################################################
# Tc procena




for Tc in [4.511,3.511,5.511]:
  listfortc=binder_u4
  for d in listfortc:
      d.x -= Tc
      d.x = d.x/Tc  
      plt.figure()
      pyalps.plot.plot(binder_u4)
      plt.xlabel('$t=(T-T_c)/T_c, T_c=%.3f$'%(Tc)) 
      plt.ylabel('Binderov kumulant U4 $g$')
      plt.title('3D Izingov model')
      plt.legend(loc='best')
      plt.savefig('figure%d.eps'%(numeratorfigs),dpi=300)
  numeratorfigs+=1


#################################################################################3
# a_nu procena

Tc=4.511 

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
      plt.savefig("figura_a_nu_BCC%d.eps"%(a*10-14),dpi=300)


#############################################################################33

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

np.savetxt('susceptibilnost_podaci_fitovanje_3d_Ising.txt',listafit,delimiter=' ')

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
plt.title('2D Izingov model, $\gamma=$ %.4s' % gamma_nu)
plt.savefig("figure35.eps",dpi=300)



#############################################
## SCIPY FIT
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
plt.savefig("figure43.eps",dpi=300)


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
plt.title(r'3D Izingov model, $\alpha=$ %.4s' % alpha_nu)
plt.savefig("figure36.eps",dpi=300)





#####
# SCIPY FIT

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


##

#################################################################################################
#
# NUMERICAL DERIVATIVE
#
#
#
#
def derivative(f,a,method='central',h=0.1):
    if method == 'central':
        return (f(a + h) - f(a - h))/(2*h)
    elif method == 'forward':
        return (f(a + h) - f(a))/h
    elif method == 'backward':
        return (f(a) - f(a - h))/h
    else:
        raise ValueError("Method must be 'central', 'forward' or 'backward'.")


file=open('binderdata_rounded_t_redosled_3d_Ising.txt')

file_data=np.loadtxt(file,usecols=(0,1))


x=file_data[:,0]
y=file_data[:,1]




llista = [2,4,6,8,10,12]
n=60




for i in range(0,len(llista)):
    exec("x%d = x[i*n:i*n+n]" % (llista[i]));

for j in range(0,len(llista)):
    exec("y%d = y[j*n:j*n+n]" % (llista[j]));

#funk = interpolate.interp1d(x, y)

for k in range(0,len(llista)):
    exec("funk%d = interpolate.interp1d(x%d, y%d)" % (llista[k],llista[k],llista[k]));

for l in range(0,len(llista)):
    exec("xizv%d = np.arange(0.1,4,0.1)" % (llista[l]));



for m in range(0,len(llista)):
    exec("yizv%d = derivative(funk%d,xizv%d)" % (llista[m],llista[m],llista[m]));
    




lista=[]
tc=4.511

for p in range(0,len(llista)):
    exec("rez = derivative(funk%d,tc)" % (llista[p]));
    lista.append(rez)





plt.figure()
plt.plot(x12,y12,label='$U_{4}$',color='b')
plt.plot(xizv12,yizv12,label='$dU_{4}/dT$',color='r')
plt.legend(loc='best')
plt.savefig("figure45.eps",dpi=300)


from scipy import optimize

def test_funk(x,a,b):
    return a*x**b

params,params_covariance=optimize.curve_fit(test_funk,llista,lista)


plt.figure()
plt.scatter(llista,lista,label='Podaci',color='b')
plt.plot(llista,test_funk(llista,params[0],params[1]),label='Fit',color='r')
plt.xlabel('$L$')
plt.ylabel(r'$dU_{4}/dT|T_{C}\approx L^{1/\nu}$')
plt.title(r'$1/\nu=$ %.13s,$\nu=$ %.13s' % (params[1],1/params[1]))
plt.legend(loc='upper left')
plt.savefig("figure46.eps",dpi=300)



##############################################################################################
#
#
#
#
#
#
Tc=4.511
a=1.6

#make a data collapse of the connected susceptibility as a function of (T-Tc)/Tc:
for d in connected_susc:
    d.x -= Tc
    d.x = d.x/Tc
    l = d.props['L']
    d.x = d.x * pow(float(l),a)

two_minus_eta=1.96 #your estimate
for d in connected_susc:
    l = d.props['L']
    d.y = d.y/pow(float(l),two_minus_eta)

plt.figure()
pyalps.plot.plot(connected_susc)
plt.xlabel('$L^a(T-T_c)/T_c, Tc=4.511, a=1.6$')
plt.ylabel(r'$L^{\gamma/\nu}\chi_c,\gamma/\nu=$ %.4s' % two_minus_eta)
plt.title('3D Ising model')
plt.savefig("figure37.eps",dpi=300)



#############################################################################################
#
#
#
#
#
#
Tc=4.511
a=1.6

#make a data collapse of the |magnetization| as a function of (T-Tc)/Tc
for d in magnetization_abs:
    d.x -= Tc
    d.x = d.x/Tc
    l = d.props['L']
    d.x = d.x * pow(float(l),a)


beta_over_nu=0.522 #your estimate    

for d in magnetization_abs:
    l = d.props['L']
    d.y = d.y / pow(float(l),-beta_over_nu)
#

plt.figure()
pyalps.plot.plot(magnetization_abs)
plt.xlabel('$L^a T-Tc/Tc')
plt.ylabel(r'Magnetizacija $|m|L^\beta/\nu, \beta/\nu=$ %.4s' % beta_over_nu)
plt.title('3D Izingov model')
plt.savefig("figure40.eps",dpi=300)
