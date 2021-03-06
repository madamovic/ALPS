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
from scipy import optimize
from scipy import interpolate
import pyalps.fit_wrapper as fw



numeratorfigs=1

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
              'SWEEPS'         : 100000,
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
plt.savefig("figure_SC%d.eps"%(numeratorfigs),dpi=300)

numeratorfigs+=1

plt.figure()
pyalps.plot.plot(connected_susc)
plt.xlabel('Temperatura $T$')
plt.ylabel('Susceptibilnost $\chi$')
plt.title('3D Izingov model')
plt.legend(loc='best')
plt.savefig("figure_SC%d.eps"%(numeratorfigs),dpi=300)

numeratorfigs+=1

plt.figure()
pyalps.plot.plot(spec_heat)
plt.xlabel('Temperatura $T$')
plt.ylabel('Specificna toplota $c_v$')
plt.title('3D Izingov model')
plt.legend(loc='best')
plt.savefig("figure_SC%d.eps"%(numeratorfigs),dpi=300)


numeratorfigs+=1

plt.figure()
pyalps.plot.plot(binder_u4)
plt.xlabel('Temperatura $T$')
plt.ylabel('Binderov kumulant U4 $g$')
plt.title('3D Izingov model')
plt.legend(loc='best')
plt.savefig("figure_SC%d.eps"%(numeratorfigs),dpi=300)

numeratorfigs+=1

plt.figure()
pyalps.plot.plot(binder_u2)
plt.xlabel('Temperatura $T$')
plt.ylabel('Binderov kumulant U2 $g$')
plt.title('3D Izingov model')
plt.legend(loc='best')
plt.savefig("figure_SC%d.eps"%(numeratorfigs),dpi=300)

numeratorfigs+=1

# OK OK OK OK OK OK 
###################################################################


f = open('binderdata_3d_Ising_SC.txt','w')
f.write(pyalps.plot.convertToText(binder_u4))
f.close()


# OK OK OK OK OK OK
#############################################################
# ROUND
r=binder_u4

for d in r:
    d.x = np.around(d.x,1)


fg = open('binderdata_rounded_t_3d_Ising_SC.txt','w')
fg.write(pyalps.plot.convertToText(r))
fg.close()

# OK OK OK OK OK OK 
####################################################################
# REDOSLED

red=binder_u4

for d in red:
    d.x=np.around(d.x,1)


fh=open('binderdata_rounded_t_3d_Ising_SC.txt','w')
fh.write(pyalps.plot.convertToText(red))
fh.close()


lvrednost=np.array([q.props['L'] for q in red])
print(lvrednost)


sel=np.argsort(lvrednost)
print(sel)

red=np.array(red)

red=red[sel]


s=open('binderdata_rounded_t_redosled_3d_Ising_SC.txt','w')
s.write(pyalps.plot.convertToText(red))
s.close()



# OK OK OK OK OK OK OK OK 
###################################################################################################################################################################
# Tc procena




for Tc in np.linspace(4.0,5.0,100):
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
      plt.savefig('figure_SC%d.eps'%(numeratorfigs),dpi=300)
  numeratorfigs+=1


# OK OK OK OK OK OK OK 
#################################################################################3
# a_nu procena

Tc=4.51

for a in np.linspace(1.0,4.0,40):
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
      plt.savefig("figure_SC%d.eps"%(numeratorfigs),dpi=300)
  numeratorfigs+=1

# OK OK OK OK OK OK OK OK OK 
#############################################################################33



def test_funk(x,a,b):
    return a*x**b

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


file=open('binderdata_rounded_t_redosled_3d_Ising_SC.txt')

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
plt.savefig("figure_SC%d.eps"%(numeratorfigs),dpi=300)

numeratorfigs+=1


params,params_covariance=optimize.curve_fit(test_funk,llista,lista)


plt.figure()
plt.scatter(llista,lista,label='Podaci',color='b')
plt.plot(llista,test_funk(llista,params[0],params[1]),label='Fit',color='r')
plt.xlabel('$L$')
plt.ylabel(r'$dU_{4}/dT|T_{C}\approx L^{1/\nu}$')
plt.title(r'$1/\nu=$ %.13s,$\nu=$ %.13s' % (params[1],1/params[1]))
plt.legend(loc='upper left')
plt.savefig("figure_SC%d.eps"%(numeratorfigs),dpi=300)


numeratorfigs+=1

# OK OK OK OK OK OK OK OK OK
##############################################################################################
#
#
#
#
#
#
Tc=4.51
a=1.6

#make a data collapse of the connected susceptibility as a function of (T-Tc)/Tc:
suscol=connected_susc

for d in suscol:
    d.x -= Tc
    d.x = d.x/Tc
    l = d.props['L']
    d.x = d.x * pow(float(l),a)

g=suscol

for two_minus_eta in np.linspace(1.0,3.0,30):
  suscol=g
  for d in suscol:
      l = d.props['L']
      d.y = d.y/pow(float(l),two_minus_eta)
      plt.figure()
      pyalps.plot.plot(suscol)
      plt.xlabel('$L^a(T-T_c)/T_c, Tc=%.3f, a=%.2f$'%(Tc,a))
      plt.ylabel(r'$L^{\gamma/\nu}\chi_c,\gamma/\nu=$ %.4s' % two_minus_eta)
      plt.title('3D Ising model')
      plt.savefig("figure_SC%d.eps"%(numeratorfigs),dpi=300)
  numeratorfigs+=1


# OK OK OK OK OK OK 
#############################################################################################
#
#
#
#
#
#
Tc=4.51
a=1.6

#make a data collapse of the |magnetization| as a function of (T-Tc)/Tc
magcol=magnetization_abs

for d in magcol:
    d.x -= Tc
    d.x = d.x/Tc
    l = d.props['L']
    d.x = d.x * pow(float(l),a)

h=magcol

for beta_over_nu in np.linspace(0.5,0.6,100):
  magcol=h  
  for d in magcol:
      l = d.props['L']
      d.y = d.y / pow(float(l),-beta_over_nu)
      plt.figure()
      pyalps.plot.plot(magnetization_abs)
      plt.xlabel('$L^a(T-T_c)/T_c, Tc=%.3f, a=%.2f$'%(Tc,a))
      plt.ylabel(r'Magnetizacija $|m|L^\beta/\nu, \beta/\nu=$ %.4s' % beta_over_nu)
      plt.title('3D Izingov model')
      plt.savefig("figure_SC%d.eps"%(numeratorfigs),dpi=300)
  numeratorfigs+=1


# OK OK OK OK OK OK
#########################
