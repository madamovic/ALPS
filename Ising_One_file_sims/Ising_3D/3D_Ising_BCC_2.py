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
plt.rcParams.update({'figure.max_open_warning': 0})



numeratorfigs=1

#prepare the input parameters
parms = []
for l in [8,12,14,16,24]: 
    for t in np.linspace(0.01,9.0,90):
        parms.append(
            { 
              'LATTICE'        : "bcc",
              'LATTICE_LIBRARY'        : "bcc.xml", 
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
input_file = pyalps.writeInputFiles('parm7d',parms)
pyalps.runApplication('spinmc',input_file,Tmin=5)
# use the following instead if you have MPI
#pyalps.runApplication('spinmc',input_file,Tmin=5,MPI=2)

pyalps.evaluateSpinMC(pyalps.getResultFiles(prefix='parm7d'))

#load the susceptibility and collect it as function of temperature T
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm7d'),['|Magnetization|', 'Connected Susceptibility', 'Specific Heat', 'Binder Cumulant', 'Binder Cumulant U2'])
magnetization_abs = pyalps.collectXY(data,x='T',y='|Magnetization|',foreach=['L'])
connected_susc = pyalps.collectXY(data,x='T',y='Connected Susceptibility',foreach=['L'])
spec_heat = pyalps.collectXY(data,x='T',y='Specific Heat',foreach=['L'])
binder_u4 = pyalps.collectXY(data,x='T',y='Binder Cumulant',foreach=['L'])
binder_u2 = pyalps.collectXY(data,x='T',y='Binder Cumulant U2',foreach=['L'])





#a_nu procena

Tc=6.4 #bestTc

n=90
llista=[8,12,14,16,24]



delta=[]
alista=[]


for a in np.linspace(1.0,4.0,40):
    binder_u4 = pyalps.collectXY(data,x='T',y='Binder Cumulant',foreach=['L'])
    listadelta=[]
    for d in binder_u4:
      d.x -= Tc
      d.x = d.x/Tc
      l = d.props['L']
      d.x = d.x * pow(float(l),a)
    lvrednost=np.array([q.props['L'] for q in binder_u4])
    sel=np.argsort(lvrednost)
    connected_susc=np.array(binder_u4)
    connected_susc=connected_susc[sel]
    s=open('binder_u4_BCC_redosled.txt','w')
    s.write(pyalps.plot.convertToText(binder_u4))
    s.close()
    file=open('binder_u4_BCC_redosled.txt')
    file_data=np.loadtxt(file,usecols=(0,1))
    x=file_data[:,0]
    y=file_data[:,1]
    for i in range(0,len(llista)):
        exec("x%d = x[i*n:i*n+n]" % (llista[i]));
    for j in range(0,len(llista)):
        exec("y%d = y[j*n:j*n+n]" % (llista[j]));
    for k in range(0,len(y8)):
        listadelta.append(abs(y8[k]-y12[k]))
    plt.figure()
    pyalps.plot.plot(s)
    plt.xlabel('$L^a(T-T_c)/T_c, a=%.2f$' % a)
    plt.ylabel('Binderov kumulant U4 $g$')
    plt.title('3D Izingov model, BCC')
    plt.legend(loc='best')
    plt.savefig("figure_2_BCC%d.eps"%(numeratorfigs),dpi=300)
    srednjadelta=sum(listadelta)/len(listadelta)
    delta.append(srednjadelta)
    alista.append(a)
    numeratorfigs+=1


abest=alista[np.argmin(delta)]

o=open('a_BCC.txt','w')
o.write(pyalps.plot.convertToText(alista[np.argmin(delta)])
o.close()
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


file=open('binderdata_rounded_t_redosled_3d_Ising_BCC.txt')

file_data=np.loadtxt(file,usecols=(0,1))


x=file_data[:,0]
y=file_data[:,1]




llista = [8,12,14,16,24]
n=90




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
tc=6.4

for p in range(0,len(llista)):
    exec("rez = derivative(funk%d,tc)" % (llista[p]));
    lista.append(rez)





plt.figure()
plt.plot(x12,y12,label='$U_{4}$',color='b')
plt.plot(xizv12,yizv12,label='$dU_{4}/dT$',color='r')
plt.legend(loc='best')
plt.savefig("figure_2_BCC%d.eps"%(numeratorfigs),dpi=300)

numeratorfigs+=1


params,params_covariance=optimize.curve_fit(test_funk,llista,lista)


plt.figure()
plt.scatter(llista,lista,label='Podaci',color='b')
plt.plot(llista,test_funk(llista,params[0],params[1]),label='Fit',color='r')
plt.xlabel('$L$')
plt.ylabel(r'$dU_{4}/dT|T_{C}\approx L^{1/\nu}$')
plt.title(r'$1/\nu=$ %.13s,$\nu=$ %.13s' % (params[1],1/params[1]))
plt.legend(loc='upper left')
plt.savefig("figure_2_BCC%d.eps"%(numeratorfigs),dpi=300)


numeratorfigs+=1

# OK OK OK OK OK OK OK OK OK
##############################################################################################
#
#
#

  
delta=[]
etalista=[]



for two_minus_eta in np.linspace(1.0,3.0,30):
    connected_susc = pyalps.collectXY(data,x='T',y='Connected Susceptibility',foreach=['L'])
    listadelta=[]
    for d in connected_susc:
        d.x -= Tc
        d.x = d.x/Tc
        l = d.props['L']
        d.x = d.x * pow(float(l),abest)
        d.y = d.y/pow(float(l),two_minus_eta)
    lvrednost=np.array([q.props['L'] for q in connected_susc])
    sel=np.argsort(lvrednost)
    connected_susc=np.array(connected_susc)
    connected_susc=connected_susc[sel]
    s=open('connectedsusc_BCC_redosled.txt','w')
    s.write(pyalps.plot.convertToText(connected_susc))
    s.close()
    file=open('connectedsusc_BCC_redosled.txt')
    file_data=np.loadtxt(file,usecols=(0,1))
    x=file_data[:,0]
    y=file_data[:,1]
    for i in range(0,len(llista)):
        exec("x%d = x[i*n:i*n+n]" % (llista[i]));
    for j in range(0,len(llista)):
        exec("y%d = y[j*n:j*n+n]" % (llista[j]));
    for k in range(0,len(y8)):
        listadelta.append(abs(y8[k]-y12[k]))
    plt.figure()
    pyalps.plot.plot(connected_susc)
    plt.xlabel('$L^a(T-T_c)/T_c, Tc=%.2f, a=%.2f$'%(Tc,abest))
    plt.ylabel(r'$L^{\gamma/\nu}\chi_c,\gamma/\nu=$ %.4s' % two_minus_eta)
    plt.title('3D Ising model, BCC')
    plt.legend(loc='best')
    plt.savefig("figure_2_BCC%d.eps"%(numeratorfigs),dpi=300)
    srednjadelta=sum(listadelta)/len(listadelta)
    delta.append(srednjadelta)
    etalista.append(two_minus_eta)
    numeratorfigs+=1


print(etalista[np.argmin(delta)])

o=open('eta_BCC.txt','w')
o.write(pyalps.plot.convertToText(etalista[np.argmin(delta)])
o.close()
# OK OK OK OK OK OK 
#############################################################################################
#
#
#
#
#
#

#make a data collapse of the |magnetization| as a function of (T-Tc)/Tc


# OK OK OK OK OK OK
#########################


delta=[]
betalista=[]



for beta_over_nu in np.linspace(0.3,0.7,40):
    magnetization_abs = pyalps.collectXY(data,x='T',y='|Magnetization|',foreach=['L'])
    listadelta=[]
    for d in magnetization_abs:
        d.x -= Tc
		d.x = d.x/Tc
		l = d.props['L']
		d.x = d.x * pow(float(l),abest)
        d.y = d.y / pow(float(l),-beta_over_nu)
    lvrednost=np.array([q.props['L'] for q in magnetization_abs])
    sel=np.argsort(lvrednost)
    magnetization_abs=np.array(magnetization_abs)
    magnetization_abs=magnetization_abs[sel]
    s=open('magnetization_abs_BCC_redosled.txt','w')
    s.write(pyalps.plot.convertToText(magnetization_abs))
    s.close()
    file=open('magnetization_abs_BCC_redosled.txt')
    file_data=np.loadtxt(file,usecols=(0,1))
    x=file_data[:,0]
    y=file_data[:,1]
    for i in range(0,len(llista)):
        exec("x%d = x[i*n:i*n+n]" % (llista[i]));
    for j in range(0,len(llista)):
        exec("y%d = y[j*n:j*n+n]" % (llista[j]));
    for k in range(0,len(y8)):
        listadelta.append(abs(y8[k]-y12[k]))
    plt.figure()
    pyalps.plot.plot(magnetization_abs)
    plt.xlabel('$L^a(T-T_c)/T_c, Tc=%.3f, a=%.2f$'%(Tc,abest))
    plt.ylabel(r'Magnetizacija $|m|L^\beta/\nu, \beta/\nu=$ %.4s' % beta_over_nu)
    plt.title('3D Izingov model, BCC')
    plt.savefig("figure_2_BCC%d.eps"%(numeratorfigs),dpi=300)
    srednjadelta=sum(listadelta)/len(listadelta)
    delta.append(srednjadelta)
    betalista.append(beta_over_nu)
    numeratorfigs+=1


print(etalista[np.argmin(delta)])

o=open('beta_BCC.txt','w')
o.write(pyalps.plot.convertToText(betalista[np.argmin(delta)])
o.close()
