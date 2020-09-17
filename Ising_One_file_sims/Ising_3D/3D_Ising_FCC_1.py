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
for l in [2,4,6,8,10,12]: 
    for t in np.linspace(0.01,6.0,60):
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



sel=np.argsort(lvrednost)


red=np.array(red)

red=red[sel]


s=open('binderdata_rounded_t_redosled_3d_Ising_SC.txt','w')
s.write(pyalps.plot.convertToText(red))
s.close()



# OK OK OK OK OK OK OK OK 
###################################################################################################################################################################
# Tc procena




for Tc in np.linspace(4.0,5.0,100):
  binder_u4 = pyalps.collectXY(data,x='T',y='Binder Cumulant',foreach=['L'])
  for d in binder_u4:
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
