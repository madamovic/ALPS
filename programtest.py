import pyalps
import matplotlib.pyplot as plt
import pyalps.plot

parms = []
for l in [16]:
    for t in [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8,
8.5, 9, 9.5, 10, 11, 12, 13, 14, 15]:
        parms.append(
            {
              'LATTICE'        : "square lattice",
              'T'              : t,
              'S'              : 2.62,
              'J0'              : 0.52 ,
              'J1'              : -0.09 ,
              'J2'              : 0.86 ,
              'THERMALIZATION' : 1000,
              'SWEEPS'         : 400000,
              'UPDATE'         : "local",
              'MODEL'          : "Ising",
              'L'              : l
            }
    )

input_file = pyalps.writeInputFiles('Square_isg',parms)
pyalps.runApplication('spinmc',input_file,Tmin=10)
pyalps.evaluateSpinMC(pyalps.getResultFiles(prefix='Square_isg'))
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='Square_isg'),['|Magnetization|',
'Specific Heat'])
magnetization_abs = pyalps.collectXY(data,x='T',y='|Magnetization|',foreach=['L'])
spec_heat = pyalps.collectXY(data,x='T',y='Specific Heat',foreach=['L'])

#make plots
plt.figure()
pyalps.plot.plot(magnetization_abs)
plt.xlabel('Temperature $T$')
plt.ylabel('Magnetization $|m|$')
plt.legend()
plt.title('2D Ising model')

plt.figure()
pyalps.plot.plot(spec_heat)
plt.xlabel('Temperature $T$')
plt.ylabel('Specific Heat $c_v$')
plt.legend()
plt.title('2D Ising model')
plt.show()