import pyalps
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import pyalps.plot 
import numpy as np 

parms=[]

for l in [4,8,16,24,48,64]:
	for t in np.linspace(0.01,2.0,200):
		parms.append(
			{
			'LATTICE': "chain lattice",
			'MODEL': "spin",
			'local_S': 1.0,
			'T': t,
			'J': 1,
			'THERMALIZATION': 10000,
			'SWEEPS': 50000,
			'L':l,
			'ALGORITHM': "loop"
			}
		)

input_file=pyalps.writeInputFiles('parm2a',parms)
pyalps.runApplication('loop',input_file)

data=pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm2a'),'Susceptibility')
susceptibility=pyalps.collectXY(data,x='T',y='Susceptibility',foreach=['L'])


plt.figure()
pyalps.plot.plot(susceptibility)
plt.xlabel('Temperature $T/J$')
plt.ylabel(r'Susceptibility $\chi J$')
plt.title('Quantum Heisenberg chain, $S=1$')
plt.legend(loc='best')
plt.savefig('figure_quantum_chain_susc_3.eps',dpi=400)
#plt.show()
