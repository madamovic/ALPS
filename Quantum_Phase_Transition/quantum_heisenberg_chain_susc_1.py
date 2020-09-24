import pyalps
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import pyalps.plot 
import numpy as np 

parms=[]

for t in np.linspace(0.01,2.0,200):
	parms.append(
		{
		'LATTICE': "chain lattice",
		'MODEL': "spin",
		'local_S': 0.5,
		'T': t,
		'J': 1,
		'THERMALIZATION': 10000,
		'SWEEPS': 50000,
		'L':60,
		'ALGORITHM': "loop"
		}
	)

input_file=pyalps.writeInputFiles('parm2a',parms)
pyalps.runApplication('loop',input_file)

data=pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm2a'),'Susceptibility')
susceptibility=pyalps.collectXY(data,x='T',y='Susceptibility')

plt.figure()
pyalps.plot.plot(susceptibility)
plt.xlabel('Temperature $T/J$')
plt.ylabel(r'Susceptibility $\chi J$')
plt.title('Quantum Heisenberg chain, $L=60$, $S=1/2$')
plt.savefig('figure_quantum_chain_susc_1.eps',dpi=400)
#plt.show()
