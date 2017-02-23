import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt('t_vs_np.dat')
x=np.linspace(1,200,200)
y=5.55/x+0.45
plt.plot(x,y) 
plt.plot(data[:,0],data[:,1],'ro')
plt.xscale('log')
plt.xlabel(r'Number of Processors', fontsize=17)
plt.ylabel(r'Execution Time (s)', fontsize=17)
plt.axis([0.7,300,0,6])
plt.title(r'Execution Time vs. Number of Processors')
plt.savefig('t_vs_np.png')
plt.clf()


