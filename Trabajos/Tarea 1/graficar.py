import numpy as np
import matplotlib.pyplot as plt
data=np.loadtxt('data.txt')
t=data[:,0]
e1=data[:,1]
e2=data[:,2]
e3=data[:,3]
plt.scatter(t,e1, s=4, c='b', lw=0, label='$E_1$')
plt.scatter(t,e2,s=4, c='red', lw=0, label='$E_2$')
plt.scatter(t,e3,s=4, c='green', lw=0, label='$E_3$')
plt.xlim(0,50000)
plt.ylim(0,0.04)
plt.xlabel("$Tiempo\ (seg)$")
plt.ylabel("$E_k$")
plt.legend(loc=0)
plt.savefig('evst.pdf')
