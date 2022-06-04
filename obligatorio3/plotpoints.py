#Programa para plotear puntos

import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

data=np.loadtxt('evolucionnorma.txt')

t=data[:,0]
solap1=data[:,1]
#error1=data[:,2]
#solap2=data[:,3]
#error2=data[:,4]
#solap3=data[:,5]
#error3=data[:,6]
#solap4=data[:,7]
#error4=data[:,8]


fig=plt.figure()
ax=fig.add_subplot(111)
ax.grid(True)


y_min = -0.2
y_max = 1.05

ax.set_ylim(y_min, y_max)

ax.set_ylabel("Tiempo", fontsize=14, fontname="Times New Roman")
ax.set_xlabel("Norma", fontsize=14, fontname="Times New Roman")

ax.plot(t,solap1, linestyle='dashed', marker='o',markersize=4, linewidth=2, color='orange')
#ax.plot(t,solap2, linestyle='dashed', linewidth=2, color='green', label='Patron 2')
#ax.plot(t,solap3, linestyle='dashed', linewidth=2, color='blue', label='Patron 3')
#ax.plot(t,solap4, linestyle='dashed', linewidth=2, color='orange', label='Patron 4')

ax.legend()

ax.hlines(y=1,xmin=0,xmax=1000, color='purple')

#plt.errorbar(t, solap1, yerr=error1,fmt=".", elinewidth=1,capthick=0,ecolor = 'lightgreen',color="orange")
#plt.errorbar(t, solap2, yerr=error2,fmt=".", elinewidth=1,capthick=0,ecolor = 'lightgreen',color="green")
#plt.errorbar(t, solap3, yerr=error3,fmt=".", elinewidth=1,capthick=0,ecolor = 'lightgreen',color="blue")
#plt.errorbar(t, solap4, yerr=error4,fmt=".", elinewidth=1,capthick=0,ecolor = 'lightgreen',color="orange")

fig.savefig('evolucionnorma.png',dpi=300)
plt.show()