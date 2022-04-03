#Programa para plotear puntos

import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

data=np.loadtxt('energiasplanetas.txt')

t=data[:,0]
mercurio=data[:,2]
venus=data[:,3]
tierra=data[:,4]
marte=data[:,5]
jupiter=data[:,5]
saturno=data[:,6]
urano=data[:,7]
neptuno=data[:,8]
total=data[:,9]

fig=plt.figure()
ax=fig.add_subplot(111)
ax.grid(True)


y_min = 0
y_max = 0.000005

ax.set_ylim(y_min, y_max)

ax.set_ylabel("Energia", fontsize=14, fontname="Times New Roman")
ax.set_xlabel("Tiempo", fontsize=14, fontname="Times New Roman")
#ax.plot(t,mercurio, linestyle='dashed', linewidth=2, color='purple')
#ax.plot(t,venus, linestyle='dashed', linewidth=2, color='brown')
#ax.plot(t,tierra, linestyle='dashed', linewidth=2, color='green')
#ax.plot(t,marte, linestyle='dashed', linewidth=2, color='red')
#ax.plot(t,jupiter, linestyle='dashed', linewidth=2, color='orange')
#ax.plot(t,saturno, linestyle='dashed', linewidth=2, color='olive')
#ax.plot(t,urano, linestyle='dashed', linewidth=2, color='blue')
#ax.plot(t,neptuno, linestyle='dashed', linewidth=2, color='cyan')
ax.plot(t,total, linestyle='dashed', linewidth=2, color='black')
fig.savefig('energias.png',dpi=300)
plt.show()