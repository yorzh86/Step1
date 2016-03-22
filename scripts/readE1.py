import matplotlib
import matplotlib.pyplot as plt
#from matplotlib import style
import numpy as np

#style.use('ggplot') #need 1.4 matplotlib

def readE(fn):
    f = open(fn)
    lines = f.readlines()
    f.close()
    t = []
    ts = []
    ke1 = []
    ke2 = []

    for k in range(2,len(lines)):
        l = map(float,lines[k].split())
        t.append(l[0])
        ts.append(l[1])
        ke1.append(l[2])
        ke2.append(l[3])

    t = np.array(t)
    ts = np.array(ts)
    ke1 = np.array(ke1)
    ke2 = np.array(ke2)

    return t, ts, abs(ke1-ke2)

def readT(fn):
    f = open(fn)
    lines = f.readlines()
    f.close()
    temps=[]
    t = []
    ts = []
    
    for k in range (2, len(lines)):
        l = map(float,lines[k].split())
        t.append(l[0])
        ts.append(l[1])
        temps.append(l[2:])

    t = np.array(t)
    ts = np.array(ts)
    temps = np.array(temps)

    return t, ts, temps

t, ts, kdif = readE('mark1.energies')
t1, ts1, temps = readT('mark1.temps')

A = 5.0*5.0 # [Angstroms]

x = t
y = kdif/0.01/A #ts = 10 fs = 0.01ps

fig = plt.figure()

ax1 = plt.subplot2grid((9,1), (0,0), rowspan=3)
ax2 = plt.subplot2grid((9,1), (3,0), rowspan=3)
ax3 = plt.subplot2grid((9,1), (6,0), rowspan=3)

plt.xlabel('time [ps]')
plt.ylabel('Flux [eV/ps/A^2]')
ax1.plot(x,y, label='Flux vs time', color = 'darkblue')
ax1.legend()
#ax1.set_xticks(np.arange(0,6,0.5))

x1=t1
xx1=0.0
y2=[]
y3=[]
for i in range(len(temps)):
    y2.append(temps[i][5]-temps[i][0])
    y3.append(temps[i][5]-temps[i][9])                      
    
ax2.plot(x1, y2, label='Gradiant #1', color='darkred')
ax2.plot(x1, y3, label='Gradiant #2', color='darkgreen')
ax2.legend()

x2 = t1
y4 = []
y5 = []
for i in range(len(temps)):
    y4.append(temps[i][5])
    y5.append(temps[i][0])

ax3.scatter(x1, y4, label='Hot', color='red',s=1)
ax3.scatter(x1, y5, label='Cold', color='blue', s=1)
ax3.legend()

ax1.grid()
ax2.grid()
ax3.grid()

fig.tight_layout()
plt.show()       
