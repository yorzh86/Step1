#!/usr/bin/env python

#todo: smoothing lines, dash line where we calculate things,
#      mirroring for plot2, and calculate conductivity.

import pylab as pl
from scipy.signal import savgol_filter
from matplotlib import style

def readE(fn):
	t, ts, ke1, ke2 = pl.loadtxt(fn,unpack=True,skiprows=2)
	return t,ts,abs(ke1-ke2)

def readT(fn):
	D = pl.loadtxt(fn,skiprows=2)
	t = D[:,0]
	ts = D[:,1]
	temps = D[:,2:]
	TC = temps[:,0]
	TH = temps[:,5]
	return t,ts,TC,TH,temps

def plotFlux(ax,t,q):
	ax.axvline(3, color='k', linestyle='--')
	ax.plot(t/1000, q,'m.', alpha = 0.6)
	ax.plot(t/1000,savgol_filter(q,11,3), color='darkviolet', linewidth=1.5)
	ax.set_xlabel('Time $t$ [ns]')
	ax.set_ylabel('Flux $q$ [eV/ps]')

def plotHC(ax,t,H,C):
	ax.axvline(3, color='k', linestyle='--')
	ax.plot(t/1000, H,'r.', label='Hot',  alpha = 0.1)
	ax.plot(t/1000, C,'c.', label='Cold', alpha = 0.1)
	ax.plot(t/1000,savgol_filter(H,101,3), color='darkred', linewidth=1.5)
	ax.plot(t/1000,savgol_filter(C,101,3), color='darkblue', linewidth=1.5)
	ax.set_xlabel('Time $t$ [ns]')
	ax.set_ylabel('Temperature $T$ [K]')
	ax.legend(loc='best',fancybox=True)

def plotG(ax,T):
	#v = bu[N:,:].sum(0)/bu[N:,:].shape[0]
	#e = bu[N:,:].std(0)
	#X,L,R = mirror(x,v)
	#X,Le,Re = mirror(x,e)
	#ax.errorbar(X,L,Le)
	#ax.errorbar(X,R,Re)
	
	
	
	
	
	a = 5.4
	z = pl.linspace(0,150*a,T.shape[1])
	#ax.errorbar(z/10,T.mean(0),T.std(0)/pl.sqrt(T.shape[0]),color='k')
	ax.plot(z/10, T.mean(0), 'darkblue')
	#ax.plot(z/10, T.mean(0) 'k')
	
	ax.set_xlim(z.min()/10,z.max()/10)
	ax.set_xlabel('Position $z$ [nm]')
	ax.set_ylabel('Temperature $T$ [K]')	

def plotC(ax,t,Ts):
	a = 5.4
	z = pl.linspace(0,150*a,Ts.shape[1])
	#pl.pcolormesh(z/10,t,Ts)
	pl.contourf(z/10,t,Ts,19)
	pl.contourf(z/10,t,Ts,20)
	#pl.colorbar()
	ax.set_xlim(z.min()/10,z.max()/10)
	ax.set_ylim(t.min(),t.max())
	ax.set_xlabel('Position $z$ [nm]')
	ax.set_ylabel('Time $t$ [ps]')

def doFigure(fnE,fnT):
	print 'Reading...'
	tE,tsE,q = readE(fnE)  #q is not flux, it's ke1-ke2
	tT,tsT,TC,TH,T = readT(fnT)
	print 'Done!'
	
	s = 1000
		
	fig = pl.figure()
	
	axC = fig.add_subplot(221)
	plotC(axC,tT[::s],T[::s,:])
	
	axG = fig.add_subplot(222)
	plotG(axG,T)
	
	axHC = fig.add_subplot(223)
	plotHC(axHC,tT[::s],TH[::s],TC[::s])
	
	axQ = fig.add_subplot(224)
	plotFlux(axQ,tE[::s],q[::s]/0.01) #10 fs = 0.01 ps
	
	fig.tight_layout()
	pl.subplots_adjust(top=0.91, hspace=0.46)
	fig.suptitle('Reverse non equilibrium method', fontsize=14, fontweight='bold')
	fig.savefig('plot.pdf')
	return q, T 

def mirror(x,u):
	N = x.size/2+1
	l = u[0:N]
	r = list(u[-1:-N:-1])
	r.insert(0,u[0])
	r = array(r)
	return x[0:N],l,r	


J2eV = 1.60217646E-19
s2ps = 1.0E-12
m2A  = 1.0E-10

# getting KE difference and temperatures:
f, T = doFigure('mark2.energies','mark2.temps')
# average flux after 3ns in [eV/ps ]:
flux = pl.mean(f[30000:]/0.01)
# mean of two gradients:
grad = ((T.mean(0)[5]-T.mean(0)[0]) + (T.mean(0)[5] - T.mean(0)[9]))/2.0
# Area [A^2]:
A = (5*5.4)*(5*5.4)
k = flux/A/grad
kSI = k*J2eV/(s2ps*m2A)

pl.title( 'kSI=%f W/m2K'%kSI)

pl.show()
