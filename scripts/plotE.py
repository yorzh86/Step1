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

def plotFlux(ax,t,q,S):
	i = S[0]/1000
	s = S[1]
	
	sfp = 11 # FIXME
	
	ax.axvline(3, color='k', linestyle='--')
	ax.plot(t[::s]/1000, q[::s],'m.', alpha = 0.6)
	ax.plot(t[::s]/1000,savgol_filter(q[::s],sfp,3), color='darkviolet', linewidth=1.5)
	ax.set_xlabel('Time $t$ [ns]')
	ax.set_ylabel('Flux $q$ [eV/ps]')
	
	dt = t[1]-t[0]
	qm = q[i:].mean()
	
	return qm/dt # eV/ps

def plotHC(ax,t,H,C,S):
	s = S[1]
	
	sfp = 101 # FIXME
	
	ax.axvline(3, color='k', linestyle='--')
	ax.plot(t[::s]/1000, H[::s],'r.', label='Hot',  alpha = 0.1)
	ax.plot(t[::s]/1000, C[::s],'c.', label='Cold', alpha = 0.1)
	ax.plot(t[::s]/1000,savgol_filter(H[::s],sfp,3), color='darkred', linewidth=1.5)
	ax.plot(t[::s]/1000,savgol_filter(C[::s],sfp,3), color='darkblue', linewidth=1.5)
	ax.set_xlabel('Time $t$ [ns]')
	ax.set_ylabel('Temperature $T$ [K]')
	ax.legend(loc='best',fancybox=True)

def plotG(ax,Ti,S):
	i = S[0]
	T = Ti[i:,:]
	
	a = 5.4
	z = pl.linspace(0,150*a,T.shape[1])
	
	Tl = pl.empty( (T.shape[0],T.shape[1]/2+1) ) # K
	Tr = pl.empty( (T.shape[0],T.shape[1]/2+1) ) # K
	zl = pl.empty( z.size/2+1 ) # A
	zr = pl.empty( z.size/2+1 ) # A
	
	Tl = T[:,:T.shape[1]/2+1]
	zl = z[:z.size/2+1]
	for k in xrange(z.size/2+1):
			Tr[:,k] = T[:,-k]
			zr[k] = z[k]

	Cl = pl.polyfit(zl,Tl.mean(0),1)
	Cr = pl.polyfit(zr,Tr.mean(0),1)
	
	grad = (Cl[0]+Cr[0])/2
	
	ax.errorbar(zl/10,Tl.mean(0),Tl.std(0)/pl.sqrt(Tl.shape[0]),color='r')
	ax.errorbar(zr/10,Tr.mean(0),Tr.std(0)/pl.sqrt(Tr.shape[0]),color='b')
	
	ax.plot(zl/10,pl.polyval(Cl,zl),'r--')
	ax.plot(zr/10,pl.polyval(Cr,zr),'b--')
		
	ax.set_xlim(zl.min()/10,zl.max()/10)
	ax.set_xlabel('Position $z$ [nm]')
	ax.set_ylabel('Temperature $T$ [K]')
	
	return grad

def plotC(ax,t,Ts,S):
	s = S[1]
	
	a = 5.4
	z = pl.linspace(0,150*a,Ts.shape[1])
	#pl.pcolormesh(z/10,t,Ts)
	pl.contourf(z/10,t[::s],Ts[::s,:],19)
	pl.contourf(z/10,t[::s],Ts[::s,:],20)
	#pl.colorbar()
	ax.set_xlim(z.min()/10,z.max()/10)
	ax.set_ylim(t[::s].min(),t[::s].max())
	ax.set_xlabel('Position $z$ [nm]')
	ax.set_ylabel('Time $t$ [ps]')

def doFigure(fnE,fnT,box):
	print 'Reading...'
	tE,tsE,q = readE(fnE)
	tT,tsT,TC,TH,T = readT(fnT)
	print 'Done!'
	
	S = (300000,1000)
		
	fig = pl.figure()
	
	axC = fig.add_subplot(221)
	plotC(axC,tT,T,S)
	
	axG = fig.add_subplot(222)
	grad = plotG(axG,T,S)
	
	axHC = fig.add_subplot(223)
	plotHC(axHC,tT,TH,TC,S)
	
	axQ = fig.add_subplot(224)
	heatRate = plotFlux(axQ,tE,q,S)
	
	A = box[0]*box[1]
	k = heatRate/(A*grad)
	
	fig.tight_layout()
	pl.subplots_adjust(top=0.91, hspace=0.46)
	fig.suptitle('Reverse non equilibrium method', fontsize=14, fontweight='bold')
	#fig.savefig('plot.pdf')
	return q, T, k

J2eV = 1.60217646E-19
s2ps = 1.0E-12
m2A  = 1.0E-10


a = 5.4
f, T, k2 = doFigure('mark2.energies','mark2.temps',a*pl.array([5,5,150]))
A = 5.4*5*5.4*5
flux = pl.mean(f[30000:]/0.1/A)
dz = 150*5.4/2
grad = ((T.mean(0)[5]-T.mean(0)[0]) + (T.mean(0)[5] - T.mean(0)[9]))/2.0/dz
k = flux/grad

kSI = k*J2eV/(s2ps*m2A)
k2SI = k2*J2eV/(s2ps*m2A)

pl.title( 'k=%f [W/m.K]; k=%f [W/m.K]'%(kSI,k2SI))
print 'k=', k2SI, '[W/m.K]'

pl.show()
