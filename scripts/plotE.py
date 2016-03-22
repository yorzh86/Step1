#!/usr/bin/env python

import pylab as pl

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
	ax.plot(t/1000,q,'m.')
	ax.set_xlabel('Time $t$ [ns]')
	ax.set_ylabel('Flux $q$ [eV/ps]')

def plotHC(ax,t,H,C):
	ax.plot(t/1000,H,'r.',label='Hot')
	ax.plot(t/1000,C,'b.',label='Cold')
	ax.set_xlabel('Time $t$ [ns]')
	ax.set_ylabel('Temperature $T$ [K]')
	ax.legend(loc='best',fancybox=True)

def plotG(ax,T):
	a = 5.4
	z = pl.linspace(0,150*a,T.shape[1])
	ax.errorbar(z/10,T.mean(0),T.std(0)/pl.sqrt(T.shape[0]),color='k')
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
	tE,tsE,q = readE(fnE)
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
	plotFlux(axQ,tE[::s],q[::s])
	
	fig.tight_layout()
	fig.savefig('plot.pdf')

doFigure('mark1.energies','mark1.temps')
pl.show()