#!/usr/bin/env python

from __future__ import print_function
from pylab import *

def getData(fn,t0):
	k,t,T,TE,KE,PE,P = loadtxt(fn,unpack=True)
	s = argmin(abs(t-t0))
	D = {}
	D['k']  = k[s:]
	D['t']  = t[s:]
	D['T']  = T[s:]
	D['TE'] = TE[s:]
	D['KE'] = KE[s:]
	D['PE'] = PE[s:]
	D['P']  = P[s:]
	return D

def plotVar(vt,vn,u,M):
	print(vt)
	cm = polyfit(M['t'],M[vn],1)
	pm = poly1d(cm)
	print('Mark1: ',['m','b'],'=',cm,M[vn].std())
	print('')
	
	figure(figsize=(5,4))
	subplot(111)
	
	plot(M['t'],M[vn],'b-',alpha=0.3,label='Sim')
	plot(M['t'],pm(M['t']),'b--',label='Fit')
	
	xlim( M['t'].min() , M['t'].max() )
	
	xlabel(r'Time $t$ [ps]')
	ylabel(r'%s $%s$ [%s]'%(vt,vn,u))
	legend(loc='best',fancybox=True,ncol=2)
	tight_layout()

def compareVar(vt,vn,u,M,L):
	print(vt)
	cm = polyfit(M['t'],M[vn],1)
	pm = poly1d(cm)
	print('Mark1: ',['m','b'],'=',cm,M[vn].std())
	cl = polyfit(L['t'],L[vn],1)
	pl = poly1d(cl)
	print('Lammps:',['m','b'],'=',cl,L[vn].std())
	re = (cl[1]-cm[1])/abs(cl[1])
	print('Rerr: ',re*100.0,'%')
	print('')
	
	figure(figsize=(5,4))
	subplot(111)
	
	plot(M['t'],M[vn],'b-',alpha=0.3,label='Mark1')
	plot(L['t'],L[vn],'r-',alpha=0.3,label='Lammps')
	
	plot(M['t'],pm(M['t']),'b--',label='Mark1 Fit')
	plot(L['t'],pl(L['t']),'r--',label='Lammps Fit')
	
	xlim( min(M['t'].min(),L['t'].min()) , max(M['t'].max(),L['t'].max()) )
	
	xlabel(r'Time $t$ [ps]')
	ylabel(r'%s $%s$ [%s]'%(vt,vn,u))
	legend(loc='best',fancybox=True,ncol=2)
	tight_layout()

M = getData('mark1.thermo',40.0)
#L = getData('lammps.thermo',40.0)

plotVar('Total Energy','TE','eV',M)
plotVar('Kinetic Energy','KE','eV',M)
plotVar('Potential Energy','PE','eV',M)
plotVar('Temperature','T','K',M)
plotVar('Pressure','P','bar',M)

show()