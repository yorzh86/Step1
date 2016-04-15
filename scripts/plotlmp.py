#!/usr/bin/env python

#As described below, the total kinetic energy transferred by these swaps is 
#computed by the fix and can be output. Dividing this quantity by time and the 
#cross-sectional area of the simulation box yields a heat flux.
#The ratio of heat flux to the slope of the temperature profile is proportional
#to the thermal conductivity of the fluid, in appropriate units.

import pylab as pl
from scipy.signal import savgol_filter

def readT(fn):
    f = open(fn)
    lines = f.readlines()
    f.close()
    hs = []
    for k in range(len(lines)):
        if lines[k].startswith('#'):
            continue
        if not lines[k].startswith(' '):
            hs.append(k)
    n_slabs = hs[1]-hs[0]-1
    
    S = []
    for h in hs:
        L = []
        for k in range(n_slabs):
            L.append(map(float,lines[h+k+1].split())[1:])
        S.append(L)
    S = pl.array(S)
    #x = S[0,:,0]
    # Nelement:rows:column]
    T = []
    for k in range(S.shape[0]):
        T.append(list(S[k,:,2]))
    T = pl.array(T)
    
    dz = 150*5.4/2
    grad = ((T.mean(0)[5]-T.mean(0)[0]) + (T.mean(0)[5] - T.mean(0)[9]))/2.0/dz
    
    return T, grad

def readKE(fn):
    time,ke = pl.loadtxt(fn, unpack=True, skiprows=1)
    flux = pl.mean(ke/time)
    
    return time, ke, flux

def plotFlux(ax, t, ke):
    
    sfp = 31 #(FIXME)
    skip = 5
    #count from (FIXME):
    i = len(t)*0.428
    
    ax.axvline(3, color='k', linestyle='--')
    ax.plot(t[::skip]/1000, (ke[::skip]/t[::skip]),'m.', alpha = 0.6)
    ax.plot(t[::skip]/1000, savgol_filter(ke[::skip]/t[::skip],sfp,3),
        color='blue', linewidth=1.0)
    ax.set_xlabel('Time $t$ [ns]')
    ax.set_ylabel('Flux $q$ [eV/ps]')
    
    q1 = ke[i:]/t[i:]
    qm = q1.mean()
    return qm #eV/ps

def plotHC(ax, t, T):
    sfp = 111 #(FIXME)
    skip = 2
    #count from (FIXME):
    c1 = len(t)*0.428
    H = T[:,5]
    C = T[:,0]
    ax.axvline(3, color='k', linestyle='--')
    ax.plot(t[::skip]/1000, H[::skip],'r.', label='Hot',  alpha = 0.2)
    ax.plot(t[::skip]/1000, C[::skip],'c.', label='Cold', alpha = 0.2)
    ax.plot(t[::skip]/1000,savgol_filter(H[::skip],sfp,3), color='darkred',
        linewidth=1.5)
    ax.plot(t[::skip]/1000,savgol_filter(C[::skip],sfp,3), color='darkblue',
        linewidth=1.5)
    ax.set_xlabel('Time $t$ [ns]')
    ax.set_ylabel('Temperature $T$ [K]')
    ax.legend(loc='best',fancybox=True)

def plotContour(ax, t, T):
    skip = 2
    z = pl.linspace(0,150*5.4,T.shape[1])
    #pl.pcolormesh(z/10,t,T)
    pl.contourf(z/10,t[::skip]/1000,T[::skip,:],19)
    pl.colorbar()
    ax.set_xlim(z.min()/10,z.max()/10)
    ax.set_ylim(t[::skip].min()/1000,t[::skip].max()/1000)
    ax.set_xlabel('Position $z$ [nm]')
    ax.set_ylabel('Time $t$ [ns]')

def plotGradient(ax, Ti):
    
    i = len(Ti)*0.428
    T = Ti[i:,:]
    
    z = pl.linspace(0,150*5.4,T.shape[1])
    
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


def doFigure(fnKE, fnT):
    temps,grad = readT(fnKE)
    time,ke,flux = readKE(fnT)
    
    fig = pl.figure()    
    
    axC = fig.add_subplot(221)
    plotContour(axC,time,temps)

    axHC = fig.add_subplot(223)
    plotHC(axHC,time,temps)

    axGrad = fig.add_subplot(222)
    gr = plotGradient(axGrad, temps)
    
    axQ = fig.add_subplot(224)
    heatrate = plotFlux(axQ,time,ke)

    fig.tight_layout()
    pl.subplots_adjust(top=0.91, hspace=0.46)
    fig.suptitle('Reverse non equilibrium method', fontsize=14, fontweight='bold') 
    
    A = (5.4*5)*(5.4*5)
    k2 = heatrate/A/gr
    
    return flux, grad, A, k2

flux, grad, A, k2 = doFigure('Tregions.dat', 'MP.dat')

J2eV = 1.60217646E-19
s2ps = 1.0E-12
m2A  = 1.0E-10

k = flux/A/grad
kSI  = k*J2eV/(s2ps*m2A)
k2SI = k2*J2eV/(s2ps*m2A)
pl.title( 'k=%f [W/m.K]; k=%f [W/m.K]'%(kSI,k2SI))
print 'k=', kSI,'[W/m.K]'
print 'k2=',k2SI,'[W/m.K]'

pl.show()
