#!/usr/bin/env python

from pylab import *

t,x,y,u,v,ax,ay = loadtxt('trajectory.dat',unpack=True)

r = sqrt(x**2+y**2)
k = u**2+v**2

s = '.' if t.size < 100 else ''

figure('Trajectory',figsize=(5,4))
subplot(111,aspect=1)
plot(x,y,'b%s-'%s,lw=1)
xl,xh = (x.min(),x.max())
xb = 0.1*(xh-xl)
xlim(xl-xb,xh+xb)
yl,yh = (y.min(),y.max())
yb = 0.1*(yh-yl)
ylim(yl-yb,yh+yb)
xlabel(r'$x$-coordinate [m]')
ylabel(r'$y$-coordinate [m]')
tight_layout()

figure('',figsize=(8,8))

subplot(221)
#~ figure('Decay',figsize=(5,4))
plot(t,r,'r%s-'%s)
yl,yh = ylim()
yb = 0.1*(yh-0)
ylim(0-yb,yh+5*yb)
xlabel(r'Time $t$ [s]')
ylabel(r'Radius $r$ [m]')
#~ tight_layout()

subplot(222)
#~ figure('Kinetic Energy',figsize=(5,4))
plot(t,k,'r%s-'%s)
yl,yh = ylim()
yb = 0.1*(yh-0)
ylim(0-yb,yh+5*yb)
xlabel(r'Time $t$ [s]')
ylabel(r'Kinetic Energy $KE$ [J]')
#~ tight_layout()

subplot(223)
#~ figure('Velocities',figsize=(5,4))
plot(t,u,'r%s-'%s,label=r'$\vec{v}\cdot\hat{e}_x$')
plot(t,v,'b%s-'%s,label=r'$\vec{v}\cdot\hat{e}_y$')
yl,yh = ylim()
yb = 0.1*(yh-yl)
ylim(yl-yb,yh+5*yb)
xlabel(r'Time $t$ [s]')
ylabel(r'Velocity $\vec{v}\cdot\hat{e}_n$ [m/s]')
legend(loc='best',fancybox=True,ncol=2)
#~ tight_layout()

subplot(224)
#~ figure('Acceleration',figsize=(5,4))
plot(t,ax,'r%s-'%s,label=r'$\vec{a}\cdot\hat{e}_x$')
plot(t,ay,'b%s-'%s,label=r'$\vec{a}\cdot\hat{e}_y$')
yl,yh = ylim()
yb = 0.1*(yh-yl)
ylim(yl-yb,yh+5*yb)
xlabel(r'Time $t$ [s]')
ylabel(r'Acceleration $\vec{a}\cdot\hat{e}_n$ [m/s$^2$]')
legend(loc='best',fancybox=True,ncol=2)
#~ tight_layout()

tight_layout()

show()
