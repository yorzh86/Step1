#!/usr/bin/env python

import pylab as pl

def readLammps(fn):
	fh = open(fn)
	lines = fh.readlines()
	fh.close()
	
	steps = []
	for kl in range(len(lines)):
		if lines[kl].startswith('ITEM: TIMESTEP'):
			step = []
			Na = int(lines[kl+3])
			for ka in range(Na):
				v = map(float,lines[kl+9+ka].split())
				step.append(v)
			steps.append(step)
			kl += 9+Na
	return pl.array(steps)

def readMark1(fn):
	fh = open(fn)
	lines = fh.readlines()
	fh.close()
	
	Na = int(lines[0])
	steps = []
	for kl in range(len(lines)):
		if lines[kl].startswith('Atoms. Timestep:'):
			step = []
			for ka in range(Na):
				v = map(float,lines[kl+1+ka].split()[1:])
				step.append(v)
			steps.append(step)
			kl += 1+Na
	return pl.array(steps)

L = readLammps('lammps.vel')
M = readMark1('mark1.vel')

U0 = pl.sqrt((L**2).sum(2)).mean()

pl.figure()

for kx in range(3):
	pl.subplot(3,1,kx+1)
	for ka in range(L.shape[1]):
		rel_diff = (L[:,ka,kx]-M[:,ka,kx])/U0
		pl.plot( rel_diff,label='lammps-mark1')
	yl = pl.ylim()
	pl.ylim([min(yl[0],0.0),max(yl[1],0.0)])
	pl.title( '%s-velocity'%('xyz'[kx]) )
pl.tight_layout()
pl.show()