load mark1.xyz
load lammps.xyz

hide everything
show spheres

run ./scripts/drawBox.py
drawBoundingBox 

color red, lammps
color cyan, mark1
