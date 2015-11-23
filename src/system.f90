module system_mod
	use kinds_mod
	implicit none
	
	!==============!
	!= Parameters =!
	!==============!
	
	real(wp),parameter::kB = 1.3806488E-23_wp
		!! Boltzmann constant in SI units
	!real(wp),parameter::E0 = 1.04233e-2
	real(wp),parameter::E0 = kB*125.7_wp
		!! Lennard Jones epsilon in SI units
	real(wp),parameter::S0 = 3.345E-10_wp
		!! Lennard Jones sigma in SI units
	real(wp),parameter::cutoff = 2.0_wp*S0
	real(wp),parameter::neighborRadius = cutoff + 1.0E-10
	
	!=========!
	!= Types =!
	!=========!
	
	type::type_t
		real(wp)::m
			!! Atomic mass
		character(2)::atom_name
			!! Atom symbol (use standard periodic tables)
	end type
	
	type::atom_t
		real(wp),dimension(3)::r
			!! Atomic position
		real(wp),dimension(3)::v
			!! Atomic velocity
		real(wp),dimension(3)::a
			!! Atomic acceleration
		integer::t
			!! Atom type
		integer,dimension(:),allocatable::neighbors
			!! Nearest neighbors
	end type
	
	!===============!
	!= Module Data =!
	!===============!
	
	logical::enableLennardJones = .false.
		!! Run Lennard-Jones potential
	
	type(type_t),dimension(:),allocatable::types
		!! Types for all atoms in system
	type(atom_t),dimension(:),allocatable::atoms
		!! All atoms in system
	
	real(wp),dimension(3)::box !3D
		!! Bounds of the simulation box
	
	integer::ts
		!! Time step counter
	real(wp)::t
		!! System time
	
contains

	subroutine buildSystem(a,N,Ti)
		implicit none
		real(wp),intent(in)::a
			!! Lattice constant
		integer,dimension(3),intent(in)::N
			!! Number of unit cells
		real(wp),intent(in)::Ti
			!! Initial temperature
		integer::i,j,k,l,idx
		real(wp), dimension(3,4), parameter::rcell= &
		reshape([0.0_wp, 0.0_wp, 0.0_wp, &
		         0.5_wp, 0.5_wp, 0.0_wp, &
		         0.0_wp, 0.5_wp, 0.5_wp, &
		         0.5_wp, 0.0_wp, 0.5_wp  ], [3,4])
		
		box = a*real(N,wp)
		
		allocate(types(1))
		allocate(atoms(size(rcell,2)*product(N)))
		
		types%m = 39.948_wp*1.6605402E-27
		types%atom_name = 'Ar'
		atoms(:)%t = 1
		
		idx = 1
		do k=1,N(3)
			do j=1,N(2)
				do i=1,N(1)
					do l=1,size(rcell,2)
						atoms(idx)%r = a*(real([i,j,k]-1,wp)+rcell(1:3,l))
						idx = idx+1
					end do
				end do
			end do
		end do
		
		do k=1,size(atoms)
			!! Create random direction of velocities
			call random_number(atoms(k)%v)
			atoms(k)%v = 2.0_wp*atoms(k)%v-1.0_wp
			do while(norm2(atoms(k)%v)>1.0_wp .and. norm2(atoms(k)%v)<0.1_wp)
				call random_number(atoms(k)%v)
				atoms(k)%v = 2.0_wp*atoms(k)%v-1.0_wp
			end do
			atoms(k)%v = atoms(k)%v/norm2(atoms(k)%v)
			!! Set velocity magnitude
			atoms(k)%v = atoms(k)%v*sqrt(2.0_wp*kB*Ti/types(atoms(k)%t)%m)*abs(randomNormal()+1.0_wp)
		end do
		forall(k=1:3) atoms(:)%v(k) = atoms(:)%v(k)-sum(atoms(:)%v(k))/real(size(atoms),wp)
		
		do i=1,size(atoms)
			call updateNeighbors(i)
		end do
		
		do k=1,size(atoms)
			atoms(k)%a = -delV(k)/types(atoms(k)%t)%m
		end do
		
		ts = 0
		t  = 0.0_wp
	end subroutine buildSystem

	subroutine writeLammpsData(fn)
		character(*),intent(in)::fn
		
		integer::N,k,iou
		
		open(file=fn,newunit=iou)
		
		write(iou,'(1A)') '# Input geometry for lammps'
		write(iou,'(1I9,1X,1A)') size(atoms),'atoms'
		write(iou,'(1I9,1X,1A)') size(types),'atom types'
		write(iou,'(2F13.6,1X,1A,1X,1A)') 0.0_wp,box(1),'xlo','xhi'
		write(iou,'(2F13.6,1X,1A,1X,1A)') 0.0_wp,box(2),'ylo','yhi'
		write(iou,'(2F13.6,1X,1A,1X,1A)') -0.5_wp,0.5_wp,'zlo','zhi'
		write(iou,'(1A)') ''
		write(iou,'(1A)') 'Masses'
		write(iou,'(1A)') ''
		do k=1,size(types)
			write(iou,'(1I4,1X,1F13.6)') k,types(k)%m
		end do
		write(iou,'(1A)') ''
		write(iou,'(1A)') 'Atoms'
		write(iou,'(1A)') ''
		do k=1,size(atoms)
			write(iou,'(1I9,1X,1I2,1X,3F13.6)') k,atoms(k)%t,atoms(k)%r,0.0_wp
		end do
		
		close(iou)
	end subroutine writeLammpsData

	pure function V(i) result(o)
		integer,intent(in)::i
		real(wp)::o
		
		o = 0.0_wp
		if(enableLennardJones) call doLennardJones(o)
		
	contains
	
		pure subroutine doLennardJones(o)
			!real(wp),intent(inout)::o
			real(wp),intent(out)::o
			real(wp),dimension(3)::d !3D
			real(wp)::l
			integer::k
			
			do k=1,size(atoms)
				if(k==i) cycle
				d = deltaR(atoms(k),atoms(i))
				l = S0/norm2(d)
				o = o+4.0_wp*E0*(l**12.0_wp-l**6.0_wp)
			end do
		end subroutine doLennardJones
	
	end function V

	pure function delV(i) result(o)
		integer,intent(in)::i
		real(wp),dimension(3)::o !3D
		real(wp),dimension(3)::r !3D
		
		integer::j,aj
		
		o = 0.0_wp
		do j=1,size(atoms(i)%neighbors)
			aj = atoms(i)%neighbors(j)
			r  = deltaR(atoms(i),atoms(aj))
			if(norm2(r)>cutoff) cycle
			o = o+delVij(i,aj,r)
		end do
	end function delV

	pure function delVij(i,j,d) result (o)
		integer,intent(in)::i,j
		real(wp),dimension(3),intent(in)::d !3D
		real(wp),dimension(3)::o			!3D
		
		o = 0.0_wp
		if(i==j) return
		
		if(enableLennardJones) call doLennardJones(o)
		
	contains
		
		pure subroutine doLennardJones(o)
			real(wp),dimension(3),intent(out)::o !3D
			real(wp)::l
			
			l = S0/norm2(d)
			!o = o+24.0_wp*E0/sum(d*d)*(l**6-2.0_wp*l**12)*d
			o = o-(48.0_wp*E0/sum(d*d))*(l**12.0_wp-0.5_wp*l**6.0_wp)*d
		end subroutine doLennardJones
		
	end function delVij

	pure function deltaR(a1,a2) result(o)
		type(atom_t),intent(in)::a1,a2
		real(wp),dimension(3)::o !3D
		
		real(wp),dimension(3)::d !3D
		d = a1%r-a2%r
		o = d-box*anint(d/box)
		
	end function deltaR

	pure function temperature() result(o)
		real(wp)::o
		
		real(wp),dimension(3)::vBulk !3D
		integer::k
		
		forall(k=1:3) vBulk(k) = sum(atoms%v(k))/real(size(atoms),wp) !3D forall (1:3)
		
		o = 0.0_wp
		do k=1,size(atoms)
			o = o+KEi(k,vBulk)
		end do
		o = o/(2.0_wp*real(size(atoms),wp)*kB)
	end function temperature

	pure function E() result(o)
	!! we never use it
		real(wp)::o
		integer::k
		
		o = 0.0_wp
		do k=1,size(atoms)
			o = o+Ei(k)
		end do
	end function E

	pure function Ei(i) result(o)
		integer,intent(in)::i
		real(wp)::o
		
		o = KEi(i) + PEi(i)
	end function Ei

	pure function KE() result (o) 
		real(wp)::o
		integer::k
		
		o = 0.0_wp
		do k=1,size(atoms)
			o = o+KEi(k)
		end do
	end function KE
	
	pure function KEi(i,vBulk) result (o)
		integer,intent(in)::i
		real(wp),dimension(3),intent(in),optional::vBulk !3D
		real(wp)::o
		
		real(wp),dimension(3)::v0 !3D
		
		v0 = 0.0_wp
		if(present(vBulk)) v0 = vBulk
		o = 0.5_wp*types(atoms(i)%t)%m*norm2(atoms(i)%v-v0)**2.0_wp
	end function KEi
	
	pure function PE() result (o)
		real(wp)::o
		integer::k
		
		o = 0.0_wp
		do k=1,size(atoms)
			o = o+PEi(k)
		end do
	end function PE

	pure function PEi(i) result (o)
		integer,intent(in)::i
		real(wp)::o
		
		o = V(i)
	end function PEi

	pure function heatflux() result(o)
		real(wp),dimension(3)::o !3D
		
		real(wp),dimension(3)::fij,vj,rij !3D
		integer::i,j
		
		o = 0.0_wp
		
		do i=1,size(atoms)
			o = o+Ei(i)*atoms(i)%v
		end do
		
		do j=1,size(atoms)
			do i=1,j-1
				rij = deltaR(atoms(i),atoms(j))
				fij = delVij(i,j,rij)
				vj  = atoms(j)%v
				o = o+dot_product(fij,vj)*rij
			end do
		end do
	end function heatflux

	subroutine updateNeighbors(i)
		integer,intent(in)::i
		
		real(wp),dimension(3)::r !3D
		integer::k
		
		atoms(i)%neighbors = [integer::]
		do k=1,size(atoms)
			if(i==k) cycle
			r  = deltaR(atoms(i),atoms(k))
			if(norm2(r)>neighborRadius) cycle
			atoms(i)%neighbors = [atoms(i)%neighbors,k]
		end do
	end subroutine updateNeighbors

	function averageNeighbors() result(o)
		real(wp)::o
		integer::k
		
		o = sum([( size(atoms(k)%neighbors) , k=1,size(atoms) )]) / real(size(atoms),wp)
	end function averageNeighbors

end module system_mod
