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
	real(wp),parameter::cutoff = 2.5_wp*S0
	real(wp),parameter::neighborRadius = cutoff*1.3_wp
	
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
		real(wp),dimension(2)::r
			!! Atomic position
		real(wp),dimension(2)::v
			!! Atomic velocity
		real(wp),dimension(2)::a
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
	
	real(wp),dimension(2)::box
		!! Bounds of the simulation box
	
	integer::ts
		!! Time step counter
	real(wp)::t
		!! Sytem time
	
contains

	subroutine buildSystem(a,N,Ti)
		real(wp),intent(in)::a
			!! Lattice constant
		integer,intent(in)::N
			!! Number of unit cells
		real(wp),intent(in)::Ti
			!! Initial temperature
		
		integer::i,j,k
		
		box = a*real([N,N],wp)!*[1.0_wp/sqrt(2.0_wp),1.0_wp/sqrt(2.0_wp)]
		!! if we change box size (increase) - we lose Potential energy and atoms behave as liquid
		
		allocate(types(1))
		allocate(atoms(N**2))
		
		!types%m = 39.948_wp
		types%m = 6.6335209E-26_wp
		types%atom_name = 'Ar'
		atoms(:)%t = 1
		
		do i=1,N
			do j=1,N
				!atoms(i+N*(j-1))%r = a*(real([i,j],wp)*[1.0_wp,1.0_wp/sqrt(2.0_wp)]+[mod(j,2)/2.0_wp,0.0_wp]) &
				!& -a/2.0_wp+a/2.0_wp*N
				atoms(i+N*(j-1))%r = a*real([i,j],wp)
			end do
		end do
		
		do k=1,N**2
			!! Create random direct
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
		
		do i=1,size(atoms)
			call updateNeighbors(i)
		end do
		
		do i=1,N
			do j=1,N
				atoms(i+N*(j-1))%a = -delV(i+N*(j-1))/types(atoms(i+N*(j-1))%t)%m
			end do
		end do
		
		forall(k=1:2) atoms(:)%v(k) = atoms(:)%v(k)-sum(atoms(:)%v(k))/real(N*N,wp)
		
		ts = 0
		t  = 0.0_wp
	end subroutine buildSystem

	pure function V(i) result(o)
		integer,intent(in)::i
		real(wp)::o
		
		o = 0.0_wp
		if(enableLennardJones) call doLennardJones(o)
		
	contains
	
		pure subroutine doLennardJones(o)
			real(wp),intent(inout)::o
			
			real(wp),dimension(2)::d
			real(wp)::l
			integer::k
			
			do k=1,size(atoms)
				if(k==i) cycle
				d = deltaR(atoms(k),atoms(i))
				l = S0/norm2(d)
				o = o+4.0_wp*E0*(l**12-l**6)
			end do
		end subroutine doLennardJones
	
	end function V

	pure function delV(i) result(o)
		integer,intent(in)::i
		real(wp),dimension(2)::o
		real(wp),dimension(2)::r
		
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
		real(wp),dimension(2),intent(in)::d
		real(wp),dimension(2)::o
		
		o = 0.0_wp
		if(i==j) return
		
		if(enableLennardJones) call doLennardJones(o)
		
	contains
		
		pure subroutine doLennardJones(o)
			real(wp),dimension(2),intent(inout)::o
			real(wp)::l
			
			l = S0/norm2(d)
			o = o+24.0_wp*E0/sum(d*d)*(l**6-2.0_wp*l**12)*d
		end subroutine doLennardJones
		
	end function delVij

	pure function deltaR(a1,a2) result(o)
		type(atom_t),intent(in)::a1,a2
		real(wp),dimension(2)::o
		
		real(wp),dimension(2)::d
		d = a1%r-a2%r
		o = d-box*nint(d/box)
		
	end function deltaR

	pure function temperature() result(o)
		real(wp)::o
		
		real(wp),dimension(2)::vBulk
		integer::k
		
		forall(k=1:2) vBulk(k) = sum(atoms%v(k))/real(size(atoms),wp)
		
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
		real(wp),dimension(2),intent(in),optional::vBulk
		real(wp)::o
		
		real(wp),dimension(2)::v0
		
		v0 = 0.0_wp
		if(present(vBulk)) v0 = vBulk
		o = 0.5_wp*types(atoms(i)%t)%m*norm2(atoms(i)%v-v0)**2
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
		real(wp),dimension(2)::o
		
		real(wp),dimension(2)::fij,vj,rij
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
		
		real(wp),dimension(2)::r
		integer::k
		
		atoms(i)%neighbors = [integer::]
		do k=1,size(atoms)
			if(i==k) cycle
			r  = deltaR(atoms(i),atoms(k))
			if(norm2(r)>neighborRadius) cycle
			atoms(i)%neighbors = [atoms(i)%neighbors,k]
		end do
	end subroutine updateNeighbors

end module system_mod
