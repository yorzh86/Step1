module system_mod
	use kinds_mod
	implicit none
	
	!==============!
	!= Parameters =!
	!==============!
	
	real(wp),parameter::kB = 8.617332478E-5_wp
		!! Boltzmann constant in metal units
	real(wp),parameter::E0 = 119.8*kB ! = 1.0298490416E-2_wp
		!! Lennard Jones epsilon
	real(wp),parameter::S0 = 3.4_wp 
		!! Lennard Jones sigma
	
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
		
		box = a*real([N,N],wp)*[1.0_wp,1.0_wp/sqrt(2.0_wp)]
		
		allocate(types(1))
		allocate(atoms(N**2))
		
		types%m = 39.9_wp
		types%atom_name = 'Ar'
		atoms(:)%t = 1
		
		forall(i=1:N,j=1:N)
			atoms(i+N*(j-1))%r = a*(real([i,j],wp)*[1.0_wp,1.0_wp/sqrt(2.0_wp)]+[mod(j,2)/2.0_wp,0.0_wp]) &
			& -a/2.0_wp+a/2.0_wp*N
			atoms(i+N*(j-1))%a = -delV(i+N*(j-1))/types(atoms(i+N*(j-1))%t)%m
		end forall
		
		do k=1,N**2
			call random_number(atoms(k)%v)
			atoms(k)%v = 2.0_wp*atoms(k)%v-1.0_wp
			do while(norm2(atoms(k)%v)>1.0_wp .and. norm2(atoms(k)%v)<0.1_wp)
				call random_number(atoms(k)%v)
				atoms(k)%v = 2.0_wp*atoms(k)%v-1.0_wp
			end do
			atoms(k)%v = atoms(k)%v/norm2(atoms(k)%v)
			atoms(k)%v = atoms(k)%v*sqrt(2.0_wp*kB*Ti/types(atoms(k)%t)%m)*abs(randomNormal()+1.0_wp)
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
		integer::j
		
		o = 0.0_wp
		do j=1,size(atoms)
			o = o+delVij(i,j)
		end do
		
	end function delV

	pure function delVij(i,j) result (o)
		integer,intent(in)::i,j
		real(wp),dimension(2)::o
		
		o = 0.0_wp
		if(i==j) return
		
		if(enableLennardJones) call doLennardJones(o)
		
	contains
		
		pure subroutine doLennardJones(o)
			real(wp),dimension(2),intent(inout)::o
			
			real(wp),dimension(2)::d
			real(wp)::l
			
			d = deltaR(atoms(i),atoms(j))
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
		integer::k
		
		! TODO: Remove bulk velocity
		o = 0.0_wp
		do k=1,size(atoms)
			o = o+KEi(k)
		end do
		o = o/(2.0_wp*real(size(atoms),wp)*kB)
	end function temperature

	pure function E() result(o)
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
		
		o = KEi(i)+PEi(i)
	end function Ei

	pure function KE() result (o) 
		real(wp)::o
		integer::k
		o = 0.0_wp
		do k=1,size(atoms)
			o = o+KEi(k)
		end do
	end function KE
	
	pure function KEi(i) result (o)
		integer,intent(in)::i
		real(wp)::o
		
		o = 0.5_wp*types(atoms(i)%t)%m*norm2(atoms(i)%v)**2
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
				fij = delVij(i,j)
				vj  = atoms(j)%v
				rij = deltaR(atoms(i),atoms(j))
				o = o+dot_product(fij,vj)*rij
			end do
		end do
	end function heatflux

end module system_mod
