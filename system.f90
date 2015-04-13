module system_mod
	use kinds_mod
	implicit none
	
	!==============!
	!= Parameters =!
	!==============!
	

	real(wp),parameter::kB = 8.617332478E-5_wp
	
	real(wp),parameter::E0 = 1.0298490416E-2_wp
	real(wp),parameter::S0 = 3.4_wp 
	
	!=========!
	!= Types =!
	!=========!
	
	type::type_t
		real(wp)::m
		character(2)::atom_name
	end type
	
	type::atom_t
		real(wp),dimension(2)::r
		real(wp),dimension(2)::v
		real(wp),dimension(2)::a
		integer::t
	end type
	
	!===============!
	!= Module Data =!
	!===============!
	
	logical::enableLennardJones = .false.
	
	type(type_t),dimension(:),allocatable::types
	type(atom_t),dimension(:),allocatable::atoms
	
	real(wp),dimension(2)::box
	
	integer::ts
	real(wp)::t
	
contains

	subroutine buildSystem(a,N,Ti)
		real(wp),intent(in)::a
		integer,intent(in)::N
		real(wp),intent(in)::Ti
		
		integer::i,j,k
		
		box = 2.0_wp*a*real([N,N],wp)
		
		allocate(types(1))
		allocate(atoms(N**2))
		
		types%m = 39.9_wp
		types%atom_name = 'Ar'
		
		atoms(:)%t = 1
		forall(i=1:N,j=1:N)
			atoms(i+N*(j-1))%r = a*real([i,j],wp)-a/2.0_wp+a/2.0_wp*N
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
				o = o
			end do
		end subroutine doLennardJones
	
	end function V

	pure function delV(i) result(o)
		integer,intent(in)::i
		real(wp),dimension(2)::o
		
		o = 0.0_wp
		if(enableLennardJones) call doLennardJones(o)
		
	contains

		pure subroutine doLennardJones(o)
			real(wp),dimension(2),intent(inout)::o
			
			real(wp),dimension(2)::d
			real(wp)::l
			integer::k
			
			do k=1,size(atoms)
				if(k==i) cycle
				d = deltaR(atoms(k),atoms(i))
				l = S0/norm2(d)
				o = o+24.0_wp*E0/sum(d*d)*(2.0_wp*l**12-l**6)*d
			end do
		end subroutine doLennardJones
	
	end function delV

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
		
		o = 0.0_wp
		do k=1,size(atoms)
			o = o+0.5_wp*types(atoms(k)%t)%m*norm2(atoms(k)%v)**2
		end do
		o = o/(2.0_wp*real(size(atoms),wp)*kB)
	end function temperature

	pure function kineticE() result (o)
		real(wp)::o
		integer::k
		o = 0.0_wp
		do k=1,size(atoms)
			o = o + 0.5_wp*types(atoms(k)%t)%m*norm2(atoms(k)%v)**2
		end do
	end function kineticE
	
	pure function potentialE(i) result (o)
		integer,intent(in)::i
		real(wp)::o
		real(wp),dimension(2)::d
		real(wp)::l
		integer::k
		
		o = 0.0_wp
		do k=1,size(atoms)
			o = o+V(k)
		end do
	end function potentialE

end module system_mod
