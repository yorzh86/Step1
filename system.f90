module system_mod
	use kinds_mod
	implicit none
	
	!==============!
	!= Parameters =!
	!==============!
	
	real(wp),parameter::G0 = 6.67E-11_wp
	real(wp),parameter::kB = 1.3806488E-23
	
	real(wp),parameter::E0 = 1.0_wp
	real(wp),parameter::S0 = 1.0_wp
	
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

	subroutine buildOrbit(T0)
		real(wp),intent(out)::T0
		
		real(wp),parameter::M0 = 1.0E100_wp
		real(wp),parameter::R0 = 2.0_wp
		
		real(wp),dimension(2)::a0
		integer::k
		
		allocate(types(2))
		allocate(atoms(2))
		
		box = 4.0_wp*R0
		
		types%m = [M0,1.0_wp]
		types%atom_name = ['Au','H ']
		
		atoms(1)%r = [2.0_wp,2.0_wp]*R0
		atoms(2)%r = [3.0_wp,2.0_wp]*R0
		
		atoms(1)%t = 1
		atoms(2)%t = 2
		
		a0 = -delV(2)/types(atoms(2)%t)%m
		
		atoms(1)%v = [0.0_wp,0.0_wp]
		atoms(2)%v = [0.0_wp,sqrt(norm2(a0)*R0)]
		
		do k=1,size(atoms)
			atoms(k)%a = -delV(k)/types(atoms(k)%t)%m
		end do
		
		t  = 0.0_wp
		ts = 0
		T0 = 2.0_wp*PI*R0/atoms(2)%v(2)
	end subroutine buildOrbit

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
			o = o+types(atoms(k)%t)%m*norm2(atoms(k)%v)**2
		end do
		o = o/(3.0_wp*real(size(atoms),wp)*kB)
	end function temperature

end module system_mod
