module lmpIntegrator_mod
	use kinds_mod
	use units_mod
	use settings_mod
	use system_mod
	use integrate_mod
	use output_mod
	
	implicit none
	private
	
	public::integrateLammps
	!integer,parameter::dp = selected_real_kind(15)
	!integer,parameter::ep = selected_real_kind(18)
	!integer,parameter::wp = dp
	!type::type_t
	!	real(wp)::m
	!		!! Atomic mass
	!	character(2)::atom_name
	!		!! Atom symbol (use standard periodic tables)
	!end type
		
	!type::atom_t
	!	real(wp),dimension(3)::r
	!		!! Atomic position
	!	real(wp),dimension(3)::v
	!		!! Atomic velocity
	!	real(wp),dimension(3)::a
	!		!! Atomic acceleration
	!	real(wp),dimension(3)::f
	!		!! Atomic force
	!	integer::t
	!		!! Atom type
	!	integer,dimension(:),allocatable::neighbors
	!		!! Nearest neighbors
	!end type
	
	!call readLammps(0,5)

contains
	subroutine readLammps(start,finish)
		integer, intent(in)::start, finish
		
		real(wp), dimension(3)::posit, velocity, force
		integer:: i, j, k, atom_id
		character(len=100)::frmt

		open(1, file='Mark1/lammps/lammps.all', status= 'old')
		frmt = '(1I4, 3F5.1, 6F13.9)'

		do i=start, finish
			do j=1, 9
				read(1,*)
			end do
			
			write(*,*) 'Timestep #', j
			do k=1, 10
				read(1, frmt) atom_id, posit, velocity, force
				write(*,frmt) atom_id, posit, velocity, force
			end do
			
			do k=1, 22
				read(1,*)
			end do
			
		end do 
		
		close(1)
		write (*,*) posit
		write (*,*) velocity
		write (*,*) force
	
	end subroutine readLammps

	subroutine integrateLammps(dt)
		real(wp),intent(in)::dt	
		integer::i, j, k, atom_id
		real(wp), dimension(3)::posit, velocity, force
		character(len=100)::frmt
		
		!open(1, file='Mark1/lammps/lammps.all', status= 'old')
		frmt = '(1I4, 3F5.1, 6F13.9)'
		
		do j=1, 9
			read(1,*)
		end do
			
		do k=1, size(atoms)
			read(1, frmt) atoms(k)%atom_id, posit, velocity, force
			atoms(k)%r = posit/1.0E10_wp
			atoms(k)%v = velocity/1.0E-2_wp
			atoms(k)%f = force/6.24150636309E8_wp
			atoms(k)%a = -delV(k)/types(atoms(k)%t)%m
		end do
		
		t  = t+dt
		ts = ts+1
		
		!close(1)
		
		!atoms(k)%r =  atoms(k)%r+atoms(k)%v*dt+0.5_wp*atoms(k)%a*dt**2
		!atoms(k)%v =  atoms(k)%v+atoms(k)%a*0.5_wp*dt
		!atoms(k)%a = -delV(k)/types(atoms(k)%t)%m-(Teta+Pepsilon)*atoms(k)%v
		!atoms(k)%f = -delV(k)
		!atoms(k)%v =  atoms(k)%v+atoms(k)%a*0.5_wp*dt
	end subroutine integrateLammps
	
end module lmpIntegrator_mod
