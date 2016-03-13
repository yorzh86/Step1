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

contains
	subroutine integrateLammps(dt)
		real(wp),intent(in)::dt	
		integer::i, j, k
		real(wp), dimension(3)::posit, velocity, force
		character(len=100)::abc

		do j=1, 9
			read(1,*)
		end do
			
		do k=1, size(atoms)
			read(1, '(1I4, 3F5.1, 6F13.9)') atoms(k)%atom_id, posit, velocity, force
			atoms(k)%r = posit/1.0E10_wp
			atoms(k)%v = velocity/1.0E-2_wp
			atoms(k)%f = force/6.24150636309E8_wp
			atoms(k)%a = -delV(k)/types(atoms(k)%t)%m
		end do
		
		t  = t+dt
		ts = ts+1
		
	end subroutine integrateLammps
	
end module lmpIntegrator_mod
