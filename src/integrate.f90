module integrate_mod
	use kinds_mod
	use system_mod
	implicit none
	private
	
	logical::doThermostat = .false.
	
	public::setThermostat
	public::velocityVerlet
	public::leapFrog
	public::doBox
	
contains

	subroutine setThermostat(state,T,tau)
		!! Turns on/off damping parameter "eta" in the Integrator
		logical,intent(in)::state
		real(wp),intent(in),optional::T,tau
		
		if(state .and. present(T) .and. present(tau) ) then
			thermostat%set = T
			thermostat%tau = tau
			thermostat%eta = 0.0_wp
		else if(state) then
			write(stdout,'(1A)') colorize('Error: Called thermostate(true) without T and tau.',[5,0,0])
			stop 1
		else
			thermostat%eta = 0.0_wp
		end if
		
		doThermostat = state
	end subroutine setThermostat

	subroutine velocityVerlet(dt)
	! Velocity-Verlet integration
		real(wp),intent(in)::dt
		
		real(wp),dimension(3)::d !3D
		integer::k
		
		do k=1,size(atoms)
			d = atoms(k)%v*dt+0.5_wp*atoms(k)%a*dt**2
			atoms(k)%r =  atoms(k)%r+d
		end do
		
		do k=1,size(atoms)
			atoms(k)%v =  atoms(k)%v+atoms(k)%a*0.5_wp*dt
		end do
		
		do k=1,size(atoms)
			atoms(k)%a = -delV(k)/types(atoms(k)%t)%m-thermostat%eta*atoms(k)%v
		end do
		
		do k=1,size(atoms)
			atoms(k)%v =  atoms(k)%v+atoms(k)%a*0.5_wp*dt
		end do
		
		if(doThermostat) thermostat%eta = thermostat%eta+DetaDt()*dt
		t  = t+dt
		ts = ts+1
	end subroutine velocityVerlet

	subroutine leapFrog(dt)
	! Leap frog integration
		real(wp),intent(in)::dt
		
		real(wp),dimension(3)::d,ao !3D
		integer::k
		
		do k=1,size(atoms)
			d = atoms(k)%v*dt+0.5_wp*atoms(k)%a*dt**2
			atoms(k)%r = atoms(k)%r+d
		end do
		do k=1,size(atoms)
			ao = atoms(k)%a
			atoms(k)%a = -delV(k)/types(atoms(k)%t)%m-thermostat%eta*atoms(k)%v
			atoms(k)%v = atoms(k)%v+0.5_wp*(ao+atoms(k)%a)*dt
		end do
		
		if(doThermostat) thermostat%eta = thermostat%eta+DetaDt()*dt
		t  = t+dt
		ts = ts+1
	end subroutine leapFrog

	subroutine doBox
	! returns moving atoms into the simulation box
		integer::k
		
		forall(k=1:3)
			where(atoms(:)%r(k)>box(k)) atoms(:)%r(k) = atoms(:)%r(k)-box(k)
			where(atoms(:)%r(k)<0.0_wp) atoms(:)%r(k) = atoms(:)%r(k)+box(k)
		end forall
	end subroutine doBox

	pure function DetaDt() result(o)
	! calculates damping parameter("eta") change over time.
		real(wp)::o
		
		o = (temperature()/thermostat%set-1.0_wp)/thermostat%tau**2
	end function DetaDt

end module integrate_mod
