module integrate_mod
	use kinds_mod
	use system_mod
	use autodiff_mod
	implicit none
	private
	
	logical::doThermostat = .false.
	logical::doBarostat = .false.
	integer::hot, cold
	
	public::setThermostat
	public::setBarostat
	public::hot, cold
   
	public::velocityVerlet
	public::doBox
	public::swapAtoms
	public::rnem
	
contains

	subroutine setThermostat(state,T,tau)
		!! Turns on/off damping parameter "eta" in the Integrator
		logical,intent(in)::state
		!real(wp),intent(in),optional::T,tau
		type(ad_t), intent(in), optional::T, tau
		
		if(state .and. present(T) .and. present(tau) ) then
			thermostat%set = T
			thermostat%tau = tau
			Teta = 0.0_wp
		else if(state) then
			write(stdout,'(1A)') colorize('Error: Called thermostate(true) without T and tau.',[5,0,0])
			stop 1
		else
			Teta = 0.0_wp
		end if
		
		doThermostat = state
	end subroutine setThermostat
	
	subroutine setBarostat(state,P,tau)
		!! Turns on/off damping parameter "eta" in the Integrator
		logical,intent(in)::state
		!real(wp),intent(in),optional::P,tau
		type(ad_t),intent(in),optional::P,tau
		
		if(state .and. present(P) .and. present(tau) ) then
			barostat%set = P
			barostat%tau = tau
			Pepsilon = 0.0_wp
		else if(state) then
			write(stdout,'(1A)') colorize('Error: Called barostate(true) without P and tau.',[5,0,0])
			stop 1
		else
			Pepsilon = 0.0_wp
		end if
		
		doBarostat = state
	end subroutine setBarostat

	subroutine velocityVerlet(dt)
		!! Velocity-Verlet integration
		!real(wp),intent(in)::dt
		type(ad_t), intent(in)::dt
		
		!real(wp),dimension(3)::d
		type(ad_t), dimension(3)::d
		integer::k
		
		do k=1,size(atoms)
			d = atoms(k)%v*dt+0.5_wp*atoms(k)%a*dt**2
			atoms(k)%r =  atoms(k)%r+d+Pepsilon*atoms(k)%r
		end do
		
		do k=1,size(atoms)
			atoms(k)%v =  atoms(k)%v+atoms(k)%a*0.5_wp*dt
		end do
		
		do k=1,size(atoms)
			atoms(k)%a = -delV(k)/types(atoms(k)%t)%m-(Teta+Pepsilon)*atoms(k)%v
			atoms(k)%f = -delV(k)
		end do
		
		do k=1,size(atoms)
			atoms(k)%v =  atoms(k)%v+atoms(k)%a*0.5_wp*dt
		end do
		
		if(doThermostat) Teta = Teta+DetaDt()*dt
		if(doBarostat) Pepsilon = Pepsilon+DepsilonDt()*dt
		box = box+Pepsilon*box
		t  = t+dt
		ts = ts+1
	end subroutine velocityVerlet

	subroutine doBox
		!! Returns moving atoms into the simulation box
		integer::k
		
		forall(k=1:3)
			where(real(atoms(:)%r(k))>real(box(k))) atoms(:)%r(k) = atoms(:)%r(k)-box(k)
			where(real(atoms(:)%r(k))<0.0_wp) atoms(:)%r(k) = atoms(:)%r(k)+box(k)
		end forall
	end subroutine doBox

	function DetaDt() result(o)
		!! Calculates damping parameter("eta") change over time
		!real(wp)::o
		type(ad_t)::o
		
		o = (1.0_wp/thermostat%tau**2)*(temperature()/thermostat%set-1.0_wp)
	end function DetaDt

	function DepsilonDt() result(o)
		!! Calculates damping parameter("eta") change over time.
		!real(wp)::o
		type(ad_t)::o
		o = (1.0_wp/barostat%tau**2)*(pressure()/barostat%set-1.0_wp)
	end function DepsilonDt
	
	subroutine rnem(k)
		integer, intent(in)::k
		
		!real(wp),dimension(N_slabs)::temperatures
		type(ad_t),dimension(N_slabs)::temperatures
		integer::j
		real(wp)::abc
		integer,dimension(:),allocatable::l
		
		abc = real(latM(3),wp)*lattice_const/real(N_slabs,wp)
		do j=1, N_slabs
			!l = regionList(j*abc - abc, j*abc)
			listofRegions = regionList(j*abc - abc, j*abc)
			temperatures(j) = listTemp(listofRegions)
			if (j==1) hot = selectHot(listofRegions)
			if (j==N_slabs/2+1) cold = selectCold(listofRegions)
		end do
		
		
		do j=1, N_slabs
			regions(j)%temps(k)=temperatures(j)
		end do
!		if(allocated(regions(k)%temps)) deallocate(regions(k)%temps)
!		allocate(regions(k)%temps(2+N_slabs))
!		regions(k)%temps(1)  = convert(real(k,wp)*dt,'s','ps')
!		regions(k)%temps(2)  = real(k,wp)
!		regions(k)%temps(3:) = temperatures
	
	end subroutine rnem

	
	subroutine swapAtoms(k)
	  !! Swapping atoms' velocities
	  !! to swap masses (add when needed):
	  !! - either change their types, or swap atoms of same type
		real(wp)::t
		integer,  intent(in)::k
		type(ad_t), dimension(3)::swapv
		t = convert((real(k,wp)*dt),'s','ps')
		
		regions(1)%energies(k/skip_swap) = KEi(cold)
		regions(2)%energies(k/skip_swap) = KEi(hot)
		
!		if(allocated(regions(k)%energies)) deallocate(regions(k)%energies)
!		allocate(regions(k)%energies(4))
!		regions(k)%energies(1) = t
!		regions(k)%energies(2) = real(k,wp)
!		regions(k)%energies(3) = KEi(hot)
!		regions(k)%energies(4) = KEi(cold)

		swapv = atoms(hot)%v
		atoms(hot)%v = atoms(cold)%v
		atoms(cold)%v = swapv
	end subroutine swapAtoms
	
end module integrate_mod
