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
		type(ad_t), intent(in)::dt
		type(ad_t), dimension(3)::d
		integer::k
		
		do k=1,size(atoms)
			d = atoms(k)%v*dt+0.5_wp*atoms(k)%a*dt**2
			atoms(k)%r =  atoms(k)%r+d+Pepsilon*atoms(k)%r*dt
			!! d = A/ps * ps + A/ps^2 * ps^2 = A + A = [A] distance
			!! atoms%r = [A] + [A] + eV / (ps * A^3 * bar)*[A] = ????  WEIRD!!!
		end do
		
		do k=1,size(atoms)
			atoms(k)%v =  atoms(k)%v+atoms(k)%a*0.5_wp*dt
			!! v = self. + A/ps^2 * ps = [A]/[ps] - correct.
		end do
		
		do k=1,size(atoms)
			atoms(k)%a = -delV(k)/types(atoms(k)%t)%m-(Teta+Pepsilon)*atoms(k)%v
			!! a = eV/A / gram/mole - A/ps^2  !first term has to be coefficient
		end do
		
		do k=1,size(atoms)
			atoms(k)%v =  atoms(k)%v+atoms(k)%a*0.5_wp*dt
			!! v = self. + A/ps^2 * ps = [A]/[ps] - correct.
		end do
		
		if(doThermostat) Teta = Teta+DetaDt()*dt
		!! teta = self. + detaDt * [ps] = 1/ps
		if(doBarostat) Pepsilon = Pepsilon+DepsilonDt()*dt
		!! pepsilon = self. + DepsilonDt * [ps]
		
		
		!! Pepsilon = eV / (ps^2 * A^3 * bar)*[ps] = eV / (ps * A^3 * bar) = 1/ps
		!! DepsilonDT = eV / (ps^2 * A^3 * bar)
		!! detaDT  = [1/ps^2]
		!! Teta = 1/ps
		
		box = box+Pepsilon*box*dt
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
		type(ad_t)::o
		
		o = (1.0_wp/thermostat%tau**2)*(temperature()/thermostat%set-1.0_wp)
		!! o = 1/ ps^2  *  K/ (K - int) = 1/ps^2
		!! DetaDT = [1/ps^2]
	end function DetaDt

	function DepsilonDt() result(o)
		!! Calculates damping parameter("eta") change over time.
		type(ad_t)::o
		o = (1.0_wp/barostat%tau**2)*(pressure()/barostat%set-1.0_wp)
		!! o = 1 / ps^2 * eV/A^3 / (bar -1) = eV / (ps^2 * A^3 * bar) !!! Weird
		
	end function DepsilonDt
	
	subroutine rnem(k)
		integer, intent(in)::k
		type(ad_t),dimension(N_slabs)::temperatures
		integer::j
		real(wp)::abc
		integer,dimension(:),allocatable::l
		
		abc = real(latM(3),wp)*lattice_const/real(N_slabs,wp)
		do j=1, N_slabs
			listofRegions = regionList(j*abc - abc, j*abc)
			temperatures(j) = listTemp(listofRegions)
			if (j==1) hot = selectHot(listofRegions)
			if (j==N_slabs/2+1) cold = selectCold(listofRegions)
		end do
		
		Totenergies(k) = E()
		
		do j=1, N_slabs
			regions(j)%temps(k)=temperatures(j)
		end do

	end subroutine rnem

	
	subroutine swapAtoms(k)
	  !! Swapping atoms' velocities
	  !! to swap masses (add when needed):
	  !! - either change their types, or swap atoms of same type
		integer,  intent(in)::k
		type(ad_t), dimension(3)::swapv
		
		regions(1)%Kenergies(k/skip_swap) = KEi(cold)
		regions(2)%Kenergies(k/skip_swap) = KEi(hot)

		swapv = atoms(hot)%v
		atoms(hot)%v = atoms(cold)%v
		atoms(cold)%v = swapv
	end subroutine swapAtoms
	
end module integrate_mod
