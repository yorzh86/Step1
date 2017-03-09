module properties_mod
	use kinds_mod
	use integrate_mod
	use utilities_mod
	use units_mod
	use settings_mod
	use system_mod
	implicit none
	
	contains

	subroutine thermalConductivity()
		type(ad_t)::grad, flux, kSI
		!http://www.knowledgedoor.com/2/elements_handbook/thermal_conductivity.html#argon

		flux = calculateFlux()
		grad = calculateGrad()

		kSI = flux/grad
		print *, 
		write(*,*)'kmetal:', kSI, '[eV/ps/A.K]'!!! 40K ~0.56 W/(m.K)
		print *, 
		
	contains
	
		function calculateFlux() result(o)
			type(ad_t)::o
			
			type(ad_t),dimension(:),allocatable::dE
			type(ad_t)::mdE, A, mdE1
			integer::i,j,skip1,k

			allocate(dE(0:N_steps/skip_swap))

			do k=0,N_steps/skip_swap
				dE(k) = regions(1)%Kenergies(k) - regions(2)%Kenergies(k)
			end do
			!write (*,'(1X,1A19, 10ES20.9 )')  "Array of KE1 - KE2:", dE

			!mdE = sum(dE)/real(N_steps/skip_swap+1, wp)
			mdE1 = sum(dE)/real(size(dE),wp) !what is the difference???

			A = (lattice_const**2)*real(latM(1)*latM(2),wp)
			
			o = mdE1/(dt*real(skip_swap,wp)*A)
			
			write(*,*)"flux:", o

		end function calculateFlux

		function calculateGrad() result(o)
			type(ad_t)::o,  mL, mR
			integer::i,j,k,skip2
			type(ad_t), dimension(:),allocatable::ar1, ar5, ar9
			type(ad_t)::m1, m5, m9, dz, slope1, slope2
			skip2 = 0
			allocate(ar1(0:N_steps))
			allocate(ar5(0:N_steps))
			allocate(ar9(0:N_steps))
			
			do k=0, N_steps
				ar1(k)= regions(2)%temps(k)
				ar5(k)= regions(6)%temps(k)
				ar9(k)= regions(10)%temps(k)
			end do
			
			m1 =  sum(ar1)/real(size(ar1),wp)
			m5 =  sum(ar5)/real(size(ar5),wp)
			m9 =  sum(ar9)/real(size(ar9),wp)
			
			dz = real(latM(3),wp)*lattice_const*0.4_wp
			slope1 = ((m5-m1)+(m5-m9))/2.0_wp/dz
			
			o = slope1

			write(*,*)"Slope: ", slope1
			
		end function calculateGrad
		

	end subroutine thermalConductivity
	
	subroutine specificHeat()
	!! subroutine used to calculate specific heat
	!! Cp = (1/M) * dE/dT   (p=const)
	!! Do not turn off barostate.
		
		type(ad_t):: mass, totEnergy, gradT
		type(ad_t):: Cp_metal, Cp_SI, metalToSI
		
		totEnergy = sum(Totenergies)/real(size(Totenergies), wp)
		
		gradT = calculateGrad()
		
		mass = real(size(atoms),wp)*types(atoms(1)%t)%m
		
		metalToSI = convert(1.0_wp, 'eV', 'J')/convert(1.0_wp, 'u', 'g')
		
		Cp_metal = (1.0_wp/mass)*(totEnergy/gradT)
		
		Cp_SI = Cp_metal*metalToSI 
		
		! should be about:
		! 34K 1bar 4.88 cal/mol.K = 0.51 J/g.K
		
		write (*,*) "mass:", mass, " [u]"
		write (*,*) "gradT:", gradT, " [K]"
		write (*,*) "total energy average:", totEnergy, " [eV]"
		
		write(*,*) "Cp_metal: ", Cp_metal, "[eV/u.K]"
		write(*,*) "Cp_SI: ", Cp_SI, "[J/g.K]"
		
		contains 
		
			
			function calculateGrad() result(o)
				type(ad_t)::o,  mL, mR
				integer::i,j,k,skip2
				type(ad_t), dimension(:),allocatable::ar1, ar5, ar9
				type(ad_t)::m1, m5, m9, dz, gradT
				skip2 = 0
				allocate(ar1(0:N_steps))
				allocate(ar5(0:N_steps))
				allocate(ar9(0:N_steps))
				
				do k=0, N_steps
					ar1(k)= regions(2)%temps(k)
					ar5(k)= regions(6)%temps(k)
					ar9(k)= regions(10)%temps(k)
				end do
				
				m1 =  sum(ar1)/real(size(ar1),wp)
				m5 =  sum(ar5)/real(size(ar5),wp)
				m9 =  sum(ar9)/real(size(ar9),wp)
				
				gradT = ((m5-m1)+(m5-m9))/2.0_wp
				
				o = gradT
				
			end function calculateGrad
	
	end subroutine specificHeat
	
	
end module properties_mod
