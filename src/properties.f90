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

		flux = calculateFlux()
		grad = calculateGrad()

		kSI = flux/grad
		print *, 
		write(*,*)'kSI:', kSI, '[W/m.K]'
		print *, 
		
	contains
	
		function calculateFlux() result(o)
			type(ad_t)::o
			
			type(ad_t),dimension(:),allocatable::dE
			type(ad_t)::mdE, A
			integer::i,j,skip1,k

			allocate(dE(0:N_steps/skip_swap))

			do k=0,N_steps/skip_swap
				dE(k) = regions(1)%energies(k) - regions(2)%energies(k)
			end do
			write (*,'(1X,1A19, 10ES20.9 )')  "Array of KE1 - KE2:", real(dE)

			mdE = sum(dE)/real(N_steps/skip_swap+1, wp)
			!mdE = sum(dE)/real(size(dE),wp)
			A = (lattice_const**2)*real(latM(1)*latM(2),wp)
				
			o = mdE/(dt*real(skip_swap,wp)*A)

			write (*,*)
		end function calculateFlux

		function calculateGrad() result(o)
			type(ad_t)::o,  mL, mR
			integer::i,j,k,skip2
			type(ad_t), dimension(:),allocatable::ar1, ar5, ar9
			type(ad_t)::m1,m5,m9, dz, slope1, slope2
			skip2 = 0
			allocate(ar1(0:N_steps))
			allocate(ar5(0:N_steps))
			allocate(ar9(0:N_steps))
			
			do k=0, N_steps
				ar1(k)= regions(2)%temps(k)%x
				ar5(k)= regions(6)%temps(k)%x
				ar9(k)= regions(10)%temps(k)%x
			end do

			m1 =  sum(ar1)/real(size(ar1),wp)
			m5 =  sum(ar5)/real(size(ar5),wp)
			m9 =  sum(ar9)/real(size(ar9),wp)
			
			dz = real(latM(3),wp)*lattice_const*0.4_wp
			slope1 = ((m5-m1)+(m5-m9))/2.0_wp/dz
			
			!mL,mR - slopes of left and right curves
			mL = slope(1,6)
			mR = slope(6,10)
			slope2 = (mL + mR)/2.0_wp
			
			o = slope2
			write(*,*)"Slope manually:", slope1
			write(*,*)"Slope least squares method:", slope2

		end function calculateGrad
		
		function slope(left, right) result (o)
			integer, intent(in)::left, right
			type(ad_t):: o
			type(ad_t)::avX, avY, sXY,sXX, sX, sY
			integer::i, j, k
			
			!sX - sum of x values (postions)
			!sXX - sum of squares
			!sXY - sum of products of corresponding x and av y (temperatures)
			!avX - mean of x values
			!avY - mean of y values
			sX  = 0.0_wp
			sY  = 0.0_wp
			sXX = 0.0_wp
			sXY = 0.0_wp
			avX = 0.0_wp
			avY = 0.0_wp
			
			do i=left, right
				sY = sY+sum(regions(i)%temps(0:N_steps))
			end do
			avY = sY/(real(N_steps,wp)*real(right-left+1,wp))
			
			do i=left,right
				sX  =  sX+real(i,wp)*lattice_const*real(latM(3),wp)/real(N_slabs,wp)
				sXX = sXX+(real(i,wp)*lattice_const*real(latM(3),wp)/real(N_slabs,wp))**2
			end do
			avX = sX/real(right-left+1,wp)

			do i=left, right
					sXY = sXY + sum(regions(i)%temps(0:N_steps))/real(N_steps,wp) * real(i,wp)*(lattice_const*real(latM(3),wp)/real(N_slabs,wp))
			end do

			!equation for slope:
			o = (sXY-sX*avY)/(sXX-sX*avX)
		end function slope

	end subroutine thermalConductivity
	
end module properties_mod
