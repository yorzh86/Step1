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
		!real(wp)::grad, flux, kSI
		type(ad_t)::grad, flux, kSI

		flux = calculateFlux()
		grad = calculateGrad()

		kSI = flux/grad
		print *, 
		write(*,'(1X, 1A4, 1F8.4, 1A8 )')'kSI:', real(kSI), '[W/m.K]'
		print *, 
		
	contains
	
		function calculateFlux() result(o)
			type(ad_t)::o
			
			type(ad_t),dimension(:),allocatable::dE
			type(ad_t)::mdE, A
			integer::i,j,skip1

			skip1 = 0
			j = 0
				
			do i=0, N_steps
				if(allocated(regions(i)%energies)) j=j+1
			end do
			allocate(dE(j))

			j=0
			do i=0, N_steps
				if(allocated(regions(i)%energies)) then
					j = j+1
					dE(j) = regions(i)%energies(3) - regions(i)%energies(4)
				end if
			end do

			mdE = sum(dE)/real(size(dE),wp)
			A = (lattice_const**2)*real(latM(1)*latM(2),wp)
				
			o = mdE/(dt*real(skip_swap,wp)*A)
		end function calculateFlux

		function calculateGrad() result(o)
			type(ad_t)::o,  mL, mR
			!type(ad_t)::avX, avY, sXY,sXX, sX, mL, mR
			integer::i,j,k,skip2
			type(ad_t), dimension(:),allocatable::ar1, ar5, ar9
			type(ad_t)::m1,m5,m9, dz, slope1
			skip2 = 0
			allocate(ar1(0:N_steps))
			allocate(ar5(0:N_steps))
			allocate(ar9(0:N_steps))
			do k=0, N_steps
				ar1(k)= regions(k)%temps(4)
				ar5(k)= regions(k)%temps(8)
				ar9(k)= regions(k)%temps(12)
			end do
			
			m1 =  sum(ar1)/real(size(ar1),wp)
			m5 =  sum(ar5)/real(size(ar5),wp)
			m9 =  sum(ar9)/real(size(ar9),wp)
			dz = real(latM(3),wp)*lattice_const*0.4_wp
			slope1 = ((m5-m1)+(m5-m9))/2.0_wp/dz
			print *, 
			print *, 'mean max', m5
			print *, 'mean min 1', m1
			print *, 'mean min 9', m9
			print *, "slope1", slope1

			print *, 
			print *, 'using least squares method:'
			!mL,mR - slopes of left and right curves
			mL = slope(3,8)
			mR = slope(8,12)
			print *, 'mL', mL
			print *, 'mR', mR
			
			o = (mL + mR)/2.0_wp
			print *, 'total slope', o
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
			
			do i=0,N_steps
				sY = sY+sum(regions(i)%temps(left:right))
			end do
			avY = sY/(real(N_steps,wp)*real(right-left,wp))
				
			do i=1,N_slabs/2+1
				sX  =  sX+lattice_const*real(latM(3),wp)/real(N_slabs*i,wp)
				sXX = sXX+(lattice_const*real(latM(3),wp)/real(N_slabs*i,wp))**2
			end do
			avX = sX/real(N_slabs/2+1,wp)
				
			do i=0,N_steps
				do j=left, right
					sXY = sXY + (regions(i)%temps(j)/real(N_steps,wp)) * (lattice_const*real(latM(3),wp)/real(N_slabs*(j-2),wp))
				end do
			end do
			
			!equation for slope:
			o = (sXY-sX*avY)/(sXX-sX*avX)
			end function slope

	end subroutine thermalConductivity

end module properties_mod
