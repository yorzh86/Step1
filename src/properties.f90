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
			write(*,'(1X, 1A4, 1F8.4, 1A8 )')'kSI:', kSI, '[W/m.K]'
			print *, 
			
			contains

				function calculateFlux() result(o)
					type(ad_t),dimension(:),allocatable::dE
					type(ad_t)::mdE, A
					integer::i,j, skip1
					type(ad_t):o
					
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

					mdE = sum(dE)/size(dE)
					A = (lattice_const**2)*latM(1)*latM(2)

					o = mdE/dt/skip_swap/A
				end function calculateFlux
				
				function calculateGrad() result(o)
					type(ad_t)::o
					type(ad_t), dimension(:),allocatable::ar0, ar5, ar9
					type(ad_t)::m0, m5, m9, dz  
					integer::i,j,skip2

					skip2 = 0

					allocate(ar0(0:N_steps))
					allocate(ar5(0:N_steps))
					allocate(ar9(0:N_steps))
					
					do i=0, N_steps
						ar0(i)= regions(i)%temps(3)
						ar5(i)= regions(i)%temps(8)
						ar9(i)= regions(i)%temps(12)
					end do
					
					m0 =  sum(ar0)/size(ar0)
					m5 =  sum(ar5)/size(ar5)
					m9 =  sum(ar9)/size(ar9)
					dz = latM(3)*lattice_const/2.0_wp

					o = ((m5-m0)+(m5-m9))/2.0_wp/dz !!!!INCORRECT!
				end function calculateGrad
				
		end subroutine thermalConductivity
end module properties_mod
