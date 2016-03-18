module system_mod
	use kinds_mod
	use utilities_mod
	use units_mod
	use settings_mod
	implicit none
	
	!=========!
	!= Types =!
	!=========!
	
	type::type_t
		real(wp)::m
			!! Atomic mass
		character(2)::atom_name
			!! Atom symbol (use standard periodic tables)
	end type
	
	type::atom_t
		real(wp),dimension(3)::r
			!! Atomic position
		real(wp),dimension(3)::v
			!! Atomic velocity
		real(wp),dimension(3)::a
			!! Atomic acceleration
		real(wp),dimension(3)::f
			!! Atomic force
		integer::atom_id
		integer::t
			!! Atom type
		integer,dimension(:),allocatable::neighbors
			!! Nearest neighbors
	end type
	
	
	!===============!
	!= Module Data =!
	!===============!
	
	logical::enableLennardJones = .false.
		!! Run Lennard-Jones potential
	
	type(type_t),dimension(:),allocatable::types
		!! Types for all atoms in system
	type(atom_t),dimension(:),allocatable::atoms
		!! All atoms in system
	
	real(wp)::Teta = 0.0_wp
		!! Thermostat DOF
	real(wp)::Pepsilon = 0.0_wp !was used as 0.01
		!! Barostat DOF
	real(wp),dimension(3)::box
		!! Bounds of the simulation box
		
	integer::ts
		!! Time step counter
	real(wp)::t
		!! System time
		
		
	
contains

	subroutine SimpleSystem (delta, N, Ti)
	implicit none
	integer, intent(in):: N
		!! N - number of atoms
	real(wp), intent(in):: delta, Ti
	
	integer::i, j, k

	box = delta*N
		
	allocate(types(1))
	allocate(atoms(N))
	
	types%m = convert(39.948_wp,'u','kg')
	types%atom_name = 'Ar'
	atoms(:)%t = 1
	
	do k=1, N
		atoms(k)%atom_id = k
		atoms(k)%r = atoms(k)%r + k*[1.0,0.0,0.0]*delta - [1.0,0.0,0.0]*delta/2.0_wp
	end do
	
	
	do k=1,N
 		!! Same treatment of velocities as in buildSystem
 		call random_number(atoms(k)%v)
 		atoms(k)%v = 2.0_wp*atoms(k)%v-1.0_wp
 		do while(norm2(atoms(k)%v)>1.0_wp .and. norm2(atoms(k)%v)<0.1_wp)
 			call random_number(atoms(k)%v)
 			atoms(k)%v = 2.0_wp*atoms(k)%v-1.0_wp
 		end do
 		atoms(k)%v = atoms(k)%v/norm2(atoms(k)%v)
 		atoms(k)%v = atoms(k)%v*sqrt(2.0_wp*kB*Ti/types(atoms(k)%t)%m)*abs(randomNormal()+1.0_wp)
 	end do
 	forall(k=1:3) atoms(:)%v(k) = atoms(:)%v(k)-sum(atoms(:)%v(k))/real(size(atoms),wp)
 	call updateAllNeighbors()
 	do k=1,N
 		atoms(k)%a = -delV(k)/types(atoms(k)%t)%m
 		atoms(k)%f = -delV(k)
 	end do
 	ts = 0
 	t  = 0.0_wp
 	
	end subroutine SimpleSystem
	
	

	subroutine buildSystem(a,N,Ti)
		implicit none
		real(wp), dimension(3)::posit, velocity, force
		real(wp),intent(in)::a
			!! Lattice constant
		integer,dimension(3),intent(in)::N
			!! Number of unit cells
		real(wp),intent(in)::Ti
			!! Initial temperature
		real(wp), dimension(3,4), parameter::rcell= &
		reshape([0.0_wp, 0.0_wp, 0.0_wp, &
		         0.5_wp, 0.5_wp, 0.0_wp, &
		         0.0_wp, 0.5_wp, 0.5_wp, &
		         0.5_wp, 0.0_wp, 0.5_wp  ], [3,4])
		integer::i,j,k,l,idx
		integer:: iou_lammps=1
		
		box = a*real(N,wp)
		
		allocate(types(1))
		allocate(atoms(size(rcell,2)*product(N)))
		
		types%m = convert(39.948_wp,'u','kg')
		types%atom_name = 'Ar'
		atoms(:)%t = 1
		
 		idx = 1
 		do k=1,N(3)
 			do j=1,N(2)
 				do i=1,N(1)
 					do l=1,size(rcell,2)
 						atoms(idx)%r = a*(real([i,j,k]-1,wp)+rcell(1:3,l))
 						idx = idx+1
 					end do
 					
 				end do
 			end do
 		end do
 		
 		do k=1,size(atoms)
 			!! Create random direction of velocities
 			call random_number(atoms(k)%v)
 			atoms(k)%v = 2.0_wp*atoms(k)%v-1.0_wp
 			do while(norm2(atoms(k)%v)>1.0_wp .and. norm2(atoms(k)%v)<0.1_wp)
 				call random_number(atoms(k)%v)
 				atoms(k)%v = 2.0_wp*atoms(k)%v-1.0_wp
 			end do
 			atoms(k)%v = atoms(k)%v/norm2(atoms(k)%v)
 			!! Set velocity magnitude
 			atoms(k)%v = atoms(k)%v*sqrt(2.0_wp*kB*Ti/types(atoms(k)%t)%m)*abs(randomNormal()+1.0_wp)
 		end do
 		forall(k=1:3) atoms(:)%v(k) = atoms(:)%v(k)-sum(atoms(:)%v(k))/real(size(atoms),wp)
 		
 		call updateAllNeighbors()
 		
 		do k=1,size(atoms)
 			atoms(k)%a = -delV(k)/types(atoms(k)%t)%m
 			atoms(k)%f = -delV(k)
 		end do
 		
! 		do k=1, 9
!			read(iou_lammps,*)
!		end do
			
!		do k=1, size(atoms)
!			read(1, '(1I4, 3F5.1, 6F13.9)') atoms(k)%atom_id, posit, velocity, force
!			atoms(k)%r = posit/1.0E10_wp
!			atoms(k)%v = velocity/1.0E-2_wp
!			atoms(k)%f = force/6.24150636309E8_wp
!			atoms(k)%a = -delV(k)/types(atoms(k)%t)%m
!		end do
 		ts = -0 !! CHANGE BACK TO ZERO!!
 		t  = 0.0_wp
	end subroutine buildSystem

	subroutine writeLammpsData(fn)
		character(*),intent(in)::fn
		
		integer::i,k,iou
		real(wp)::E0,S0
		
		E0 = lj%coeffs(1)
		S0 = lj%coeffs(2)
		
		open(file=fn,newunit=iou)
		
		write(iou,'(1A)') '# Input geometry for lammps'
		write(iou,'(1I9,1X,1A)') size(atoms),'atoms'
		write(iou,'(1I9,1X,1A)') size(types),'atom types'
		write(iou,'(2F13.6,1X,1A,1X,1A)') 0.0_wp,convert(box(1),'m','A'),'xlo','xhi'
		write(iou,'(2F13.6,1X,1A,1X,1A)') 0.0_wp,convert(box(2),'m','A'),'ylo','yhi'
		write(iou,'(2F13.6,1X,1A,1X,1A)') 0.0_wp,convert(box(3),'m','A'),'zlo','zhi'
		write(iou,'(1A)') ''
		write(iou,'(1A)') 'Masses'
		write(iou,'(1A)') ''
		do k=1,size(types)
			write(iou,'(1I4,1X,1F13.6)') k,convert(types(k)%m,'kg','u')
		end do
		write(iou,'(1A)') ''
		write(iou,'(1A)') 'PairIJ Coeffs # lj/cut'
		write(iou,'(1A)') ''
		write(iou,'(2I3,3F13.6)') 1, 1 , convert(E0,'J','eV'), convert(S0,'m','A'), convert(lj%cutoff,'m','A')
		write(iou,'(1A)') ''
		write(iou,'(1A)') 'Atoms'
		write(iou,'(1A)') ''
		do k=1,size(atoms)
			write(iou,'(1I9,1X,1I2,1X,3E25.15)') k,atoms(k)%t,[( convert(atoms(k)%r(i),'m','A') , i=1,3 )]
		end do
		write(iou,'(1A)') ''
		write(iou,'(1A)') 'Velocities'
		write(iou,'(1A)') ''
		do k=1,size(atoms)
			write(iou,'(1I9,1X,1X,3E25.15)') k,[( convert(atoms(k)%v(i),'m/s','A/ps') , i=1,3 )]
		end do
		
		write(iou,'(1A)') ''
		write(iou,'(1A)') 'Forces'
		write(iou,'(1A)') ''
		do k=1,size(atoms)
			write(iou,'(1I9,1X,1X,3E25.15)') k,[( convert(atoms(k)%f(i),'N','eV/A') , i=1,3 )]
		end do
		
		close(iou)
	end subroutine writeLammpsData

	pure function V(i) result(o)
		integer,intent(in)::i
		real(wp)::o
		
		o = 0.0_wp
		if(enableLennardJones) call doLennardJones(o)
		
	contains
	
		pure subroutine doLennardJones(o)
			!real(wp),intent(inout)::o
			real(wp),intent(out)::o
			
			real(wp),dimension(3)::d
			real(wp)::l,E0,S0
			integer::j,aj
			
			E0 = lj%coeffs(1)
			S0 = lj%coeffs(2)
			
			do j=1,size(atoms(i)%neighbors)
				aj = atoms(i)%neighbors(j)
				d = deltaR(atoms(aj),atoms(i))
				if( norm2(d)>lj%cutoff ) cycle
				l = S0/norm2(d)
				o = o+( 0.5_wp ) * ( 4.0_wp*E0*l**6*(l**6-1.0_wp) )
					!! Only include half of the bond energy for each atom in the bond
			end do
			!do j=1,size(atoms)
			!	if(j==i) cycle
			!	d = deltaR(atoms(j),atoms(i))
			!	l = S0/norm2(d)
			!	o = o+(0.5_wp) * 4.0_wp*E0*(l**12.0_wp-l**6.0_wp)
			!end do
			
		end subroutine doLennardJones
	
	end function V

	pure function delV(i) result(o)
	!! Total force on atom
		integer,intent(in)::i
		real(wp),dimension(3)::o
		real(wp),dimension(3)::r
		
		integer::j,aj
		
		o = 0.0_wp
		do j=1,size(atoms(i)%neighbors)
			aj = atoms(i)%neighbors(j)
			r  = deltaR(atoms(i),atoms(aj))
			if( norm2(r)>lj%cutoff ) cycle
			o = o+delVij(i,aj,r)
			!o = o + atoms(i)%f
		end do
	end function delV
	
	pure function fnorm() result(o)
		real(wp)::o
		integer::k
		
		o = 0.0_wp
		do k=1,size(atoms)
			o = o + sum(atoms(k)%f**2)
		end do
		o = sqrt(o)
	end function fnorm
	


	pure function delVij(i,j,d) result (o)
	!! Force between two atoms
		integer,intent(in)::i,j
		real(wp),dimension(3),intent(in)::d
		real(wp),dimension(3)::o
		
		o = 0.0_wp
		if(i==j) return
		
		if(enableLennardJones) call doLennardJones(o)
		
	contains
		
		pure subroutine doLennardJones(o)
			real(wp),dimension(3),intent(out)::o
			real(wp)::l,E0,S0
			
			E0 = lj%coeffs(1)
			S0 = lj%coeffs(2)
			
			l = S0/norm2(d)
			o = o+24.0_wp*E0/sum(d*d)*(l**6)*(1.0_wp-2.0_wp*l**6)*d
			!o = atoms(i)%f
		end subroutine doLennardJones
		
	end function delVij

	pure function deltaR(a1,a2) result(o)
		type(atom_t),intent(in)::a1,a2
		real(wp),dimension(3)::o
		
		real(wp),dimension(3)::d
		d = a1%r-a2%r
		o = d-box*real(nint(d/box),wp)
		
	end function deltaR

	pure function temperature() result(o)
		real(wp)::o
		
		real(wp),dimension(3)::vBulk
		real(wp)::SKE
		integer::k
		
		forall(k=1:3) vBulk(k) = sum(atoms%v(k))/real(size(atoms),wp)
		
		SKE = 0.0_wp
		do k=1,size(atoms)
			SKE = SKE+KEi(k,vBulk)
		end do
		o = (2.0_wp*SKE)/(3.0_wp*real(size(atoms),wp)*kB)
	end function temperature

	pure function E() result(o)
		real(wp)::o
		integer::k
		
		o = 0.0_wp
		do k=1,size(atoms)
			o = o+Ei(k)
		end do
	end function E

	pure function Ei(i) result(o)
		integer,intent(in)::i
		real(wp)::o
		
		o = KEi(i) + PEi(i)
	end function Ei

	pure function KE() result (o) 
		real(wp)::o
		integer::k
		
		o = 0.0_wp
		do k=1,size(atoms)
			o = o+KEi(k)
		end do
	end function KE
	
	pure function KEi(i,vBulk) result (o)
		integer,intent(in)::i
		real(wp),dimension(3),intent(in),optional::vBulk
		real(wp)::o
		
		real(wp),dimension(3)::v0
		
		v0 = 0.0_wp
		if(present(vBulk)) v0 = vBulk
		o = 0.5_wp*types(atoms(i)%t)%m*norm2(atoms(i)%v-v0)**2
	end function KEi
	
	pure function PE() result (o)
		real(wp)::o
		integer::k
		
		o = 0.0_wp
		do k=1,size(atoms)
			o = o+PEi(k)
		end do
	end function PE

	pure function PEi(i) result (o)
		integer,intent(in)::i
		real(wp)::o
		
		o = V(i)
	end function PEi

	pure function Si(i) result(o)
		integer,intent(in)::i
		real(wp),dimension(3,3)::o
		
		real(wp),dimension(3)::v,r,F
		integer::j,aj
		
		v = atoms(i)%v
		o = -types(atoms(i)%t)%m*matmul(asCol(v),asRow(v))
		
		do j=1,size(atoms(i)%neighbors)
			aj = atoms(i)%neighbors(j)
			r  = deltaR(atoms(i),atoms(aj))
			if( norm2(r)>lj%cutoff ) cycle
			F = delVij(i,aj,r)
			!F = atoms(i)%f
			o = o-0.5_wp*(matmul(asCol(r),asRow(F))+matmul(asCol(F),asRow(r)))
		end do
	
	contains
	
		pure function asCol(v) result(o)
			real(wp),dimension(:),intent(in)::v
			real(wp),dimension(size(v),1)::o
			
			o(:,1) = v(:)
		end function asCol

		pure function asRow(v) result(o)
			real(wp),dimension(:),intent(in)::v
			real(wp),dimension(1,size(v))::o
			
			o(1,:) = v(:)
		end function asRow
		
	end function Si

	pure function heatflux() result(o)
		real(wp),dimension(3)::o
		
		integer::i,j,aj
		real(wp),dimension(3)::Fij,rij
		
		o = 0.0_wp

! 		do i=1,size(atoms)
! 			o = o+ ( Ei(i)*atoms(i)%v-matmul(Si(i),atoms(i)%v) )
! 		end do

		do i=1,size(atoms)
			o = o+Ei(i)*atoms(i)%v
			do j=1,size(atoms(i)%neighbors)
				aj = atoms(i)%neighbors(j)
				rij  = deltaR(atoms(i),atoms(aj))
				if( norm2(rij)>lj%cutoff ) cycle
				Fij = delVij(i,aj,rij)
				!Fij = atoms(i)%f
				o = o+dot_product(Fij,atoms(i)%v)*rij
			end do
		end do
	end function heatflux

	subroutine updateAllNeighbors()
		integer::k
		
		do k=1,size(atoms)
			call updateNeighbors(k)
		end do
	end subroutine updateAllNeighbors

	subroutine updateNeighbors(i)
		integer,intent(in)::i
		
		real(wp),dimension(3)::r
		integer::k
		
		atoms(i)%neighbors = [integer::]
		do k=1,size(atoms)
			if(i==k) cycle
			r  = deltaR(atoms(i),atoms(k))
 			if( norm2(r)>(lj%cutoff+lj%skin) ) cycle
			atoms(i)%neighbors = [atoms(i)%neighbors,k]
		end do
	end subroutine updateNeighbors

	function averageNeighbors() result(o)
		real(wp)::o
		integer::k
		
		o = sum([( size(atoms(k)%neighbors) , k=1,size(atoms) )]) / real(size(atoms),wp)/2.0_wp
	end function averageNeighbors

	pure function virial() result(o)
		real(wp)::o
		
		real(wp),dimension(3)::F,r
		integer::i,j,aj
		
		o = 0.0_wp
		do i=1,size(atoms)
			do j=1,size(atoms(i)%neighbors)
				aj = atoms(i)%neighbors(j)
				r  = deltaR(atoms(i),atoms(aj))
				if( norm2(r)>lj%cutoff ) cycle
				F = delVij(i,aj,r)
				o = o+0.5_wp*dot_product(F,r)
			end do
		end do
	end function virial

	pure function pressure() result(o)
		real(wp)::o
		
		o = (real(size(atoms),wp)*kB*temperature()-virial()/3.0_wp)/product(box)
	end function pressure

end module system_mod
