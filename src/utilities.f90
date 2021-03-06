module utilities_mod
	!! Utility module containing miscellaneous tools that don't
	!! quite fit anywhere else.
	use iso_fortran_env
	use kinds_mod
	implicit none
	private
	
	integer,parameter::stdin  = INPUT_UNIT
	integer,parameter::stdout = OUTPUT_UNIT
	
	interface mixval
		!! Return a 2-vector comprising the minimum and maximum values of an array
		module procedure mixval_1
		module procedure mixval_2
		module procedure mixval_3
	end interface
	
	interface span
		!! Return a the maximum-minumum values of an array
		module procedure span_1
		module procedure span_2
		module procedure span_3
	end interface
	
	interface flatten
		!! Reduce an array to one dimension
		module procedure flatten_2
		module procedure flatten_3
	end interface
	
	public::stdin
	public::stdout
	
	public::mixval
	public::span
	public::linspace
	
	public::startsWith
	public::endsWith
	
	public::meshGridX
	public::meshGridY
	
	public::setRandomSeed
	public::randomNormal
	public::randomUniform
	public::mean
	public::standardDeviation
	
	public::flatten
	
	public::colorize
	public::int2char
	public::real2char
	
	public::showProgress
	
contains

	function mixval_1(x) result(b)
		!! Return [hi,low] for an array
		real(wp),dimension(:),intent(in)::x
			!! Array to find extrema in
		real(wp),dimension(2)::b
		
		b = [minval(x),maxval(x)]
	end function mixval_1

	function mixval_2(x) result(b)
		!! Return [hi,low] for an array
		real(wp),dimension(:,:),intent(in)::x
			!! Array to find extrema in
		real(wp),dimension(2)::b
		
		b = [minval(x),maxval(x)]
	end function mixval_2

	function mixval_3(x) result(b)
		!! Return [hi,low] for an array
		real(wp),dimension(:,:,:),intent(in)::x
			!! Array to find extrema in
		real(wp),dimension(2)::b
		
		b = [minval(x),maxval(x)]
	end function mixval_3

	function span_1(x) result(o)
		real(wp),dimension(:),intent(in)::x
		real(wp)::o
		
		o = maxval(x)-minval(x)
	end function span_1

	function span_2(x) result(o)
		real(wp),dimension(:,:),intent(in)::x
		real(wp)::o
		
		o = maxval(x)-minval(x)
	end function span_2

	function span_3(x) result(o)
		real(wp),dimension(:,:,:),intent(in)::x
		real(wp)::o
		
		o = maxval(x)-minval(x)
	end function span_3

	function linspace(l,h,N) result(o)
		!! Return an array of evenly-spaced values
		real(wp),intent(in)::l
			!! Low-bound for values
		real(wp),intent(in)::h
			!! High-bound for values
		integer,intent(in),optional::N
			!! Number of values (default 20)
		real(wp),dimension(:),allocatable::o
		
		integer::Nl,i
		
		Nl = 20
		if(present(N)) Nl = N
		
		o = [( (h-l)*real(i-1,wp)/real(Nl-1,wp)+l , i=1 , Nl )]
	end function linspace

	function startsWith(text,str) result(o)
		!! Test if text starts with str
		character(*),intent(in)::text
			!! Text to search
		character(*),intent(in)::str
			!! String to look for
		logical::o
		integer::k
		
		k = len(str)
		o = text(1:k)==str
	end function startsWith

	function endsWith(text,str) result(o)
		!! Test if text ends with str
		character(*),intent(in)::text
			!! Text to search
		character(*),intent(in)::str
			!! String to look for
		logical::o
		integer::k
		
		k = len(text)
		o = text(k-len(str)+1:k)==str
	end function endsWith

	subroutine setRandomSeed(S)
		integer::S
		integer::k,N

		call random_seed(size=N)
		call random_seed(put=[(k-1,k=1,N)]*S)
	end subroutine setRandomSeed
	
	function randomNormal() result(o)
		!! Return a sample from an approximate normal distribution
		!! with a mean of \(\mu=0\) and a standard deviation of
		!! \(\sigma=1\). In this approximate distribution, \(x\in[-6,6]\).
		real(wp)::o
		real(wp),dimension(12)::x

		call random_number(x)
		o = sum(x)-6.0_wp
	end function randomNormal

	function randomUniform() result(o)
		!! Return a sample from a uniform distribution
		!! in the range \(x\in[-1,1]\).
		real(wp)::o

		call random_number(o)
		o = o*2.0_wp-1.0_wp
	end function randomUniform

	function flatten_2(A) result(o)
		real(wp),dimension(:,:),intent(in)::A
		real(wp),dimension(:),allocatable::o
		
		o = reshape(A,[size(A)])
	end function flatten_2

	function flatten_3(A) result(o)
		real(wp),dimension(:,:,:),intent(in)::A
		real(wp),dimension(:),allocatable::o
		
		o = reshape(A,[size(A)])
	end function flatten_3

	function meshGridX(x,y) result(o)
		real(wp),dimension(:),intent(in)::x,y
		real(wp),dimension(:,:),allocatable::o
		
		integer::Nx,Ny
		integer::i,j
		
		Nx = size(x)
		Ny = size(y)
		
		allocate(o(Nx,Ny))
		
		forall(i=1:Nx,j=1:Ny) o(i,j) = x(i)
	end function meshGridX

	function meshGridY(x,y) result(o)
		real(wp),dimension(:),intent(in)::x,y
		real(wp),dimension(:,:),allocatable::o
		
		integer::Nx,Ny
		integer::i,j
		
		Nx = size(x)
		Ny = size(y)
		
		allocate(o(Nx,Ny))
		
		forall(i=1:Nx,j=1:Ny) o(i,j) = y(j)
	end function meshGridY

	function colorize(s,c) result(o)
		character(*),intent(in)::s
		integer,dimension(3)::c ! c in [0,5]
		character(:),allocatable::o
		
		character(1),parameter::CR  = achar(13)
		character(1),parameter::ESC = achar(27)
		
		character(20)::pre
		character(3)::cb
		
		write(cb,'(1I3)') 6**2*c(1)+6*c(2)+c(3)+16
		pre = ESC//'[38;5;'//trim(adjustl(cb))//'m'
		o = trim(pre)//s//ESC//'[0m'
	end function colorize

	elemental function real2char(a,f,l) result(o)
		real(wp),intent(in)::a
		character(*),optional,intent(in)::f
		integer,optional,intent(in)::l
		character(:),allocatable::o
		
		character(128)::buf
		
		if(present(l)) then
			allocate(character(l)::o)
			if(present(f)) then
				write(o,'('//f//')') a
			else
				write(o,*) a
			end if
		else
			if(present(f)) then
				write(buf,'('//f//')') a
			else
				write(buf,*) a
			end if
			o = trim(adjustl(buf))
		end if
	end function real2char

	elemental function real2time(a) result(o)
		real(wp),intent(in)::a
		character(:),allocatable::o
		
		integer::d,r,h,m,s
		
		r = nint(a)
		
		d = r/(3600*24)
		r = mod(r,3600*24)
		
		h = r/3600
		r = mod(r,3600)
		
		m = r/60
		r = mod(r,60)
		
		s = r
		
		o = ''
		if(d>0) o = o//int2char(d)//'d '
		if(h>0.or.d>0) o = o//int2char(h)//'h '
		if(m>0.or.h>0.or.d>0) o = o//int2char(m)//'m '
		o = o//int2char(s)//'s'
	end function real2time

	elemental function int2char(a,f,l) result(o)
		integer,intent(in)::a
		character(*),optional,intent(in)::f
		integer,optional,intent(in)::l
		character(:),allocatable::o
		
		character(128)::buf
		
		if(present(l)) then
			allocate(character(l)::o)
			if(present(f)) then
				write(o,'('//f//')') a
			else
				write(o,*) a
			end if
		else
			if(present(f)) then
				write(buf,'('//f//')') a
			else
				write(buf,*) a
			end if
			o = trim(adjustl(buf))
		end if
	end function int2char

	subroutine showProgress(m,p,ml)
		character(*),intent(in)::m
		real(wp),intent(in)::p
		integer,intent(in),optional::ml
		
		real(wp)::r
		real(wp),save::po
		real(wp),save::tStart
		real(wp)::tNow
		integer::mld
		integer::N,k
		
		N = 40
		mld = 40
		if(present(ml)) mld = ml
		
		if(p<=0.0_wp) then
!~ 			write(stdout,'(1A)') ''
			po = p
			tStart = wallTime()
		else if(p-po<=0.01 .and. p<1.0_wp) then
			return
		else
			po = p
		end if
		tNow = wallTime()
		
		write(stdout,'(1A)',advance='no') achar(13)//colorize(m//repeat(' ',mld-len(m))//' [',[5,5,0])
		do k=1,N
			r = real(k-1,wp)/real(N-1,wp)
			if(r<=p) then
				write(stdout,'(1A)',advance='no') colorize('=',cmap(r,[0.0_wp,1.0_wp]))
			else
				write(stdout,'(1A)',advance='no') colorize(' ',[0,0,0])
			end if
		end do
		write(stdout,'(1A,1A,1X,1A,1A,1A,1A,1A,1A)',advance='no') colorize('] ',[5,5,0]), &
		& colorize(real2char(100.0_wp*p,'1F5.1'),cmap(p,[0.0_wp,1.0_wp])), &
		& colorize('%',[5,5,0]), colorize(' (',[5,5,0]), real2time(tNow-tStart), &
		& colorize(' / ',[5,5,0]), real2time( (tNow-tStart)/(p+0.0001_wp) ), colorize(')',[5,5,0])
		if(p>=1.0_wp) write(stdout,'(1A)') ''
		flush(stdout)
	end subroutine showProgress

	function cmap(v,r) result(c)
		real(wp),intent(in)::v
		real(wp),dimension(2),intent(in)::r
		integer,dimension(3)::c
		
		integer::s
		
		if(v<sum(r)/2.0_wp) then
			s = nint((v-r(1))/(sum(r)/2.0_wp-r(1))*5.0_wp)
			c = [s,s,5]
		else
			s = 5-nint((v-sum(r)/2.0_wp)/(r(2)-sum(r)/2.0_wp)*5.0_wp)
			c = [5,s,s]
		end if
	end function cmap

	function mean(d) result(o)
		real(wp),dimension(:),intent(in)::d
		real(wp)::o
		
		o = sum(d)/real(size(d),wp)
	end function mean

	function standardDeviation(d) result(o)
		real(wp),dimension(:),intent(in)::d
		real(wp)::o
		
		o = sqrt(sum((d-mean(d))**2)/real(size(d)-1,wp))
	end function standardDeviation

	function wallTime() result(o)
		real(wp)::o
		
		integer,parameter::ip = selected_int_kind(15)
		integer(ip)::ticks,tickRate,r
		
		call system_clock(ticks,tickRate)
		
		r = mod(ticks,tickRate)
		
		o = real(ticks/tickRate,wp)+real(r,wp)/real(tickRate,wp)
	end function wallTime
	
end module utilities_mod
