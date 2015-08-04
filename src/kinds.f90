module kinds_mod
	implicit none
	private
	
	!==============!
	!= Real Kinds =!
	!==============!
	
	integer,parameter::sp = selected_real_kind(6)
	integer,parameter::dp = selected_real_kind(15)
	integer,parameter::ep = selected_real_kind(18)
	integer,parameter::qp = selected_real_kind(32)
	integer,parameter::wp = ep
	
	!==================!
	!= Math Constants =!
	!==================!
	
	real(wp),parameter::PI = 4.0_wp*atan(1.0_wp)
!~ 	real(wp),parameter::E  = exp(1.0_wp)
	
	!===========!
	!= Exports =!
	!===========!
	
	public::wp
	public::PI
	
	public::printTypes
	public::randomNormal
	
contains

	subroutine printTypes
		write(*,*) 'sp: ',sp
		write(*,*) 'dp: ',dp
		write(*,*) 'ep: ',ep
		write(*,*) 'qp: ',qp
		write(*,*) 'wp: ',wp
	end subroutine printTypes

	function randomNormal() result(o)
		real(wp)::o
		real(wp),dimension(12)::x
		
		call random_number(x)
		o = sum(x)-6.0_wp
	end function randomNormal

end module kinds_mod
