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
	real(wp),parameter::E  = exp(1.0_wp)
	
	!===========!
	!= Exports =!
	!===========!
	
	public::wp
	public::PI,E
	
	public::printTypes
	
contains

	subroutine printTypes
		write(*,*) 'sp: ',sp
		write(*,*) 'dp: ',dp
		write(*,*) 'ep: ',ep
		write(*,*) 'qp: ',qp
		write(*,*) 'wp: ',wp
	end subroutine printTypes

end module kinds_mod
