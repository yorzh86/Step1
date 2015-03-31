program main_prg
	use kinds_mod
	use system_mod
	use output_mod
	implicit none
	
	integer::iou
	
	open(file='out.xyz',newunit=iou)
	
	call buildSystem(1.0_wp,3,20.0_wp)
	call writeStepXYZ(iou)
	
	close(iou)
	
contains



end program main_prg 
