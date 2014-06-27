!###################################################################
!#                                                                 #
!#                 Ines Flores (L. M. 05/Jun/2007)                 #
!#                                                                 #
!# Reduced version of the program read_snapshot_single.f90         #
!# in order to read only the header of the simulation.             #
!# Check it for further information on its origin                  #
!#                                                                 #
!###################################################################

! NOTE 1: Snapshot is written with BIG ENDIAN FORMAT => Should be
!         compiled with -convert big_endian option or afterwards
!         executing the command setenv F_UFMTENDIAN "big"

! NOTE 2: Reading the simulations requires lots of RAM memory (~1GB) so
!         before executing this program check that the stacksize status
!         is unlimited (type limit in the shell). If it's not unlimited,
!         type: limit stacksize unlimited
!         (for tcsh shell)

PROGRAM read_header

IMPLICIT NONE

! Declaration of variables
CHARACTER(LEN=128)  :: fname, flnm_out
CHARACTER*4         :: block_id
INTEGER*4           :: npart(6), flag_sfr, flag_feedback, npartot(6), N, N_gas
INTEGER*4           :: Cool, n_files
REAL*8              :: massarr(6), time, redshift, L, h, Omega_m, Omega_L, Box
INTEGER             :: u_in, u_out, I

u_in  = 20
u_out = 30

PRINT*, 'Write the filename or path of the simulation'
READ(*,*) fname

OPEN(u_in, file=fname, form='unformatted')
PRINT*, 'Dentro de la lectura'
READ(u_in) block_id
PRINT*, block_id
READ(u_in) npart, massarr, time, redshift, flag_sfr, flag_feedback, npartot, cool, n_files, L, Omega_m, Omega_L, h

CLOSE(u_in)

! Conversion of units for the masses
massarr=massarr*1e10
! Idem for the box size
Box=L/1000

PRINT*, 'Write the output filename'
READ(*,*) flnm_out

OPEN(u_out, file=TRIM(flnm_out), status='unknown', action='write')

WRITE(u_out,*) 'Simulation:', TRIM(fname)
WRITE(u_out,*) '  '
WRITE(u_out,*) 'HEADER CONTENTS'
WRITE(u_out,*) '==============='
WRITE(u_out,*) '     N part       Masses(Msun)'
DO I=1,6
  WRITE(u_out,*) npart(i), massarr(i)
ENDDO
WRITE(u_out,*) '  '
WRITE(u_out,*) 'z=            ', redshift
WRITE(u_out,*) 'time=         ', time
WRITE(u_out,*) 'h0=           ', h
WRITE(u_out,*) 'Omega_m=      ', Omega_m
WRITE(u_out,*) 'Omega_L=      ', Omega_L
WRITE(u_out,*) 'Box(h^-1 Mpc)=', Box
WRITE(u_out,*) 'Flags=        ', flag_sfr, flag_feedback, cool
WRITE(u_out,*) 'a=            ', 1/(redshift+1)
WRITE(u_out,*) 'files=        ', n_files
CLOSE(u_out)

END PROGRAM read_header
