!     Last change:  EY   15 September 2018    8:43 am
module para
!
! Notes:
!
! value of cell_index =
!   1  : head (entire)
!   2  : neck (entire)
!   3  : upper chest (including shoulders and down to just below the nipples)
!   4  : lower abdomen (including stomach and down to the waist)
!   5  : arms (both left and right)
!   6  : hands (left and right)
!   7  : legs (left and right, front and back)
!   8  : feet (left and right)
!   9  : back
!  10  : buttocks
!  11  : crotch area (frontal area below the waist)
!
! Note - cell_index is more aptly referred to as body part index

implicit none
! challengeGV4a enforces double precision for all real variables
integer, parameter :: ep = kind(0.0d0)

integer, parameter :: nBodyParts = 11                   ! # body parts in division of the human form

integer, parameter :: nvertices = 89360                 ! # nodes (vertices) used in tesellation of human form
integer, parameter :: ncells = 178716                   ! # triangular cells (elements) used in tesselation of human form

integer, parameter :: nBins = 13                        ! # particle size bins used to represent droplet size distribution

! Parameters defining lognormal particle size distribution
! (depracated - use JEM to determine challenge levels instead)
!real(ep), parameter :: dgv = 15.0                           ! mass median diameter (um)
!real(ep), parameter :: sigG = 1.5                           ! geometric standard deviation

! Evaporation constant (D2 evaporation law)
! (depracated - use JEM to determine challenge levels instead )

!real(ep), parameter :: Kp = 0.013                           ! evaporation constant (um^2/s)
                                                        ! set Kp = 0, if one wants an invariant droplet size distribution with time

! Reference inlet wind speed used to obtain deposition probability in CFD simulations
! (depracated in current version)
!real(ep), parameter :: Uref = 2.5                           ! inlet wind speed (m/s) used for CFD simulation
                                                        ! used as dimensionless correction factor to reduce actual wind speed U used in determination
                                                        ! of deposition density on the body (probably not required, so not implemented)

! Probability of adhension

real(ep), parameter :: pAdhesion = 1.0_ep

!____________________________________________________________________________________________________________________________
! Some physical constants required for the computation of the normalized deposition velocity

real(ep), parameter :: kb = 1.38E-16_ep                        ! Boltzmann's constant in erg/K
real(ep), parameter :: lambda = 0.653E-5_ep                    ! mean free path of molecules at atmospheric pressure in cm
real(ep), parameter :: eta = 1.83e-4_ep                        ! dynamic viscosity of air in poise (g/(cm s))
real(ep), parameter :: rho = 0.001225_ep                       ! density of air in g/cm^3
real(ep), parameter :: rhoP = 1.0_ep                           ! density of particles in g/cm^3
real(ep), parameter :: nu = 0.149388_ep                        ! kinematic viscosity of air in cm^2/s = eta/rho
real(ep), parameter :: SRatio = 816.327_ep                     ! particle-to-air density ratio [-] = rhoP/rho
real(ep), parameter :: g = 980.665_ep                          ! acceleration due to gravity in cm/s^2

real(ep), parameter :: absT = 293.0_ep                         ! air temperature in K

!_____________________________________________________________________________________________________________________________

!____________________
!____________________

integer :: iProtect                                     ! indicator: iProtect = 0 (no protection - penetration probability = 1 for all
                                                        ! particle size bins); iProtect /= 0 (protection with clothing ensemble - penetration
                                                        ! probability encoded in routine pene
                                                        ! For protected human in IPE, the material systems available are as follows:
                                                        !
                                                        ! 1. iProtect = 1, Combat uniform
                                                        ! 2. iProtect = 2, Horizon 1
                                                        ! 3. iProtect = 3, CPCU MkIII
                                                        ! 4. iProtect = 4, CPCU MkIV
                                                        !
integer :: iWrite                                       ! indicator: iWrite = 0, do not write out deposition density files (ASCII and TecPlot);
                                                        ! iWrite /= 0, write out deposition density files (ASCII and TecPlot)
                                                        ! Note - iWrite = 10, write TecPlot 360 file, any other non-zero value of
                                                        ! iWrite, write TecPlot 10 file

integer :: kIndx                                        ! index denoting receptor location

integer :: nTimeSteps                                   ! # time steps in challenge level file (containing concentration-time history)
integer :: nFiles                                       ! # files to process (each file corresponds to receptor at a given
                                                        ! distance downwind of the release) => nFiles = # receptor locations
real(ep) :: dt                                          ! time increment for each time step (s) used in challenge level calculation
real(ep) :: U                                           ! advection (wind) velocity (m/s) - used for calculation of mass flux
!___________________________________________________________________________________________________________________________________
real(ep) :: Umodel                                          ! model advection velocity (m/s) used in determination of u* distribution on human form
                                                        ! Umodel is used in calculation of deposition probability (using u* to get deposition velocity vd)
!___________________________________________________________________________________________________________________________________

integer, dimension(ncells) :: cell_index                ! body part index (1 - 11) - see note above
integer, dimension(ncells,3) :: vert                    ! connectivity information (which nodes are connected
                                                        ! to each other to form a triangular cell or element)

real(ep) :: totMass                                     ! total mass deposited on the human body, normalized by Q

real(ep), dimension(nvertices) :: xv, yv, zv            ! coordinates of nodes: x, y, and z - these nodes form
                                                        ! the vertices of the triangular cells
!___________________________________________________________________________________________________________________________________
real(ep), dimension(ncells) :: ustarV                   ! u* (m/s) on each cell used for tesslation of human form
real(ep), dimension(ncells) :: Vf                       ! face velocity Vf (m/s) on each cell used for tesselation of (protected) human form
!___________________________________________________________________________________________________________________________________
real(ep), dimension(ncells) :: cell_area                ! area of triangular cells used in tesselation of human form (m^2)
real(ep), dimension(ncells,nBins) :: dep_prob           ! deposition probability (associated with each triangular cell)
                                                        ! as function of particle size (binned into nBins)
!___________________________________________________________________________________________________________________________________
real(ep), dimension(ncells,nBins) :: pen_prob           ! pentration probability through protective clothing ensemble as
                                                        ! function of location on body surface (cell) and particle size (binned into nBins)
!___________________________________________________________________________________________________________________________________

real(ep), allocatable, dimension(:,:) :: massFrac       ! mass fraction for each particle size bin at every time step
                                                        ! massFrac(i,j) - mass fraction in bin i at time step j
                                                        ! i \in [1, nBins] and j \in [1, nTimeSteps]

real(ep), allocatable, dimension(:) :: t                ! time (s) from initial release
real(ep), allocatable, dimension(:) :: C                ! concentration time history (C/Q with units 1/m^3)

real(ep), allocatable, dimension(:) :: cumMass          ! cumulative mass deposited on the body normalized by Q [-]

real(ep), dimension(ncells) :: dep_den                  ! deposition density normalized by Q (1/m^2)

real(ep), dimension(nBodyParts) :: massBP               ! mass deposited on body part normalized by Q [-]

real(ep), dimension(nBins+1) :: d                       ! bin boundaries for particle sizes (there are nBins class intervals)
                                                        ! d(i) and d(i+1) bracket class interval i (i = 1,2, ..., nBins)

real(ep), allocatable, dimension(:,:) :: massFracMax    ! mass fraction distribution corresponding to maximum liquid concentration at receptor location
                                                        ! massFracMax(i,j) : i - bin index, j = receptor index
real(ep), allocatable, dimension(:) :: LiqConcMax       ! maximum liquid concentration (per unit mass Q released) at receptor location

! Input file names

!___________________________________________________________________________________________________________________________________________
character(len=60) :: fnameGeom                              ! input file name containing information of grid geometry over human body
character(len=60) :: fnameUstar                             ! input file name containing distribution of u* (shear velocity) over the human body
character(len=60) :: fnameVf                                ! input file name containing distribution of Vf (face velocity) over human form
                                                            ! normal velocity to cell: (+ve value of Vf implies velocity points AWAY from body,
                                                            ! whereas -ve value of Vf implies velocity points TOWARDS body) - unit outward normal
                                                            ! to each cell contained in file whose name is given by fnameGeom
!____________________________________________________________________________________________________________________________________________
character(len=60), allocatable, dimension(:) :: fnameCL     ! input file containing challenge level (concentration
                                                            ! time history)

! Output file names

character(len=60), allocatable, dimension(:) :: fnameTecPlot     ! name of TecPlot file to output results (deposition density)
character(len=60), allocatable, dimension(:) :: fnameOut         ! name of ASCII file to output results (deposition density)
character(len=60), allocatable, dimension(:) :: fnameOutMass     ! name of ASCII file to output results (mass on each body part)
character(len=60), allocatable, dimension(:) :: fnameOutMassHistory   ! name of ASCII file to output results (mass deposited on entire body as function of time)

character(len=60) :: fnameMassDistributionMax                    ! name of file containing mass fraction distribution
                                                                 !   associated with the time of maximum liquid aerosol concentration at receptors
character(len=60) :: fnameMaxLiquidConcentration                 ! name of file containing maximum liquid aerosol concentration
                                                                 !   (normalized by release mass Q) at receptors

save

end module para


program challengeG
!
! Program to determine deposition densities on the human form (in normal and
! protected postures) for a given challenge level (concentration time history)
!
! Input file required: WhatToProcessV4.txt
! Refer to comments in subroutine ioFileNamesV2 for the required contents of the
! input file
!
use para

implicit none

!_____________________________________________________________________________________


! Get names of all input files and output files

call ioFileNamesV2

! Set up bin boundaries for particle size bins

call bb

! Read in information contained in common input files
! This information is common for processing of challenge level
! to determine deposition densities at every receptor location.
! After relevant information is read in, this information is
! used to compute deposition probability for each particle size bin.
! In addition, the face velocity for each patch (cell) on human form
! is read in.

call readInputCommonV2

! Specify protection afforded by protective clothing ensemble (if
! it exists)

call peneV2

! For each receptor location compute the deposition density on human form

do kIndx=1,nFiles

! Read in challenge level information

   call readInput

! Determine the mass distribution associated with the maximum liquid
! aerosol concentration at receptor location, as well as the maximum
! liquid aerosol concentration

   call extractMassDistribution

! Compute deposition density on the surface of human form

   call ddV2

! Compute the mass on each body part comprising human form

   call mass

! Write output files containing deposition density

   call writeOutput

end do     ! next receptor location

! Write out information pertaining to mass distribution associated with
! maximum liquid aerosol concentration at the various receptor points,
! as well as the value of the maximum liquid aerosol concentration at
! these receptor points

!call writeMassDistribution

stop

end program challengeG


!______________________________________________________________________________________________________________________________

subroutine ioFileNamesV2
!
! Set names of files required as input/output
!
use para
implicit none

integer :: i
integer :: condition

!------------------------------------------------------------------
! Open file to read in necessary information
! Required (ASCII) input file name: WhatToProcessV4.txt

write(*,*) ' Reading file WhatToProcessV4.txt'

open(unit=21,file='WhatToProcessV4.txt',status='old',action='read',form='formatted',iostat=condition)

! Read in some parameter values

read(21,*) iProtect                     ! indicator (iProtect = 0, no protection; iProtect /= 0, protection as specified in
                                        ! routine pene
                                        ! For protected human in IPE, the material systems available are as follows:
                                        !
                                        ! 1. iProtect = 1, Combat uniform
                                        ! 2. iProtect = 2, Horizon 1
                                        ! 3. iProtect = 3, CPCU MkIII
                                        ! 4. iProtect = 4, CPCU MkIV
read(21,*) iWrite                       ! indicator (iWrite = 0, do not write out deposition density files (ASCII and TecPlot);
                                        ! iWrite /= 0, write out deposition density files (ASCII and TecPlot)
read(21,*) nTimeSteps                   ! # time steps in challenge level file (containing concentration-time history)
read(21,*) dt                           ! time increment of each time step(s)
read(21,*) U                            ! advection velocity (m/s) - used for calculation of mass flux
read(21,*) Umodel                       ! model advection velocity (m/s) used in calculation of u* on human body surface
                                        !  (this velocity is used in calculation of the deposition probability)
read(21,*) nFiles                       ! # files to process ( = # receptors downwind of release)

! Allocate some arrays

allocate( t(nTimeSteps) )
allocate( C(nTimeSteps) )
allocate( cumMass(nTimeSteps) )

allocate( massFrac(nBins,nTimeSteps) )

allocate( massFracMax(nBins,nFiles) )
allocate( LiqConcMax(nFiles) )

!---------------------------------------------
! Input files

! Read in name of input file containing geometry information on human form

read(21,*) fnameGeom

! Read in name of input file containing information on the distribution of u* on human form

read(21,*) fnameUstar

! Read in name of input file containing information on the distribution of face velocity Vf on human form

read(21,*) fnameVf

! Read in names of input files containing challenge level (concentration time history) for
! each receptor location

allocate( fnameCL(nFiles) )

do i=1,nFiles
   read(21,*) fnameCL(i)
end do

!----------------------------------------------------------
! Output files
! NB: Do not include extensions to the file name for all output files
! These extensions will be added in depending on whether the files correspond
! to unprotected or protected personnel


! Read in names of output files that will contain deposition density for TecPlot plotting
! (one file associated with each input file containing challenge level)

allocate( fnameTecPlot(nFiles) )

do i=1,nFiles
   read(21,*) fnameTecPlot(i)
end do

! Read in names of output files that will contain deposition density
! (one file associated with each input file containing challenge level)

allocate( fnameOut(nFiles) )

do i=1,nFiles
   read(21,*) fnameOut(i)
end do

! Read in names of output files that will contain the total mass deposited
! on body surface as function of time
! (one file associated with each input file containing challenge level)

allocate( fnameOutMassHistory(nFiles) )

do i=1,nFiles
   read(21,*) fnameOutMassHistory(i)
end do

! Read in names of output files that will contain mass deposited on each body part
! (one file associated with each input file containing challenge level)

allocate( fnameOutMass(nFiles) )

do i=1,nFiles
   read(21,*) fnameOutMass(i)
end do

! Read in names of output files that will contain the mass distribution associated
! with the time of maximum liquid aerosol concentration (normalized by release mass Q)
! and that will contain the maximum liquid aerosol concentration at each of the
! receptor locations

read(21,*) fnameMassDistributionMax
read(21,*) fnameMaxLiquidConcentration

close(unit=21)
write(*,*) '    .... done'


! Append the proper extension to the output file names depending on whether
! one is dealing with unprotected (-U) or protected (-P) personnel

select case (iProtect)

   case (0)      ! unprotected personnel
      do i=1,nFiles
         fnameTecPlot(i) = trim(fnameTecPlot(i))//'-U.plt'
         fnameOut(i) = trim(fnameOut(i))//'-U.dat'
         fnameOutMassHistory(i) = trim(fnameOutMassHistory(i))//'-U.dat'
         fnameOutMass(i) = trim(fnameOutMass(i))//'-U.dat'
      end do
   case (1)      ! Combat uniform
      do i=1,nFiles
         fnameTecPlot(i) = trim(fnameTecPlot(i))//'-P1.plt'
         fnameOut(i) = trim(fnameOut(i))//'-P1.dat'
         fnameOutMassHistory(i) = trim(fnameOutMassHistory(i))//'-P1.dat'
         fnameOutMass(i) = trim(fnameOutMass(i))//'-P1.dat'
      end do
   case (2)      ! Horizon 1
      do i=1,nFiles
         fnameTecPlot(i) = trim(fnameTecPlot(i))//'-P2.plt'
         fnameOut(i) = trim(fnameOut(i))//'-P2.dat'
         fnameOutMassHistory(i) = trim(fnameOutMassHistory(i))//'-P2.dat'
         fnameOutMass(i) = trim(fnameOutMass(i))//'-P2.dat'
      end do
   case (3)      ! CPCU MkIII
      do i=1,nFiles
         fnameTecPlot(i) = trim(fnameTecPlot(i))//'-P3.plt'
         fnameOut(i) = trim(fnameOut(i))//'-P3.dat'
         fnameOutMassHistory(i) = trim(fnameOutMassHistory(i))//'-P3.dat'
         fnameOutMass(i) = trim(fnameOutMass(i))//'-P3.dat'
      end do
   case (4)     ! CPCU MkIV
      do i=1,nFiles
         fnameTecPlot(i) = trim(fnameTecPlot(i))//'-P4.plt'
         fnameOut(i) = trim(fnameOut(i))//'-P4.dat'
         fnameOutMassHistory(i) = trim(fnameOutMassHistory(i))//'-P4.dat'
         fnameOutMass(i) = trim(fnameOutMass(i))//'-P4.dat'
      end do

end select

return

end subroutine ioFileNamesV2


subroutine bb
!
! Set up bin boundaries for particle size bins
! ** This subroutine should be changed if different bin boundaries
!     are used for discretization of droplet size distribution
!
use para
implicit none

integer :: i

do i=0,nBins

   d(i+1) = -0.9 + 0.2*real(i,ep)    ! log10(d_i)
   d(i+1) = 10**(d(i+1))             ! d (um)

end do

return

end subroutine bb


!_______________________________________________________________________________________________________________
subroutine readInputCommonV2
!
! Read in all input data from various files. Determine the deposition probability.
!
use para
implicit none

integer :: condition
integer :: i, j

real(ep) :: nx, ny, nz

real(ep) :: vdPlus         ! Declare function (for computation of normalized deposition velocity)

real(ep) :: dp
real(ep) :: depVel
real(ep) :: us

character(len=128) :: header

real(ep) :: bodyPartArea(11), totCellArea

real(ep), parameter :: check = 0              ! turn on (1) or off (0) check on area of body part

! Read in geometry and connectivity information used for
! tesselation of the body surface of the human form; also,
! read in information on cell areas and body part index associated with each
! cell (element) used in representation of huamn form.
! Let nv = nvertices, nc = ncells.
! File format:
! xv(i) : i=1,nv  , x-coordinates of nodes
! yv(i) : i=1,nv  , y-coordinates of nodes
! zv(i) : i=1,nv  , z-coordinates of nodes
! nx(j) : j=1,nc  , x-component of unit outward normal to element
! ny(j) : j=1,nc  , y-component of unit outward normal to element
! nz(j) : j=1,nc  , z-component of unit outward normal to element
! cell_area(j) : j=1,nc     , cell area of element (m^2)
! cell_index(j) : j=1,nc    , body part index (1:11)
! vert(j,1), vert(j,2), vert(j,3) : j=1,nc    , connectivity of nodes used to form
!                                               triangular elements used for representation
!                                               of human form

open(unit=21,file=trim(fnameGeom),status='old',action='read',form='formatted',iostat=condition)

write(*,*) ' Reading file ',trim(fnameGeom)

! Skip 10 header lines

do i=1,10
   read(21,*) header
end do

!    ... read in x-, y-, and z-coordinates of nodes used to represent human form
do i=1,nvertices
   read(21,*) xv(i)
end do
do i=1,nvertices
   read(21,*) yv(i)
end do
do i=1,nvertices
   read(21,*) zv(i)
end do
!   ... read in x-, y-, z-components of unit outward normal to elements of human form
do j=1,ncells
   read(21,*) nx
end do
do j=1,ncells
   read(21,*) ny
end do
do j=1,ncells
   read(21,*) nz
end do
!   ... read in cell areas (m^2), viz. the area of elements used to represent human form
do j=1,ncells
   read(21,*) cell_area(j)
end do
!   ... read in body part index associated with each element used to represent human form
do j=1,ncells
   read(21,*) cell_index(j)
end do
!   ... read in element connectivity information
do j=1,ncells
   read(21,*) vert(j,1),vert(j,2),vert(j,3)
end do

close(unit=21)
write(*,*) '    .... done'

! As a check, calculate body part area as well as total body area
if (check == 1) then

   bodyPartArea = 0.0

   do j=1,ncells
      bodyPartArea( cell_index(j) ) = bodyPartArea( cell_index(j) ) + cell_area(j)
   end do

   totCellArea = sum( bodyPartArea(1:nBodyParts) )

   write(*,*) 'Body part index           Area (m^2)'
   do i=1,nBodyParts
      write(*,*) '       ',i,'               ',bodyPartArea(i)
   end do
   write(*,*)
   write(*,*) 'Area of body surface (m^2): ',totCellArea
   write(*,*)

end if

! Read in distribution of u* on the human body surface
! File format:
!    nc rows - u* values for each cell or patch on human form (m/s)

open(unit=21,file=trim(fnameUstar),status='old',action='read',form='formatted',iostat=condition)

write(*,*) ' Reading file ',trim(fnameUstar)

! ... read in u* values (m/s)
do j=1,ncells
   read(21,*) ustarV(j)
end do

close(unit=21)
write(*,*) '    .... done'


! Compute the deposition probability as a function of particle size

do i=1,nBins                   ! for each particle size bin
   dp = 10**( 0.5*( log10( d(i) ) + log10( d(i+1) ) ) )  ! particle diameter in um (midpoint of bin boundaries on common
                                                                ! logarithmic scale)
   dp = dp*1.0E-4                                        ! particle diameter in cm
   do j=1,ncells

      us = ustarV(j)*100.0                               ! convert u* in m/s to cm/s
      depVel = vdPlus(dp,us)                             ! normalized deposition velocity vd+
      depVel = depVel*us                                 ! unnormalize vd+ (gives vd in cm/s)

      dep_prob(j,i) = depVel/(Umodel*100.0)              ! normalize vd by Umodel (in cm/s) to get deposition probability

   end do
end do

!
! Read in distribution of Vf (face velocity) on the protected human body surface
! File format:
!    nc rows - Vf values for each cell or patch on human form (m/s)

open(unit=21,file=trim(fnameVf),status='old',action='read',form='formatted',iostat=condition)

write(*,*) ' Reading file ',trim(fnameVf)

! ... read in Vf values (m/s) for each element
do j=1,ncells
   read(21,*) Vf(j)
end do

close(unit=21)
write(*,*) '    .... done'

return

end subroutine readInputCommonV2
!___________________________________________________________________________________________________________________________

subroutine peneV2
!
! Specify the penetration probability through protective clothing
! ensemble as function of face velocity and particle size bins
! Note: setting pen_prob = 1.0 implies there is no protection (viz.,
! the computation is for unprotected individuals)
!
! For unprotected human, use iProtect = 0
!
! For protected human in IPE, the material systems available are as follows:
!
! 1. iProtect = 1, Combat uniform
! 2. iProtect = 2, Horizon 1
! 3. iProtect = 3, CPCU MkIII
! 4. iProtect = 4, CPCU MkIV
!
! The models for penetration probability as function of face velocity and particle
! size have been obtained using classical filtration theory. The parameters of the
! model have been obtained through nonlinear least squares fitting of limited penetration
! data for material system (summarized above) obtained through swatch-level testing.
!
use para
implicit none

integer :: i, j

real(ep) :: Pe        ! Declare function (Peclet number)
real(ep) :: Stk       ! Declare function (Stokes number)

real(ep) :: dp

real(ep) :: aPar,bPar,cPar,dPar         ! parameters that define penetration model for each material system
real(ep) :: dfPar, tPar                 ! effective fibre diameter (in cm) and fabric thickness (in cm)

real(ep) :: ymod

real(ep) :: velF                        ! face velocity (cm/s)

select case (iProtect)
    case (0)                    ! unprotected
       pen_prob = 1.0
    case (1)                    ! combat uniform
!
! Define values for parameters that determine penetration for material system
       tPar = 0.04              ! thickness of material system (in cm)
       dfPar = 1.0e-4           ! effective fibre diameter of material system (in cm)
       aPar = 0.101359
       bPar = -0.519135
       cPar = 0.000444615
       dPar = 1.05104

       do i=1,nBins

          dp = 10**( 0.5*( log10( d(i) ) + log10( d(i+1) ) ) )  ! particle diameter in um (midpoint of bin boundaries on common
                                                                ! logarithmic scale)
          dp = dp*1.0E-4                                        ! particle diameter in cm

          do j=1,ncells

            velF = Vf(j)*100.0                     ! face velocity for cell in cm/s
!            velF = MIN(0.0,cvelF)
            if (velF < 0.0) then                   ! penetration if face velocity is towards the body

               ymod = aPar*(Stk(abs(velF),dp,dfPar))**bPar + cPar*(Pe(abs(velF),dp,dfPar))**dPar       ! [ymod] = cm^(-1)

               pen_prob(j,i) = exp( -ymod*tPar )

            else                                   ! otherwise, no penetration

               pen_prob(j,i) = 0.0

            end if

          end do

       end do
    case (2)                   ! Horizon 1
!
! Define values for parameters that determine penetration for material system
       tPar = 0.13              ! thickness of material system (in cm)
       dfPar = 1.0e-4           ! effective fibre diameter of material system (in cm)
       aPar = 4.52095
       bPar = 0.792578
       cPar = 9.64952
       dPar = -0.384131

       do i=1,nBins

          dp = 10**( 0.5*( log10( d(i) ) + log10( d(i+1) ) ) )  ! particle diameter in um (midpoint of bin boundaries on common
                                                                ! logarithmic scale)
          dp = dp*1.0E-4                                        ! particle diameter in cm

          do j=1,ncells

            velF = Vf(j)*100.0                     ! face velocity for cell in cm/s
!            velF = MIN(0.0,velF)
            if (velF < 0.0) then                   ! penetration if face velocity is towards the body

               ymod = aPar*(Stk(abs(velF),dp,dfPar))**bPar + cPar*(Pe(abs(velF),dp,dfPar))**dPar       ! [ymod] = cm^(-1)

               pen_prob(j,i) = exp( -ymod*tPar )

            else

              pen_prob(j,i) = 0.0                  ! otherwise, no penetration

            end if

          end do

       end do
    case (3)                    ! CPCU MkIII
!
! Define values for parameters that determine penetration for material system
       tPar = 0.122             ! thickness of material system (in cm)
       dfPar = 0.1e-4           ! effective fibre diameter of material system (in cm)
       aPar = 0.413873
       bPar = -0.559138
       cPar = 14.472
       dPar = 0.274692

       do i=1,nBins

          dp = 10**( 0.5*( log10( d(i) ) + log10( d(i+1) ) ) )  ! particle diameter in um (midpoint of bin boundaries on common
                                                                ! logarithmic scale)
          dp = dp*1.0E-4                                        ! particle diameter in cm

          do j=1,ncells

            velF = Vf(j)*100.0                     ! face velocity for cell in cm/s
!            velF = MIN(0.0,velF)
            if (velF < 0.0) then                   ! penetration if face velocity is towards the body

               ymod = aPar*(Stk(abs(velF),dp,dfPar))**bPar + cPar*(Pe(abs(velF),dp,dfPar))**dPar      ! [ymod] = cm^(-1)

               pen_prob(j,i) = exp( -ymod*tPar )

            else

               pen_prob(j,i) = 0.0                ! otherwise, no penetration

            end if

          end do

       end do
     case (4)                    ! CPCU MkIV
!
! Define values for parameters that determine penetration for material system
       tPar = 0.133             ! thickness of material system (in cm)
       dfPar = 0.25e-4          ! effective fibre diameter of material system (in cm)
       aPar = 67.3268
       bPar = 0.229111
       cPar = 10.278
       dPar = -0.894001

       do i=1,nBins

          dp = 10**( 0.5*( log10( d(i) ) + log10( d(i+1) ) ) )  ! particle diameter in um (midpoint of bin boundaries on common
                                                                ! logarithmic scale)
          dp = dp*1.0E-4                                        ! particle diameter in cm

          do j=1,ncells

            velF = Vf(j)*100.0                     ! face velocity for cell in cm/s
!            velF = MIN(0.0,velF)
            if (velF < 0.0) then                   ! penetration if face velocity is towards the body

               ymod = aPar*(Stk(abs(velF),dp,dfPar))**bPar + cPar*(Pe(abs(velF),dp,dfPar))**dPar     ! [ymod] = cm^(-1)

               pen_prob(j,i) = exp( -ymod*tPar )

            else

               pen_prob(j,i) = 0.0                 ! otherwise, no penetration

            end if

          end do

       end do

end select

!do i=1,nBins
!   dp = 10**( 0.5*( log10( d(i) ) + log10( d(i+1) ) ) )
!   write(*,*) dp, pen_prob(i)
!end do

!stop

return

end subroutine peneV2
!__________________________________________________________________________________________________________________________________


subroutine readInput
!
! Read in challenge level information for particular receptor
! indexed by kIndx
!
use para
implicit none

integer :: condition
integer :: i, j, k

integer :: ntme

real(ep), dimension(3) :: concentration

! Read in challenge level information
! File format:
! first row of file gives the number of time points
! each successive row contains 17 columns:
!   - col 1 = time from start of release (s)
!   - col 2 - 14 = mass fraction of total liquid aerosol concentration in particle size bins (implied 13 bins)
!   - col 15 = sum of liquid aerosol concentration obtained by summing concentration over all particle size bins
!   - col 16 = total liquid aerosol concentration (per unit release mass Q) [C/Q with units 1/m^3]
!   - col 17 = total vapor concentration (per unit release mass Q) [C/Q with units 1/m^3]

open(unit=21,file=trim(fnameCL(kIndx)),status='old',action='read',form='formatted',iostat=condition)

write(*,*)' Reading file ',trim(fnameCL(kIndx)),' corresponding to receptor #: ',kIndx

read(21,*) ntme     ! # time points (steps) in file

do i=1,nTimeSteps

   read(21,*) t(i), (massFrac(j,i),j=1,13), (concentration(k),k=1,3)
   C(i) = concentration(1)                       ! extract total liquid aerosol concentration

end do

close(unit=21)

write(*,*) '    .... done'

return

end subroutine readInput


!____________________________________________________________________________________________________________________________________
subroutine ddV2
!
! Compute the deposition density on surface of the human form
! Also compute the cumulative total mass deposited on the body
!
use para
implicit none

integer :: i, j, k

real(ep) :: sumI

dep_den = 0.0
cumMass = 0.0

do k=1,ncells    ! for each triangular element forming human body surface

   do j=1,nTimeSteps      ! for each time step associated with the challenge

      if( C(j) < 1.0e-30 ) then   ! skip time step if concentration is effectively zero
         cumMass(j) = cumMass(j-1)
         cycle
      end if

      sumI = 0.0

      do i=1,nBins           ! for each particle size bin

         sumI = sumI + pAdhesion*dep_prob(k,i)*pen_prob(k,i)*massFrac(i,j)

      end do  ! next bin

      dep_den(k) = dep_den(k) + U*C(j)*dt*sumI

      cumMass(j) = cumMass(j) + dep_den(k)*cell_area(k)

   end do  ! next time step

end do  ! next triangular cell

return

end subroutine ddV2
!__________________________________________________________________________________________________________________________________________


subroutine mass
!
! Compute mass of material deposited on each body part
!
use para
implicit none

integer :: i

massBP = 0.0

do i=1,ncells

   massBP( cell_index(i) ) = massBP( cell_index(i) ) + dep_den(i)*cell_area(i)

end do

totMass = sum( massBP(1:nBodyParts) )

return

end subroutine mass


subroutine extractMassDistribution
!
! Extract the mass distribution corresponding to the maximum of the liquid aerosol
! concentration (normalized by Q) at a receptor location (indexed by kIndx)
!
use para
implicit none

! Local variables

integer :: iMax
integer :: i

!------------------

iMax = maxloc(C,dim=1)     ! array index of maximum liquid aerosol concentration
                           ! iMax = time step associated with maximum liquid aerosol concentration

LiqConcMax(kIndx) = C(iMax)

do i=1,nBins
   massFracMax(i,kIndx) = massFrac(i,iMax)
end do

return

end subroutine extractMassDistribution

!____________________________________________________________________________________________________________________

real(ep) function vdPlus(dp,ustar)
!
! This function computes the normalized deposition velocity [-] (defined as the deposition velocity normalized by the
! shear (friction) velocity. The inputs to the routine are the particle diameter dp (in cm) and the shear (friction)
! velocity of the surface (in cm/s)
!
! NB: This form of normalized deposition velocity is appropriate for vertical surface with gravity in the flow direction.
! One needs to add the term ut for horizontal surface with gravity perpendicular to the flow.
!
use para
implicit none

! Declare calling parameters
real(ep), intent(in) :: dp
real(ep), intent(in) :: ustar

! Declare local parameters
real(ep) :: Sc        ! Declare function
real(ep) :: tauPlus   ! Declare function
real(ep) :: ut        ! Declare function

! Evaluate normalized deposition velocity

vdPlus = 0.057*Sc(dp)**(-2./3.) + 4.5E-4*tauPlus(dp,ustar)**2 + ut(dp,ustar)

end function vdPlus


real(ep) function Cc(dp)
!
! This function computes the Cunningham slip correction factor. The input to the routine is the particle
! diameter dp in cm.
!
use para
implicit none

! Declare calling parameters
real(ep), intent(in) :: dp

! Evaluate Cunningham slip correction factor

Cc = 1.0 + 1.246*2.0*lambda/dp + 0.42*2.0*(lambda/dp)*exp(-0.87*dp/(2.0*lambda) )

end function Cc

real(ep) function diffCoef(dp)
!
! This function computes the diffusion coefficient (diffusivity) of particle in cm^2/s. The input to the routine is
! the particle diameter in cm.
!
use para
implicit none

! Declare calling parameters
real(ep), intent(in) :: dp

! Declare local parameters
real(ep) :: Cc       ! Declare function
real(ep), PARAMETER :: Pi = 3.141592654

! Evaluate diffusivity

diffCoef = (kb*absT*Cc(dp))/(3.0*Pi*eta*dp)

end function diffCoef

real(ep) function Sc(dp)
!
! This function computes the Schmidt number for the particle. The input to the routine is the
! particle diameter in cm.
!
use para
implicit none

! Declare calling parameters
real(ep), intent(in) :: dp

! Declare local parameters
real(ep) :: diffCoef    ! Declare function

! Evaluate Schmidt number

Sc = nu/diffCoef(dp)       ! nu is kinematic viscosity of air in cm^2/s

end function Sc

real(ep) function Pe(velF,dp,df)
!
! This function computes the Peclet number for the (Brownian) diffusion. The inputs to the routine
! are the face velocity velF in cm/s, the particle diameter in cm, and the equivalent (effective) fibre
! diameter for the fabric in cm.
!
use para
implicit none

! Declare calling parameters
real(ep), intent(in) :: velF
real(ep), intent(in) :: dp
real(ep), intent(in) :: df

! Declare local parameters
real(ep) :: diffCoef    ! Declare function

! Evaluate the Peclet number

Pe = velF*df/diffCoef(dp)

end function Pe

real(ep) function Stk(velF,dp,df)
!
! This function evaluates the Stokes number (measure of the importance of inertial impaction in the
! filtration process). The inputs to the routine are the face velocity velF in cm/s, the particle diameter
! in cm, and the equivalent (effective) fibre diameter for the fabric in cm.
!
use para
implicit none

! Declare calling parameters
real(ep), intent(in) :: velF
real(ep), intent(in) :: dp
real(ep), intent(in) :: df

! Declare local parameters
real(ep) :: Cc       ! Declare function

! Evaluate Stokes number

Stk = ( rhop*(dp**2)*velF*Cc(dp) )/( 18.0*eta*df )    ! rhop is particle density in g/cm^3, eta is dynamic viscosity of air (poise)


end function Stk


real(ep) function ut(dp,ustar)
!
! This function computes the contribution to the deposition velocity due to gravitational sedimentation
! on horizontal surface. The input to the routine are the particle diameter dp in cm and the shear (or
! friction) velocity of the surface where deposition is to occur in cm/s.
!
use para
implicit none

! Declare calling parameters
real(ep), intent(in) :: dp
real(ep), intent(in) :: ustar

! Declare local parameters
real(ep) :: tauPlus    ! Declare function

! Evaluate contribution

if( ustar > 0.0 ) then
   ut = tauPlus(dp,ustar)*nu*g/(ustar**3)    ! nu is kinematic viscosity of air (cm^2/s)
                                             ! g is acceleration due to gravity (cm/s^2)
else
   ut = 0.0
end if

end function ut


real(ep) function tauPLus(dp,ustar)
!
! This function computes the non-dimensional particle relaxation time. The inputs to the routine are
! the particle diameter dp in cm and the shear (or friction) velocity of the surface where deposition
! is to occur in cm/s
!
use para
implicit none

! Declare calling parameters
real(ep), intent(in) :: dp
real(ep), intent(in) :: ustar

! Declare local parameters
real(ep) :: Cc    ! Declare function

! Evaluate the non-dimensional particle relaxation time

tauPLus = ((SRatio*(dp**2)*(ustar**2))/(18.0*nu**2))*Cc(dp)

end function tauPlus

!_________________________________________________________________________________________________________



subroutine writeOutput
!
! Write output files containing information on deposition density
!
use para
implicit none

integer :: i, j

!-------------------------------------------------
!  ... write out mass of material deposited on each body part
!  Note - this is normalized by Q, so is dimensionless; alternatively,
!  the value here can be interpreted as mass deposited on human body part
!  per unit of release mass at the source

OPEN(unit=22,file=trim(fnameOutMass(kIndx)),status='REPLACE',action='write',form='FORMATTED')

write(*,*) ' Writing file ',trim(fnameOutMass(kIndx)),' corresponding to receptor #: ',kIndx

write(22,*) ' Body Part i = '
write(22,*) '   1  : head (entire)'
write(22,*) '   2  : neck (entire)'
write(22,*) '   3  : upper chest (including shoulders and down to just below the nipples)'
write(22,*) '   4  : lower abdomen (including stomach and down to the waist)'
write(22,*) '   5  : arms (both left and right)'
write(22,*) '   6  : hands (left and right)'
write(22,*) '   7  : legs (left and right, front and back)'
write(22,*) '   8  : feet (left and right)'
write(22,*) '   9  : back'
write(22,*) '  10  : buttocks'
write(22,*) '  11  : crotch area (frontal area below the waist)'
write(22,*)
write(22,*)
write(22,*) ' Body index        Mass (per unit release mass)               Fractional contribution'
write(22,*)

do i=1,nBodyParts

   write(22,*) i,massBP(i),massBP(i)/totMass

end do

write(22,*)
write(22,*)
if (iProtect == 0) then   ! for unprotected human
   write(22,*) 'Mass deposited on head, neck, and hands (per unit release mass): ',massBP(1)+massBP(2)+massBP(6)
   write(22,*)
end if
write(22,*) '--------------------------------------------------------------------'
write(22,*)
write(22,*) 'Total mass deposited (per unit release mass): ',totMass

close(unit=22)

write(*,*) '    .... done'

!__________________________________________________

!--------------------------------------------------
!  ... write out total mass deposited on the human body as a function of time (since release of material
!      from source

OPEN(unit=22,file=trim(fnameOutMassHistory(kIndx)),status='REPLACE',action='write',form='FORMATTED')

write(*,*) ' Writing file ',trim(fnameOutMassHistory(kIndx)),' corresponding to receptor #: ',kIndx

write(22,*) ' Time [s]         Mass (per unit release mass)'

do j=1,nTimeSteps

   write(22,*) real(j,ep)*dt,cumMass(j)

end do

close(unit=22)

write(*,*) '    .... done'

!_________________________________________________

if (iWrite /= 0) then   ! only write out these very large files if requested

!-------------------------------------------------
!   ... write out deposition density data
!   Note - dep_den has units of 1/m^2 (viz., it is normalized by Q)

   OPEN(unit=22,file=trim(fnameOut(kIndx)),status='REPLACE',action='write',form='FORMATTED')

   write(*,*) ' Writing file ',trim(fnameOut(kIndx)),' corresponding to receptor #: ',kIndx

   write(22,*) ncells
   do i=1,ncells
      write(22,*) i, cell_area(i), cell_index(i), dep_den(i)
   end do

   close(unit=22)

   write(*,*) '    .... done'

!__________________________________________________

!--------------------------------------------------
!   ... write out TecPlot file containing node coordinates and connectivity information
!   for tessellation of the surface of the human form, as well as deposition density for each triangular
!   element

   OPEN(unit=22,file=trim(fnameTecPlot(kIndx)),status='REPLACE',action='write',form='FORMATTED')

   write(*,*) ' Writing file ',trim(fnameTecPlot(kIndx)),' corresponding to receptor #: ',kIndx

   write(22,*) 'TITLE = "Deposition density"'
   write(22,*) 'VARIABLES = "X"'
   write(22,*) '"Y"'
   write(22,*) '"Z"'
   write(22,*) '"Deposition Density"'
   if (iWrite == 10) then
    write(22,"(A108)") 'ZONE N=89360, E=178716, DATAPACKING=BLOCK, ZONETYPE=FETRIANGLE, VARLOCATION=([1-3]=NODAL, [4]=CellCentered)'
   else
    write(22,"(A118)") 'ZONE NODES=89360, ELEMENTS=178716, DATAPACKING=BLOCK, ZONETYPE=FEtriangle,VARLOCATION=([1-3]=NODAL, [4]=CellCentered)'
   end if

   do i=1,nvertices
      write(22,*) xv(i)
   end do
   do i=1,nvertices
      write(22,*) yv(i)
   end do
   do i=1,nvertices
     write(22,*) zv(i)
   end do
   do i=1,ncells
      if(dep_den(i) <= 0.0) then
         write(22,*) -32.0
      else
         write(22,*) log10( dep_den(i) )
      end if
   end do
   do i=1,ncells
     write(22,*) vert(i,1), vert(i,2), vert(i,3)
   end do

   close(unit=22)

   write(*,*) '    .... done'

end if

!________________________________________________

return

end subroutine writeOutput


subroutine writeMassDistribution
!
! Write output file containing information on mass distribution
! of the challenge levels at each receptor location. The mass distribution
! corresponds to that time at each receptor location associated with the
! maximum liquid aerosol concentration (normalized by release mass Q).
! In addition, another output file is written containing the maximum
! liquid aerosol concentration (normalized by release mass Q) at each
! receptor.
!
use para
implicit none

! Local variables

integer :: i,j

real(ep) :: dp

!__________________________________________________

!--------------------------------------------------
!  ... write out mass distribution at each receptor location
!      each row contains: particle diameter (um), mass fraction (receptor 1),
!      mass fraction (receptor 2), ... , mass fraction (receptor nFiles)
!      so, first column contain particle size for bin
!      each successive column contains the mass distribution for a receptor location

OPEN(unit=22,file=trim(fnameMassDistributionMax),status='REPLACE',action='write',form='FORMATTED')

write(*,*) ' Writing file ',trim(fnameMassDistributionMax)

do i=1,nBins

   dp = 10**( 0.5*( log10( d(i) ) + log10( d(i+1) ) ) )  ! particle diameter in um (midpoint of bin boundaries on common
                                                         ! logarithmic scale)

   write(22,*) dp, (massFracMax(i,j), j=1,nFiles)

end do

close(unit=22)

write(*,*) '    .... done'

!__________________________________________________

!--------------------------------------------------
!  ... write out the maximum liquid aerosol concentration at each receptor location
!      one column of data containing maximum liquid aerosol concentration (normalized by
!      release mass Q) [1/m^3] as function of receptor location


OPEN(unit=22,file=trim(fnameMaxLiquidConcentration),status='REPLACE',action='write',form='FORMATTED')

write(*,*) ' Writing file ',trim(fnameMaxLiquidConcentration)

do j=1,nFiles

   write(22,*) LiqConcMax(j)

end do

close(unit=22)

write(*,*) '    .... done'

return

end subroutine writeMassDistribution
