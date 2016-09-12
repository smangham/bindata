program main
	IMPLICIT NONE

	!Program parameters
	integer, parameter 			:: iFileIn=100,iFileOut=101
	integer, parameter 			:: iKindDP=selected_real_kind(15,300)

	!Physical parameters
	real(iKindDP), parameter 	:: rcPi = DACOS(-1.D0), rSecsToDays=86400.0, rcC=29979245800.0 !cm/s
	real(iKindDP), parameter 	:: rcG = 6.67408e-8 !cgs
	real(iKindDP), parameter 	:: rcMSol = 1.98855e33 !g

	!Internals
	logical				:: bFound, bMessy=.FALSE.
	logical				:: bAllScat=.FALSE., bNoLog=.FALSE., bLineMalformed=.FALSE.
	integer				:: iDimX=100, iDimY=100, iDimR=100
	integer				:: i,j,k,iErr=0, iEOF=0, iDummy, iBinX, iBinY, iBinR, iPhot=0,iPhotR=0
	integer				:: iArg=1, iArg_iter
	integer 			:: iObserver, iObserverMin=0, iObserverMax=0, iObservers=1, iObs
	integer, allocatable			:: aiMap(:,:,:),aiMapX(:,:),aiMapY(:,:),aiMapR(:)
	real(iKindDP), allocatable		:: arMap(:,:,:),arMapX(:,:),arMapY(:,:),arBinX(:),arBinY(:)
	real(iKindDP)		:: rDummy, rLambda, rWeight, rDelay, rPosX,rPosY,rPosZ,rErr
	real(iKindDP)		:: rMinX=-1., rMaxX=-1., rMinY=-1., rMaxY=-1., rRngX=-1., rRngY=-1., rMinR=-1.,rMaxR=-1., rRngR=-1.
	real(iKindDP)		:: rMinP=1e300_iKindDP, rMaxP=0.,rMinI=1e300_iKindDP,rMaxI=0.,rRad=-1, rTemp=0.0
	character(len=512) 	:: cFileIn="",cFileOut="", cDummy, cArg, cTicks='"%g"'
	character(len=512)  :: cBuffer

	!Variables for reading in dump file
	integer				:: iNScat, iLine=0

	!Variables for scatter selection
	integer 			:: iNRScat, iNRScatMin=1, iNRScatMax=999
	integer				:: iNCScat, iNCScatMin=0, iNCScatMax=999

	!Variables for gnuplot
	logical 			:: bChangeCB=.FALSE., bNoKey=.FALSE., bNoTicks=.FALSE.
	real(iKindDP)		:: rMinCB, rMaxCB

	!Reweighting variables
	logical 			:: bReweight=.FALSE., bReweightBinLog=.FALSE., bReweightBinGeom=.FALSE., bLookupY=.TRUE.
	integer				:: iReweightGeom
	real(iKindDP)		:: rReweightPow, rReweightBase
	real(iKindDP), allocatable	:: arMapR(:),arBinR(:),arPosR(:),arReweightMult(:)

	!Line mode variables
	logical 			:: bLineMode=.FALSE., bLineFound=.FALSE., bLineBHEstimate=.FALSE., bLineBHUseRMS=.FALSE.
	integer				:: iLines=0, iNRes
	integer, allocatable 		:: aiLine(:)
	real(iKindDP)		:: rLineLambda, rLineLambdaUpper, rLineLambdaLower, rLineVelUpper, rLineVelLower, rLineVelMax, rLineVelMin
	logical 			:: bLineVel=.FALSE.

	!Error tracking variables
	integer				:: iErrWeight=0, iErrMalf=0, iErrLog=0

	!Pointwise mode variables
	logical 			:: bPointwise=.FALSE.,bPointwiseOnly=.FALSE.
	real(iKindDP)		:: rPathPeak, rPathFWHMlower, rPathFWHMupper, rPeakFlux, rPathCent, rFluxCent, rPathCentU, rPathCentL,rFluxTemp
	integer 			:: iPathPeak
	real(iKindDP)		:: rPosXL, rPosXU

	!Variables for specifying origin
	logical 			:: bOriginMode=.FALSE., bOriginFound=.FALSE.
	integer				:: iOrigins=0, iOrigin
	integer, allocatable		:: aiOrigin(:)

	!Variables for outputting CDF
	logical 			:: bCDF=.FALSE.
	!Centroid delay variables
	real(iKindDP)		:: rDelayCent, rDelayL, rDelayU

	if(command_argument_count().EQ.0)then
		print *,"DESCRIPTION:"
		print *,"The bindata utility is intended to bin up delay dump outputs produced by PYTHON."
		print *,""
		print *,"ARGUMENTS:"
		print *,"	-i FILE"
		print *,"Input file, no suffix. "//&
		"If no -i argument is provided, the first un-flagged argument is taken as the input."
		print *,""
		print *,"	-o FILE"
		print *,"Output file base, no suffix. "//&
		"If no -o argument is provided, the input file is used as a base."
		print *,""
		print *,"	-d N M"
		print *,"Bin on a grid with N by M dimensions. Default 100 by 100."
		print *,""
		print *,"	-r VAL VAL"
		print *,"Minimum & maximum wind radius in file. Default is to find from input file."
		print *,""
		print *,"	-s VAL VAL"
		print *,"Minimum  & maximum number of scatters. Default is 1-999."
		print *,""
		print *,"	-sc VAL VAL"
		print *,"Minimum  & maximum number of continuum scatters. Default is 0-999."
		print *,""
		print *,"	-p VAL VAL"
		print *,"Minimum & maximum path distances to plot. Default is .9-.3.25x wind radius."
		print *,""
		print *,"	-v VAL VAL"
		print *,"Wavelength range to bin. Default is to use spectrum_wavemin-max from input file."
		print *,""
		print *,"	-c VAL VAL"
		print *,"Intensity range for colour plot. Default is to use gnuplot's automatic mode."
		print *,""
		print *,"	-l VAL [VAL] [VAL] [...]"
		print *,"Macro-atom lines to plot. May be arbitrarily long."
		print *,""
		print *,"	-or VAL [VAL] [VAL] [...]"
		print *,"Origins to plot. May be arbitrarily long. Disk = 2, Wind = 3, AGN = 4."		
		print *,""
		print *,"	-rwp VAL"
		print *,"Reweight mode: take data and reweight to r^VAL power law surface brightness."
		print *,""
		print *,"   -rwb VAL"
		print *,"Reweight mode: reweight using VAL bins (default 100)."
		print *,""
		print *,"   -obs VAL"
		print *,"Select points from observer VAL (default 0)."
		print *,""
		print *,"	-lw VAL"
		print *,"When in line mode, if investigating single line for BH mass estimate, specify its wavelength."
		print *,""
		print *,"	-rms"
		print *,"When making virial mass estimate, use line RMS error instead of FWHM."
		print *,""
		print *,"	-bl"
		print *,"Reweight mode: Bin logarithmically, default is linear."
		print *,""
		print *,"	-bg"
		print *,"Reweight mode: Bin geometrically, default is linear."
		print *,""
		print *,"	-m"
		print *,"Messy, do not delete intermediate files."
		print *,""		
		print *,"	-nk"
		print *,"No key, remove colour key from plots."
		print *,""
		print *,"	-nl"
		print *,"No log, plot linear colour scale."
		print *,""
		print *,"	-nt"
		print *,"No ticks, remove tick numbers from plots."
		print *,""
		print *,"	-x"
		print *,"Pointwise mode."
		print *,""
		print *,"	-xo"
		print *,"Pointwise only mode."
		print *,""		
		print *,"	-vel"
		print *,"Plot using velocity, not frequency. Requires -lw."
		STOP
	endif

	do while(iArg.LE.command_argument_count())
		call get_command_argument(iArg, cArg)
		if(cArg.EQ."-d".OR.cArg.EQ."-D")then
			call get_command_argument(iArg+1, cArg)
			read(cArg,*,iostat=iErr)iDimX
			if(iErr.NE.0.OR.iDimX.LT.1)then
				print *,"ERROR: X dimension argument '"//trim(cArg)//"' invalid!"
				STOP
			endif
			call get_command_argument(iArg+2, cArg)
			read(cArg,*,iostat=iErr)iDimY
			if(iErr.NE.0.OR.iDimY.LT.1)then
				print *,"ERROR: Y dimension argument '"//trim(cArg)//"' invalid!"
				STOP
			endif
			iArg=iArg+3

		else if(cArg.EQ."-r".OR.cArg.EQ."-R")then
			call get_command_argument(iArg+1, cArg)
			read(cArg,*,iostat=iErr)rMinR
			if(iErr.NE.0 .OR.rMinR.LT.0)then
				print *,"ERROR: Minimum radius argument '"//trim(cArg)//"' invalid!"
				STOP
			endif
			call get_command_argument(iArg+2, cArg)
			read(cArg,*,iostat=iErr)rMaxR
			if(iErr.NE.0 .OR.rMaxR.LT.rMinR)then
				print *,"ERROR: Maximum radius argument '"//trim(cArg)//"' invalid!"
				STOP
			endif
			iArg=iArg+3
			rMinY=rMaxR*0.90
			rMaxY=rMaxR*3.25
				
		else if(cArg.EQ."-v".OR.cArg.EQ."-V")then
			call get_command_argument(iArg+1, cArg)
			read(cArg,*,iostat=iErr)rMinX
			if(iErr.NE.0 .OR.rMinX.LT.0)then
				print *,"ERROR: Minimum wavelength argument '"//trim(cArg)//"' invalid!"
				STOP
			endif
			call get_command_argument(iArg+2, cArg)
			read(cArg,*,iostat=iErr)rMaxX
			if(iErr.NE.0 .OR.rMaxX.LT.rMinX)then
				print *,"ERROR: Maximum wavelength argument '"//trim(cArg)//"' invalid!"
				STOP
			endif
			iArg=iArg+3

		else if(cArg.EQ."-p".OR.cArg.EQ."-P")then
			call get_command_argument(iArg+1, cArg)
			read(cArg,*,iostat=iErr)rMinY
			if(iErr.NE.0)then
				print *,"ERROR: Minimum path argument '"//trim(cArg)//"' invalid!"
				STOP
			endif
			call get_command_argument(iArg+2, cArg)
			read(cArg,*,iostat=iErr)rMaxY
			if(iErr.NE.0 .OR. rMaxY.LT.rMinY)then
				print *,"ERROR: Maximum path argument '"//trim(cArg)//"' invalid!"
				STOP
			endif
			bLookupY=.FALSE.
			iArg=iArg+3

		else if(cArg.EQ."-c".OR.cArg.EQ."-C")then
			call get_command_argument(iArg+1, cArg)
			read(cArg,*,iostat=iErr)rMinCB
			if(iErr.NE.0)then
				print *,"ERROR: Minimum intensity argument '"//trim(cArg)//"' invalid!"
				STOP
			endif
			call get_command_argument(iArg+2, cArg)
			read(cArg,*,iostat=iErr)rMaxCB
			if(iErr.NE.0 .OR. rMaxY.LT.rMinY)then
				print *,"ERROR: Maximum intensity argument '"//trim(cArg)//"' invalid!"
				STOP
			endif
			bChangeCB=.TRUE.
			iArg=iArg+3

		else if(cArg.EQ."-s".OR.cArg.EQ."-S")then
			call get_command_argument(iArg+1, cArg)
			read(cArg,*,iostat=iErr)iNRScatMin
			if(iErr.NE.0 .OR.iNRScatMin.LT.0)then
				print *,"ERROR: Minimum resonant scatters argument '"//trim(cArg)//"' invalid!"
				STOP
			endif
			call get_command_argument(iArg+2, cArg)
			read(cArg,*,iostat=iErr)iNRScatMax
			if(iErr.NE.0 .OR.iNRScatMax.LT.iNRScatMin)then
				print *,"ERROR: Maximum resonant scatters argument '"//trim(cArg)//"' invalid!"
				STOP
			endif
			iArg=iArg+3

		else if(cArg.EQ."-sc".OR.cArg.EQ."-SC")then
			call get_command_argument(iArg+1, cArg)
			read(cArg,*,iostat=iErr)iNCScatMin
			if(iErr.NE.0 .OR.iNCScatMin.LT.0)then
				print *,"ERROR: Minimum continuum scatters argument '"//trim(cArg)//"' invalid!"
				STOP
			endif
			call get_command_argument(iArg+2, cArg)
			read(cArg,*,iostat=iErr)iNCScatMax
			if(iErr.NE.0 .OR.iNCScatMax.LT.iNCScatMin)then
				print *,"ERROR: Maximum continuum scatters argument '"//trim(cArg)//"' invalid!"
				STOP
			endif
			iArg=iArg+3

		else if(cArg.EQ."-l".OR.cArg.EQ."-L")then
			iArg_iter=iArg+1
			iLines=0
			bLineMode=.TRUE.
			do while(iArg_iter.LE.command_argument_count())
				call get_command_argument(iArg_iter, cArg)
				read(cArg,*,iostat=iErr)iLine
				if(iErr.EQ.0 .AND.iLine.GE.0) then
					iLines = iLines +1
				else
					EXIT
				endif
				iArg_iter = iArg_iter +1
			end do

			if(iLines.EQ.0)then
				print *,"ERROR: No valid lines listed!"
				STOP
			endif
			allocate(aiLine(iLines))
			write(*,'(A)',advance='no')'Line mode- tracking line(s):'
			do iArg_iter=1,iLines
				call get_command_argument(iArg+iArg_iter, cArg)
				read(cArg,*,iostat=iErr)aiLine(iArg_iter)
				write(*,'(X,I0)',advance='no')aiLine(iArg_iter)
			end do
			write(*,*)''
			iArg=iArg+iLines+1

		else if(cArg.EQ."-or".OR.cArg.EQ."-OR")then
			iArg_iter=iArg+1
			iOrigins=0
			bOriginMode=.TRUE.
			do while(iArg_iter.LE.command_argument_count())
				call get_command_argument(iArg_iter, cArg)
				read(cArg,*,iostat=iErr)iOrigin
				if(iErr.EQ.0 .AND.iOrigin.GE.-1) then
					iOrigins = iOrigins +1
				else
					EXIT
				endif
				iArg_iter = iArg_iter +1
			end do

			if(iOrigins.EQ.0)then
				print *,"ERROR: No valid origins listed!"
				STOP
			endif
			allocate(aiOrigin(iOrigins))
			write(*,'(A)',advance='no')'Origin mode- tracking origin(s):'
			do iArg_iter=1,iOrigins
				call get_command_argument(iArg+iArg_iter, cArg)
				read(cArg,*,iostat=iErr)aiOrigin(iArg_iter)
				write(*,'(X,I0)',advance='no')aiOrigin(iArg_iter)
			end do
			write(*,*)''
			iArg=iArg+iOrigins+1			

		else if(cArg.EQ."-rwp".OR.cArg.EQ."-RWP")then
			call get_command_argument(iArg+1, cArg)
			read(cArg,*,iostat=iErr)rReweightPow
			if(iErr.NE.0 )then
				print *,"ERROR: Reweighting power argument '"//trim(cArg)//"' invalid!"
				STOP
			endif
			iArg=iArg+2
			bReweight=.TRUE.
		else if(cArg.EQ."-rwb".OR.cArg.EQ."-RWB")then
			call get_command_argument(iArg+1, cArg)
			read(cArg,*,iostat=iErr)iDimR
			if(iErr.NE.0 .OR.iDimR.LT.0)then
				print *,"ERROR: Reweighting bins argument '"//trim(cArg)//"' invalid!"
				STOP
			endif
			iArg=iArg+2

		else if(cArg.EQ."-obs".OR.cArg.EQ."-OBS")then
			call get_command_argument(iArg+1, cArg)
			read(cArg,*,iostat=iErr)iObserverMin
			if(iErr.NE.0 .OR.iObserverMin.LT.0)then
				print *,"ERROR: Observer argument '"//trim(cArg)//"' invalid!"
				STOP
			endif
			iArg=iArg+2

			call get_command_argument(iArg, cArg)
			read(cArg,*,iostat=iErr)iObserverMax
			if(iErr.NE.0)then
				!If there is no valid upper bound
				iObserverMax = iObserverMin
			else if(iObserverMax.LT.iObserverMin)then
				print *,"ERROR: Observer range must be low to high!"
				STOP					
			else
				iArg=iArg+1				
			endif

		else if(cArg.EQ."-lw".OR.cArg.EQ."-LW")then
			call get_command_argument(iArg+1, cArg)
			read(cArg,*,iostat=iErr)rLineLambda
			if(iErr.NE.0 )then
				print *,"ERROR: Line wavelength argument '"//trim(cArg)//"' invalid!"
				STOP
			endif
			iArg=iArg+2
			bLineBHEstimate=.TRUE.

		else if(cArg.EQ."-bg".OR.cArg.EQ."-BG")then
			iArg=iArg+1
			bReweightBinGeom=.TRUE.
		else if(cArg.EQ."-bl".OR.cArg.EQ."-BL")then
			iArg=iArg+1
			bReweightBinLog=.TRUE.

		else if(cArg.EQ."-m".OR.cArg.EQ."-M")then
			bMessy=.TRUE.
			iArg=iArg+1
		else if(cArg.EQ."-nk".OR.cArg.EQ."-NK")then
			bNoKey=.TRUE.
			iArg=iArg+1
		else if(cArg.EQ."-nl".OR.cArg.EQ."-NL")then
			bNoLog=.TRUE.
			iArg=iArg+1
		else if(cArg.EQ."-nt".OR.cArg.EQ."-NT")then
			bNoTicks=.TRUE.
			cTicks='" "'
			iArg=iArg+1
		else if(cArg.EQ."-x".OR.cArg.EQ."-X")then
			bPointwise=.TRUE.
			iArg=iArg+1
		else if(cArg.EQ."-xo".OR.cArg.EQ."-XO")then
			bPointwise=.TRUE.
			bPointwiseOnly=.TRUE.
			iArg=iArg+1
		else if(cArg.EQ."-vel".OR.cArg.EQ."-VEL")then
			bLineVel=.TRUE.
			iArg=iArg+1

		else if(cArg.EQ."-i".OR.cArg.EQ."-I")then
			call get_command_argument(iArg+1, cFileIn)
			iArg=iArg+2
		else if(cArg.EQ."-o".OR.cArg.EQ."-O")then
			call get_command_argument(iArg+1, cFileOut)
			iArg=iArg+2
		else
			call get_command_argument(iArg, cFileIn)
			iArg=iArg+1			
		endif
	end do

	if(cFileOut.EQ."")cFileOut=cFileIn

	inquire(file=trim(cFileIn)//".delay_dump", exist=bFound)
	if(.NOT.bFound)then
		print *,"ERROR: Input file '"//trim(cFileIn)//".delay_dump' does not exist!"
		STOP
	endif
	if(bLineBHEstimate.AND..NOT.bLineMode)then
		print *,'Trying to estimate BH mass from a line, but no line is specified!'
		STOP
	endif
	if(bLineBHEstimate.AND.iLines.GT.1)then
		print *,'Trying to estimate BH mass from a line, but too many lines specified!'
		STOP
	endif
	if(bLineVel.AND..NOT.bLineBHEstimate)then
		print *,'Trying to plot line velocity without knowing line wavelength!'
		STOP
	endif

	if(bLookupY .AND. rMinR.LT.0 .AND. bReweight)then
		rRad  = rfGetKeyword(trim(cFileIn)//".pf","rstar")							!Units = cm
		rMinR = rfGetKeyword(trim(cFileIn)//".pf","wind_keplerian.diskmin")*rRad 	!Units = rstar
		rMaxR = rfGetKeyword(trim(cFileIn)//".pf","wind_keplerian.diskmax")*rRad 	!Units = rstar
	endif
	if(rMinX.LT.0)then
		rMinX = rfGetKeyword(trim(cFileIn)//".pf","spectrum_wavemin")
		rMaxX = rfGetKeyword(trim(cFileIn)//".pf","spectrum_wavemax")
		print '(X,A,ES8.2,A,ES8.2,A,I0,A)','Binning wavelengths from ',rMinX,' to ',rMaxX,' in ',iDimX,' steps'
	endif
	if(bLookupY)then
		rRad = rfGetKeyword(trim(cFileIn)//".pf","wind.radmax")
		rMinY= rRad*0.50/(rSecsToDays*rcC) !Find delay in days from travel time in light-seconds, then seconds->days
		rMaxY= rRad*4.50/(rSecsToDays*rcC)
		print '(X,A,ES8.2,A,ES8.2,A,I0,A)','Binning paths from ',rMinY,' to ',rMaxY,'cm in ',iDimY,' steps'
	endif
	rRngX=rMaxX-rMinX
	rRngY=rMaxY-rMinY
	rRngR=rMaxR-rMinR

	if(bPointwise)then
		print '(A)','Plotting pointwise mode'
	endif
	iObservers = 1 +(iObserverMax - iObserverMin)
	if(iObservers > 1)then
		print '(A,I0,A,I0)','Plotting observers ',iObserverMin,' to ',iObserverMax
	else
		print '(A,I0)','Plotting observer ',iObserverMin
	endif


!	============================================================================
!	ALLOCATION AND ZEROING SECTION
!	----------------------------------------------------------------------------
	allocate(aiMapX(iDimX,iObserverMin:iObserverMax), aiMapY(iDimY,iObserverMin:iObserverMax))
	allocate(aiMapR(iDimR), aiMap(iDimX,iDimY,iObserverMin:iObserverMax))
	allocate(arMapX(iDimX,iObserverMin:iObserverMax), arMapY(iDimY,iObserverMin:iObserverMax))
	allocate(arMapR(iDimR), arMap(iDimX,iDimY,iObserverMin:iObserverMax))
	allocate(arBinX(iDimX+1), arBinY(iDimY+1))

	do i=0,iDimX
		arBinX(i+1)=rMinX+i*(rMaxX-rMinX)/real(iDimX)
	end do
	do i=0,iDimY
		arBinY(i+1)=rMinY+i*(rMaxY-rMinY)/real(iDimY)
	end do

	aiMap=0
	aiMapX=0
	aiMapY=0

	arMap=0.0
	arMapX=0.0
	arMapY=0.0

!	============================================================================
!	REWEIGHTING SECTION
!	----------------------------------------------------------------------------
	if(bReweight)then
		! Allocate and zero arrays.
		allocate(arBinR(iDimR+1), arPosR(iDimR), arReweightMult(iDimR))
		if(bReweightBinLog)then
			rReweightBase = (log(rMaxR) - log(rMinR)) / iDimR
			arBinR(1) = rMinR
			do i=1,iDimR
				arBinR(i+1)= exp(log(rMinR) + rReweightBase*i)
			end do

		elseif(bReweightBinGeom)then
			iReweightGeom=0
			do i=1,iDimR
				iReweightGeom = iReweightGeom+i
			end do
			rReweightBase = (rMaxR-rMinR) / iReweightGeom

			arBinR(1) = rMinR
			do i=1,iDimR
				arBinR(i+1) = arBinR(i) + rReweightBase*i
			end do

		else
			do i=0,iDimR
				arBinR(i+1)=rMinR+i*(rMaxR-rMinR)/real(iDimR)
			end do	
		endif

		aiMapR=0
		arMapR=0.0
		arReweightMult=1.0

		print '(X,A,ES9.3,A,ES9.3,A,I0,A,ES9.2,A)',&
				'Reweighting range from ',rMinR,' to ',rMaxR,' in ',iDimR,' bins to r^',rReweightPow,'.'
		iErr=0
		iEOF=0
		open(iFileIn,file=trim(cFileIn)//".delay_dump",status="OLD",action="READ")
		read(iFileIn,*) cDummy
		do while(iErr.EQ.0.AND.iEOF.EQ.0)
			read(iFileIn,'(A512)',iostat=iEOF) cBuffer
			if(cBuffer(1:1).NE."#")then
				read(cBuffer,*,iostat=iErr) rDummy, rLambda, rWeight, rPosX, rPosY, rPosZ, &
											iNScat, iNRScat, rDelay, iObserver, iOrigin, iNRes

				if(iErr.GT.0)then
					iErr=0
					iEOF=0
				elseif(iNScat.EQ.1) then
					iBinR = ifLookupIndex(arBinR, sqrt(rPosX**2+rPosY**2) )
					if(iBinR.GT.0)then
						aiMapR(iBinR) = aiMapR(iBinR) + 1
						arMapR(iBinR) = arMapR(iBinR) + rWeight
					endif
				endif
			endif
		end do
		close(iFileIn)

		bReweight=.FALSE.												!Have we successfully reweighted?
		open(iFileOut,file=trim(cFileIn)//".brightness",status="REPLACE",action="WRITE")
		do i=1,iDimR
			arPosR(i) = (arBinR(i)+arBinR(i+1)) / 2.0
			if(aiMapR(i).EQ.0)then
				print '(X,A,ES10.4,A,ES10.4,A)','WARNING: Unsampled radial bin from ',&
				arBinR(i),' to ',arBinR(i+1),'cm.'
				write(iFileOut,'(ES12.5,X)') arPosR(i)
			else
				bReweight=.TRUE.
				arMapR(i) = arMapR(i) / (rcPi * (arBinR(i+1)**2-arBinR(i)**2))	!Area correction
				if(i>1) arReweightMult(i) = arMapR(1)/arMapR(i)					!Build reweight map
				write(iFileOut,'(3(ES12.5,X),I0)') arPosR(i), arMapR(i),arReweightMult(i),aiMapR(i)	!Write to brightness file
			endif
		enddo
		close(iFileOut)
		if(.NOT.bReweight)then
			print '(X,A)','ERROR: No radial bins sampled, could not reweight.'
			STOP
		endif
	endif
!	============================================================================

	print *,"Reading in '"//trim(cFileIn)//".delay_dump'..."
	iPhotR =0
	iErr=0
	iEOF=0
	iLine=0
	open(iFileIn,file=trim(cFileIn)//".delay_dump",status="OLD",action="READ")
	read(iFileIn,*) cDummy
	iLine=iLine+1
	do while(iErr.EQ.0 .AND.iEOF.EQ.0)
		read(iFileIn,'(A512)',iostat=iEOF) cBuffer
		iLine=iLine+1
		
		if(cBuffer(1:1).NE."#")then
			iPhot=iPhot+1
			read(cBuffer,*,iostat=iErr) rDummy, rLambda, rWeight, rPosX, rPosY, rPosZ, &
										iNScat, iNRScat, rDelay, iObserver, iOrigin, iNRes
			iNCScat = iNScat - iNRScat

			!If we're in line mode, check to see if this photon's origin line is in the list of tracked lines
			if(bLineMode)then
				bLineFound=.FALSE.
				do i=1,iLines
					if(iNRes.EQ.aiLine(i)) bLineFound=.TRUE.
				enddo
			endif
			!If we're in origin mode, check to see if this photon's origin is in the list of tracked origins
			if(bOriginMode)then
				bOriginFound=.FALSE.
				do i=1,iOrigins
					if(iOrigin.EQ.aiOrigin(i)) bOriginFound=.TRUE.
				enddo
			endif

			if(iErr.GT.0)then
				iErr=0
				iEOF=0
				iErrMalf=iErrMalf+1
				if(iErrMalf.LT.10)then
					print '(X,A,I0)','WARNING: Malformed entry on line ',iLine
				elseif(iErrMalf.EQ.10)then
					print '(X,A,I0,A)','WARNING: Malformed entry on line ',iLine,'. Suppressing further warnings.'
				endif

			elseif(iObserver.LT.iObserverMin.OR.iObserver.GT.iObserverMax)then
				!Do nothing
			elseif(bLineMode.AND..NOT.bLineFound)then
				!Do nothing
			elseif(bOriginMode.AND..NOT.bOriginFound)then
				!Do nothing
			elseif(iNCScat.LT.iNCScatMin.OR.iNCScat.GT.iNCScatMax)then
				!Do nothing
			elseif(iNRScat.LT.iNRScatMin.OR.iNRScat.GT.iNRScatMax)then
				!Do nothing
			else 

				iPhotR= iPhotR+1
				iBinX = ifLookupIndex(arBinX,rLambda)
				iBinY = ifLookupIndex(arBinY,rDelay/rSecsToDays)
				if(iBinX.LT.1 .OR. iBinY.LT.1)then
					iErrLog=iErrLog+1
					if(iErrLog.LT.10) then
						print '(X,A,ES10.4,A,ES10.4,A,I0)','WARNING: Point at frequency ',rLambda,&
								', path ',rDelay/rSecsToDays,' lies outside range #',iErrLog
					elseif(iErrLog.EQ.10)then
						print '(X,A,ES10.4,A,ES10.4,A,I0,A)','WARNING: Point at frequency ',rLambda,&
								', path ',rDelay/rSecsToDays,' lies outside range #',iErrLog,". Suppressing further warnings."
					endif
				else
					!If using surface brightness correction, then adjust
					if(bReweight)then
						rRad    = sqrt(rPosX**2+rPosY**2)
						iBinR	= ifLookupIndex(arBinR, rRad)
						if(iBinR .GT. -1)then
							rWeight = arReweightMult(iBinR) * rWeight * rRad**(-1.5)
							aiMap(iBinX,iBinY,iObserver) 	= aiMap(iBinX,iBinY,iObserver) + 1
							arMap(iBinX,iBinY,iObserver) 	= arMap(iBinX,iBinY,iObserver) + rWeight
							aiMapX(iBinX,iObserver) 		= aiMapX(iBinX,iObserver) + 1
							arMapX(iBinX,iObserver) 		= arMapX(iBinX,iObserver) + rWeight
							aiMapY(iBinY,iObserver) 		= aiMapY(iBinY,iObserver) + 1
							arMapY(iBinY,iObserver) 		= arMapY(iBinY,iObserver) + rWeight
						else
							iErrWeight=iErrWeight+1
							if(iErrWeight.LT.10)then
								print '(X,A,ES9.2,A,I0)','WARNING: Point at radius ',rRad,' could not be reweighted #',iErrWeight
							elseif(iErrWeight.EQ.10)then
								print '(X,A,ES9.2,A,I0,A)','WARNING: Point at radius ',rRad,' could not be reweighted #',&
														iErrWeight,'. Suppressing further warnings.'
							endif
						endif
					else
						aiMap(iBinX,iBinY,iObserver) = 	aiMap(iBinX,iBinY,iObserver) + 1
						arMap(iBinX,iBinY,iObserver) = 	arMap(iBinX,iBinY,iObserver) + rWeight
						aiMapX(iBinX,iObserver) = 		aiMapX(iBinX,iObserver) + 1
						arMapX(iBinX,iObserver) = 		arMapX(iBinX,iObserver) + rWeight
						aiMapY(iBinY,iObserver) = 		aiMapY(iBinY,iObserver) + 1
						arMapY(iBinY,iObserver) = 		arMapY(iBinY,iObserver) + rWeight
					endif

				endif
			endif
		endif
	end do
	close(iFileIn)

	if(iErrLog.GT.0) print '(X,A,I0,A)',"WARNING: ",iErrLog," points lay outside the map range."
	if(iErrMalf.GT.0)print '(X,A,I0,A)',"WARNING: ",iErrMalf," lines in the file were malformed."

	print '(X,I0,A,I0,A,I0,A)',iPhotR,"/",iPhot," photons scattered."

	if(iPhotR.EQ.0)then
		print *,'ERROR: No photons fit scattering parameters'
		STOP
	endif

	!For each observer, output a plot
	do iObs=iObserverMin,iObserverMax
		!If there is only one observer, don't bother adding the name
		if(iObservers.GT.1)then
			print '(A,I0,A)','Writing to "'//trim(cFileOut)//'.',iObs,'.eps"'
		else
			print '(A)','Writing to "'//trim(cFileOut)//'.eps"'
		endif

		!If we are tracking a single line to find the mass
		if(bLineMode.AND.iLines.EQ.1)then
			call sLineAnalysis(arMapX(:,iObs),arMapY(:,iObs),arBinX,arBinY, rLineLambda)
		endif
		
		if(.not.bPointwiseOnly)then
			open(iFileOut,file=trim(cFileOut)//".bin_XY",status="REPLACE",action="WRITE")
			do j=1,iDimY
				rPosY = rMinY+(j-.5)*(rMaxY-rMinY)/real(iDimY)
				do i=1,iDimX
					rPosX = rMinX+(i-.5)*(rMaxX-rMinX)/real(iDimX)
					if(aiMap(i,j,iObs).GT.0 .AND. arMap(i,j,iObs).GT.0) then 
						rErr = sqrt(REAL(aiMap(i,j,iObs)))/aiMap(i,j,iObs)
						write(iFileOut,'(4(ES12.5,1X))') rPosX, rPosY, arMap(i,j,iObs), rErr
					else
						write(iFileOut,'(4(ES12.5,1X))') rPosX, rPosY, 0.0, 1.0
					endif
				end do
			end do
			close(iFileOut)
		endif
		if(bPointwise)then
			!If we're doing this pointwise, then for each
			open(iFileOut,file=trim(cFileOut)//".bin_XYp",status="REPLACE",action="WRITE")
			do i=1,iDimX
				rPathPeak = rfFindPeak(arMap(i,:,iObs), arBinY, iPathPeak)
				rPathCent = rfFindCentroid(arMap(i,:,iObs), arBinY, rPathCentL, rPathCentU, 0.8*arBinY(iPathPeak))

				if(rPathCent > 0)then
					write(iFileOut,'(6(ES12.5,1X))') (arBinX(i)+arBinX(i+1))/2, arBinX(i), arBinX(i+1),&
													 rPathCent, rPathCentL, rPathCentU
				endif
				
			end do
			close(iFileOut)
		endif
	
		open(iFileOut,file=trim(cFileOut)//".bin_X",status="REPLACE",action="WRITE")
		do i=1,iDimX
			rPosX = rMinX+(i-.5)*(rMaxX-rMinX)/real(iDimX)
			if(aiMapX(i,iObs).GT.0 .AND. arMapX(i,iObs).GT.0) then 
				rErr = sqrt(REAL(aiMapX(i,iObs)))/aiMapX(i,iObs)	
				write(iFileOut,'(3(ES12.5,1X))') rPosX, arMapX(i,iObs), rErr
			else
				write(iFileOut,'(3(ES12.5,1X))') rPosX, 0.0, 1.0
			endif
		end do
		close(iFileOut)

		open(iFileOut,file=trim(cFileOut)//".bin_Y",status="REPLACE",action="WRITE")
		do i=1,iDimY
			rPosY = rMinY+(i-.5)*(rMaxY-rMinY)/real(iDimY)
			if(aiMapY(i,iObs).GT.0 .AND. arMapY(i,iObs).GT.0) then 
				rErr = sqrt(REAL(aiMapY(i,iObs)))/aiMapY(i,iObs)
				write(iFileOut,'(3(ES12.5,1X))') rPosY, arMapY(i,iObs), rErr
			else
				write(iFileOut,'(3(ES12.5,1X))') rPosY, 0.0, 1.0
			endif
		end do
		close(iFileOut)

		open(iFileOut,file=trim(cFileOut)//".plot",status="REPLACE",action="WRITE")
		write(iFileOut,'(A)')'set term postscript eps color enhanced'
		if(iObservers.GT.1)then
			write(iFileOut,'(A,I0,A)')'set output "'//trim(cFileOut)//'.',iObs,'.eps"'
		else
			write(iFileOut,'(A)')'set output "'//trim(cFileOut)//'.eps"'
		endif
		write(iFileOut,'(A)')'set multiplot'
		write(iFileOut,'(A)')'set bmargin 0; set tmargin 0; set lmargin 0; set rmargin 0'
		if(.NOT.bNoLog.AND..NOT.bPointwiseOnly)then
			write(iFileOut,'(A)')'set log cb'
		endif
		write(iFileOut,'(A)')'set colorbox user origin 0.65,0.05 size .05,0.3'
		write(iFileOut,'(A)')'set cblabel "Luminosity (ergs s^{-1})"'
		if(bNoKey.OR.bPointwiseOnly)then
			write(iFileOut,'(A)')'unset colorbox'		
		endif

		write(iFileOut,'(A)')'set origin 0.1,0.4'
		write(iFileOut,'(A)')'set size 0.5,0.5'
		write(iFileOut,'(A)')'set xrange ['//trim(r2c(rMinX))//':'//trim(r2c(rMaxX))//']'
		write(iFileOut,'(A)')'set yrange ['//trim(r2c(rMinY))//':'//trim(r2c(rMaxY))//']'

		if(bLineVel)then
			rLineLambdaUpper = rMaxX - modulo(rMaxX, 100.0)
			rLineLambdaLower = rMinX + (100- modulo(rMinX, 100.0))
			rLineVelMax = 100.0 * rcC * (rMaxX - rLineLambda)/ rLineLambda
			rLineVelMin = 100.0 * rcC * (rLineLambda - rMinX)/ rLineLambda
			rLineVelUpper = rLineVelMax - modulo(rLineVelMax, 100.0)
			if(modulo(rLineVelMax,100.0)>0)then
				rLineVelLower = rLineVelMin + (100 - modulo(rLineVelMax, 100.0))
			else
				rLineVelLower = rLineVelMin
			endif
			rLineLambdaUpper = rLineLambda + (rLineLambda * rLineVelUpper / (rcC*100.0))
			rLineLambdaLower = rLineLambda + (rLineLambda * rLineVelLower / (rcC*100.0))

			write(iFileOut,'(A)')'set xtics ('//trim(r2c(rMinX))//&
							', "'//trim(r2c(rLineVelLower))//'" '//trim(r2c(rLineLambdaLower))//&
							', "0" '//trim(r2c(rLineLambda))//&
							', "'//trim(r2c(rLineVelUpper))//'" '//trim(r2c(rLineLambdaUpper))//&
							', '//trim(r2c(rMaxX))//') mirror format ""'	
		else
			write(iFileOut,'(A)')'set xtics ('//trim(r2c(rMinX))//', '//trim(r2c(rMinX+.25*rRngX))//&
							', '//trim(r2c(rMinX+.5*rRngX))//', '//trim(r2c(rMinX+.75*rRngX))//&
							', '//trim(r2c(rMaxX))//') mirror format ""'
		endif

		write(iFileOut,'(A)')'unset xlabel'
		write(iFileOut,'(A)')'set ytics ('//trim(r2c(rMinY))//', '//trim(r2c(rMinY+.25*rRngY))//&
						', '//trim(r2c(rMinY+.5*rRngY))//', '//trim(r2c(rMinY+.75*rRngY))//&
						', '//trim(r2c(rMaxY))//') mirror format ""'
		write(iFileOut,'(A)')'unset ylabel'
		write(iFileOut,'(A)')'set tics front'
		write(iFileOut,'(A)')'set palette rgb -34,-35,-36'
		if(bChangeCB)then
			write(iFileOut,'(A)')'set cbrange ['//trim(r2c(rMinCB))//':'//trim(r2c(rMaxCB))//']'
		endif
		if(bPointwiseOnly)then
			write(iFileOut,'(A)')'plot "'//trim(cFileOut)//'.bin_XYp" u 1:4:2:3:5:6 w xyerrorbars notitle'
		elseif(bPointwise)then
			write(iFileOut,'(A)')'plot "'//trim(cFileOut)//'.bin_XY" u 1:2:3 w ima notitle, \'
			write(iFileOut,'(A)')' "'//trim(cFileOut)//'.bin_XYp" u 1:4:2:3:5:6 w xyerrorbars notitle linewidth 3.0 lc rgb "light-magenta"'
		else
			write(iFileOut,'(A)')'plot "'//trim(cFileOut)//'.bin_XY" u 1:2:3 w ima notitle'
		endif

		write(iFileOut,'(A)')'set origin 0.6,0.4'
		write(iFileOut,'(A)')'set size 0.15,0.5'
		write(iFileOut,'(A)')'set xrange [*:*] reverse'
		write(iFileOut,'(A)')'set yrange ['//trim(r2c(rMinY))//':'//trim(r2c(rMaxY))//']'

		write(iFileOut,'(A)')'set xtics autofreq format ""'
		write(iFileOut,'(A)')'unset xlabel'
		write(iFileOut,'(A)')'unset ytics'
		write(iFileOut,'(A)')'unset ylabel'
		write(iFileOut,'(A)')'unset colorbox'
		write(iFileOut,'(A)')'unset log x'
		write(iFileOut,'(A)')'unset log y'
		write(iFileOut,'(A)')'unset log cb'


		write(iFileOut,'(A)')'set y2tics ('//trim(r2c(rMinY))//', '//trim(r2c(rMinY+.25*rRngY))//&
						', '//trim(r2c(rMinY+.5*rRngY))//', '//trim(r2c(rMinY+.75*rRngY))//&
						', '//trim(r2c(rMaxY))//') mirror format '//trim(cTicks)
		if(.NOT.bNoTicks)write(iFileOut,'(A)')'set y2label "Delay (days)"'

		write(iFileOut,'(A)')'plot "'//trim(cFileOut)//'.bin_Y" u 2:1 w l notitle'

		write(iFileOut,'(A)')'set origin 0.1,0.25'
		write(iFileOut,'(A)')'set size 0.5,0.15'
		write(iFileOut,'(A)')'set xrange ['//trim(r2c(rMinX))//':'//trim(r2c(rMaxX))//'] noreverse'
		write(iFileOut,'(A)')'set yrange [*:*]'

		if(bNoTicks)then
			continue
		elseif(bLineVel)then
			write(iFileOut,'(A)')'set xtics ('//trim(r2c(rMinX))//&
							', "'//trim(r2c(rLineVelLower))//'" '//trim(r2c(rLineLambdaLower))//&
							', "0" '//trim(r2c(rLineLambda))//&
							', "'//trim(r2c(rLineVelUpper))//'" '//trim(r2c(rLineLambdaUpper))//&
							', '//trim(r2c(rMaxX))//') mirror format '//trim(cTicks)	
		else
			write(iFileOut,'(A)')'set xlabel "Wavelength (10^{-10}cm)"'
			write(iFileOut,'(A)')'set xtics ('//trim(r2c(rMinX))//', '//trim(r2c(rMinX+.25*rRngX))//&
							', '//trim(r2c(rMinX+.5*rRngX))//', '//trim(r2c(rMinX+.75*rRngX))//&
							', '//trim(r2c(rMaxX))//') mirror format '//trim(cTicks)
		endif

		write(iFileOut,'(A)')'set ytics autofreq format ""'
		write(iFileOut,'(A)')'unset ylabel'
		write(iFileOut,'(A)')'unset y2tics'
		write(iFileOut,'(A)')'unset y2label'
		write(iFileOut,'(A)')'unset log x'
		write(iFileOut,'(A)')'unset log y'

		write(iFileOut,'(A)')'plot "'//trim(cFileOut)//'.bin_X" u 1:2 w l notitle'
		close(iFileOut)

		call  EXECUTE_COMMAND_LINE('gnuplot '//trim(cFileOut)//'.plot')
		if(.NOT.bMessy)then
			call EXECUTE_COMMAND_LINE('rm '//trim(cFileOut)//'.bin_* '//trim(cFileOut)//'.plot ')
		endif
	end do

	print *,"Finished"

contains
	Subroutine sLineAnalysis(arMapWave, arMapPath, arBinWave, arBinPath, rLineWave)
		Real(iKindDP), intent(in), dimension(:)	:: arMapWave, arMapPath, arBinWave, arBinPath
		Real(iKindDP), intent(in)				:: rLineWave
		Real(iKindDP)	:: rWaveWidth, rFreqPeak, rPathPeak, rPathCent, rPathCentU, rPathCentL
		Real(iKindDP)	:: rVelocity, rMass, rRadiusMin, rRadiusMax, rMassMin, rMassMax, rRadius
		Real(iKindDP)	:: rFormFactor, rFormFactorError
		Integer 		:: iPathPeak

		rPathPeak = rfFindPeak(arMapPath, arBinPath, iPathPeak)
		rPathCent = rfFindCentroid(arMapPath, arBinPath, rPathCentL, rPathCentU, 0.8*arBinPath(iPathPeak))
		if(bLineBHUseRMS)then
			rWaveWidth = rfFindRMS(arMapWave, arBinWave, rPathCent)
		else
			rWaveWidth = rfFindFWHM(arMapWave, arBinWave)
		endif
		print *,'Path centroid for line: '//trim(r2c(rPathCent))//' days (1σ: '//trim(r2c(rPathCentL))//' - '//trim(r2c(rPathCentU))//' days)'
		print *,'Width for line:'//trim(r2c(rWaveWidth))
		
		rFormFactor = 5.0
		rFormFactorError = 1.0

		rVelocity 	= rcC 		  * rWaveWidth  / rLineWave
		rRadius 	= rPathCent   * rSecsToDays * rcC
		rRadiusMin 	= rPathCentL  * rSecsToDays * rcC
		rRadiusMax 	= rPathCentU  * rSecsToDays * rcC
		rMass 		= rFormFactor * rRadius 	* rVelocity**2 / (rcG * rcMSol)
		rMassMin 	= (rFormFactor - rFormFactorError) * rRadiusMin	* rVelocity**2 / (rcG * rcMSol)
		rMassMax 	= (rFormFactor + rFormFactorError) * rRadiusMax * rVelocity**2 / (rcG * rcMSol)

		print *,'Mass for line: '//trim(r2c(rMass))//' ('//trim(r2c(rMassMin))//' - '//trim(r2c(rMassMax))//')'
	End Subroutine

	Real(iKindDP) Function rfFindFWHM(arVal, arBin)
		Real(iKindDP), intent(in), dimension(:)		:: arVal, arBin
		Real(iKindDP)	:: rPeak, rPeakVal, rHalfL, rHalfU
		Integer 		:: i, iPeak

		rPeak = rfFindPeak(arVal, arBin, iPeak, rPeakVal)

		do i=iPeak,size(arVal)
			if(arVal(i).LE.rPeakVal/2.0)then
				rHalfU = arBin(i)
				EXIT
			endif
		enddo
		do i=iPeak,1,-1
			if(arVal(i).LE.rPeakVal/2.0)then
				rHalfL = arBin(i+1)
				EXIT
			endif			 
		enddo
		rfFindFWHM = rHalfU - rHalfL
	End Function

	Real(iKindDP) Function rfFindRMS(arVal, arBin, rPeak)
		Real(iKindDP), intent(in), dimension(:)		:: arVal, arBin
		Real(iKindDP), intent(in)					:: rPeak
		Real(iKindDP)	:: rRMS, rValTot
		Integer 		:: i, iPeak

		rValTot = 0.0
		do i=1,size(arVal)
			rValTot = rValTot + arVal(i)
		enddo

		do i=1,size(arVal)
			rRMS = rRMS + (arVal(i) * (((arBin(i)+arBin(i+1))/2) - rPeak) ** 2.0)
		enddo

		rfFindRMS = sqrt(rRMS/rValTot)
	End Function

	Real(iKindDP) Function rfFindCentroid(arVal, arBin, rCentLOpt, rCentUOpt, rThresholdOpt)
		Real(iKindDP), intent(in), dimension(:) 	:: arVal, arBin
		Real(iKindDP), intent(in), optional 		:: rThresholdOpt
		Real(iKindDP), intent(out),optional			:: rCentLOpt, rCentUOpt
		Real(iKindDP)	:: rVal, rValCent, rCentL, rCentU, rCent, rThreshold
		Integer 		:: i
		Logical 		:: bFoundCentL, bFoundCentU, bFoundCent

		rValCent 	= 0.0
		rThreshold 	= 0.0
		if(present(rThresholdOpt)) rThreshold = rThresholdOpt
	
		do i=1,size(arVal)
			if(arBin(i) .GE. rThreshold) then
				rValCent	= rValCent + arVal(i)
			endif
		enddo

		if(rValCent.GT.0.0)then
			rCentL  	= rThreshold
			rVal 		= 0.0
			bFoundCentL = .FALSE.
			bFoundCentU = .FALSE.
			bFoundCent  = .FALSE.

			do i=1,size(arVal)
				if(arBin(i) .GE. rThreshold) then
					rVal = rVal + arVal(i)
					if(.not.bFoundCentL.AND.rVal.GT.rValCent*0.15865)then
						rCentL = arBin(i)
						bFoundCentL = .TRUE.
					endif
					if(.not.bFoundCent .AND.rVal.GT.rValCent*0.5)then
						rCent = (arBin(i)+arBin(i+1))/2
						bFoundCent  = .TRUE.
					endif
					if(.not.bFoundCentU.AND.rVal.GT.rValCent*0.84135)then
						rCentU = arBin(i+1)
						bFoundCentU = .TRUE.
					endif
				endif
			enddo
		else
			rCent 	= 0.0
			rCentL	= 0.0
			rCentU	= 0.0
		endif
		if(present(rCentLOpt)) rCentLOpt = rCentL
		if(present(rCentUOpt)) rCentUOpt = rCentU
		rfFindCentroid = rCent
	End Function

	Real(iKindDP) Function rfFindPeak(arVal, arBin, iPeakOpt, rPeakValOpt)
		Real(iKindDP), intent(in),	dimension(:)	:: arVal, arBin
		Real(iKindDP), intent(out), optional		:: rPeakValOpt
		Integer, intent(out), optional				:: iPeakOpt
		Real(iKindDP)	:: rPeakVal, rPeak
		Integer 		:: i, iPeak
		rPeakVal = 0.0

		do i=1, size(arVal) 
			if(arVal(i).GT.rPeakVal) then
				rPeakVal = arVal(i)
				rPeak = (arBin(i)+arBin(i+1))/2
				iPeak = i
			endif
		enddo
		if(present(iPeakOpt)) iPeakOpt = iPeak
		if(present(rPeakValOpt)) rPeakValOpt = rPeakVal
		rfFindPeak = rPeak
	End Function

	character(len=32) function r2c(rIn)
		real(iKindDP), intent(in) :: rIn
		write(r2c,'(ES12.5)')rIn
		r2c=adjustl(r2c)
	end function

	integer function ifLookupIndex(arBin, rVal)
		real(iKindDP), intent(in) :: arBin(:), rVal
		integer :: i

		ifLookupIndex = -1

		do i=1,size(arBin)-1
			if(rVal.GE.arBin(i).AND.rVal.LE.arBin(i+1)) ifLookupIndex = i
		end do
	end function
	
	real(iKindDP) function rfGetKeyword(cFile, cKey)
		character(len=*), intent(in) :: cFile, cKey
		
		logical :: bFound
		integer :: iFile = 102, iEOF, iErr,i
		character(len=512) :: cBuffer
		
		bFound=.FALSE.
		iErr=0
		iEOF=0
				
		open(iFile,file=cFile,status="OLD",action="READ",iostat=iErr)
		if(iErr.NE.0)then
			print *,"ERROR: Could not open file '"//&
					trim(cFile)//"' to read keyword "//trim(cKey)
			STOP
		endif
		
		do while(iEOF.EQ.0 .AND..NOT.bFound)
			read(iFile,'(A512)',iostat=iEOF)cBuffer
			cBuffer=adjustl(cBuffer)
			if(cBuffer(1:1).NE."#")then
				do i=1,len_trim(cBuffer)
					if(cBuffer(i:i).EQ."(".OR.cBuffer(i:i).EQ." ")then
						if(trim(cBuffer(:i-1)).EQ.trim(cKey))then
							bFound=.TRUE.
						endif
						EXIT
					endif
				enddo
				if(bFound)then		
					do i=i,len_trim(cBuffer)
						if(cBuffer(i:i).EQ.")".OR.cBuffer(i:i).EQ." ")then
							read(cBuffer(i+1:),*,iostat=iErr)rfGetKeyword
							if(iErr.NE.0)then
								print *,"ERROR: Could not read keyword "//trim(cKey)//&
										" value from '"//trim(adjustl(cBuffer(i+1:)))//"'"
								STOP
							endif					
						endif
					enddo
				endif	
			endif	
		end do
		close(iFile)
		
		if(.NOT.bFound)then
			print *,"ERROR: Could not find keyword '"//trim(cKey)//&
					"' in file '"//trim(cFile)//"'"
			STOP
		endif
	end function
	
	logical function fbIsLetter(cIn)
		character, intent(in) :: cIn
		if(cIn.GE.'a'.AND.cIn.LE.'z'.OR.cIn.GE.'A'.AND.cIn.LE.'Z')then
			fbIsLetter = .TRUE.
		else
			fbIsLetter = .FALSE.
		endif
	end function
	
	logical function fbIsKey(cIn)
		character, intent(in) :: cIn
		if(cIn.EQ."_".OR.cIn.EQ.".".OR.fbIsLetter(cIn))then
			fbIsKey = .TRUE.
		else
			fbIsKey = .FALSE.
		endif
	end function

end program main
