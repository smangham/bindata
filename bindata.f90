program main
	IMPLICIT NONE

	integer, parameter 				:: iFileIn=100,iFileOut=101
	integer, parameter 				:: iKindDP=selected_real_kind(15,300)
	real(iKindDP), parameter 		:: rcPi = DACOS(-1.D0), rSecsToDays=86400

	logical				:: bFound, bMessy=.FALSE., bNoKey=.FALSE., bNoTicks=.FALSE.
	logical				:: bAllScat=.FALSE., bNoLog=.FALSE., bLineMalformed=.FALSE., bUseExtracted=.FALSE.
	integer				:: iDimX=100, iDimY=100, iDimR=100
	integer				:: i,j,iErr=0, iEOF=0, iDummy, iBinX, iBinY, iBinR, iPhot=0,iPhotR=0, iPhotRE=0,iExtracted
	integer				:: iNScat, iNRScat, iNScatMin=1, iNScatMax=999, iScat, iArg=1,iLine=0
	integer 			:: iObserver, iObserverMin=0, iObserverMax=0, iObservers=1, iObs
	integer, allocatable			:: aiMap(:,:,:),aiMapX(:,:),aiMapY(:,:),aiMapR(:)
	real(iKindDP), allocatable		:: arMap(:,:,:),arMapX(:,:),arMapY(:,:),arBinX(:),arBinY(:)
	real(iKindDP)		:: rDummy, rLambda, rWeight, rDelay, rPosX,rPosY,rPosZ,rErr
	real(iKindDP)		:: rMinX=-1., rMaxX=-1., rMinY=-1., rMaxY=-1., rRngX=-1., rRngY=-1., rMinR=-1.,rMaxR=-1., rRngR=-1.
	real(iKindDP)		:: rMinP=1e300_iKindDP, rMaxP=0.,rMinI=1e300_iKindDP,rMaxI=0.,rRad=-1, rTemp=0.0
	character(len=512) 	:: cFileIn="",cFileOut="", cDummy, cArg, cTicks='"%g"'
	character(len=512)  :: cBuffer

	logical 			:: bReweight=.FALSE., bReweightBinLog=.FALSE., bReweightBinGeom=.FALSE., bLookupY=.TRUE.
	integer				:: iReweightGeom
	real(iKindDP)		:: rReweightPow, rReweightBase
	real(iKindDP), allocatable	:: arMapR(:),arBinR(:),arPosR(:),arReweightMult(:)

	integer				:: iErrWeight=0, iErrMalf=0, iErrLog=0
	real(iKindDP)		:: rPathMax=0;

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
		print *,"	-p VAL VAL"
		print *,"Minimum & maximum path distances to plot. Default is .9-.3.25x wind radius."
		print *,""
		print *,"	-v VAL VAL"
		print *,"Wavelength range to bin. Default is to use spectrum_wavemin-max from input file."
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
		print *,"	-bl"
		print *,"Reweight mode: Bin logarithmically, default is linear."
		print *,""
		print *,"	-bg"
		print *,"Reweight mode: Bin geometrically, default is linear."
		print *,""
		print *,"	-a"
		print *,"All scatters, bin both resonant and non-resonant."
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
		print *,"	-e"
		print *,"Extract mode, use only extracted photons."
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

		else if(cArg.EQ."-s".OR.cArg.EQ."-S")then
			call get_command_argument(iArg+1, cArg)
			read(cArg,*,iostat=iErr)iNScatMin
			if(iErr.NE.0 .OR.iNScatMin.LT.0)then
				print *,"ERROR: Minimum scatters argument '"//trim(cArg)//"' invalid!"
				STOP
			endif
			call get_command_argument(iArg+2, cArg)
			read(cArg,*,iostat=iErr)iNScatMax
			if(iErr.NE.0 .OR.iNScatMax.LT.iNScatMin)then
				print *,"ERROR: Maximum scatters argument '"//trim(cArg)//"' invalid!"
				STOP
			endif
			iArg=iArg+3

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
				!Nothing
			else if(iObserverMax.LT.iObserverMin)then
				iObserver = iObserverMin
				iObserverMin = iObserverMax
				iObserverMax = iObserver	
				iArg=iArg+1							
			else
				iArg=iArg+1				
			endif

		else if(cArg.EQ."-bg".OR.cArg.EQ."-BG")then
			iArg=iArg+1
			bReweightBinGeom=.TRUE.
		else if(cArg.EQ."-bl".OR.cArg.EQ."-BL")then
			iArg=iArg+1
			bReweightBinLog=.TRUE.

		else if(cArg.EQ."-e".OR.cArg.EQ."-E")then
			bUseExtracted=.TRUE.
			iArg=iArg+1
		else if(cArg.EQ."-a".OR.cArg.EQ."-A")then
			bAllScat=.TRUE.
			iArg=iArg+1
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
		rMinY= rRad*0.50
		rMaxY= rRad*4.50
		print '(X,A,ES8.2,A,ES8.2,A,I0,A)','Binning paths from ',rMinY,' to ',rMaxY,'cm in ',iDimY,' steps'
	endif
	rRngX=rMaxX-rMinX
	rRngY=rMaxY-rMinY
	rRngR=rMaxR-rMinR

	iObservers = 1 +(iObserverMax - iObserverMin)
	print '(A,I0,A,I0)','Plotting observers ',iObserverMin,' to ',iObserverMax

	allocate(aiMapX(iDimX,iObservers), aiMapY(iDimY,iObservers), aiMapR(iDimR), aiMap(iDimX,iDimY,iObservers))
	allocate(arMapX(iDimX,iObservers), arMapY(iDimY,iObservers), arMapR(iDimR), arMap(iDimX,iDimY,iObservers))
	allocate(arBinX(iDimX+1), arBinY(iDimY+1))

	do i=0,iDimX
		arBinX(i+1)=rMinX+i*(rMaxX-rMinX)/real(iDimX)
	end do
	do i=0,iDimY
		arBinY(i+1)=(rMinY+i*(rMaxY-rMinY)/real(iDimY)) / rSecsToDays
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
											iNScat, iNRScat, rDelay, iExtracted

				if(iErr.GT.0)then
					iErr=0
					iEOF=0
				elseif(iNScat.EQ.1 .AND. iExtracted.EQ.0) then
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
	iPhotRE=0
	iErr=0
	iEOF=0
	iLine=0
	open(iFileIn,file=trim(cFileIn)//".delay_dump",status="OLD",action="READ")
	read(iFileIn,*) cDummy
	iLine=iLine+1
	do while(iErr.EQ.0.AND.iEOF.EQ.0)
		read(iFileIn,'(A512)',iostat=iEOF) cBuffer
		iLine=iLine+1
		if(cBuffer(1:1).NE."#")then
			iPhot=iPhot+1
			read(cBuffer,*,iostat=iErr) rDummy, rLambda, rWeight, rPosX, rPosY, rPosZ, &
										iNScat, iNRScat, rDelay, iExtracted, iObserver
			if(bAllScat)iNRScat=iNScat

			if(iErr.GT.0)then
				iErr=0
				iEOF=0
				iErrMalf=iErrMalf+1
				if(iErrMalf.LT.10)then
					print '(X,A,I0)','WARNING: Malformed entry on line ',iLine
				elseif(iErrMalf.EQ.10)then
					print '(X,A,I0,A)','WARNING: Malformed entry on line ',iLine,'. Suppressing further warnings.'
				endif

			elseif(iExtracted.EQ.1 .AND..NOT.bUseExtracted)then
				!Do nothing
			elseif(iExtracted.EQ.0 .AND.bUseExtracted)then
				!Do nothing
			elseif(iObserver.LT.iObserverMin.OR.iObserver.GT.iObserverMax)then
				!Do nothing

			else if(iNRScat.GE.iNScatMin.AND.iNRScat.LE.iNScatMax.AND.iNScat.LE.iNScatMax)then
				if(rDelay.GT.rPathMax)rPathMax=rDelay !DEBUG

				iPhotR= iPhotR+1
				iPhotRE=iPhotRE+iExtracted
				iBinX = ifLookupIndex(arBinX,rLambda)
				iBinY = ifLookupIndex(arBinY,rDelay/rSecsToDays)
				if(iBinX.LT.1 .OR. iBinY.LT.1)then
					iErrLog=iErrLog+1
					if(iErrLog.LT.10) then
						print '(X,A,ES10.4,A,ES10.4,A,I0)','WARNING: Point at frequency ',rLambda,&
								', path ',rDelay,' lies outside range #',iErrLog
					elseif(iErrLog.EQ.10)then
						print '(X,A,ES10.4,A,ES10.4,A,I0,A)','WARNING: Point at frequency ',rLambda,&
								', path ',rDelay,' lies outside range #',iErrLog,". Suppressing further warnings."
					endif
				else
					!If using surface brightness correction, then adjust
					if(bReweight)then
						rRad    = sqrt(rPosX**2+rPosY**2)
						iBinR	= ifLookupIndex(arBinR, rRad)
						if(iBinR .GT. -1)then
							rWeight = arReweightMult(iBinR) * rWeight * rRad**(-1.5)
							aiMap(iBinX,iBinY,iObserver+1) 	= aiMap(iBinX,iBinY,iObserver+1) + 1
							arMap(iBinX,iBinY,iObserver+1) 	= arMap(iBinX,iBinY,iObserver+1) + rWeight
							aiMapX(iBinX,iObserver+1) 		= aiMapX(iBinX,iObserver+1) + 1
							arMapX(iBinX,iObserver+1) 		= arMapX(iBinX,iObserver+1) + rWeight
							aiMapY(iBinY,iObserver+1) 		= aiMapY(iBinY,iObserver+1) + 1
							arMapY(iBinY,iObserver+1) 		= arMapY(iBinY,iObserver+1) + rWeight
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
						aiMap(iBinX,iBinY,iObserver+1) = 	aiMap(iBinX,iBinY,iObserver+1) + 1
						arMap(iBinX,iBinY,iObserver+1) = 	arMap(iBinX,iBinY,iObserver+1) + rWeight
						aiMapX(iBinX,iObserver+1) = 		aiMapX(iBinX,iObserver+1) + 1
						arMapX(iBinX,iObserver+1) = 		arMapX(iBinX,iObserver+1) + rWeight
						aiMapY(iBinY,iObserver+1) = 		aiMapY(iBinY,iObserver+1) + 1
						arMapY(iBinY,iObserver+1) = 		arMapY(iBinY,iObserver+1) + rWeight
					endif

				endif
			endif
		endif
	end do
	!print '(1X,A,ES9.3,X,ES9.3)','Maximum path traveled [raw/radii]: ',rPathMax,rPathMax/rMaxR

	close(iFileIn)

	if(iErrLog.GT.0) print '(X,A,I0,A)',"WARNING: ",iErrLog," points lay outside the map range."
	if(iErrMalf.GT.0)print '(X,A,I0,A)',"WARNING: ",iErrMalf," lines in the file were malformed."

	print '(X,I0,A,I0,A,I0,A)',iPhotR,"/",iPhot," photons resonance scattered, of which ",iPhotRE," were extracted."
	if(iPhotR.EQ.0)then
		print *,'ERROR: No photons fit scattering parameters'
		STOP
	endif

	do iObs=iObserverMin+1,iObserverMax+1
		if(iObservers.GT.1)then
			print '(A,I0,A)','Writing to "'//trim(cFileOut)//'.',iObs-1,'.eps"'
		else
			print '(A)','Writing to "'//trim(cFileOut)//'.eps"'
		endif

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
			write(iFileOut,'(A,I0,A)')'set output "'//trim(cFileOut)//'.',iObs-1,'.eps"'
		else
			write(iFileOut,'(A)')'set output "'//trim(cFileOut)//'.eps"'
		endif
		write(iFileOut,'(A)')'set multiplot'
		write(iFileOut,'(A)')'set bmargin 0; set tmargin 0; set lmargin 0; set rmargin 0'
		if(.NOT.bNoLog) write(iFileOut,'(A)')'set log cb'
		write(iFileOut,'(A)')'set colorbox user origin 0.65,0.05 size .05,0.3'
		write(iFileOut,'(A)')'set cblabel "Flux (erg s^{-1} cm^{2})"'
		if(bNoKey)then
			write(iFileOut,'(A)')'unset colorbox'		
		endif

		write(iFileOut,'(A)')'set origin 0.1,0.4'
		write(iFileOut,'(A)')'set size 0.5,0.5'
		write(iFileOut,'(A)')'set xrange ['//trim(r2c(rMinX))//':'//trim(r2c(rMaxX))//']'
		write(iFileOut,'(A)')'set yrange ['//trim(r2c(rMinY))//':'//trim(r2c(rMaxY))//']'

		write(iFileOut,'(A)')'set xtics ('//trim(r2c(rMinX))//', '//trim(r2c(rMinX+.25*rRngX))//&
						', '//trim(r2c(rMinX+.5*rRngX))//', '//trim(r2c(rMinX+.75*rRngX))//&
						', '//trim(r2c(rMaxX))//') mirror format ""'
		write(iFileOut,'(A)')'unset xlabel'
		write(iFileOut,'(A)')'set ytics ('//trim(r2c(rMinY))//', '//trim(r2c(rMinY+.25*rRngY))//&
						', '//trim(r2c(rMinY+.5*rRngY))//', '//trim(r2c(rMinY+.75*rRngY))//&
						', '//trim(r2c(rMaxY))//') mirror format ""'
		write(iFileOut,'(A)')'unset ylabel'
		write(iFileOut,'(A)')'set tics front'
		write(iFileOut,'(A)')'set palette rgb -34,-35,-36'
		write(iFileOut,'(A)')'plot "'//trim(cFileOut)//'.bin_XY" u 1:2:3 w ima notitle'

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

		if(.NOT.bNoTicks)write(iFileOut,'(A)')'set xlabel "Wavelength (10^{-10}cm)"'
		write(iFileOut,'(A)')'set xtics ('//trim(r2c(rMinX))//', '//trim(r2c(rMinX+.25*rRngX))//&
						', '//trim(r2c(rMinX+.5*rRngX))//', '//trim(r2c(rMinX+.75*rRngX))//&
						', '//trim(r2c(rMaxX))//') mirror format '//trim(cTicks)
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
