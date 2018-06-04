!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     AlphaMLP.f90
!
! DESCRIPTION:
!> @brief    A new implimentation of gene prob
!
!> @details  A revised version of the Kerr and Kinghorn (1995) algorithm designed
! to handle multiple snps and hopefully simplify the notation.
! General TODO:
! Include better updates of MAF estimates
! Swap multiplication across the board to log addition.
! Check accuracy.
!> @author   Andrew Whalen, awhalen@roslin.ed.ac.uk
!
!> @date     December 7, 2016
!
!> @version  0.0.1 (alpha)
!
! REVISION HISTORY:
! 2016-12-7 AWhalen - Initial Version
! 2016-12-7 AWhalen - Change to log scale
!-------------------------------------------------------------------------------


module AlphaMLPModule
    use globalGP
    use blas95
    use iso_fortran_env
    use f95_precision
    use IFPORT
    use AlphaMLPInputModule

    interface runAlphaMLP
        ! module procedure runAlphaMLPAlphaImpute
        module procedure runAlphaMLPIndependently
    end interface runAlphaMLP
    contains 





    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      run the AlphaMLP program as a function from alphaimpute
    !
    !> @details   This function contains paralellization with openMp. This is done with the idea that it can be wrapped up in another parallel framework around it
    !
    !> @author    David Wilson, david.wilson@roslin.ed.ac.uk
    !
    !> @date       Febuary 7, 2016
    !---------------------------------------------------------------------------
    subroutine runAlphaMLPAlphaImpute(startSnp, endSnp, ped, AlphaMLPOutput, Maf)

        use globalGP, only :pedigree, nsnps, inputParams


        integer, intent(in) :: startSnp, endSnp
        integer :: i, j
        type(PedigreeHolder) ,target:: ped
        real(kind=real64), dimension(:,:,:), allocatable, intent(out) :: AlphaMLPOutput !< output 3 dimensional array as requeuired by alphaimpute. 1:pedigree%pedigreeSize, nSnps, nHaplotypes
        real(kind=real64),allocatable,dimension (:), intent(out) :: Maf !< double vector containing MaF for each Snp
        integer :: nHaplotypes
        pedigree => ped
        nHaplotypes = 4

       

        ! if (.not. allocated(AlphaMLPOutput)) then
        !     allocate(AlphaMLPOutput(pedigree%pedigreesize, startSnp:endSnp, 4))
        ! endif
        call pedigree%getMatePairsAndOffspring(offspringList, listOfParents, nMatingPairs)
        inputParams = AlphaMLPInput(startSnp,endSnp,"single",.true.)
        
        call setupPhaseChildOfFounders()
        call setupTraceTensor
        call setupGenerations()
        nSnps = inputParams%endSnp-inputParams%startSnp+1

        nSnpsAll = inputParams%nSnp
        nAnimals = pedigree%pedigreeSize
        founders = pedigree%founders%convertToArrayIDs()



        if(inputParams%runType == "single") call runSingleLocus(AlphaMLPOutput, maf,output=.false. )
        if(inputParams%runType == "multi") call runIndependentMultiLocus()
        ! do i = startSnp, endSnp
        !     do j = 1, 4
        !         AlphaMLPOutput(:, i, j) =  currentPeelingEstimates(i)%haplotypeEstimates(j,:)
        !     enddo
        !     maf(i-inputParams%startSnp+1) = currentPeelingEstimates(i-inputParams%startSnp+1)%maf
        ! enddo

        ! deallocate(currentPeelingEstimates)

    end subroutine runAlphaMLPAlphaImpute




    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      run the AlphaMLP program independently
    !
    !> @details   This function contains paralellization with openMp.
    !>   Outputs are written to a file
    !
    !> @author    David Wilson, david.wilson@roslin.ed.ac.uk
    !
    !> @date       Febuary 7, 2016
    !---------------------------------------------------------------------------
    subroutine runAlphaMLPIndependently()
        use globalGP, only: nSnps, sequenceData
        use alphaFullChromModule
        implicit none
        character(len=4096) :: SpecFile
        type(AlphaMLPInput) :: inputParamsIn

        if (Command_Argument_Count() > 0) then
            call Get_Command_Argument(1,SpecFile)
        else
            SpecFile="AlphaMLPSpec.txt"
        end if


        inputParamsIn = AlphaMLPInput(SpecFile)

        if (inputParamsIn%plinkinputfile /= "") then
            call runPlink(inputParamsIn%plinkinputfile, inputParamsIn, runAlphaMLPBasedOnSpec)

        else
            call runAlphaMLPBasedOnSpec(inputParamsIn)
        endif


    end subroutine runAlphaMLPIndependently

    subroutine runAlphaMLPBasedOnSpec(inputParamsIn, pedigreeIn)         

        use globalGP, only: nSnps, sequenceData, inputParams, nAnimals, founders, pedigree
        implicit none
        class(baseSpecFile), target :: inputParamsIn
        type(pedigreeHolder), target, optional :: pedigreeIn

        real(kind=real64), dimension(:,:,:), allocatable:: outputHaplotypes
        real(kind=real64), dimension(:), allocatable :: maf


        select type(inputParamsIn)

            type is (AlphaMLPInput)
            inputParams = inputParamsIn
            class default
            write(error_unit, *) "ERROR: AlphaImpute given correct object type as input"
            call abort()
        end select

        nSnps = inputParams%endSnp-inputParams%startSnp+1
        nSnpsAll = inputParams%nSnp


        inputParams%prefix = trim(inputParams%resultFolderPath)//trim(inputParams%basePrefix)
        if (inputParams%plinkinputfile /= "") then
            inputParams%endSnp = inputParamsIn%nSnp
            inputParams%startSnp = 1
            pedigree => pedigreeIn
        else
            call setupPedigree(inputParams)
        endif   

        nAnimals = pedigree%pedigreeSize
        founders = pedigree%founders%convertToArrayIDs()

        call pedigree%getMatePairsAndOffspring(offspringList, listOfParents, nMatingPairs)
        call setupPhaseChildOfFounders()
        call setupTraceTensor
        call setupGenerations()


        if(inputParams%runType == "single") call runSingleLocus(outputHaplotypes, maf ,output=.true.)
        if(inputParams%runType == "multi") call runIndependentMultiLocus()

    end subroutine

    subroutine runIndependentMultiLocus()
        use globalGP, only: nSnps, sequenceData, inputParams
        implicit none
        integer :: i
        type(peelingEstimates), dimension(:), pointer :: currentPeelingEstimates

        call runMultiLocusAlphaMLP(currentPeelingEstimates)
        
        call writeOutputsToFileMultiLocus(currentPeelingEstimates)
        
        if(associated(currentPeelingEstimates)) then
            do i = 1, nSnps
                call currentPeelingEstimates(i)%deallocatePeelingEstimates()
            enddo
            currentPeelingEstimates => null()
        endif

    end subroutine



    

    subroutine runSingleLocus(outputHaplotypes, maf, output)
        use globalGP, only: nSnps, inputParams, mapIndexes, mapDistance
        implicit none
        integer :: i, snpID
        type(peelingEstimates), pointer :: markerEstimates
        logical, optional :: output
        integer(kind=1), dimension(:), allocatable :: tmpGenotypes
        real(kind=real64), dimension(:,:), allocatable:: bayesianProduct
        real(kind=real64), dimension(:,:,:), allocatable:: outputHaplotypes
        real(kind=real64), dimension(:,:), pointer:: outputDosages
        real(kind=real64), dimension(nSnps) :: markerError
        real(kind=real64), dimension(:), allocatable :: maf
        real(kind=real64), dimension(:,:,:), allocatable:: segregationEstimates
        real(kind=real64), dimension(:,:), allocatable:: markerSegregation
        integer :: segregationOffset, prevSnpSegID, nextSnpSegID
        
        allocate(outputHaplotypes(nHaplotypes, nSnps, nAnimals))
        allocate(outputDosages(nSnps, nAnimals))

        allocate(mapIndexes(2, nSnpsAll))
        allocate(mapDistance(nSnpsAll))
        allocate(maf(nsnps))
        !print *, "Running in single locus peeling mode"
        if(inputParams%segFile .ne. "No seg") then 
            call readSegregationFile(inputParams, segregationEstimates, segregationOffset, mapIndexes, mapDistance, pedigree)
        endif      
        

        !!$OMP parallel do &
        !!$OMP default(shared) & 
        !!$OMP private(i, snpID, markerEstimates,prevSnpSegID,nextSnpSegID,markerSegregation)
        do i = 1, nSnps
            print *, i
            allocate(markerEstimates)
            snpID = inputParams%startSnp + i - 1
            !Note: We want the SNP by order in which it is loaded in, not the snp id
            call markerEstimates%initializePeelingEstimates(nHaplotypes, nAnimals, nMatingPairs, nSnpsAll, i) 
            if(inputParams%segFile .ne. "No seg") then 
                !Get the estimate midway between the two markers. 

                prevSnpSegID = mapIndexes(1, snpID) - segregationOffset 
                nextSnpSegID = mapIndexes(2, snpID) - segregationOffset
                markerSegregation = (1-mapDistance(i))*segregationEstimates(:,prevSnpSegID,:) + &
                                            mapDistance(i)*segregationEstimates(:,nextSnpSegID,:) 

                markerEstimates%fullSegregation = markerSegregation
            else
                markerEstimates%fullSegregation = .25
            endif
            call runSingleIndex(markerEstimates)
            outputHaplotypes(:, i, :) = markerEstimates%haplotypeEstimates
            outputDosages(i, :) = markerEstimates%genotypeEstimates

            markerError(i) = markerEstimates%genotypingErrorRate
            maf(i) = markerEstimates%maf

            call markerEstimates%deallocateMarkerVariables()
            call markerEstimates%deallocatePeelingEstimates()
            deallocate(markerEstimates)
        enddo
        !!$OMP end parallel do

        if (output) then
            call writeOutputsToFileSingleLocus(outputHaplotypes, outputDosages, markerError, maf)
        endif


        ! deallocate(outputHaplotypes)
        ! deallocate(outputDosages)
    end subroutine runSingleLocus

    ! Two subroutines below for running the HMM.

    subroutine runMultiLocusAlphaMLP(currentPeelingEstimates)
        use globalGP, only: nAnimals, nHaplotypes, nMatingPairs, nSnps, inputParams, nSnpsAll
        implicit none
        type(peelingEstimates), dimension(:), pointer, intent(out) :: currentPeelingEstimates
        integer :: i, nCycles, cycleIndex
        logical :: converged

        !print *, "Running in multi locus peeling mode"

        allocate(currentPeelingEstimates(nSnps))
        do i = 1, nSnps
            !print *, "Allocating ", i
            call currentPeelingEstimates(i)%initializePeelingEstimates(nHaplotypes, nAnimals, nMatingPairs, nSnpsAll, inputParams%startSnp + i - 1)
        enddo


        converged = .false.

        nCycles = inputParams%nCycles
        !!print *, nCycles

        !Handle the Multilocus Peeler
        cycleIndex = 1
        converged = .false.
        do while(cycleIndex < nCycles .and. .not. converged)
            ! Forward Pass
            print *, "cycle ", cycleIndex, ", Forward "
            do i = 2, nSnps
                call runIndex(pedigree%getAllGenotypesAtPositionWithUngenotypedAnimals(i), i, currentPeelingEstimates, 1)
                if(mod(i, 100) .eq. 0) print *, "cycle ", cycleIndex, ", Forward ", i
                 ! print *, "cycle ", cycleIndex, ", Forward ", i
            enddo
            
            ! Backward Pass
            !!print *, "cycle ", cycleIndex, ", Backward "
            do i = nSnps-1, 1, -1
                call runIndex(pedigree%getAllGenotypesAtPositionWithUngenotypedAnimals(i), i, currentPeelingEstimates, 2)
                if(mod(i, 100) .eq. 0) print *, "cycle ", cycleIndex, ", Backward ", i
            enddo

            ! Join Pass                
            !!print *, "cycle ", cycleIndex, ", Join "
            do i = nSnps, 1, -1
                call runIndex(pedigree%getAllGenotypesAtPositionWithUngenotypedAnimals(i), i, currentPeelingEstimates, 3, .true.)
                if(mod(i, 100) .eq. 0) print *, "cycle ", cycleIndex, ", Join ", i
            enddo

            if(cycleIndex > 1) converged = checkConvergence(currentPeelingEstimates)
            cycleIndex = cycleIndex + 1
            ! call updateAllRecombinationRates(currentPeelingEstimates)

        enddo
    
    end subroutine

    subroutine runSingleIndex(markerEstimates)
        use globalGP, only: nAnimals, nHaplotypes, nMatingPairs, nSnps, inputParams, nSnpsAll
        implicit none
        real(kind=real64), dimension(nHaplotypes,0:9) :: genotypesToHaplotypes

        type(peelingEstimates), pointer, intent(inout) :: markerEstimates
        integer :: i, nCycles, cycleIndex
        logical :: converged
        real(kind=real64), dimension(:,:), allocatable :: haplotypesOut
        real(kind=real64), dimension(:), allocatable :: previousGenotypeEstimate
        real(kind=real64) :: error, p, q, pf
        integer, dimension(:), allocatable :: ref, alt
        integer :: iteration, maxIteration

        call markerEstimates%allocateMarkerVariables(nHaplotypes, nAnimals, nMatingPairs)
        markerEstimates%currentSegregationEstimate = markerEstimates%fullSegregation
        !Setup Penetrance

        ! !print *, markerEstimates%penetrance
        ! !print *, markerEstimates%penetrance(:, 1:5)
        !Make sure this is set from a file.
        call buildSegregationTraceTensor(markerEstimates)
        maxIteration = inputParams%nCycles
        iteration = 0
        converged = .false.
        do while (.not. converged .and. iteration < maxIteration)
            ! !print *, "Running iteration", iteration
            iteration = iteration + 1
            ! print *, iteration
            !Create Anterior
            !Reset penetrance to account for new error rate.
            if (.not. inputParams%isSequence) then 
                call markerEstimates%setPenetranceFromGenotypes()
            else 
                call markerEstimates%setPenetranceFromSequence()
            endif


            call markerEstimates%setAnterior(markerEstimates%maf)

            call performPeeling(markerEstimates, doUpdateSegregation = .false.)

            
            previousGenotypeEstimate = markerEstimates%genotypeEstimates
            markerEstimates%genotypeEstimates = markerEstimates%getGenotypeEstimates()

            markerEstimates%estimatedError = sum(abs(previousGenotypeEstimate - markerEstimates%genotypeEstimates))/size(markerEstimates%genotypeEstimates)

            if(.not. inputParams%isSequence) call markerEstimates%updateGenotypeErrorRates()
            if(inputParams%isSequence) call markerEstimates%updateSequenceErrorRates()
            call updateMafEstimates(markerEstimates%genotypeEstimates, markerEstimates)
            converged = markerEstimates%estimatedError < .000001

        enddo
        ! !print *, "Finished after", iteration


    end subroutine
    
    function checkConvergence(currentPeelingEstimates) result(res)
        type(peelingEstimates), dimension(:), pointer, intent(in) :: currentPeelingEstimates
        real(kind=real64), dimension(nSnps) :: markerError
        logical :: res
        integer :: i
        do i = 1, nSnps
            markerError(i) = currentPeelingEstimates(i)%estimatedError
        enddo
        !print *, "MaxMarkerError", maxval(markerError) 
        res = maxval(markerError) < .0001

    end function
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Runs an index.
    !
    !> @details     Runs an index. By setting up an initial set of haplotype estimates and then peeling up and peeling down.
    ! Still need to add in updating MAF information, and revising the anterior estimate with new maf information. 
    !
    !> @author     Andrew Whalen, awhalen@roslin.ed.ac.uk
    !
    !> @date       December 7, 2016
    !
    ! PARAMETERS:
    !> 
    !---------------------------------------------------------------------------
    

    subroutine runIndex(genotypes, indexNumber, currentPeelingEstimates, runType, updateParamsTemp)
        !Variable Declarations
        use globalGP, only : nHaplotypes, nAnimals, inputParams
        use alphaStatMod, only: CorrelationReal64, cor
        implicit none
        integer(kind=1), dimension(:), intent(in) :: genotypes
        type(peelingEstimates), dimension(:), pointer, intent(inout) :: currentPeelingEstimates
        integer, intent(in) :: indexNumber, runType

        logical, optional :: updateParamsTemp
        logical :: updateParams

        real(kind=real64), dimension(nHaplotypes,nAnimals):: logSum, haplotypeEstimates, anterior
        real(kind=real64), dimension(nHaplotypes,0:9) :: genotypesToHaplotypes
        real(kind=real64), dimension(nAnimals) :: genotypeEstimates, weights, previousGenotypeEstimate
        real(kind=real64), dimension(:), allocatable :: overrideEstimate
        real(kind=real64), dimension(:,:), allocatable :: tmpVect
        real(kind=real64), dimension(:,:), pointer :: posterior

        integer, dimension(:), allocatable :: tmpFamilyList
        integer :: i, j, fam, father, mate, k, blockSize
        type(peelingEstimates), pointer :: markerEstimates
        logical :: usePhaseOverride

        if(present(updateParamsTemp)) then
            updateParams = updateParamsTemp
        else
            updateParams = .false.
        endif


        markerEstimates => currentPeelingEstimates(indexNumber)
        call markerEstimates%allocateMarkerVariables(nHaplotypes, nAnimals, nMatingPairs)
        
        !Setup Penetrance
        if (.not. inputParams%isSequence) then 
            call markerEstimates%setPenetranceFromGenotypes()
        else 
            call markerEstimates%setPenetranceFromSequence()
        endif

        !Create Anterior
        call markerEstimates%setAnterior(markerEstimates%maf)

        markerEstimates%anterior = anterior

        !Create posterior
        !Not using pointers... maybe should?
        markerEstimates%sirePosteriorMate = markerEstimates%sirePosteriorMateAll(:,:,runType)
        markerEstimates%damePosteriorMate = markerEstimates%damePosteriorMateAll(:,:,runType)
        markerEstimates%posterior = markerEstimates%posteriorAll(:,:,runType)
        posterior => markerEstimates%posterior
        !Build segregation.
        
        call getSegregationEstimate(markerEstimates, currentPeelingEstimates, indexNumber, runType)
        call buildSegregationTraceTensor(markerEstimates)


        !Run Analysis

        usePhaseOverride = .true.
        if(usePhaseOverride) then
            do j=1, size(phaseChildren)
                if(phaseChildren(j) > 1) then
                    overrideEstimate = phaseChildrenOverride(:,j) * markerEstimates%currentSegregationEstimate(:,phaseChildren(j))
                    markerEstimates%currentSegregationEstimate(:,phaseChildren(j)) = overrideEstimate
                endif
            enddo
        endif

        call performPeeling(markerEstimates, doUpdateSegregation = .true.)
        

        !Now do post segregation updating.
        !Isn't this all currentSegregationEstimate -- yes, but it's saving it to different slots.
        if(runType == 1) markerEstimates%transmitForward = markerEstimates%currentSegregationEstimate
        if(runType == 2) markerEstimates%transmitBackward = markerEstimates%currentSegregationEstimate
        if(runType == 3) markerEstimates%fullSegregation = markerEstimates%currentSegregationEstimate

        do i = 1, nAnimals            
            markerEstimates%pointSegregation(:,i) = markerEstimates%pointSegregation(:,i)/sum(markerEstimates%pointSegregation(:,i))

            if(runType == 1) markerEstimates%transmitForward(:,i) = markerEstimates%transmitForward(:,i)/sum(markerEstimates%transmitForward(:,i))
            if(runType == 2) markerEstimates%transmitBackward(:,i) = markerEstimates%transmitBackward(:,i)/sum(markerEstimates%transmitBackward(:,i))
            if(runType == 3) markerEstimates%fullSegregation(:,i) = markerEstimates%fullSegregation(:,i)/sum(markerEstimates%fullSegregation(:,i))
        enddo

        !Save values and update Params
        if(runType == 3) then
            logSum = markerEstimates%anterior + markerEstimates%posterior + markerEstimates%penetrance
            do i = 1, size(logSum, 2)
                haplotypeEstimates(:,i) = lhtp(logSum(:,i))
            enddo
            markerEstimates%haplotypeEstimates = haplotypeEstimates
            previousGenotypeEstimate = markerEstimates%genotypeEstimates

            genotypeEstimates=0
            call gemv(haplotypeEstimates, haplotypesToGenotypes, genotypeEstimates, trans = "T")
            weights = sum(haplotypeEstimates, 1)
            genotypeEstimates = genotypeEstimates/weights
            markerEstimates%genotypeEstimates = genotypeEstimates

            markerEstimates%estimatedError = sum(abs(previousGenotypeEstimate - markerEstimates%genotypeEstimates))/size(markerEstimates%genotypeEstimates)
        endif
        
        if(updateParams) then
            call updateRecombinationRate(markerEstimates, currentPeelingEstimates, indexNumber)
            ! if(markerEstimates%recombinationRate > .01) !print *, markerEstimates%recombinationRate
            if(.not. inputParams%isSequence) call markerEstimates%updateGenotypeErrorRates()
            if(inputParams%isSequence) call markerEstimates%updateSequenceErrorRates()
            call updateMafEstimates(genotypeEstimates, markerEstimates)
        endif
        
        markerEstimates%sirePosteriorMateAll(:,:,runType) = markerEstimates%sirePosteriorMate 
        markerEstimates%damePosteriorMateAll(:,:,runType) = markerEstimates%damePosteriorMate 
        markerEstimates%posteriorAll(:,:,runType) = markerEstimates%posterior

        call markerEstimates%deallocateMarkerVariables
    end subroutine runIndex

    subroutine performPeeling(markerEstimates, doUpdateSegregation) 
        use globalGP, only : nAnimals, segregationTensor, familiesInGeneration, childrenAtGeneration, listOfParents
        implicit none
        type(PeelingEstimates), pointer :: markerEstimates
        integer, dimension(:), allocatable :: tmpFamilyList
        real(kind=real64), dimension(:,:), pointer :: posterior
        real(kind=real64), dimension(:,:), allocatable :: oldPosteriorSire, oldPosteriorDame
        logical :: doUpdateSegregation
        integer :: fam, i, j, father, mate
        
        posterior => markerEstimates%posterior
        oldPosteriorSire = markerEstimates%sirePosteriorMate
        oldPosteriorDame = markerEstimates%damePosteriorMate


    ! PEEL DOWN 
        do i = 1, nGenerations
            tmpFamilyList = familiesInGeneration(i)%array
            !!$OMP PARALLEL DO DEFAULT(SHARED) &
            !!$OMP PRIVATE(fam,j)
            do j = 1, size(tmpFamilyList)
                fam = tmpFamilyList(j)        
                call peelDown(markerEstimates, fam)
            enddo
            !!$OMP END PARALLEL DO
        enddo 


! PEEL UP + SEGREGATION ESTIMATES
        do i = nGenerations, 1, -1
            tmpFamilyList = familiesInGeneration(i)%array
            if(doUpdateSegregation) then
                !!$OMP PARALLEL DO DEFAULT(SHARED) &
                !!$OMP PRIVATE(fam,j)
                do j = 1, size(tmpFamilyList)
                    fam = tmpFamilyList(j)        
                    call updateSegregation(markerEstimates, fam)
                enddo
                !!$OMP END PARALLEL DO
                markerEstimates%currentSegregationEstimate(:,childrenAtGeneration(i)%array) = &
                                        markerEstimates%currentSegregationEstimate(:, childrenAtGeneration(i)%array) *&
                                        markerEstimates%pointSegregation(:, childrenAtGeneration(i)%array)

                call buildNewSegregationTensors(markerEstimates, childrenAtGeneration(i)%array)
            endif
            

            !!$OMP PARALLEL DO DEFAULT(SHARED) &
            !!$OMP PRIVATE(fam,j)
            do j = 1, size(tmpFamilyList)
                fam = tmpFamilyList(j)        
                call peelUp(markerEstimates, fam)
            enddo
            !!$OMP END PARALLEL DO


            !!$OMP PARALLEL DO DEFAULT(SHARED) &
            !!$OMP PRIVATE(fam, father,mate, j)
            do j = 1, size(tmpFamilyList)
                fam = tmpFamilyList(j)
                father = listOfParents(1, fam)
                mate = listOfParents(2, fam)
                
                posterior(:,father) = posterior(:,father) - oldPosteriorSire(:,fam) + markerEstimates%sirePosteriorMate(:, fam)
                posterior(:,father) = posterior(:,father) - maxval(posterior(:,father))
                
                posterior(:,mate) = posterior(:,mate) - oldPosteriorDame(:,fam) + markerEstimates%damePosteriorMate(:, fam)
                posterior(:,mate) = posterior(:,mate) - maxval(posterior(:,mate))
            enddo
            !!$OMP END PARALLEL DO 
        enddo 
    end subroutine


!---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Perform a peeling down, updating the anterior probabilities.
    !
    !> @details     Perform a peeling down, updating the anterior probabilities. 
    ! This is done by iterating over mating paris and updating the anterior probabilities of their children.
    ! TODO: Include updating MAF information.
    !
    !> @author     Andrew Whalen, awhalen@roslin.ed.ac.uk
    !
    !> @date       December 7, 2016
    !
    ! PARAMETERS:
    !> 
    !---------------------------------------------------------------------------
    
    
    subroutine peelDown(markerEstimates, fam)
        use globalGP, only: nHaplotypes, listOfParents, offspringList
        use blas95
        implicit none        
        type(PeelingEstimates), pointer, intent(inout) :: markerEstimates
        integer, intent(in) :: fam
        real(kind=real64), dimension(:,:), pointer :: posterior, penetrance, anterior,sirePosteriorMate, damePosteriorMate
        real(kind=real64), dimension(nHaplotypes) :: pfather, pmate
        real(kind=real64), dimension(nHaplotypes, nHaplotypes) :: pjoint
        real(kind=real64), dimension(:,:,:), allocatable :: childEstimate
        real(kind=real64), dimension(:,:,:), pointer, contiguous :: temp 
        integer, dimension(:), allocatable :: offspring
        integer :: father, mate, nChildren, child
        integer j

        anterior => markerEstimates%anterior
        penetrance => markerEstimates%penetrance
        posterior => markerEstimates%posterior
        sirePosteriorMate => markerEstimates%sirePosteriorMate
        damePosteriorMate => markerEstimates%damePosteriorMate

        father = listOfParents(1, fam)
        mate = listOfParents(2, fam)
        pfather = anterior(:, father)+penetrance(:,father)+posterior(:,father)-sirePosteriorMate(:,fam)
        pmate = anterior(:, mate)+penetrance(:,mate)+posterior(:,mate)-damePosteriorMate(:,fam)
    
        pjoint = 0
        call additiveOuterProduct(pmate, pfather, pjoint)
        call additiveOuterProductSpread(pmate, pfather, pjoint)
        
        
        offspring = offspringList(fam)%convertToArrayIDs()
        nChildren = offspringList(fam)%length
        allocate(childEstimate(nHaplotypes,nHaplotypes, nChildren))
        
        do j=1, nChildren
            child = offspring(j)
            temp => markerEstimates%currentSegregationTensors(:,:,:,child)
            childEstimate(:,:,j) = childTraceMultiply(penetrance(:, child) + posterior(:, child), temp)
            ! childEstimate(:,:,j) = childTraceMultiplyMKL(penetrance(:, child) + posterior(:, child), temp)
        enddo

        pjoint = pjoint + sum(childEstimate, dim=3)
        !MP
        do j=1, nChildren
            child = offspring(j)
            anterior(:,child) = parentTraceMultiply(pjoint-childEstimate(:,:,j), markerEstimates%currentSegregationTensors(:,:,:,child))
            anterior(:,child) = anterior(:,child) - maxval(anterior(:,child))
            !isPseudoFounder may not be allocated otherwise.
        enddo
        deallocate(offspring, childEstimate)

    end subroutine


    subroutine updateSegregation(markerEstimates, fam)
        use globalGP, only: nHaplotypes, listOfParents, offspringList
        use blas95
        implicit none
        type(PeelingEstimates), pointer, intent(inout) :: markerEstimates
        REAL(kind=real64), dimension(:,:), pointer :: anterior, posterior, penetrance, sirePosteriorMate, damePosteriorMate
        REAL(kind=real64), dimension(:,:), pointer :: newPointEstimate
        integer :: father, mate, nChildren, child
        real(kind=real64), dimension(nHaplotypes) :: pfather, pmate
        real(kind=real64), dimension(nHaplotypes, nHaplotypes) :: pjoint
        real(kind=real64), dimension(:,:,:), allocatable :: childEstimate

        real(kind=real64), dimension(:,:,:), pointer, contiguous :: temp 

        integer, dimension(:), allocatable :: offspring
        integer fam, j


        anterior => markerEstimates%anterior
        penetrance => markerEstimates%penetrance
        posterior => markerEstimates%posterior
        sirePosteriorMate => markerEstimates%sirePosteriorMate
        damePosteriorMate => markerEstimates%damePosteriorMate
        newPointEstimate => markerEstimates%pointSegregation

        father = listOfParents(1, fam)
        mate = listOfParents(2, fam)
        pfather = anterior(:, father)+penetrance(:,father)+posterior(:,father)-sirePosteriorMate(:,fam)
        pmate = anterior(:, mate)+penetrance(:,mate)+posterior(:,mate)-damePosteriorMate(:,fam)
        pjoint = 0
        call additiveOuterProduct(pmate, pfather, pjoint)
        call additiveOuterProductSpread(pmate, pfather, pjoint)

        
        nChildren = offspringList(fam)%length

        allocate(offspring(nChildren))
        allocate(childEstimate(nHaplotypes,nHaplotypes, nChildren))
        offspring = offspringList(fam)%convertToArrayIDs()

        !MP
        do j=1, nChildren
            child = offspring(j)
            temp => markerEstimates%currentSegregationTensors(:,:,:,child)
            childEstimate(:,:,j) = childTraceMultiply(penetrance(:, child)+posterior(:, child), temp)
            ! childEstimate(:,:,j) = childTraceMultiplyMKL(penetrance(:, child)+posterior(:, child), temp)
            pjoint = pjoint + childEstimate(:,:,j)
        enddo
        !MP
        do j=1, nChildren
            child = offspring(j)
            if(.not. isPhasedChild(child)) then
                newPointEstimate(:,child) = reduceSegregationTensor(posterior(:,child)+penetrance(:,child), pjoint-childEstimate(:,:,j))
            endif
        enddo
        deallocate(offspring, childEstimate)

    end subroutine


    subroutine peelUp(markerEstimates, fam)        
        use globalGP, only: nHaplotypes, listOfParents, offspringList
        use blas95
        implicit none
        type(PeelingEstimates), pointer, intent(inout) :: markerEstimates
        real(kind=real64), dimension(:,:), pointer :: posterior, penetrance, anterior,sirePosteriorMate, damePosteriorMate
        integer :: father, mate, nChildren, child
        real(kind=real64), dimension(nHaplotypes) :: pfather, pmate
        real(kind=real64), dimension(nHaplotypes, nHaplotypes) :: pjoint, tmp
        real(kind=real64), dimension(nHaplotypes) :: tempEstimate
        integer, dimension(:), allocatable :: offspring
        integer fam, j

        anterior => markerEstimates%anterior
        penetrance => markerEstimates%penetrance
        posterior => markerEstimates%posterior
        sirePosteriorMate => markerEstimates%sirePosteriorMate
        damePosteriorMate => markerEstimates%damePosteriorMate


        father = listOfParents(1, fam)
        mate = listOfParents(2, fam)
        pfather = anterior(:, father)+penetrance(:,father)+posterior(:,father)-sirePosteriorMate(:,fam)
        pmate = anterior(:, mate)+penetrance(:,mate)+posterior(:,mate)-damePosteriorMate(:,fam)
        
        nChildren = offspringList(fam)%length
        allocate(offspring(nChildren))
        offspring = offspringList(fam)%convertToArrayIDs()

        pjoint = 0
        !MP
        do j=1, nChildren
            !May be able to pre-multiply this stuff together.
            child = offspring(j)
            ! call printTrace(markerEstimates%currentSegregationTensors(:,:,:,child))
            tmp = childTraceMultiply(penetrance(:, child)+posterior(:, child), markerEstimates%currentSegregationTensors(:,:,:,child))
            ! tmp = childTraceMultiplyMKL(penetrance(:, child) + posterior(:, child), childTrace)
        
            pjoint = pjoint + tmp
           
        enddo
        tempEstimate = jointByMateToFather(pmate, pjoint)
        sirePosteriorMate(:,fam) = tempEstimate - maxval(tempEstimate)
        
        tempEstimate = jointByFatherToMate(pfather, pjoint)
        damePosteriorMate(:,fam) = tempEstimate - maxval(tempEstimate)
        deallocate(offspring)

    end subroutine


    subroutine setupGenerations() 
        use globalGP, only: nMatingPairs, familiesInGeneration, childrenAtGeneration, offspringList, nGenerations, nAnimals
        use graphModule, only : selectIndexesBasedOnMask
        integer, dimension(nMatingPairs) :: familyGeneration
        integer, dimension(nAnimals) :: animalGeneration
        integer :: i

        if(allocated(familiesInGeneration)) deallocate(familiesInGeneration)
        if(allocated(childrenAtGeneration)) deallocate(childrenAtGeneration)
        
        do i = 1, nMatingPairs
            familyGeneration(i) = offspringList(i)%first%item%generation
        enddo 

        nGenerations = maxval(familyGeneration)
        allocate(familiesInGeneration(nGenerations))
        do i = 1, nGenerations
            familiesInGeneration(i)%array = selectIndexesBasedOnMask(familyGeneration == i)
        enddo

        allocate(childrenAtGeneration(0:nGenerations))
        do i = 1, nAnimals
            animalGeneration(i) = pedigree%pedigree(i)%generation
        enddo
        do i = 0, nGenerations
            childrenAtGeneration(i)%array = selectIndexesBasedOnMask(animalGeneration == i) 
        enddo

    end subroutine

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief     Sets up pedigree
    !
    !> @details   This function contains paralellization with openMp. This is done with the idea that it can be wrapped up in another parallel framework around it
    !
    !> @author    David Wilson, david.wilson@roslin.ed.ac.uk
    !
    !> @date       Febuary 7, 2016
    !---------------------------------------------------------------------------
    subroutine setupPedigree(inputParams)
        use globalGP, only: pedigree, nAnimals, founders, sequenceData, nSnps
        use PedigreeModule
        use IndividualModule, only: Individual
        implicit  none
        type(AlphaMLPInput),intent(inout) :: inputParams 


        !Depreciating old genotype format.
        allocate(pedigree)
        call initPedigree(pedigree,trim(inputParams%pedFile), nsnps=nSnps)
        if (.not. inputParams%isSequence) then 
            call pedigree%addGenotypeInformationFromFile(inputParams%inputFile,inputParams%nsnp, startSnp=inputParams%startSnp, endSnp=inputParams%endSnp)     
        else 
            call pedigree%addSequenceFromFile(inputParams%sequenceFile, inputParams%nsnp, startSnp=inputParams%startSnp, endSnp=inputParams%endSnp)
        endif

    end subroutine
    

    subroutine deteriminePeelingOrder() 
        use globalGP, only: nMatingPairs, listOfParents
        use graphModule
        use IndividualModule, only: Individual
        use IndividualHelperModule
        use IndividualLinkedListModule
        
        implicit none
        
        integer, dimension(nMatingPairs, nMatingPairs) :: distance
        type(Individual) :: ind
        integer, dimension(:), allocatable :: targetPair, parents
        integer, dimension(nMatingPairs) :: individualIDs
        real(kind=real64), dimension(nMatingPairs) :: individualScore

        type(graph) :: g
        integer :: i, j
        logical :: related
        
        related = .false.

        distance = 0
        !print *, "setting up graph to decompose."
        do i = 1, nMatingPairs
            do j = 1, i-1
                !See if j is the parent of either individual in i.

                related = .false.

                targetPair = listOfParents(:,j)
                
                ind = pedigree%Pedigree(listOfParents(1, i))
                parents = ind%getIntegerVectorOfRecodedIds()
                related = related .or. (targetPair(1) == parents(2) .and. targetPair(2) == parents(3))
                ind = pedigree%Pedigree(listOfParents(2, i))
                parents = ind%getIntegerVectorOfRecodedIds()
                related = related .or. (targetPair(1) == parents(2) .and. targetPair(2) == parents(3))

                targetPair = listOfParents(:,i)
                
                ind = pedigree%Pedigree(listOfParents(1, j))
                parents = ind%getIntegerVectorOfRecodedIds()
                related = related .or. (targetPair(1) == parents(2) .and. targetPair(2) == parents(3))
                
                ind = pedigree%Pedigree(listOfParents(2, j))
                parents = ind%getIntegerVectorOfRecodedIds()
                related = related .or. (targetPair(1) == parents(2) .and. targetPair(2) == parents(3))
                
                if(related) then
                    distance(i, j) = 1
                    distance(j, i) = 1
                endif

            enddo
        enddo
        individualScore = 1
        individualIDs = [(i, i=1,nMatingPairs)]
        call g%initializeGraph(individualIDs, distance, individualScore)


    end subroutine

    subroutine setupPhaseChildOfFounders()
        use globalGP, only: pedigree, founders, phaseChildren, isPhasedChild
        use PedigreeModule
        use IndividualModule, only: Individual
        use IndividualHelperModule
        use IndividualLinkedListModule
        implicit  none
        type(individualPointerContainer), dimension(:), allocatable  :: children
        type(IndividualPointerContainer) :: child
        type(Individual) :: parent, childInd
        real(kind=real64), dimension(:), allocatable :: genotypeDensity, allGenotypeDensity
        real(kind=real64), dimension(4, 2) :: newArrayOptions
        integer i, j, founder, childId
        integer, dimension(:), allocatable :: maxl
        
        newArrayOptions(:, 1) = [0, 0, 1, 1]
        newArrayOptions(:, 2) = [0, 1, 0, 1]
        if(allocated(phaseChildren)) deallocate(phaseChildren)
        if(allocated(phaseChildrenOverride)) deallocate(phaseChildrenOverride)
        if(allocated(isPhasedChild)) deallocate(isPhasedChild)

        allocate(phaseChildren(size(founders)))
        allocate(phaseChildrenOverride(4, size(founders)))
        allocate(isPhasedChild(nAnimals))
        allGenotypeDensity = pedigree%getGenotypePercentage()
        isPhasedChild = .false.
        do i = 1, size(founders)
            founder = founders(i)
            parent = pedigree%pedigree(founder)
            children = parent%Offsprings
            if(parent%nOffs == 0) then
                phaseChildren(i) = 0
            else 
                allocate(genotypeDensity(parent%nOffs))
                do j = 1, parent%nOffs
                    genotypeDensity(j) = allGenotypeDensity(children(j)%p%id)
                enddo
                maxl = maxloc(genotypeDensity)
                child = children(maxl(1))
                childInd = child%p
                childId = childInd%id
                isPhasedChild(childId) = .true.
                phaseChildren(i) = childId
                phaseChildrenOverride(:,i) = newArrayOptions(:, parent%gender)
                deallocate(genotypeDensity)
            endif
        enddo

    end subroutine
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Sets up the trace tensor array
    !
    !> @details     Sets up the trace tensor array. This is currently fixed to handle 4 possible haplptypes. More work will be needed to extend it beyond that point.
    !
    !> @author     Andrew Whalen, awhalen@roslin.ed.ac.uk
    !
    !> @date       December 7, 2016
    !
    ! PARAMETERS:
    !> globalGP: traceTensor(i,j,k) probability of observing haplotypes in individual i, with father j, and mother k.
    !---------------------------------------------------------------------------
    
    subroutine setupTraceTensor
        use globalGP, only: traceTensorExp, nHaplotypes
        implicit none
        integer :: segregation, allele
        integer, dimension(4) :: fatherAlleleCoding, motherAlleleCoding
        real(kind=real64), dimension(0:1, nHaplotypes) :: paternalTransmission, maternalTransmission, father, mother
        real(kind=real64), dimension(nHaplotypes, nHaplotypes) :: zeroArray
        real(kind = real64) :: error
        real(kind=real64), dimension(nHaplotypes,nHaplotypes,nHaplotypes) :: traceTensor
        
        if(allocated(traceTensorExp)) deallocate(traceTensorExp)
        if(associated(segregationTensor)) segregationTensor => null()
        if(associated(segregationTensorParentsFirst)) segregationTensorParentsFirst => null()


        allocate(traceTensorExp(nHaplotypes,nHaplotypes,nHaplotypes))
        !Trace tensor is of the form "child, father, mother"
        !assume the haplotypes are of the form aa, aA, Aa, AA, in order Father,Mother
        !This will need some work for other haplotypes
        
        error = 1E-7 !Mutation rate
        paternalTransmission(0,:) = [1.-error, 1.-error, error, error] !likelihood of transmitting "a" given you are transmitting the paternal haplotype
        paternalTransmission(1,:) = 1. - paternalTransmission(0,:)
        
        maternalTransmission(0,:) = [1.-error, error, 1.-error, error]
        maternalTransmission(1,:) = 1 - maternalTransmission(0,:) 
        
        fatherAlleleCoding = [0, 0, 1, 1]
        motherAlleleCoding = [0, 1, 0, 1]

        !                  fm  fm  fm  fm 
        !segregationOrder: pp, pm, mp, mm

        allocate(segregationTensor(nHaplotypes,nHaplotypes,nHaplotypes, 4))
        allocate(segregationTensorParentsFirst(nHaplotypes,nHaplotypes,nHaplotypes, 4))
        segregationTensor = 0
        do segregation = 1, 4
            if(segregation == 1) then !pp
                father = paternalTransmission(:,:)
                mother = paternalTransmission(:,:)
            endif
            if(segregation == 2) then !pm
                father = paternalTransmission(:,:)
                mother = maternalTransmission(:,:)
            endif
            if(segregation == 3) then !mp
                father = maternalTransmission(:,:)
                mother = paternalTransmission(:,:)
            endif
            if(segregation == 4) then !mm
                father = maternalTransmission(:,:)
                mother = maternalTransmission(:,:)
            endif

            !alleles: aa, aA, Aa, AA
            do allele = 1, 4
                zeroArray = 0
                !ger as a rank one update of zero array.
                !likelihood that the father transmits the requested paternal alleles is: father(fatherAlleleCoding(allele),:)
                call ger(zeroArray, father(fatherAlleleCoding(allele),:), mother(motherAlleleCoding(allele),:)) 
                segregationTensor(allele, :, :, segregation) = zeroArray
                segregationTensorParentsFirst(:, :, allele, segregation) = zeroArray
            enddo
        enddo

        traceTensor = sum(segregationTensor, dim=4)/4.
        traceTensorExp = traceTensor


    end subroutine
    
    function calculateSegregationTransmissionMatrix(recombinationRate) result(segregationTransmissionMatrix)
        implicit none
        real(kind=real64), intent(in):: recombinationRate
        real(kind=real64) :: error
        real(kind=real64), dimension(0:1,4) :: paternalTransmission, maternalTransmission
        real(kind=real64), dimension(4,4) :: segregationTransmissionMatrix

        error = recombinationRate !recombinationRate
        paternalTransmission(0,:) = [1.-error/2, 1.-error/2, error/2, error/2] !Now this is structured error for the transmission rate thingys
        paternalTransmission(1,:) = 1. - paternalTransmission(0,:)
        
        maternalTransmission(0,:) = [1.-error/2, error/2, 1.-error/2, error/2]
        maternalTransmission(1,:) = 1. - maternalTransmission(0,:) 

        !We should also do recombination here...
        segregationTransmissionMatrix(1,:) = paternalTransmission(0,:)
        segregationTransmissionMatrix(2,:) = paternalTransmission(0,:)
        segregationTransmissionMatrix(3,:) = paternalTransmission(1,:)
        segregationTransmissionMatrix(4,:) = paternalTransmission(1,:)

        segregationTransmissionMatrix(:,1) = segregationTransmissionMatrix(:,1) * maternalTransmission(0,:)
        segregationTransmissionMatrix(:,2) = segregationTransmissionMatrix(:,2) * maternalTransmission(1,:)
        segregationTransmissionMatrix(:,3) = segregationTransmissionMatrix(:,3) * maternalTransmission(0,:)
        segregationTransmissionMatrix(:,4) = segregationTransmissionMatrix(:,4) * maternalTransmission(1,:)


    end function

!---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Collapses the trace matrix over a child's haplotype probabilities
    !
    !> @details    Collapses the trace matrix over a child's haplotype probabilities
    !
    !> @author     Andrew Whalen, awhalen@roslin.ed.ac.uk
    !
    !> @date       December 7, 2016
    !
    ! PARAMETERS:
    !> @param[out] 
    !---------------------------------------------------------------------------

    subroutine buildSegregationTraceTensor(markerEstimates)
        use globalGP, only : nAnimals, segregationTensor
        implicit none
        type(PeelingEstimates), pointer :: markerEstimates
        real(kind=real64), dimension(:,:), pointer :: currentSegregationTensors, currentSegregationEstimates, matSegTensor
        
        ! if(allocated(markerEstimates%currentSegregationTensors)) deallocate(markerEstimates%currentSegregationTensors)
        
        matSegTensor(1:64, 1:4) => segregationTensor
        currentSegregationEstimates => markerEstimates%currentSegregationEstimate !(1:4, nAnimals)
        currentSegregationTensors(1:64, 1:nAnimals) => markerEstimates%currentSegregationTensors

        call gemm(matSegTensor, currentSegregationEstimates, currentSegregationTensors)

    end subroutine



    subroutine buildNewSegregationTensors(markerEstimates, subset)
        use globalGP, only : nAnimals, segregationTensor
        implicit none
        type(PeelingEstimates), pointer :: markerEstimates
        real(kind=real64), dimension(:,:,:,:), pointer, contiguous :: subsetSegregationEstimate
        real(kind=real64), dimension(:,:), pointer :: currentSegregationTensors, matSegTensor
        real(kind=real64), dimension(:,:), allocatable :: currentSegregationEstimates
        integer, dimension(:) :: subset
        integer :: subsetSize
        subsetSize = size(subset)
        
        matSegTensor(1:64, 1:4) => segregationTensor
        ! currentSegregationEstimates = markerEstimates%currentSegregationEstimate(:, subset) !(1:4, nAnimals)
        currentSegregationEstimates = markerEstimates%currentSegregationEstimate(:, subset) !(1:4, nAnimals)

        allocate(subsetSegregationEstimate(4, 4, 4, subsetSize)) 
        currentSegregationTensors(1:64, 1:subsetSize) => subsetSegregationEstimate

        call gemm(matSegTensor, currentSegregationEstimates, currentSegregationTensors)

        markerEstimates%currentSegregationTensors(:,:,:,subset) = subsetSegregationEstimate
        deallocate(subsetSegregationEstimate)
    end subroutine


    subroutine getSegregationEstimate(markerEstimates, currentPeelingEstimates, indexNumber, runType)
        use globalGP, only : nAnimals, nSnps
        implicit none
        type(PeelingEstimates), pointer :: markerEstimates
        type(PeelingEstimates), dimension(:), pointer :: currentPeelingEstimates
        real(kind=real64), dimension(4,nAnimals) :: previousSegregation, nextSegregation, nextSegTransmitted, previousSegTransmitted
        real(kind=real64), dimension(4,4) :: prevTransmissionMatrix, nextTransmissionMatrix
        integer :: runType
        integer :: indexNumber

        if(indexNumber > 1 .and. (runType == 1 .or. runType == 3)) then        
            previousSegTransmitted = 0
            previousSegregation(:, :) = currentPeelingEstimates(indexNumber-1)%transmitForward
            prevTransmissionMatrix = calculateSegregationTransmissionMatrix(currentPeelingEstimates(indexNumber-1)%recombinationRate)
            call gemm(prevTransmissionMatrix, previousSegregation, previousSegTransmitted)
        else 
            previousSegTransmitted = 1 
        endif

        if(indexNumber < nSnps .and. (runType == 2 .or. runType == 3)) then
            nextSegTransmitted = 0
            nextSegregation(:, :) = currentPeelingEstimates(indexNumber+1)%transmitBackward
            nextTransmissionMatrix = calculateSegregationTransmissionMatrix(currentPeelingEstimates(indexNumber)%recombinationRate)
            call gemm(nextTransmissionMatrix, nextSegregation, nextSegTransmitted)
        else 
            nextSegTransmitted = 1 
        endif

        markerEstimates%currentSegregationEstimate = previousSegTransmitted * nextSegTransmitted


    end subroutine
    function reduceSegregationTensorMKL(childGenotypes, parentJointGenotypes) result(collapsedEstimate)
        use globalGP, only: nHaplotypes, segregationTensorParentsFirst
        implicit none
        real(kind=real64), dimension(nHaplotypes), intent(in) :: childGenotypes
        real(kind=real64), dimension(nHaplotypes, nHaplotypes), intent(in) :: parentJointGenotypes

        real(kind=real64), dimension(:,:), pointer, contiguous :: expJoint
        real(kind=real64), dimension(nHaplotypes) :: collapsedEstimate
        real(kind=real64), dimension(:,:), pointer :: reducedSegParentFirst, reorderedReduced
        real(kind=real64), dimension(:), pointer :: expJointReduced, reducedSegTensor
        
        real(kind=real64), dimension(nHaplotypes) :: expChildGenotypes
        
        allocate(expJoint(nHaplotypes, nHaplotypes))
        allocate(reducedSegTensor(nHaplotypes * 4))
        expChildGenotypes = exp(childGenotypes-maxval(childGenotypes))
        expJoint = exp(parentJointGenotypes-maxval(parentJointGenotypes))

        !Sort this out.
        reducedSegParentFirst(1:nHaplotypes**2, 1:nHaplotypes*4) => segregationTensorParentsFirst
        expJointReduced(1:nHaplotypes**2) => expJoint
        call gemv(reducedSegParentFirst, expJointReduced, reducedSegTensor, trans='T') !(h**2 x 4h)T * (h**2) => 4h
        !then reduce again
        reorderedReduced(1:nHaplotypes, 1:4) => reducedSegTensor
        call gemv(reorderedReduced, expChildGenotypes, collapsedEstimate, trans='T') !(h x 4)T * (h) => 4
        
        collapsedEstimate = collapsedEstimate/sum(collapsedEstimate)

        deallocate(expJoint, reducedSegTensor)
      
    end function

    function reduceSegregationTensor(childGenotypes, parentJointGenotypes) result(collapsedEstimate)
        use globalGP, only: nHaplotypes
        implicit none
        real(kind=real64), dimension(nHaplotypes), intent(in) :: childGenotypes
        real(kind=real64), dimension(nHaplotypes, nHaplotypes), intent(in) :: parentJointGenotypes

        real(kind=real64), dimension(:,:), allocatable :: expJoint
        real(kind=real64), dimension(nHaplotypes) :: collapsedEstimate
        
        real(kind=real64), dimension(nHaplotypes) :: expChildGenotypes
        integer :: seg, allele

        ! for each seg option pp, pm, mp, mm
        expChildGenotypes = exp(childGenotypes-maxval(childGenotypes))
        expJoint = exp(parentJointGenotypes-maxval(parentJointGenotypes))
        collapsedEstimate = 0
        do seg=1, 4
            !Do for each child allele, aa, aA, Aa, AA
			!!!$OMP DO SIMD REDUCTION(+:collapsedEstimate)
            do allele = 1, 4
                collapsedEstimate(seg) = collapsedEstimate(seg) + sum(segregationTensor(allele, :, :, seg)*expJoint) * expChildGenotypes(allele)
            enddo
        enddo

    end function

    function childTraceMultiply(vect, traceTensorExp) result(collapsedEstimate)
        use globalGP, only: nHaplotypes
        implicit none
        real(kind=real64), dimension(nHaplotypes), intent(in) :: vect
        real(kind=real64), dimension(nHaplotypes) :: vectExp
        real(kind=real64), dimension(nHaplotypes, nHaplotypes) :: collapsedEstimate
        real(kind=real64), dimension(nHaplotypes, nHaplotypes, nHaplotypes), intent(in) :: traceTensorExp
        integer :: i, j
        real(kind=real64)  :: max
        
        max = maxval(vect)
        vectExp = vect - max
        vectExp = exp(vectExp)

        !we need to rebuild the trace tensor
        ! traceTensorExp = buildTraceTensor(childSegregation)
        do i=1,nHaplotypes
			!!!$OMP DO SIMD
           do j = 1, nHaplotypes
                collapsedEstimate(i, j) = log(sum(vectExp * traceTensorExp(:, i, j))) + max
           enddo
        enddo
        
    end function

    function childTraceMultiplyMKL(vect, traceTensorLocal) result(collapsedEstimate)
        use globalGP, only: nHaplotypes
        implicit none
        real(kind=real64), dimension(nHaplotypes), intent(in) :: vect
        real(kind=real64), dimension(nHaplotypes) :: vectExp
        real(kind=real64), dimension(:,:), allocatable :: collapsedEstimate
        ! real(kind=real64), dimension(:,:,:), intent(in) :: traceTensorExp
        real(kind=real64), dimension(:,:,:), pointer, contiguous :: traceTensorLocal

        real(kind=real64), dimension(:,:), pointer :: reducedTraceTensor
        real(kind=real64), dimension(:), allocatable :: reducedCollapsedEstimate

        real(kind=real64)  :: max
        
        allocate(collapsedEstimate(nHaplotypes, nHaplotypes))
        allocate(reducedCollapsedEstimate(1:nHaplotypes**2))

        max = maxval(vect)
        vectExp = exp(vect - max)

        reducedTraceTensor(1:nHaplotypes, 1:nHaplotypes**2) => traceTensorLocal

        call gemv(reducedTraceTensor, vectExp, reducedCollapsedEstimate, trans = "T") !(hxh**2) x h => h**2

        collapsedEstimate = reshape(reducedCollapsedEstimate, [nHaplotypes,nHaplotypes])
        collapsedEstimate = log(collapsedEstimate)

        deallocate(reducedCollapsedEstimate)
    end function
   
    function parentTraceMultiply(parentMat, traceTensorExp) result(collapsedEstimate)
        use globalGP, only: nHaplotypes
        implicit none
        real(kind=real64), dimension(nHaplotypes, nHaplotypes), intent(in) :: parentMat
        ! real(kind=real64), dimension(nHaplotypes), intent(in) :: childSegregation
        real(kind=real64), dimension(nHaplotypes) :: collapsedEstimate
        real(kind=real64), dimension(nHaplotypes, nHaplotypes) :: parentMatExp
        real(kind=real64), dimension(nHaplotypes, nHaplotypes, nHaplotypes), intent(in) :: traceTensorExp
        integer :: i
        real(kind=real64)  :: max

        max = maxval(parentMat)
        parentMatExp = parentMat - max
        parentMatExp = exp(parentMatExp)
        !!$OMP DO SIMD
        do i=1,nHaplotypes
            collapsedEstimate(i) = sum(parentMatExp * traceTensorExp(i,:,:))
        enddo

        collapsedEstimate = log(collapsedEstimate)
    end function


    function buildTraceTensor(childEstimate2) result(traceTensor)
        use globalGP, only: nHaplotypes, segregationTensor
        implicit none
        real(kind=real64), dimension(:,:,:), pointer, contiguous :: traceTensor
        ! real(kind=real64), dimension(:,:,:), pointer, contiguous :: traceTensorOut
        real(kind=real64), dimension(4), intent(in) :: childEstimate2
        real(kind=real64), dimension(4) :: childEstimate
        integer :: i
        allocate(traceTensor(nHaplotypes,nHaplotypes,nHaplotypes))
        traceTensor=0
        ! traceTensor = 0
        childEstimate = childEstimate2/sum(childEstimate2)
        do i = 1, 4
            traceTensor = traceTensor + childEstimate(i)*segregationTensor(:,:,:, i)
        enddo

    end function

    function jointByMateToFather(vect, joint) result(collapsedEstimate)
        use globalGP, only: nHaplotypes
        implicit none
        real(kind=real64), dimension(nHaplotypes), intent(in) :: vect
        real(kind=real64), dimension(nHaplotypes, nHaplotypes), intent(in) :: joint
        real(kind=real64), dimension(nHaplotypes) :: collapsedEstimate
        integer :: i
        ! joint is of the form (or should be of the form) father, mother
        !MP
        do i=1,nHaplotypes
            collapsedEstimate(i) = logAdd1t0(vect + joint(i, :))
        enddo

    end function
    function jointByFatherToMate(vect, joint) result(collapsedEstimate)
        use globalGP, only: nHaplotypes
        implicit none
        real(kind=real64), dimension(nHaplotypes), intent(in) :: vect
        real(kind=real64), dimension(nHaplotypes, nHaplotypes), intent(in) :: joint
        real(kind=real64), dimension(nHaplotypes) :: collapsedEstimate
        integer :: i
        ! joint is of the form (or should be of the form) father, mother
        !MP
		!!$OMP DO SIMD
        do i=1,nHaplotypes
            collapsedEstimate(i) = logAdd1t0(vect + joint(:, i))
        enddo

    end function

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Calculates genotype estimates
    !
    !> @details     Calculates genotype estimates
    !
    !> @author     Andrew Whalen, awhalen@roslin.ed.ac.uk
    !
    !> @date       December 7, 2016
    !
    ! PARAMETERS:
    !> @param[out] genotypeEstimate 
    !---------------------------------------------------------------------------
    
    subroutine updateMafEstimates(genotypeEstimates, markerEstimates)
        use globalGP, only: pedigree
        use fixedPointModule
        implicit none
        type(peelingEstimates), intent(inout) :: markerEstimates
        real(kind=real64), dimension(:), intent(in) :: genotypeEstimates
        real(kind=real64) :: observedMaf
        type(fixedPointEstimator), pointer :: currentMafEstimator


        observedMaf = (sum(genotypeEstimates(founders))+3)/(2*(pedigree%founders%length+3))
        currentMafEstimator => markerEstimates%mafEstimator
        call currentMafEstimator%addObservation(markerEstimates%maf, observedMaf)
        markerEstimates%maf = currentMafEstimator%secantEstimate(isLogitIn=.true.)

    end subroutine


    subroutine updateRecombinationRate(markerEstimates, currentPeelingEstimates, indexNumber)
        use globalGP, only:  nAnimals
        use fixedPointModule
        implicit none
        type(peelingEstimates) :: markerEstimates
        type(PeelingEstimates), dimension(:), pointer :: currentPeelingEstimates

        real(kind=real64), dimension(4,nAnimals) :: zi, plusZeroSegregation, plusOneSegregation
        real(kind=real64), dimension(4) :: normalizedEstimate
        real(kind=real64), dimension(4,4) :: transmissionMatrix, estimate
        real(kind=real64) :: nChanges, nObservations, observedChangeRate
        real(kind=real64) :: diagonalSum
        ! type(fixedPointEstimator), pointer :: currentRecombinationEstimator
        integer i, j, indexNumber, ind

        !We update n to n+1 here.
        transmissionMatrix = calculateSegregationTransmissionMatrix(markerEstimates%recombinationRate)


        nChanges = 1d0/nSnpsAll
        nObservations = 1.0

        if(indexNumber < nSnps) then 
            nObservations = nObservations + nAnimals
            plusZeroSegregation = markerEstimates%transmitForward
            plusOneSegregation = currentPeelingEstimates(indexNumber+1)%transmitBackward
            do ind = 1, nAnimals
                estimate = transmissionMatrix
                do j = 1, 4
                    estimate(j, :) = estimate(j,:) * plusZeroSegregation(:, ind)
                enddo
                do j = 1, 4
                    estimate(:,j) = estimate(:,j) * plusOneSegregation(:, ind)
                enddo
                estimate = estimate/sum(estimate)
                diagonalSum = 0
                do j = 1, 4
                    diagonalSum = diagonalSum + estimate(j, j)
                enddo
                nChanges = nChanges + 1-diagonalSum
            enddo

        endif
        observedChangeRate = nChanges/nObservations

        ! currentRecombinationEstimator => markerEstimates%recombinationEstimator
        ! call currentRecombinationEstimator%addObservation(markerEstimates%recombinationRate, observedChangeRate)
        ! markerEstimates%recombinationRate = currentRecombinationEstimator%secantEstimate(isLogitIn=.true.)
        markerEstimates%recombinationRate = observedChangeRate
    

    end subroutine



    ! subroutine updateAllRecombinationRates(currentPeelingEstimates)
    !     use globalGP, only:  nAnimals, nSnps
    !     use fixedPointModule
    !     implicit none
    !     type(peelingEstimates), pointer :: markerEstimates
    !     type(PeelingEstimates), dimension(:), pointer :: currentPeelingEstimates

    !     real(kind=real64), dimension(4,nAnimals) :: zi, plusZeroSegregation, plusOneSegregation
    !     real(kind=real64), dimension(4) :: normalizedEstimate
    !     real(kind=real64), dimension(4,4) :: transmissionMatrix
    !     real(kind=real64) :: nChanges, nObservations, observedChangeRate
    !     ! type(fixedPointEstimator), pointer :: currentRecombinationEstimator
    !     integer i, j, indexNumber

    !     !We update n to n+1 here.
    !     markerEstimates => currentPeelingEstimates(1)
    !     !print *, markerEstimates%recombinationRate


    !     transmissionMatrix = calculateSegregationTransmissionMatrix(markerEstimates%recombinationRate)
    !     !print *, transmissionMatrix

    !     nChanges = .01
    !     nObservations = 1.0
    !     do indexNumber = 1, nSnps-1
    !         if(indexNumber < nSnps) then 
    !             nObservations = nObservations + nAnimals
    !             plusZeroSegregation = currentPeelingEstimates(indexNumber)%fullSegregation
    !             plusOneSegregation = currentPeelingEstimates(indexNumber+1)%fullSegregation

    !             zi = 0
    !             do i = 1, nAnimals
    !                 do j= 1, 4
    !                     normalizedEstimate = plusZeroSegregation(:,i) * transmissionMatrix(:,j)
    !                     normalizedEstimate = normalizedEstimate/sum(normalizedEstimate)
    !                     zi(j, i) = zi(j, i) + plusOneSegregation(j,i) * (1-normalizedEstimate(j))
    !                 enddo
    !             enddo
    !             nChanges = nChanges + sum(zi)
    !         endif
    !     enddo
    !     ! !print *, "obs", nChanges, nObservations
    !     observedChangeRate = nChanges/nObservations
    !     !print *, indexNumber, " ", nChanges, " ", nObservations
    !     !print *, indexNumber, " ", observedChangeRate
        
    !     ! currentRecombinationEstimator => markerEstimates%recombinationEstimator
    !     ! call currentRecombinationEstimator%addObservation(markerEstimates%recombinationRate, observedChangeRate)
    !     ! observedChangeRate = currentRecombinationEstimator%secantEstimate(isLogitIn=.true.)

    !     !print *, "change rate", observedChangeRate

    !     do indexNumber = 1, nSnps
    !         currentPeelingEstimates(indexNumber)%recombinationRate = observedChangeRate
    !     enddo

    ! end subroutine


    subroutine writeOutputsToFileMultiLocus(currentPeelingEstimates)
        use globalGP
        implicit none
        type(peelingEstimates), pointer :: markerEstimates
        type(peelingEstimates), dimension(:), pointer :: currentPeelingEstimates
        real(kind=real64), dimension(:,:,:), allocatable :: combinedHaplotypes
        real(kind=real64), dimension(:,:), allocatable :: combinedGenotypes
        integer, dimension(nSnps) :: individualGenotype
        real(kind=real64), dimension(:), allocatable :: thresholds
        CHARACTER(LEN=30) :: rowfmt
        integer :: haplotypeFile, dosageFile, segregationFile, paramaterFile
        integer :: i, j, tmp


        open(newunit = haplotypeFile, FILE = trim(inputParams%prefix) // ".haps", status="replace")
        open(newunit = dosageFile, FILE = trim(inputParams%prefix) // ".dosages", status="replace")
        open(newunit = segregationFile, FILE = trim(inputParams%prefix) // ".seg", status="replace")
        open(newunit = paramaterFile, FILE = trim(inputParams%prefix) // ".params", status="replace")
    
        !print *, "Writting outputs"
        write(paramaterFile, '(a)') "maf ", "genotypeError ", "recombinationRate"
        do i = 1, nSnps
            markerEstimates => currentPeelingEstimates(i)
            write(paramaterFile, '(3f12.7)') markerEstimates%maf, markerEstimates%genotypingErrorRate, markerEstimates%recombinationRate
        enddo

        allocate(combinedHaplotypes(4, nSnps, nAnimals))


        do i = 1, nSnps
            markerEstimates => currentPeelingEstimates(i)
            combinedHaplotypes(:, i, :) = markerEstimates%haplotypeEstimates(:,:)
        enddo

        WRITE(rowfmt,'(A,I9,A)') '(a,',nSnps+10,'f10.4)'
        !print *, "row format", rowfmt
        do i = 1, nAnimals
            tmp = haplotypeFile
            write(tmp, rowfmt) pedigree%pedigree(i)%originalID, combinedHaplotypes(1,:, i)
            write(tmp, rowfmt) pedigree%pedigree(i)%originalID, combinedHaplotypes(2,:, i)
            write(tmp, rowfmt) pedigree%pedigree(i)%originalID, combinedHaplotypes(3,:, i)
            write(tmp, rowfmt) pedigree%pedigree(i)%originalID, combinedHaplotypes(4,:, i)
            ! write(outputFile(index),'(a)') " " 
        enddo
        if(allocated(inputParams%thresholds)) then
            thresholds = inputParams%thresholds
            do i = 1, size(thresholds)
                call callAlleles(thresholds(i), combinedHaplotypes)
            enddo
        endif


        deallocate(combinedHaplotypes)
        allocate(combinedHaplotypes(4, nSnps, nAnimals))

        !Output Segregation
        do i = 1, nSnps
            markerEstimates => currentPeelingEstimates(i)
            combinedHaplotypes(:, i, :) = markerEstimates%fullSegregation(:,:)
        enddo

        WRITE(rowfmt,'(A,I9,A)') '(a,',nSnps+10,'f10.4)'
        !print *, "row format", rowfmt
        do i = 1, nAnimals
            tmp = segregationFile
            write(tmp, rowfmt) pedigree%pedigree(i)%originalID, combinedHaplotypes(1,:, i)
            write(tmp, rowfmt) pedigree%pedigree(i)%originalID, combinedHaplotypes(2,:, i)
            write(tmp, rowfmt) pedigree%pedigree(i)%originalID, combinedHaplotypes(3,:, i)
            write(tmp, rowfmt) pedigree%pedigree(i)%originalID, combinedHaplotypes(4,:, i)
            ! write(outputFile(index),'(a)') " " 
        enddo
        deallocate(combinedHaplotypes)

        allocate(combinedGenotypes(nSnps, nAnimals))
        do i = 1, nSnps
            markerEstimates => currentPeelingEstimates(i)
            combinedGenotypes(i, :) = markerEstimates%genotypeEstimates(:)
        enddo

        do i = 1, nAnimals
            write(dosageFile, rowfmt) pedigree%pedigree(i)%originalID, combinedGenotypes(:, i)
        enddo

      
    end subroutine

    subroutine writeOutputsToFileSingleLocus(combinedHaplotypes, combinedGenotypes, markerError, maf)
        use globalGP
        implicit none
        real(kind=real64), dimension(:,:,:), allocatable :: combinedHaplotypes
        real(kind=real64), dimension(:,:), pointer :: combinedGenotypes
        real(kind=real64), dimension(:), intent(inout) :: markerError, maf
        real(kind=real64), dimension(:), allocatable :: thresholds
        CHARACTER(LEN=128) :: rowfmt
        integer, dimension(nSnps) :: individualGenotype
        integer :: haplotypeFile, dosageFile, paramaterFile

        integer :: i, j, tmp


        !print *, "Writting outputs"

        if(inputParams%writeHaps) open(newunit = haplotypeFile, FILE = trim(inputParams%prefix) // ".haps", status="replace")
        if(inputParams%writeDosages) open(newunit = dosageFile, FILE = trim(inputParams%prefix) // ".dosages", status="replace")
        if(inputParams%writeParams) open(newunit = paramaterFile, FILE = trim(inputParams%prefix) // ".params", status="replace")

        ! write(paramaterFile, '(a)') "maf ", "gError "
        if(inputParams%writeParams) then
            do i = 1, nSnps
                write(paramaterFile, '(3f12.7)') maf(i), markerError(i)
            enddo
        endif
        WRITE(rowfmt,'(A,I9,A)') '(a,',nSnps+10,'f10.4)'

        if(inputParams%writeHaps) then
            !print *, "row format", rowfmt
            do i = 1, nAnimals
                tmp = haplotypeFile
                write(tmp, rowfmt) pedigree%pedigree(i)%originalID, combinedHaplotypes(1,:, i)
                write(tmp, rowfmt) pedigree%pedigree(i)%originalID, combinedHaplotypes(2,:, i)
                write(tmp, rowfmt) pedigree%pedigree(i)%originalID, combinedHaplotypes(3,:, i)
                write(tmp, rowfmt) pedigree%pedigree(i)%originalID, combinedHaplotypes(4,:, i)
                ! write(outputFile(index),'(a)') " " 
            enddo
        endif
        !print *, "row format", rowfmt
        !print *, size(combinedGenotypes, 1), size(combinedGenotypes, 2)
        ! do i = 1, nAnimals
        !     !print *,  pedigree%pedigree(i)%originalID, combinedGenotypes(1:50, i)
        ! enddo
        if(inputParams%writeDosages) then
            do i = 1, nAnimals
                write(dosageFile, rowfmt) pedigree%pedigree(i)%originalID, combinedGenotypes(:, i)
            enddo
        endif      

        !Output based on threshold
        if(allocated(inputParams%thresholds)) then
            thresholds = inputParams%thresholds
            do i = 1, size(thresholds)
                call callAlleles(thresholds(i), combinedHaplotypes)
            enddo
        endif

    end subroutine

    subroutine callAlleles(threshold, combinedHaplotypes) 
        use globalGP
        use CompatibilityModule
        implicit none
        real(kind=real64), dimension(:,:,:), allocatable :: combinedHaplotypes
        real(kind=real64) :: threshold, currentVal, score
        CHARACTER(LEN=128) :: rowfmt, fileName
        integer(kind=1),allocatable,dimension (:,:) :: individualGenotype
        integer :: i, j, tmp, genoCall

        allocate(individualGenotype(nAnimals, nSnps))

        write(fileName, '(F8.5)') threshold
        fileName =  trim(inputParams%prefix) // 'pointGenotypes_' // trim(adjustl(fileName)) // '.txt'
        ! open(newunit = tmp, FILE = trim(fileName), status="replace")

        ! WRITE(rowfmt,'(A,I9,A)') '(a,',nSnps+10,'I2)'

        individualGenotype = 9
        do i = 1, nAnimals
            do j = 1, nSnps
                genoCall = 9
                currentVal = -1
                score = combinedHaplotypes(1, j, i)
                if(score > threshold .and. score > currentVal) genoCall = 0
                
                score = combinedHaplotypes(2, j, i) + combinedHaplotypes(3, j, i)
                if(score > threshold .and. score > currentVal) genoCall = 1
                
                score = combinedHaplotypes(4, j, i)
                if(score > threshold .and. score > currentVal) genoCall = 2

                individualGenotype(i, j) = genoCall
            enddo


            ! write(tmp, rowfmt) pedigree%pedigree(i)%originalID, individualGenotype
        enddo
        call pedigree%addGenotypeInformationFromArray(individualGenotype)

        if(inputParams%binaryCallOutput) then
            call writeOutPlinkBinary(pedigree, fileName, startPosIn= inputParams%startSnp)
        else
            call pedigree%writeOutGenotypes(fileName)
        endif
        ! close(tmp)
    end subroutine


    subroutine additiveOuterProduct(pfather, pmate, pjoint)
        use globalGP, only: nHaplotypes
        implicit none
        REAL(kind=real64), dimension(nHaplotypes), intent(in) :: pfather, pmate
        REAL(kind=real64), dimension(nHaplotypes, nHaplotypes) :: pjoint
        integer i
        do i = 1, nHaplotypes
            pjoint(:,i) = pfather(i) + pmate
        enddo
    end subroutine

    subroutine additiveOuterProductSpread(pfather, pmate, pjoint)
        use globalGP, only: nHaplotypes
        implicit none
        REAL(kind=real64), dimension(nHaplotypes), intent(in) :: pfather, pmate
        REAL(kind=real64), dimension(nHaplotypes, nHaplotypes) :: pjoint

        pjoint = spread(pfather, 1, nHaplotypes) + spread(pmate, 2, nHaplotypes)
    end subroutine

    function logAdd1t0(mat) result(out)
        implicit none
        REAL(kind=real64), dimension(:), intent(in) :: mat
        REAL(kind=real64) :: out, vmax
        vmax = maxval(mat)
        out = vmax + log( sum(exp(mat-vmax)) )
    end function


    subroutine ppa2(a)
        implicit none
        integer, dimension(:, :), intent(in) :: a
        integer i
        !print *, " "
        do i = 1, size(a, 2)
            write (*,"(200f35.30)",advance="no") a(:, i)
            !print *, ""
        enddo
    end subroutine

    subroutine ppa3(a)
        implicit none
        real(kind=real64), dimension(:, :,:), intent(in) :: a
        integer i,j
        do j = 1, size(a, 3)
            do i = 1, size(a, 2)
                !print *, a(:, i,j)
            enddo
            !print *, ""
        enddo
    end subroutine

    subroutine ppa4(a)
        implicit none
        real(kind=real64), dimension(:, :,:,:), intent(in) :: a
        integer i,j,k
        do k = 1, size(a, 4)
            do j = 1, size(a, 3)
                do i = 1, size(a, 2)
                    !print *, nint(a(:, i,j,k))
                enddo
                !print *, ""
            enddo
            !print *, "----"
        enddo
    end subroutine

    subroutine printTrace(a)
        implicit none
        real(kind=real64), dimension(:, :,:), intent(in) :: a
        integer i,j
        do i = 1, size(a, 1)
            do j = 1, size(a, 2)
                write (*,"(200f12.7)",advance="no") a(i, j, :)
                !print *, ""
            enddo
            !print *, ""
        enddo
    end subroutine

    subroutine printSegregation(a)
        implicit none
        real(kind=real64), dimension(:, :,:,:), allocatable, intent(in) :: a
        integer i,j,k
        do i = 1, size(a, 1)
            do j = 1, size(a, 2)
                do k = 1, size(a, 3)
                    !print *, nint(a(i, j, k, :))
                    ! !print *, (a(j, :,i,k))
                enddo
                !print *, ""
            enddo
            !print *, "----"
        enddo
    end subroutine


end module

