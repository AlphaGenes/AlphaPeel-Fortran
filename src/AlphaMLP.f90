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
        module procedure runAlphaMLPAlphaImpute
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

    use globalGP, only :pedigree, nsnps
        
        integer, intent(in) :: startSnp, endSnp
        integer :: i, j
        type(PedigreeHolder) :: ped
        real(kind=real64), dimension(:,:,:), allocatable, intent(out) :: AlphaMLPOutput !< output 3 dimensional array as requeuired by alphaimpute. 1:pedigree%pedigreeSize, nSnps, nHaplotypes
        real(kind=real64),allocatable,dimension (:), intent(out) :: Maf !< double vector containing MaF for each Snp
        type(AlphaMLPInput) :: inputParams
        integer :: nHaplotypes

        type(peelingEstimates), dimension(:), pointer:: currentPeelingEstimates

        pedigree = ped
        nHaplotypes = 4

       

        if (.not. allocated(AlphaMLPOutput)) then
            allocate(AlphaMLPOutput(pedigree%pedigreesize, startSnp:endSnp, 4))
        endif
        call pedigree%getMatePairsAndOffspring(offspringList, listOfParents, nMatingPairs)
        inputParams = AlphaMLPInput(startSnp,endSnp)
        
        call setupPhaseChildOfFounders()
        
        call setupTraceTensor
        nSnps = inputParams%endSnp-inputParams%startSnp+1
 
        nAnimals = pedigree%pedigreeSize
        founders = pedigree%founders%convertToArrayIDs()


        call runMultiLocusAlphaMLP(currentPeelingEstimates)
        do i = startSnp, endSnp
            do j = 1, 4
                AlphaMLPOutput(:, i, j) =  currentPeelingEstimates(i)%haplotypeEstimates(j,:)
            enddo
            maf(i) = currentPeelingEstimates(i)%maf
        enddo

        deallocate(currentPeelingEstimates)

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
        use globalGP, only: nSnps, sequenceData, inputParams
        implicit none
        integer :: i
        character(len=1), dimension(2) :: x
        type(peelingEstimates), dimension(:), pointer :: currentPeelingEstimates


        x = ["1", "2"]
        inputParams = AlphaMLPInput()        
        print *, "setup Pedigree"

        call setupPedigree(inputParams)
        call pedigree%getMatePairsAndOffspring(offspringList, listOfParents, nMatingPairs)

        print *, "read genotypes"

        if(inputParams%isSequence) call readSequence(inputParams, pedigree, sequenceData)


        print *, "setup founder phasing"

        call setupPhaseChildOfFounders()

        ! call deteriminePeelingOrder()

        print *, "setup trace"

        call setupTraceTensor
        !For each allele, run gene prob on that index.

        do i = 1, 2
            open(newunit = outputFile(i), FILE = "multiLocusHaplotypes" // x(i) // ".txt", status="replace", access="append")
            open(newunit = auxFile(i), FILE = "multiLocusGenotypesNoPhase" // x(i) // ".txt", status="replace", access="append")
            open(newunit = segregationFile(i), FILE = "segregationEstimates" // x(i) //".txt", status="replace", access="append")
            open(newunit = consensusFile(i), FILE = "pointGenotypes" // x(i) //".txt", status="replace", access="append")
            open(newunit = paramaterFile(i), FILE = "paramaterEstimates" // x(i) // ".txt", status="replace", access="append")
        enddo
        nSnps = inputParams%endSnp-inputParams%startSnp+1
        print *, "run AlphaMLP"
        call setupGenerations()
        call runMultiLocusAlphaMLP(currentPeelingEstimates, 1)

        deallocate(currentPeelingEstimates)
    end subroutine runAlphaMLPIndependently

    ! Two subroutines below for running the HMM.


    subroutine setupGenerations() 
        use globalGP, only: nMatingPairs, familiesInGeneration, offspringList, nGenerations
        use graphModule, only : selectIndexesBasedOnMask
        integer, dimension(nMatingPairs) :: familyGeneration
        integer :: i

        do i = 1, nMatingPairs
            familyGeneration(i) = offspringList(i)%first%item%generation
        enddo 

        nGenerations = maxval(familyGeneration)
        allocate(familiesInGeneration(nGenerations))
        do i = 1, nGenerations
            familiesInGeneration(i)%array = selectIndexesBasedOnMask(familyGeneration == i)
            print *, i, familiesInGeneration(i)%array
        enddo
        print *, nGenerations
    end subroutine

    subroutine setupHMM()
        use Global
        use AlphaImputeModule
        use informationModule
        use GlobalVariablesHmmMaCH
        use Output
        use AlphaImputeInMod
        use Imputation
        use InputMod
        use globalGP, only : nPseudoFounders
        implicit none
        character(len=4096) :: SpecFile
        specfile="AlphaImputeSpec.txt"

        allocate(defaultInput)
        call defaultInput%ReadInParameterFile(SpecFile)
        inputParams => defaultInput
        call CountInData
        call ReadInData
        print *, "Loaded Input:", defaultInput%nSnp
        allocate(inputGenotypeProbabilitiesFromGeneprob(4, nPseudoFounders, defaultInput%nSnp))
        call MakeFiles

        if(allocated(ImputeGenos)) deallocate(ImputeGenos)
        if(allocated(ImputePhase)) deallocate(ImputePhase)
        allocate(ImputeGenos(1:nPseudoFounders,inputParams%nsnp))
        allocate(ImputePhase(1:nPseudoFounders,inputParams%nsnp,2))
        ImputePhase=9
        ImputeGenos=9
        inputParams%HMMOption = RUN_HMM_GeneProb
        inputParams%nroundshmm = 20
        inputParams%hmmburninround = 5
        call ped%addGenotypeInformation(inputParams%GenotypeFile,inputParams%nsnp,NanisG)
        print *, "finished adding genotypes"
    end subroutine

    subroutine runHMMOnFounders(genotypes, hmmEstimate)
        use Global
        use GlobalVariablesHmmMaCH
        use Imputation
        implicit none
        real(kind=real64), dimension(:,:,:), intent(inout):: hmmEstimate
        integer(kind=1), dimension(:,:), allocatable, intent(in) :: genotypes
        
        print *, "Really running HMM"

        print *, size(ImputeGenos, 1), size(ImputeGenos,2)
        print *, size(genotypes, 1), size(genotypes,2)
        ImputeGenos(1:nPseudoFounders,:) = genotypes

        print *, "Really, really running HMM"
        call MaCHController(RUN_HMM_GeneProb)
    
        print *, size(realGenosCounts, 1), size(realGenosCounts, 2), size(realGenosCounts, 3)
        print *, size(hmmEstimate, 1),size(hmmEstimate, 2), size(hmmEstimate ,3)
        hmmEstimate(1,:,:) = 1. - realGenosCounts(:,:, 1) - realGenosCounts(:,:, 2)
        hmmEstimate(2,:,:) = realGenosCounts(:,:, 1)
        hmmEstimate(3,:,:) = realGenosCounts(:,:, 2)

        print *, hmmEstimate(1:3, 1:5, 1:5)

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
        use globalGP, only: pedigree, nAnimals, founders, sequenceData
        use PedigreeModule
        use IndividualModule, only: Individual
        implicit  none
        type(AlphaMLPInput),intent(in) :: inputParams        

        if (inputParams%pedFile /= "No Pedigree") then
            pedigree = PedigreeHolder(inputParams%pedFile)
            
            if (.not. inputParams%isSequence) then 
                call pedigree%addGenotypeInformationFromFile(inputParams%inputFile,inputParams%nsnp)     
            else 
                call pedigree%addSequenceFromFile(inputParams%sequenceFile, inputParams%nsnp, inputParams%nGenotypedAnimals, sequenceData)
            endif
        else
            ! Init pedigree with format of genotype file
            ! assume old pedigree file
            if (.not. inputParams%isSequence) then
            pedigree = PedigreeHolder(inputParams%inputFile)
                call readGenotypes(inputParams, pedigree)
            else 
                call readSequence(inputParams, pedigree, sequenceData)
            endif
        endif
        nAnimals = pedigree%pedigreeSize
        founders = pedigree%founders%convertToArrayIDs()
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
        print *, "setting up graph to decompose."
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

        print *, "decomposingGraph."

        print *, "decomposition order", g%decomposeGraph()


    end subroutine

    subroutine setupPseudoFounders()
        use globalGP, only: pedigree, pseudoFounders, nPseudoFounders
        use PedigreeModule
        use IndividualModule, only: Individual
        use IndividualHelperModule
        use IndividualLinkedListModule
        use graphModule
        implicit  none
        integer :: numberOfGenerations, i
        type(IndividualLinkedList) :: genotypedFounders
        integer, dimension(:), allocatable :: individualIDs

        numberOfGenerations = 1
        genotypedFounders = pedigree%getGenotypedFounders(numberOfGenerations)        
        individualIDs = genotypedFounders%convertToArrayIDs()
        nPseudoFounders = size(individualIDs)
        pseudoFounders = individualIDs

        allocate(isPseudoFounder(nAnimals))
        isPseudoFounder = .false.
        do i = 1, nPseudoFounders
            isPseudoFounder(pseudoFounders(i)) = .true.
        enddo


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
        allocate(traceTensorExp(nHaplotypes,nHaplotypes,nHaplotypes))
        !Trace tensor is of the form "child, father, mother"
        !assume the haplotypes are of the form aa, aA, Aa, AA, in order Father,Mother
        !This will need some work for other haplotypes
        
        error = .0001 !Mutation rate
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

    subroutine runMultiLocusAlphaMLP(currentPeelingEstimates,writeOutputs)
        use globalGP, only: nAnimals, nHaplotypes, nMatingPairs, nSnps, nPseudoFounders, inputParams
        implicit none
        type(peelingEstimates), dimension(:), pointer, intent(out) :: currentPeelingEstimates
        type(individual) :: tmpInd
        integer(kind=1), dimension(:,:), allocatable :: genotypes
        integer, optional :: writeOutputs
        integer :: i, roundNumber, nCycles, cycleIndex, nRounds
        logical :: converged, firstRun
        real(kind=real64), dimension(:,:,:), allocatable :: hmmEstimate

        allocate(currentPeelingEstimates(nSnps))
        do i = 1, nSnps
            call currentPeelingEstimates(i)%initializePeelingEstimates(nHaplotypes, nAnimals, nMatingPairs)
        enddo


        converged = .false.
        firstRun = .true.

        nRounds = 1
        if(inputParams%useHMM) then
            print *, "Using HMM"
            nRounds = 2 !No HMM: nRounds = 1; HMM: nRounds = 2
            call setupPseudoFounders()
        endif
        nCycles = inputParams%nCycles

        do roundNumber = 1, nRounds
            print *, "Round ", roundNumber

            if(.not. firstRun .and. inputParams%useHMM) then
                !Handle the HMM
                call setupHMM()
                allocate(hmmEstimate(3,nPseudoFounders,nSnps))
                allocate(genotypes(nPseudoFounders, nSnps))
                do i = 1, nPseudoFounders
                    tmpInd = pedigree%pedigree(i)
                    genotypes(i,:) = tmpInd%individualGenotype%toIntegerArray()
                enddo
                !Pass a subset of genotypes (only psuedoFounders) and use those genotypes to phase.
                call runHMMOnFounders(genotypes, hmmEstimate)
                do i = 1, nSnps
                    currentPeelingEstimates(i)%postHMM = .true.
                    currentPeelingEstimates(i)%hmmEstimate = hmmEstimate(:,:,i)
                enddo
            endif
            firstRun = .false.

            !Handle the Multilocus Peeler
            cycleIndex = 1
            converged = .false.
            do while(cycleIndex < nCycles .and. .not. converged)
                ! Forward Pass
                do i = 2, nSnps
                    call runIndex(pedigree%getAllGenotypesatPosition(i), i, currentPeelingEstimates, 1)
                    if(mod(i, 100) .eq. 0) print *, "Round ", roundNumber,  ", cycle ", cycleIndex, ", Forward ", i
                enddo
                ! Backward Pass
                do i = nSnps-1, 1, -1
                    call runIndex(pedigree%getAllGenotypesatPosition(i), i, currentPeelingEstimates, 2)
                    if(mod(i, 100) .eq. 0) print *, "Round ", roundNumber,  ", cycle ", cycleIndex, ", Backward ", i
                enddo
                ! Join Pass
                do i = nSnps, 1, -1
                    call runIndex(pedigree%getAllGenotypesatPosition(i), i, currentPeelingEstimates, 3, .true.)
                    if(mod(i, 100) .eq. 0) print *, "Round ", roundNumber,  ", cycle ", cycleIndex, ", Join ", i
                enddo
                if(cycleIndex > 1) converged = checkConvergence(currentPeelingEstimates)
                cycleIndex = cycleIndex + 1
                call updateAllRecombinationRates(currentPeelingEstimates)
            enddo

            if (present(writeOutputs)) then
                call writeOutputsToFile(roundNumber, currentPeelingEstimates)
            endif
        enddo

    
    end subroutine
    
    function checkConvergence(currentPeelingEstimates) result(res)
        type(peelingEstimates), dimension(:), pointer, intent(in) :: currentPeelingEstimates
        real(kind=real64), dimension(nSnps) :: markerError
        logical :: res
        integer :: i
        do i = 1, nSnps
            markerError(i) = currentPeelingEstimates(i)%estimatedError
        enddo
        print *, "MaxMarkerError", maxval(markerError) 
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
        real(kind=real64), dimension(:,:), allocatable :: oldPosteriorSire, oldPosteriorDame

        integer, dimension(:), allocatable :: ref, alt
        integer, dimension(:), allocatable :: tmpFamilyList
        ! real(kind=real64), dimension(4,4) :: segregationTransmissionMatrix
        real(kind=real64) :: error
        real(kind=real64) :: p, q, pf
        integer :: i, j, fam, father, mate
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
        error = markerEstimates%genotypingErrorRate
        genotypesToHaplotypes(:,0) = [ (1-error*2/3), error/6, error/6, error/3 ]
        genotypesToHaplotypes(:,1) = [ error/3, (1-error*2/3)/2, (1-error*2/3)/2, error/3 ]
        genotypesToHaplotypes(:,2) = [ error/3, error/6, error/6, (1-error*2/3) ]
        genotypesToHaplotypes(:,9) = [ .25, .25, .25, .25 ]

        if(.not. inputParams%isSequence) then
            markerEstimates%penetrance = log(genotypesToHaplotypes(:, genotypes))
        else
            error = markerEstimates%genotypingErrorRate
            ref = sequenceData(indexNumber, 1, :)
            alt = sequenceData(indexNumber, 2, :)

            p = log(1-error)
            q = log(error)
            pf = log(.5)

            markerEstimates%penetrance(1, :) = p*ref + q*alt
            markerEstimates%penetrance(2, :) = pf*ref + pf*alt - log(2D0)
            markerEstimates%penetrance(3, :) = pf*ref + pf*alt - log(2D0)
            markerEstimates%penetrance(4, :) = q*ref + p*alt

        endif
        if(markerEstimates%postHMM .and. allocated(markerEstimates%hmmEstimate)) then
            allocate(tmpVect(4, nPseudoFounders))
            call gemm(genotypesToHaplotypes(:,0:2), markerEstimates%hmmEstimate, tmpVect)
            markerEstimates%penetrance(:, pseudoFounders) = log(tmpVect)
        endif

        !Create Anterior
        p = log(markerEstimates%maf)
        q = log(1-markerEstimates%maf)
        anterior(1,:) = 2*q
        anterior(2,:) = p+q
        anterior(3,:) = p+q
        anterior(4,:) = 2*p

        if(markerEstimates%postHMM) then
            anterior(:, pseudoFounders) = 0
        endif
        markerEstimates%anterior = anterior

        !Create posterior
        !Not using pointers... maybe should?
        markerEstimates%sirePosteriorMate = markerEstimates%sirePosteriorMateAll(:,:,runType)
        markerEstimates%damePosteriorMate = markerEstimates%damePosteriorMateAll(:,:,runType)
        markerEstimates%posterior = markerEstimates%posteriorAll(:,:,runType)
        posterior => markerEstimates%posterior
        oldPosteriorSire = markerEstimates%sirePosteriorMate
        oldPosteriorDame = markerEstimates%damePosteriorMate

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
        
        do i = 1, nGenerations
            tmpFamilyList = familiesInGeneration(i)%array
            ! print *, i, "size", size(tmpFamilyList)
           !$omp parallel do &
           !$omp private(j, fam)  
            do j= 1,size(tmpFamilyList)      
                fam = tmpFamilyList(j)        
                call peelDown(markerEstimates, fam)
            enddo
           !$omp end parallel do
        enddo 
        
        do i = nGenerations, 1, -1
            tmpFamilyList = familiesInGeneration(i)%array
            
            !$omp parallel do &     
            !$omp private(j, fam)   
            do j= 1, size(tmpFamilyList)   
                fam = tmpFamilyList(j)        
                call updateSegregation(markerEstimates, fam)
                call peelUp(markerEstimates, fam)
            enddo
            !$omp end parallel do

            do j = 1, size(tmpFamilyList)
                fam = tmpFamilyList(j)
                father = listOfParents(1, fam)
                mate = listOfParents(2, fam)
                
                posterior(:,father) = posterior(:,father) - oldPosteriorSire(:,fam) + markerEstimates%sirePosteriorMate(:, fam)
                posterior(:,father) = posterior(:,father) - maxval(posterior(:,father))
                
                posterior(:,mate) = posterior(:,mate) - oldPosteriorDame(:,fam) + markerEstimates%damePosteriorMate(:, fam)
                posterior(:,mate) = posterior(:,mate) - maxval(posterior(:,mate))
            enddo

        enddo 

        !Now do post segregation updating.
        !Isn't this all currentSegregationEstimate -- yes, but it's saving it to different slots.
        if(runType == 1) markerEstimates%transmitForward = markerEstimates%currentSegregationEstimate * markerEstimates%pointSegregation
        if(runType == 2) markerEstimates%transmitBackward = markerEstimates%currentSegregationEstimate * markerEstimates%pointSegregation
        if(runType == 3) markerEstimates%fullSegregation = markerEstimates%currentSegregationEstimate * markerEstimates%pointSegregation

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
            ! call updateRecombinationRate(markerEstimates, currentPeelingEstimates, indexNumber)
            ! if(markerEstimates%recombinationRate > .01) print *, markerEstimates%recombinationRate
            if(.not. inputParams%isSequence) call updateGenotypeErrorRates(genotypes, haplotypeEstimates, markerEstimates)
            if(inputParams%isSequence) call updateSequenceErrorRates(ref, alt, haplotypeEstimates, markerEstimates)
            call updateMafEstimates(genotypeEstimates, markerEstimates)
        endif
        
        markerEstimates%sirePosteriorMateAll(:,:,runType) = markerEstimates%sirePosteriorMate 
        markerEstimates%damePosteriorMateAll(:,:,runType) = markerEstimates%damePosteriorMate 
        markerEstimates%posteriorAll(:,:,runType) = markerEstimates%posterior

        call markerEstimates%deallocateMarkerVariables
    end subroutine runIndex

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
        pfather = anterior(:, father)+penetrance(:,father)+posterior(:,father)-sirePosteriorMate(:,fam) !lhtp
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
            if(markerEstimates%postHMM) then
                if(isPseudoFounder(child)) then
                    anterior(:, child) = 0    
                endif 
            endif
        enddo
        deallocate(offspring, childEstimate)

    end subroutine


!---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief      Performs a series of peeling up operations.
    !
    !> @details     Peels up the pedigree updating posterior probabilities. This is done based on mating pairs in inverse order. 
    ! This is *only* done for mating pairs. If an individual is not in a mating pair then they have no children and so have a posterior estimate of 1.
    !
    !> @author     Andrew Whalen, awhalen@roslin.ed.ac.uk
    !
    !> @date       December 7, 2016
    !
    ! PARAMETERS:
    !> @param[out] genotypes integer array, read genotype file.
    !---------------------------------------------------------------------------    
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
        real(kind=real64), dimension(:,:,:), pointer, contiguous :: childTrace
        integer, dimension(:), allocatable :: offspring
        integer fam, j

        anterior => markerEstimates%anterior
        penetrance => markerEstimates%penetrance
        posterior => markerEstimates%posterior
        sirePosteriorMate => markerEstimates%sirePosteriorMate
        damePosteriorMate => markerEstimates%damePosteriorMate


        father = listOfParents(1, fam)
        mate = listOfParents(2, fam)
        pfather = anterior(:, father)+penetrance(:,father)+posterior(:,father)-sirePosteriorMate(:,fam) !lhtp
        pmate = anterior(:, mate)+penetrance(:,mate)+posterior(:,mate)-damePosteriorMate(:,fam)
        
        nChildren = offspringList(fam)%length
        allocate(offspring(nChildren))
        offspring = offspringList(fam)%convertToArrayIDs()

        pjoint = 0
        !MP
        do j=1, nChildren
            !May be able to pre-multiply this stuff together.
            child = offspring(j)
            childTrace => buildTraceTensor(markerEstimates%currentSegregationEstimate(:,child)*markerEstimates%pointSegregation(:,child))
            tmp = childTraceMultiply(penetrance(:, child)+posterior(:, child), childTrace)
            ! tmp = childTraceMultiplyMKL(penetrance(:, child) + posterior(:, child), childTrace)
        
            pjoint = pjoint + tmp
           
            deallocate(childTrace)
        enddo

        tempEstimate = jointByMateToFather(pmate, pjoint)
        ! posterior(:,father) = posterior(:,father) - sirePosteriorMate(:,fam)  + tempEstimate
        ! posterior(:,father) = posterior(:,father) - maxval(posterior(:,father))
        sirePosteriorMate(:,fam) = tempEstimate - maxval(tempEstimate)

        tempEstimate = jointByFatherToMate(pfather, pjoint)
        ! posterior(:,mate) = posterior(:,mate) - damePosteriorMate(:,fam) + tempEstimate
        ! posterior(:,mate) = posterior(:,mate) - maxval(posterior(:,mate))
        damePosteriorMate(:,fam) = tempEstimate - maxval(tempEstimate)
    
        deallocate(offspring)

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
        pfather = anterior(:, father)+penetrance(:,father)+posterior(:,father)-sirePosteriorMate(:,fam) !lhtp
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
                ! newPointEstimate(:,child) = reduceSegregationTensor(anterior(:,child)+posterior(:,child)+penetrance(:,child), pjoint-childEstimate(:,:,j))
                newPointEstimate(:,child) = reduceSegregationTensor(posterior(:,child)+penetrance(:,child), pjoint-childEstimate(:,:,j))
                ! newPointEstimate(:,child) = reduceSegregationTensorMKL(posterior(:,child)+penetrance(:,child), pjoint-childEstimate(:,:,j))
            endif
        enddo
        deallocate(offspring, childEstimate)

    end subroutine

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
        do seg=1, 4
            !Do for each child allele, aa, aA, Aa, AA
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
            ! print *, i
           do j = 1, nHaplotypes
                ! print *, i, " ", j
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
        !MP
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
        real(kind=real64), dimension(4,4) :: transmissionMatrix
        real(kind=real64) :: nChanges, nObservations, observedChangeRate
        ! type(fixedPointEstimator), pointer :: currentRecombinationEstimator
        integer i, j, indexNumber

        !We update n to n+1 here.
        transmissionMatrix = calculateSegregationTransmissionMatrix(markerEstimates%recombinationRate)


        nChanges = .01
        nObservations = 1.0

        if(indexNumber < nSnps) then 
            nObservations = nObservations + nAnimals
            plusZeroSegregation = markerEstimates%fullSegregation
            plusOneSegregation = currentPeelingEstimates(indexNumber+1)%fullSegregation
            
            zi = 0
            do i = 1, nAnimals
                do j= 1, 4
                    normalizedEstimate = plusZeroSegregation(:,i) * transmissionMatrix(:,j)
                    normalizedEstimate = normalizedEstimate/sum(normalizedEstimate)
                    zi(j, i) = zi(j, i) + plusOneSegregation(j,i) * (1-normalizedEstimate(j))
                enddo
            enddo
            nChanges = nChanges + sum(zi)
        endif

        observedChangeRate = nChanges/nObservations
        print *, "change rate", observedChangeRate
        ! print *, indexNumber, " ", nChanges, " ", nObservations
        ! print *, indexNumber, " ", observedChangeRate

        ! currentRecombinationEstimator => markerEstimates%recombinationEstimator
        ! call currentRecombinationEstimator%addObservation(markerEstimates%recombinationRate, observedChangeRate)
        ! markerEstimates%recombinationRate = currentRecombinationEstimator%secantEstimate(isLogitIn=.true.)
        markerEstimates%recombinationRate = observedChangeRate
    

    end subroutine



    subroutine updateAllRecombinationRates(currentPeelingEstimates)
        use globalGP, only:  nAnimals, nSnps
        use fixedPointModule
        implicit none
        type(peelingEstimates), pointer :: markerEstimates
        type(PeelingEstimates), dimension(:), pointer :: currentPeelingEstimates

        real(kind=real64), dimension(4,nAnimals) :: zi, plusZeroSegregation, plusOneSegregation
        real(kind=real64), dimension(4) :: normalizedEstimate
        real(kind=real64), dimension(4,4) :: transmissionMatrix
        real(kind=real64) :: nChanges, nObservations, observedChangeRate
        ! type(fixedPointEstimator), pointer :: currentRecombinationEstimator
        integer i, j, indexNumber

        !We update n to n+1 here.
        markerEstimates => currentPeelingEstimates(1)
        print *, markerEstimates%recombinationRate
        transmissionMatrix = calculateSegregationTransmissionMatrix(markerEstimates%recombinationRate)


        nChanges = .01
        nObservations = 1.0
        do indexNumber = 1, nSnps-1
            if(indexNumber < nSnps) then 
                nObservations = nObservations + nAnimals
                plusZeroSegregation = currentPeelingEstimates(indexNumber)%fullSegregation
                plusOneSegregation = currentPeelingEstimates(indexNumber+1)%fullSegregation

                zi = 0
                do i = 1, nAnimals
                    do j= 1, 4
                        normalizedEstimate = plusZeroSegregation(:,i) * transmissionMatrix(:,j)
                        normalizedEstimate = normalizedEstimate/sum(normalizedEstimate)
                        zi(j, i) = zi(j, i) + plusOneSegregation(j,i) * (1-normalizedEstimate(j))
                    enddo
                enddo
                nChanges = nChanges + sum(zi)
            endif
        enddo
        ! print *, "obs", nChanges, nObservations
        observedChangeRate = nChanges/nObservations
        print *, indexNumber, " ", nChanges, " ", nObservations
        print *, indexNumber, " ", observedChangeRate
        
        ! currentRecombinationEstimator => markerEstimates%recombinationEstimator
        ! call currentRecombinationEstimator%addObservation(markerEstimates%recombinationRate, observedChangeRate)
        ! observedChangeRate = currentRecombinationEstimator%secantEstimate(isLogitIn=.true.)

        print *, "change rate", observedChangeRate

        do indexNumber = 1, nSnps
            currentPeelingEstimates(indexNumber)%recombinationRate = observedChangeRate
        enddo

    end subroutine


    subroutine updateGenotypeErrorRates(genotypes, haplotypes, markerEstimates)
        use globalGP, only: nAnimals
        use fixedPointModule
        implicit none
        integer(kind=1), dimension(:), intent(in) :: genotypes
        real(kind=real64), dimension(:,:), intent(in) :: haplotypes        
        type(peelingEstimates), intent(inout) :: markerEstimates

        real(kind=real64), dimension(3,nAnimals) :: zi, reducedHaplotypes, recodedGenotypes
        real(kind=real64), dimension(3, 0:9) :: genotypesToHaplotypes
        real(kind=real64) :: nChanges, nObservations, observedChangeRate
        type(fixedPointEstimator), pointer :: currentErrorEstimator


        genotypesToHaplotypes(:,0) = [1, 0, 0]
        genotypesToHaplotypes(:,1) = [0, 1, 0]
        genotypesToHaplotypes(:,2) = [0, 0, 1]
        genotypesToHaplotypes(:,9) = [0, 0, 0]

        recodedGenotypes = genotypesToHaplotypes(:, genotypes)

        reducedHaplotypes(1,:) = haplotypes(1,:) 
        reducedHaplotypes(2,:) = haplotypes(2,:) + haplotypes(3,:)
        reducedHaplotypes(3,:) = haplotypes(4,:) 


        zi = recodedGenotypes * (1-reducedHaplotypes)

        nChanges = 0.05*2 + sum(zi)
        nObservations = 1.0*2 + sum(recodedGenotypes)

        observedChangeRate = nChanges/nObservations
        currentErrorEstimator => markerEstimates%genotypingErrorEstimator
        call currentErrorEstimator%addObservation(markerEstimates%genotypingErrorRate, observedChangeRate)
        markerEstimates%genotypingErrorRate = currentErrorEstimator%secantEstimate(isLogitIn=.true.)
    end subroutine


    subroutine updateSequenceErrorRates(ref, alt, haplotypes, markerEstimates)
        use globalGP, only: nAnimals
        use fixedPointModule
        implicit none
        integer, dimension(:), intent(in) :: ref, alt
        integer, dimension(:), allocatable :: totReads
        real(kind=real64), dimension(:,:), intent(in) :: haplotypes        
        type(peelingEstimates), intent(inout) :: markerEstimates

        real(kind=real64), dimension(3,nAnimals) :: reducedHaplotypes
        real(kind=real64) :: nChanges, nObservations, observedChangeRate
        type(fixedPointEstimator), pointer :: currentErrorEstimator


        totReads = ref + alt

        reducedHaplotypes(1,:) = haplotypes(1,:) 
        reducedHaplotypes(2,:) = haplotypes(2,:) + haplotypes(3,:)
        reducedHaplotypes(3,:) = haplotypes(4,:) 

        nChanges = 2 + sum(alt*reducedHaplotypes(1, :) +  ref*reducedHaplotypes(3, :) )
        nObservations = 4 + sum(totReads*reducedHaplotypes(1,:) + totReads*reducedHaplotypes(3,:))

        observedChangeRate = nChanges/nObservations
        currentErrorEstimator => markerEstimates%genotypingErrorEstimator
        call currentErrorEstimator%addObservation(markerEstimates%genotypingErrorRate, observedChangeRate)
        markerEstimates%genotypingErrorRate = currentErrorEstimator%secantEstimate(isLogitIn=.true.)
    end subroutine


    subroutine writeOutputsToFile(index, currentPeelingEstimates)
        use globalGP
        type(peelingEstimates), pointer :: markerEstimates
        type(peelingEstimates), dimension(:), pointer :: currentPeelingEstimates
        real(kind=real64), dimension(:,:,:), allocatable :: combinedHaplotypes
        integer, dimension(nSnps) :: individualGenotype
        real(kind=real64) :: threshold

        integer :: index, i, j


        print *, "Writting outputs"
        write(paramaterFile(index), '(a, a, a)') "maf ", "gError ", "tError"
        do i = 1, nSnps
            markerEstimates => currentPeelingEstimates(i)
            write(auxFile(index), '(60000f10.4)') markerEstimates%genotypeEstimates(:)
            write(paramaterFile(index), '(3f12.7)') markerEstimates%maf, markerEstimates%genotypingErrorRate, markerEstimates%recombinationRate
        enddo

        allocate(combinedHaplotypes(4, nSnps, nAnimals))


        do i = 1, nSnps
            markerEstimates => currentPeelingEstimates(i)
            combinedHaplotypes(:, i, :) = markerEstimates%haplotypeEstimates(:,:)
        enddo

        do i = 1, nAnimals
            write(outputFile(index), '(60000f10.4)') combinedHaplotypes(1,:, i)
            write(outputFile(index), '(60000f10.4)') combinedHaplotypes(2,:, i)
            write(outputFile(index), '(60000f10.4)') combinedHaplotypes(3,:, i)
            write(outputFile(index), '(60000f10.4)') combinedHaplotypes(4,:, i)
            write(outputFile(index), *) " " 
        enddo
        threshold = .9
        do i = 1, nAnimals
            individualGenotype = 9
            do j = 1, nSnps
                if(combinedHaplotypes(1, j, i) > threshold) individualGenotype(j) = 0
                if(combinedHaplotypes(2, j, i) > threshold) individualGenotype(j) = 1
                if(combinedHaplotypes(3, j, i) > threshold) individualGenotype(j) = 1
                if(combinedHaplotypes(4, j, i) > threshold) individualGenotype(j) = 2
            enddo

            write(consensusFile(index), '(60000i2)') individualGenotype
        enddo


      
    end subroutine

    function lhtp(logVector) result(prob)
        use globalGP, only: nHaplotypes
        implicit none
        real(kind=real64), dimension(nHaplotypes), intent(in) :: logVector
        real(kind=real64), dimension(nHaplotypes) :: prob
        prob = exp(logVector - maxval(logVector))
        prob = prob/sum(prob)
    end function

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
        print *, " "
        do i = 1, size(a, 2)
            write (*,"(200f35.30)",advance="no") a(:, i)
            print *, ""
        enddo
    end subroutine

    subroutine ppa3(a)
        implicit none
        real(kind=real64), dimension(:, :,:), intent(in) :: a
        integer i,j
        do j = 1, size(a, 3)
            do i = 1, size(a, 2)
                print *, a(:, i,j)
            enddo
            print *, ""
        enddo
    end subroutine

    subroutine ppa4(a)
        implicit none
        real(kind=real64), dimension(:, :,:,:), allocatable, intent(in) :: a
        integer i,j,k
        do k = 1, size(a, 4)
            do j = 1, size(a, 3)
                do i = 1, size(a, 2)
                    print *, nint(a(:, i,j,k))
                enddo
                print *, ""
            enddo
            print *, "----"
        enddo
    end subroutine

    subroutine printTrace(a)
        implicit none
        real(kind=real64), dimension(:, :,:), intent(in) :: a
        integer i,j
        do i = 1, size(a, 1)
            do j = 1, size(a, 2)
                write (*,"(200f12.7)",advance="no") a(i, j, :)
                print *, ""
            enddo
            print *, ""
        enddo
    end subroutine

    subroutine printSegregation(a)
        implicit none
        real(kind=real64), dimension(:, :,:,:), allocatable, intent(in) :: a
        integer i,j,k
        do i = 1, size(a, 1)
            do j = 1, size(a, 2)
                do k = 1, size(a, 3)
                    print *, nint(a(i, j, k, :))
                    ! print *, (a(j, :,i,k))
                enddo
                print *, ""
            enddo
            print *, "----"
        enddo
    end subroutine


end module

