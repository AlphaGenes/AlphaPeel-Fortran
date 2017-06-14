
module globalGP
    use fixedPointModule
    use iso_fortran_env
    use blas95
    use f95_precision
    use IFPORT
    use PedigreeModule
    use IndividualLinkedListModule
    use AlphaMLPInputModule, only : AlphaMLPInput
    implicit none

    type(AlphaMLPInput) :: inputParams

    integer(kind=4) :: nHaplotypes = 4, nAnimals, nMatingPairs
    integer :: nSnps, nPseudoFounders, nSnpsAll
    integer, dimension(:), allocatable :: founders, phaseChildren, generations
    integer, dimension(:, :), allocatable :: listOfParents

    real(kind=real64), dimension(:,:,:,:), pointer, CONTIGUOUS :: segregationTensor, segregationTensorParentsFirst
    real(kind=real64), dimension(:,:,:), allocatable :: traceTensorExp
    real(kind=real64) :: recombinationRate = .01

    real(kind=real64), dimension(:,:), allocatable :: phaseChildrenOverride

    real(kind=real64), dimension(4), parameter :: haplotypesToGenotypes = (/0, 1, 1, 2/)

    type(PedigreeHolder) :: pedigree
    type(IndividualLinkedList), dimension(:), allocatable :: offspringList
    ! type(peelingEstimates), dimension(:), pointer :: currentPeelingEstimates

    ! units for file
    integer :: outputFile
    integer :: auxFile
    integer :: segregationFile, consensusFile, pointEstimate
    integer :: paramaterFile
    
    integer, dimension(:,:,:), allocatable :: sequenceData

    logical, dimension(:), allocatable :: isPhasedChild

    integer, dimension(:, :), allocatable :: mapIndexes
    real(kind=real64), dimension(:), allocatable :: mapDistance

    type arrayContainer
        integer, dimension(:), allocatable :: array
    end type arrayContainer

    type(arrayContainer), dimension(:), allocatable :: familiesInGeneration 
    type(arrayContainer), dimension(:), allocatable :: childrenAtGeneration 
    integer :: nGenerations

    type peelingEstimates
        real(kind=real64), dimension(:,:,:), allocatable :: posteriorAll, sirePosteriorMateAll, damePosteriorMateAll

        real(kind=real64), dimension(:,:), allocatable :: anterior, penetrance, posterior
        real(kind=real64), dimension(:,:), allocatable :: sirePosteriorMate, damePosteriorMate

        real(kind=real64), dimension(:,:), pointer :: currentSegregationEstimate
        real(kind=real64), dimension(:,:), allocatable :: pointSegregation, fullSegregation, transmitForward, transmitBackward

        real(kind=real64), dimension(:,:), allocatable :: hmmEstimate, haplotypeEstimates
        real(kind=real64), dimension(:), allocatable :: genotypeEstimates

        real(kind=real64), dimension(:,:,:,:), pointer, CONTIGUOUS :: currentSegregationTensors
        real(kind=real64) :: estimatedError

        ! real(kind=real64), dimension(:,:), allocatable :: penetrance, anterior
        real(kind=real64) :: genotypingErrorRate, maf, recombinationRate
        type(fixedPointEstimator), pointer:: recombinationEstimator, genotypingErrorEstimator, mafEstimator
        contains
        procedure :: initializePeelingEstimates
        procedure :: deallocatePeelingEstimates
        procedure :: allocateMarkerVariables
        procedure :: deallocateMarkerVariables
        procedure :: getHaplotypeEstimates
        procedure :: getGenotypeEstimates
        procedure :: setHaplotypeEstimates
        procedure :: setAnterior
        procedure :: setPenetrance
    end type peelingEstimates

    contains 
        subroutine initializePeelingEstimates(this, nHaplotypes, nAnimals, nMatingPairs, nSnpsAll)
            implicit none
            class(peelingEstimates) :: this
            integer :: nHaplotypes, nAnimals, nMatingPairs, nSnpsAll

            allocate(this%posteriorAll(nHaplotypes, nAnimals,3))
            allocate(this%sirePosteriorMateAll(nHaplotypes, nMatingPairs,3))
            allocate(this%damePosteriorMateAll(nHaplotypes, nMatingPairs,3))

            allocate(this%haplotypeEstimates(nHaplotypes, nAnimals))
            allocate(this%genotypeEstimates(nAnimals))

            allocate(this%fullSegregation(4, nAnimals))
            allocate(this%transmitForward(4, nAnimals))
            allocate(this%transmitBackward(4, nAnimals))



            this%transmitForward = 1
            this%transmitBackward = 1
            this%fullSegregation = 1.0/4.0

            this%genotypeEstimates = 0

            this%posteriorAll = 0
            this%sirePosteriorMateAll = 0
            this%damePosteriorMateAll = 0

            this%genotypingErrorRate = .01
            this%recombinationRate = 1d0/nSnpsAll
            this%maf = .5
            allocate(this%recombinationEstimator)
            allocate(this%genotypingErrorEstimator)
            allocate(this%mafEstimator)
            call this%recombinationEstimator%initialize
            call this%genotypingErrorEstimator%initialize
            call this%mafEstimator%initialize
            
        end subroutine

        subroutine deallocatePeelingEstimates(this)
            class(peelingEstimates) :: this

            deallocate(this%posteriorAll)
            deallocate(this%sirePosteriorMateAll)
            deallocate(this%damePosteriorMateAll)
            deallocate(this%haplotypeEstimates)
            deallocate(this%genotypeEstimates)
            deallocate(this%fullSegregation)
            deallocate(this%transmitForward)
            deallocate(this%transmitBackward)
            deallocate(this%recombinationEstimator)
            deallocate(this%genotypingErrorEstimator)
            deallocate(this%mafEstimator)
        end subroutine

        subroutine allocateMarkerVariables(this, nHaplotypes, nAnimals, nMatingPairs) 
            implicit none
            class(peelingEstimates) :: this
            integer :: nHaplotypes, nAnimals, nMatingPairs

            allocate(this%pointSegregation(4, nAnimals))
            allocate(this%currentSegregationEstimate(4, nAnimals))
            this%pointSegregation = 1
            this%currentSegregationEstimate = 1


            allocate(this%penetrance(nHaplotypes, nAnimals))
            allocate(this%anterior(nHaplotypes, nAnimals))
            ! print *, "Hello world!", nHaplotypes, nMatingPairs
            allocate(this%posterior(nHaplotypes, nAnimals))
            allocate(this%sirePosteriorMate(nHaplotypes, nMatingPairs))
            allocate(this%damePosteriorMate(nHaplotypes, nMatingPairs))

            this%penetrance = 1
            this%anterior = 1
            this%posterior = 1
            this%sirePosteriorMate = 1
            this%damePosteriorMate = 1


            allocate(this%currentSegregationTensors(nHaplotypes, nHaplotypes, nHaplotypes, nAnimals))

        end subroutine
        subroutine deallocateMarkerVariables(this)
            implicit none
            class(peelingEstimates) :: this

            deallocate(this%pointSegregation)
            deallocate(this%currentSegregationEstimate)
            deallocate(this%penetrance)
            deallocate(this%anterior)
            deallocate(this%posterior)
            deallocate(this%sirePosteriorMate)
            deallocate(this%damePosteriorMate)
            deallocate(this%currentSegregationTensors)
        end subroutine


        function getGenotypeEstimates(this) result(res)
            class(peelingEstimates) :: this
            real(kind=real64), dimension(nHaplotypes, nAnimals) :: haplotypeEstimates
            real(kind=real64), dimension(nAnimals) :: res, weights

            haplotypeEstimates = this%getHaplotypeEstimates()
            res=0
            call gemv(haplotypeEstimates, haplotypesToGenotypes, res, trans = "T")
            weights = sum(haplotypeEstimates, 1)
            res = res/weights


        end function

        subroutine setHaplotypeEstimates(this)
            class(peelingEstimates) :: this
            real(kind=real64), dimension(nHaplotypes, nAnimals) :: logHaplotypes, haplotypeEstimates
            real(kind=real64), dimension(nAnimals) :: weights
            integer :: i 

            logHaplotypes = this%anterior + this%posterior + this%penetrance
            do i = 1, nAnimals
                haplotypeEstimates(:,i) = lhtp(logHaplotypes(:,i))
            enddo
            this%haplotypeEstimates = haplotypeEstimates
        end subroutine

        function getHaplotypeEstimates(this) result(res)
            class(peelingEstimates) :: this
            real(kind=real64), dimension(nHaplotypes, nAnimals) :: res
 
            call this%setHaplotypeEstimates()
            res = this%haplotypeEstimates
        end function

        function lhtp(logVector) result(prob)
            implicit none
            real(kind=real64), dimension(nHaplotypes), intent(in) :: logVector
            real(kind=real64), dimension(nHaplotypes) :: prob
            prob = exp(logVector - maxval(logVector))
            prob = prob/sum(prob)
        end function

        subroutine setPenetrance(this, genotypes)
            implicit none
            class(peelingEstimates) :: this
            integer(kind=1), dimension(:), intent(in) :: genotypes
            real(kind=real64) :: error
            real(kind=real64), dimension(nHaplotypes,0:9) :: genotypesToHaplotypes
            integer, dimension(:), allocatable :: ref, alt
            real(kind=real64) :: p, q, pf !pf is log of .5
            integer :: indexNumber 
            error = this%genotypingErrorRate
            genotypesToHaplotypes(:,0) = [ (1-error*2/3), error/6, error/6, error/3 ]
            genotypesToHaplotypes(:,1) = [ error/3, (1-error*2/3)/2, (1-error*2/3)/2, error/3 ]
            genotypesToHaplotypes(:,2) = [ error/3, error/6, error/6, (1-error*2/3) ]
            genotypesToHaplotypes(:,9) = [ .25, .25, .25, .25 ]

            if(.not. inputParams%isSequence) then
                this%penetrance = log(genotypesToHaplotypes(:, genotypes))
            else
                error = this%genotypingErrorRate
                ref = sequenceData(indexNumber, 1, :)
                alt = sequenceData(indexNumber, 2, :)

                p = log(1-error)
                q = log(error)
                pf = log(.5)

                this%penetrance(1, :) = p*ref + q*alt
                this%penetrance(2, :) = pf*ref + pf*alt - log(2D0)
                this%penetrance(3, :) = pf*ref + pf*alt - log(2D0)
                this%penetrance(4, :) = q*ref + p*alt

            endif
            
            !Maybe use this as a second error measure for sequence.

            ! if(markerEstimates%postHMM .and. allocated(markerEstimates%hmmEstimate)) then
            !     allocate(tmpVect(4, nPseudoFounders))
            !     call gemm(genotypesToHaplotypes(:,0:2), markerEstimates%hmmEstimate, tmpVect)
            !     markerEstimates%penetrance(:, pseudoFounders) = log(tmpVect)
            ! endif


        end subroutine

        subroutine setAnterior(this, maf)
            implicit none
            class(peelingEstimates) :: this

            real(kind=real64) :: p, q, maf
            real(kind=real64), dimension(nHaplotypes, nAnimals) :: res
           
            p = log(maf)
            q = log(1-maf)
            res(1,:) = 2*q
            res(2,:) = p+q
            res(3,:) = p+q
            res(4,:) = 2*p
            this%anterior = res

        end subroutine

end module


