
module globalGP
    use fixedPointModule
    use iso_fortran_env
    use PedigreeModule
    use IndividualLinkedListModule
    use AlphaMLPInputModule, only : AlphaMLPInput
    implicit none

    type(AlphaMLPInput) :: inputParams

    integer(kind=4) :: nHaplotypes = 4, nAnimals, nMatingPairs
    integer :: nSnps, nPseudoFounders, nSnpsAll
    integer, dimension(:), allocatable :: founders, phaseChildren, generations, pseudoFounders
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
    integer, dimension(2) :: outputFile
    integer, dimension(2) :: auxFile
    integer, dimension(2) :: segregationFile, consensusFile, pointEstimate
    integer, dimension(2) :: paramaterFile
    
    integer, dimension(:,:,:), allocatable :: sequenceData

    logical, dimension(:), allocatable :: isPseudoFounder, isPhasedChild


    type arrayContainer
        integer, dimension(:), allocatable :: array
    end type arrayContainer

    type(arrayContainer), dimension(:), allocatable :: familiesInGeneration 
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
        logical :: postHMM

        ! real(kind=real64), dimension(:,:), allocatable :: penetrance, anterior
        real(kind=real64) :: genotypingErrorRate, maf, recombinationRate
        type(fixedPointEstimator), pointer:: recombinationEstimator, genotypingErrorEstimator, mafEstimator
        contains
        procedure :: initializePeelingEstimates
        procedure :: allocateMarkerVariables
        procedure :: deallocateMarkerVariables
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

            this%postHMM = .false.

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

end module


