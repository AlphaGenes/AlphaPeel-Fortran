!-----------------------------------------------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-----------------------------------------------------------------------------------------------------------------------
!
! MODULE: AlphaMLPInputModule
!
!> @file        AlphaMLPInputModule.f90
!
! DESCRIPTION:
!> @brief       Module holding input parameters
!>
!> @details     This MODULE contains a class which contains all input parameters read in from a spec file.
!> It also contains the default container object for the spec file, defaultInput.
!
!> @author      David Wilson, david.wilson@roslin.ed.ac.uk
!
!> @date        Feb 07, 2016
!
!> @version     0.0.1 (alpha)
!
! REVISION HISTORY:
! 2016.11.07  DWilson - Initial Version
!
!-----------------------------------------------------------------------------------------------------------------------
 

 Module AlphaMLPInputModule


    use iso_fortran_env
    use ConstantModule, only : FILELENGTH
    

    type AlphaMLPInput
        integer :: nGenotypedAnimals
        integer :: nSnp
        character(len=FILELENGTH) :: inputFile
        character(len=FILELENGTH) :: outputFile
        integer :: startSnp
        integer :: endSnp
        integer :: nCycles
        logical :: isSequence, useHMM
        character(len=FILELENGTH) :: sequenceFile
        character(len=FILELENGTH) :: pedFile

    end type AlphaMLPInput

    integer(KIND=1), allocatable, dimension(:,:) :: defaultInputGenotypes



    interface AlphaMLPInput

        module procedure initFromFile
        module procedure initFromParams
        module procedure initFromSnps
    end interface AlphaMLPInput
    contains

        !---------------------------------------------------------------------------
        !> @brief Constructor for AlphaMLP input based on passed in parameters
        !> @author  David Wilson david.wilson@roslin.ed.ac.uk
        !> @date    Febuary 08, 2016
        !> @return AlphaMLPInput of info for spec file
        !---------------------------------------------------------------------------
        function initFromParams(nGenotypedAnimals,nSnp,inputFile,outputFile,startSnp,endSnp) result(res)
            integer,intent(in) :: nGenotypedAnimals
            integer,intent(in) :: nSnp
            character(len=*),intent(in) :: inputFile
            character(len=*),intent(in) :: outputFile
            integer,intent(in) :: startSnp
            integer,intent(in) :: endSnp

            type(AlphaMLPInput) :: res


            res%nGenotypedAnimals = nGenotypedAnimals
            res%nsnp = nsnp
            res%inputFile = inputFile
            res%outputFile = outputFile
            res%startSnp = startSnp
            res%endSnp = endSnp



        end function initFromParams

                !---------------------------------------------------------------------------
        !> @brief Constructor for AlphaMLP input based on passed in parameters
        !> @author  David Wilson david.wilson@roslin.ed.ac.uk
        !> @date    Febuary 08, 2016
        !> @return AlphaMLPInput of info for spec file
        !---------------------------------------------------------------------------
        function initFromSnps(startSnp,endSnp) result(res)

            integer,intent(in) :: startSnp
            integer,intent(in) :: endSnp

            type(AlphaMLPInput) :: res


            res%startSnp = startSnp
            res%endSnp = endSnp


        end function initFromSnps

        
                !---------------------------------------------------------------------------
        !> @brief Constructor for AlphaMLP input based on a file that will be read in
        !> @author  David Wilson david.wilson@roslin.ed.ac.uk
        !> @date    Febuary 08, 2016
        !> @return AlphaMLPInput of info for spec file
        !---------------------------------------------------------------------------
        function initFromFile(SpecFileIn) result(res)
            use AlphaHouseMod, only: parseToFirstWhitespace,splitLineIntoTwoParts,toLower

            character(len=*),optional, intent(in) :: SpecFileIn !< Spec file input
            character(len=FILELENGTH) :: SpecFile
            type(AlphaMLPInput) :: res
            integer :: unit,IOStatus
            character(len=300) :: first, line, tmp
            character(len=:), allocatable::tag
            character(len=300),dimension(:), allocatable :: second

            if (present(SpecFileIn)) then
                specFile = SpecFileIn
            else
                specFile = "AlphaMLPSpec.txt"
            endif

            ! init everything
            res%nGenotypedAnimals = 0
            res%nSnp = 0
            res%startSnp = 0
            res%endSnp = 0
            res%inputFile = "AlphaMLPGenotypes.txt"
            res%outputFile = "AlphaMLPs.txt"
            res%pedFile = "No Pedigree"

            res%isSequence = .false.
            res%sequenceFile = ""

            res%nCycles = 10
            res%useHMM = .false.

            open(newunit=unit, file=SpecFile, action="read", status="old")
            IOStatus = 0
            READFILE: do while (IOStatus==0)
                read(unit,"(A)", IOStat=IOStatus)  line
                if (len_trim(line)==0) then
                    CYCLE
                end if

                call splitLineIntoTwoParts(trim(line), first, second)
                tag = parseToFirstWhitespace(first)
                if (first(1:1)=="=" .or. len(trim(line))==0) then
                    cycle
                else
                    select case(trim(tag))

                        case("nanis")
                            read(second(1),*) res%nGenotypedAnimals

                        case("nsnp")
                            read(second(1),*) res%nsnp

                        case("inputfilepath")
                            if (.not. allocated(second)) then
                                write(*, "(A,A)") "No input file specified. Using default filename: ", res%inputFile
                            else
                                write(res%inputFile, "(A)") trim(second(1))
                            end if

                        case("outputfilepath")
                            if (.not. allocated(second)) then
                                write(*, "(A,A)") "No output file specified. Using default filename: ", res%outputFile
                            else
                                write(res%outputFile, "(A)") trim(second(1))
                            end if
                        
                        case("pedigree")
                            write(res%pedFile, "(A)") trim(second(1))
                        
                        case("startsnp")
                            read(second(1),*) res%startsnp

                        case("endsnp")
                            read(second(1),*) res%endSnp

                        case("usesequence")
                            read(second(1),*) tmp
                            if(tmp == "yes") res%isSequence = .true.

                        case("sequencefile")
                            read(second(1),*) res%sequenceFile
                        case("usehmm")
                            read(second(1),*) tmp
                            if(tmp == "yes") res%useHMM = .true.
                        case("ncycles")
                            read(second(1),*) res%nCycles
                        


                        case default
                            write(*,"(A,A)") trim(tag), " is not valid for the AlphaMLP Spec File."
                            cycle
                    end select
                endif

                
            enddo READFILE
            
        end function initFromFile



        subroutine readGenotypes(input, pedigree)
            use PedigreeModule
            use ConstantModule, only : IDLENGTH,DICT_NULL

            type(AlphaMLPInput),intent(in) :: input
            type(PedigreeHolder) , intent(inout) :: pedigree
            ! type(Pedigreeholder), intent(inout) :: genotype
            integer(KIND=1), allocatable, dimension(:) :: tmp
            integer :: unit, tmpID,i
            character(len=IDLENGTH) :: seqid, seqsire, seqdam !placeholder variables
            open(newunit=unit,FILE=trim(input%inputFile),STATUS="old") !INPUT FILE
        
            ! allocate(res(input%nGenotypedAnimals,input%endSnp-input%startSnp+1))
            allocate(tmp(input%endSnp-input%startSnp+1))
            tmp = 9
            do i=1,input%nGenotypedAnimals
                print *, i
                read (unit,*) seqid, seqsire, seqdam, tmp(:)

                tmpID = pedigree%dictionary%getValue(seqid)

                if (tmpID /= DICT_NULL) then
                    call pedigree%pedigree(tmpID)%setGenotypeArray(tmp)
                endif
            end do

            close(unit)
        
        end subroutine readGenotypes

        subroutine readSequence(input, pedigree, sequenceData)
            use PedigreeModule
            use ConstantModule, only : IDLENGTH,DICT_NULL

            type(AlphaMLPInput),intent(in) :: input
            type(PedigreeHolder) , intent(inout) :: pedigree
            ! type(Pedigreeholder), intent(inout) :: genotype
            integer(KIND=1), allocatable, dimension(:) :: tmp, ref, alt
            integer(KIND=1), allocatable, dimension(:,:) :: genoEst
            integer, allocatable, dimension(:,:,:) :: sequenceData
            integer :: unit, tmpID,i, j
            character(len=IDLENGTH) :: seqid !placeholder variables
            real(kind=real64) :: err, p, q, pf

            allocate(sequenceData(input%endSnp-input%startSnp+1, 2, input%nGenotypedAnimals))

            open(newunit=unit,FILE=trim(input%sequenceFile),STATUS="old") !INPUT FILE
        
            ! allocate(res(input%nGenotypedAnimals,input%endSnp-input%startSnp+1))
            allocate(ref(input%endSnp-input%startSnp+1))
            allocate(alt(input%endSnp-input%startSnp+1))
            allocate(genoEst(input%endSnp-input%startSnp+1, 3))
            allocate(tmp(input%endSnp-input%startSnp+1))

            err = 0.01
            p = log(err)
            q = log(1-err)
            pf = log(.5)

            ! tmp = 9
            sequenceData = 0
            do i=1,input%nGenotypedAnimals
                read (unit,*) seqid, ref(:)
                read (unit,*) seqid, alt(:)

                tmpID = pedigree%dictionary%getValue(seqid)
                print *, tmpID
                sequenceData(:, 1, tmpId) = ref(:)
                sequenceData(:, 2, tmpId) = alt(:)

                genoEst(:, 1) = p*ref + q*alt
                genoEst(:, 2) = pf*ref + pf*alt
                genoEst(:, 3) = q*ref + p*alt

                tmp = maxloc(genoEst, dim=2) - 1
                where(ref+alt < 15) tmp = 9

                if (tmpID /= DICT_NULL) then
                    call pedigree%pedigree(tmpID)%setGenotypeArray(tmp)
                endif
            end do

            close(unit)
        end subroutine readSequence


end module AlphaMLPInputModule