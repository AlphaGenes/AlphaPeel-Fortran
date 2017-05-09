
module fixedPointModule

    use iso_fortran_env    
    implicit none


    type fixedPointEstimator
        real(kind=real64), dimension(:), allocatable :: values, responses
        integer :: maxIndex, currentIndex
        contains
            procedure :: secantEstimate
            procedure :: addObservation
            procedure :: initialize
            procedure :: grow
    end type fixedPointEstimator

    contains
        subroutine initialize(this)
            implicit none
            class(fixedPointEstimator) :: this

            this%maxIndex = 20
            this%currentIndex = 0

            allocate(this%values(this%maxIndex))
            allocate(this%responses(this%maxIndex))

        end subroutine
        
        subroutine addObservation(this, value, response)
            implicit none
            class(fixedPointEstimator) :: this
            real(kind=real64), intent(in) :: value, response
            integer :: index
            if(this%currentIndex == this%maxIndex) then
                call this%grow(this%maxIndex * 2)
            endif

            index = this%currentIndex + 1
            this%values(index) = value
            this%responses(index) = response
            this%currentIndex = index

        end subroutine

        subroutine grow(this, newLength)
            !This is probably a really bad way to do this.
            implicit none
            class(fixedPointEstimator) :: this
            integer :: newLength, oldLength
            real(kind=real64), dimension(newLength) :: values, responses
            values = 0
            responses = 0

            oldLength = this%maxIndex
            this%maxIndex = newLength
            values(1:oldLength) = this%values
            responses(1:oldLength) = this%responses
            this%values = values
            this%responses = responses
        end subroutine

        function secantEstimate(this, indexIn, isLogitIn) result(prediction)
            implicit none
            class(fixedPointEstimator) :: this
            real(kind=real64) :: prediction
            real(kind=real64) :: x1, x2, y1, y2, a, b
            real(kind=real64) :: maxJump
            logical, optional, intent(in) :: isLogitIn
            integer, optional, intent(in) :: indexIn
            logical :: isLogit
            integer :: index
            !Set optional values
            if(.not. present(isLogitIn)) then
                isLogit = .false.
            else
                isLogit = isLogitIn
            endif

            if(.not. present(indexIn)) then
                index = this%currentIndex
            else
                index = indexIn
            endif
            if(index == 0) then
                print *, "no data to estimate off of, returning 0/.5"
                prediction = 0
                if(isLogit) prediction = .5
                RETURN
            endif

            !Calculate estimate
            if(index == 1) then
                prediction = this%responses(index)
                RETURN
            endif
            !index > 1 so index-1 and index are both defined.
            x1 = this%values(index)
            x2 = this%values(index-1)
            y1 = this%responses(index)
            y2 = this%responses(index-1)
            if(abs(x1 - x2) < .0000000001) then
                prediction = x2
                RETURN
            endif
            if(isLogit) then 
                x1 = logit(x1)
                x2 = logit(x2)
                y1 = logit(y1)
                y2 = logit(y2)
            endif

            y1 = y1-x1
            y2 = y2-x2

            a = (y1-y2)/(x1-x2)
            b = y2 - a*x2
            prediction = -b/a

            maxJump = .5
            if(abs(prediction - x1) > maxJump) then
                if(prediction < x1) prediction = x1- maxJump
                if(prediction > x1) prediction = x1 + maxJump
            endif
            if(isLogit) then
                prediction = invlogit(prediction)
            endif
        end function


        function logit(a) result(res)
            implicit none
            real(kind=real64), intent(in) :: a
            real(kind=real64) res

            res = log(a/(1-a))
        end function

        function invlogit(a) result(res)
            implicit none
            real(kind=real64), intent(in) :: a
            real(kind=real64) res

            res = exp(a)/(1+exp(a))
        end function



end module

