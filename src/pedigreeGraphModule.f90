module pedigreeGraphModule

use ConstantModule
use PedigreeModule
    implicit none

contains
subroutine calculateCommonAndDistance(commonAnDistance, ped)

    integer, dimension(:,:), allocatable, intent(out) :: commonAnDistance
    type(IndividualLinkedListNode) ,pointer ::indNode,indNode2
    type(Individual), pointer :: tmpSire, tmpDam
    integer :: sireVal, damVal
    integer :: i,j,k,h
    type(pedigreeHolder), intent(in) :: ped

    allocate(commonAnDistance(ped%pedigreeSize,ped%pedigreeSize))
    commonAnDistance = -9999
        indNode => ped%generations(0)%first
        do i=1,ped%generations(0)%length
            indNode2 => ped%generations(0)%first
            do h=1,ped%generations(0)%length
                if (h == i) then
                    commonAnDistance(indNode%item%id, indNode2%item%id) = 0
                    commonAnDistance(indNode2%item%id, indNode%item%id) = 0
                else 
                    commonAnDistance(indNode%item%id, indNode2%item%id) = -1
                    commonAnDistance(indNode2%item%id, indNode%item%id) = -1
                endif
                indNode2 => indNode2%next
            enddo
            indNode => indNode%next

        end do    

        ! THIS has to be performed in order
        do i =1, ped%maxGeneration
            indNode => ped%generations(i)%first
            do h =1, ped%generations(i)%length

                do j=1, i-1
                    indNode2 => ped%generations(j)%first
                    do k=1, ped%generations(j)%length
                        tmpSire =>indNode%item%sirePointer
                        tmpDam => indNode%item%dampointer

                        sireVal =  commonAnDistance(tmpSire%id,indNode2%item%id)
                        damVal = commonAnDistance(tmpDam%id,indNode2%item%id)

                        if (sireVal == -1) then
                            commonAnDistance(indNode%item%id, indNode2%item%id) = damVal+1
                            commonAnDistance(indNode2%item%id, indNode%item%id) = damVal+1
                        else if (damVal == -1) then
                            commonAnDistance(indNode%item%id, indNode2%item%id) = sireVal+1
                            commonAnDistance(indNode2%item%id, indNode%item%id) = sireVal+1
                        else
                            commonAnDistance(indNode%item%id, indNode2%item%id) =min(sireVal,damVal)+1
                            commonAnDistance(indNode2%item%id, indNode%item%id) = min(sireVal,damVal)+1
                        endif
                        indNode2 => indNode2%next
                    enddo

                    
                enddo
                indNode => indNode%next
            enddo
            indNode => ped%generations(i)%first
            do h=1,ped%generations(i)%length
                indNode2 = ped%generations(i)%first
                do k = 1, ped%generations(i)%length
                    if(h == k) then
                        commonAnDistance(indNode%item%id, indNode2%item%id) = 0
                    else
                        tmpSire =>indNode%item%sirePointer
                        tmpDam => indNode%item%dampointer

                        sireVal =  commonAnDistance(tmpSire%id,indNode2%item%id)
                        damVal = commonAnDistance(tmpDam%id,indNode2%item%id)

                        if (sireVal == -1) then
                            commonAnDistance(indNode%item%id, indNode2%item%id) = damVal+1
                            commonAnDistance(indNode2%item%id, indNode%item%id) = damVal+1
                        else if (damVal == -1) then
                            commonAnDistance(indNode%item%id, indNode2%item%id) = sireVal+1
                            commonAnDistance(indNode2%item%id, indNode%item%id) = sireVal+1
                        else
                            commonAnDistance(indNode%item%id, indNode2%item%id) =min(sireVal,damVal)+1
                            commonAnDistance(indNode2%item%id, indNode%item%id) = min(sireVal,damVal)+1
                        endif
                    endif
                    indNode2 => indNode2%next
                enddo
                indNode=> indNode%next
            enddo 
        enddo
    end subroutine calculateCommonAndDistance

    
    

        


end module pedigreeGraphModule
