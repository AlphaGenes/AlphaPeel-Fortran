
module graphModule
    use iso_fortran_env
    implicit none

    type graph
        integer, dimension(:), allocatable :: id
        integer, dimension(:,:), allocatable :: network
        real(kind=real64), dimension(:,:), allocatable :: weightedNetwork
        integer, dimension(:), allocatable :: degree
        real(kind=real64)  :: score
        real(kind=real64), dimension(:), allocatable  :: individualScore
        integer :: size

        contains
        procedure :: removeIndex
        procedure :: removeIndexes
        procedure :: initializeGraph
        procedure :: cleanInitialize
        procedure :: deconstructGraph
        procedure :: brute_minimalVertexCover
        procedure :: removeEdge_minimalVertexCover
        procedure :: getSubGraphs
        procedure :: getMinimalCover
        procedure :: decomposeGraph

    end type graph

    contains 
        subroutine initializeGraph(this, id, network, individualScore)
            implicit none
            class(graph) :: this
            integer, dimension(:), intent(in) :: id
            integer, dimension(:,:), intent(in) :: network
            integer, dimension(:,:), allocatable :: localNetwork
            real(kind=real64), dimension(:,:), allocatable :: weightedNetwork
            real(kind=real64), dimension(:), intent(in) :: individualScore
            integer :: i,j

            if(size(network, 1) .ne. size(network, 2)) then
                print *, "Non square matrix as a network!"
                stop
            endif
            localNetwork = min(1, transpose(network) + network) !This ensures this is symetric with zero diagonal
            
            weightedNetwork = localNetwork
            do i = 1, size(localNetwork, 1)
                localNetwork(i,i) = 0
            enddo
            this%degree = sum(localNetwork, dim = 1)

            do i = 1, size(network, 1)
                do j = 1, size(network, 2)
                    weightedNetwork(i, j) = localNetwork(i, j) * 1./individualScore(i) * 1./individualScore(j)
                enddo
            enddo

            this%id = id
            this%network = localNetwork
            this%weightedNetwork = weightedNetwork
            this%score = 0
            ! this%degree = sum(network, dim = 1)
            this%size = size(id)
            this%individualScore = individualScore

        end subroutine

        subroutine cleanInitialize(this, id, network, weightedNetwork, individualScore, degree)
            implicit none
            class(graph) :: this
            integer, dimension(:), intent(in) :: id, degree
            integer, dimension(:,:), intent(in) :: network
            real(kind=real64), dimension(:), intent(in) :: individualScore
            real(kind=real64), dimension(:,:), intent(in) :: weightedNetwork

            this%degree = degree
            this%id = id
            this%network = network
            this%weightedNetwork = weightedNetwork
            this%size = size(id)
            this%individualScore = individualScore

            this%score = 0

        end subroutine

        subroutine deconstructGraph(this)
            implicit none
            class(graph) :: this
            deallocate(this%id)
            deallocate(this%network)
            deallocate(this%degree)

        end subroutine

        function decomposeGraph(this) result(res)
            implicit none
            class(graph) :: this
            type(graph) :: g
            integer, dimension(:), allocatable :: res, locA, indexes
            integer :: loc, currentIndex, count
            allocate(res(this%size))
            print *, this%size

            res = 0
            g = this
            currentIndex = 1
            count = 0
            do while(g%size > 0 .and. count < 100)
                count = count + 1
                locA = minloc(g%degree)
                loc = locA(1)
                !If 0 or 1, we can remove them all at the same time.
                if(g%degree(loc) .le. 1) then
                    indexes = selectIndexesBasedOnMask(g%degree == g%degree(loc))
                else
                    indexes = [loc]
                endif
                
                print *, "degree", g%degree(loc), "removing", indexes

                res(currentIndex:(currentIndex+size(indexes)-1)) = g%id(indexes)
                currentIndex = currentIndex + size(indexes)
                g = g%removeIndexes(indexes)
            enddo

        end function

        function removeIndex(this, toRemove) result(res)
            implicit none
            class(graph) :: this
            type(graph) :: res
            integer toRemove

            res = this%removeIndexes([toRemove])

        end function
        function removeIndexes(this, toRemove) result(res)
            implicit none
            class(graph) :: this
            type(graph) :: res
            integer, dimension(:), allocatable :: otherIndexes
            integer, dimension(:), intent(in) :: toRemove
            integer, dimension(:), allocatable :: newId, newDegree
            integer, dimension(:,:), allocatable :: newNetwork
            real(kind=real64), dimension(:), allocatable :: newScores
            real(kind=real64), dimension(:,:), allocatable :: newWeightedNetwork
            
            otherIndexes = selectOthers(toRemove, this%size)

            newId = this%id(otherIndexes)
            newScores = this%individualScore(otherIndexes)

            newNetwork = this%network(otherIndexes, otherIndexes)
            newWeightedNetwork = this%weightedNetwork(otherIndexes, otherIndexes)

            newDegree = this%degree - sum(this%network(:, toRemove), 2)
            newDegree = newDegree(otherIndexes)


            call res%cleanInitialize(newId, newNetwork, newWeightedNetwork, newScores, newDegree)

        end function

        recursive function brute_minimalVertexCover(this) result(res)
            implicit none
            class(graph) :: this
            type(graph) :: res, g1, g2, r1, r2
            integer, dimension(2) :: loc

            if(sum(this%degree) == 0) then
                res = this
                return
            endif

            loc = maxloc(this%weightedNetwork)

            g1 = this%removeIndex(loc(1))
            g2 = this%removeIndex(loc(2))
            r1 = g1%brute_minimalVertexCover()
            r2 = g2%brute_minimalVertexCover()

            ! if(r1%size > r2%size) then
            if(sum(r1%individualScore) > sum(r1%individualScore)) then
                res = r1
                return
            endif
            res = r2
        end function

        recursive function removeEdge_minimalVertexCover(this) result(res)
            implicit none
            class(graph) :: this
            type(graph) :: res, g1
            integer, dimension(2) :: loc

            if(all(this%degree == 0)) then
                res = this
                return
            endif

            g1 = this
            do while(.not. all(g1%degree == 0)) 
                loc = maxloc(g1%weightedNetwork)
                g1 = g1%removeIndexes([loc(1), loc(2)])
            enddo

            res = g1
        end function

        function getMinimalCover(this) result(res)
            class(graph) :: this
            type(graph), dimension(:), allocatable :: graphs, covGraphs
            type(graph) :: res
            integer :: i

            graphs = this%getSubGraphs()

            allocate(covGraphs(size(graphs)))
            do i = 1, size(graphs)
                if(graphs(i)%size < 15) then
                    covGraphs(i) = brute_minimalVertexCover(graphs(i))
                else
                    covGraphs(i) = removeEdge_minimalVertexCover(graphs(i))
                endif
            enddo
            res = mergeGraphs(covGraphs)


        end function

        function mergeGraphs(graphs) result(res)
            type(graph), dimension(:) :: graphs
            type(graph) :: res
            integer, dimension(:), allocatable :: id
            real(kind=real64), dimension(:), allocatable :: individualScore
            integer, dimension(:,:), allocatable :: network
            integer :: currentIndex, i, totalSize
            currentIndex = 1
            totalSize = 0

            do i = 1, size(graphs)
                totalSize = totalSize + graphs(i)%size
            enddo

            allocate(id(totalSize))
            allocate(individualScore(totalSize))
            allocate(network(totalSize, totalSize))
            individualScore = 0
            id = 0
            network = 0
            do i = 1, size(graphs)
                id(currentIndex:(currentIndex + graphs(i)%size -1)) = graphs(i)%id
                individualScore(currentIndex:(currentIndex + graphs(i)%size -1)) = graphs(i)%individualScore
                network(currentIndex:(currentIndex + graphs(i)%size -1), currentIndex:(currentIndex + graphs(i)%size -1)) = graphs(i)%network(:,:)
                currentIndex = currentIndex + graphs(i)%size
            enddo
            call res%initializeGraph(id, network, individualScore)

        end function

        function getSubGraphs(this) result(resArray)
            implicit none
            class(graph) :: this
            type(graph), dimension(:), allocatable :: resArray
            integer, dimension(:,:), allocatable :: network
            integer, dimension(this%size) :: assignment 
            integer, dimension(:), allocatable :: indexes, tempIndex

            integer :: currentValue, index, nGraphs, i
            
            assignment = 0
            currentValue = 0
            
            network = this%network

            do while (minval(assignment) == 0)
                currentValue = currentValue + 1
                tempIndex = minloc(assignment)
                index = tempIndex(1)
                assignment(index) = currentValue
                assignment = addNeighborsToAssignment(assignment, index, network)
            enddo   
            nGraphs = currentValue
            allocate(resArray(nGraphs))

            do i = 1, nGraphs
                indexes = selectIndexesBasedOnMask(assignment == i)
                call resArray(i)%initializeGraph(this%id(indexes), this%network(indexes, indexes), this%individualScore(indexes))
            enddo

        end function
        
        recursive function addNeighborsToAssignment(assignment, index, network) result(newAssignement)
            implicit none
            integer, dimension(:,:), allocatable :: network
            integer, dimension(:), intent(inout) :: assignment
            integer, dimension(:), allocatable :: newAssignement, neighbors 
            integer, dimension(:), allocatable :: newNeighbors
            integer :: index, i, currentValue
            
            currentValue = assignment(index)

            neighbors = selectIndexesBasedOnMask(network(index,:) > 0)

            newNeighbors = selectBasedOnMask(neighbors, assignment(neighbors) == 0)
            
            if(size(newNeighbors) == 0) then
                newAssignement = assignment
                return
            endif
            do i = 1, size(newNeighbors)
                assignment(newNeighbors(i)) = currentValue
                assignment = addNeighborsToAssignment(assignment, newNeighbors(i), network)
            enddo

            newAssignement = assignment
        end function

        function selectIndexesBasedOnMask(mask) result(subArray)
            implicit none
            logical, dimension(:) :: mask
            integer, dimension(:), allocatable :: subArray
            integer :: nElements, i, index
            
            nElements = count(mask)
            allocate(subArray(nElements))
            index = 1
            do i = 1, size(mask)
                if(mask(i)) then
                    subArray(index) = i
                    index = index + 1
                endif
            enddo
        end function

        function selectBasedOnMask(array, mask) result(subArray)
            implicit none
            integer, dimension(:) :: array
            logical, dimension(:) :: mask
            integer, dimension(:), allocatable :: subArray
            integer :: nElements, i, index
            
            nElements = count(mask)
            allocate(subArray(nElements))
            index = 1
            do i = 1, size(array)
                if(mask(i)) then
                    subArray(index) = array(i)
                    index = index + 1
                endif
            enddo
        end function

        function selectOthers(notArray, max) result(res)
            implicit none
            integer, dimension(:), intent(in) :: notArray
            integer, intent(in) :: max
            integer, dimension(:), allocatable :: res 

            integer :: i, numOut, index
            numOut = max - size(notArray)
            allocate(res(numOut))

            index = 1
            do i = 1, max 
                if(.not. any(notArray == i)) then
                    res(index) = i
                    index = index + 1
                endif
            enddo
        end function
end module


