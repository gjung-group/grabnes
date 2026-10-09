module sparse_matrix_module
    implicit none
contains
    subroutine sparse_matvec(nn, row_ptr, col_ind, values, x, y)
        implicit none
        integer, intent(in) :: nn, row_ptr(:), col_ind(:)
        complex*16, intent(in) :: values(:), x(:)
        complex*16, intent(out) :: y(nn)
        integer :: i, j

        y = (0.0, 0.0)
        do i = 1, nn
            print *, 'Row ', i, ':'
            do j = row_ptr(i), row_ptr(i+1) - 1
                if (col_ind(j) < 1 .or. col_ind(j) > nn) then
                    print *, 'Error: col_ind out of bounds: ', col_ind(j)
                    stop
                end if
                print *, '  Adding value ', values(j), ' * x(', col_ind(j), ') = ', x(col_ind(j))
                y(i) = y(i) + values(j) * x(col_ind(j))
            end do
        end do

        ! Debugging prints
        print *, 'sparse_matvec input x: ', x
        print *, 'sparse_matvec output y: ', y
    end subroutine sparse_matvec

    subroutine convert_to_csr(HLoc, N, row_ptr, col_ind, values)
        implicit none
        integer, intent(in) :: N
        complex*16, intent(in) :: HLoc(N, N)
        integer, allocatable, intent(out) :: row_ptr(:), col_ind(:)
        complex*16, allocatable, intent(out) :: values(:)
        integer :: nnz, k, i, j

        ! Count the number of non-zero elements
        nnz = count(HLoc /= (0.0, 0.0))

        ! Allocate CSR arrays
        allocate(row_ptr(N+1), col_ind(nnz), values(nnz))

        row_ptr(1) = 1
        k = 1

        ! Populate CSR arrays
        do i = 1, N
            do j = 1, N
                if (HLoc(i, j) /= (0.0, 0.0)) then
                    values(k) = HLoc(i, j)
                    col_ind(k) = j
                    print *, 'Adding to CSR: ', 'row', i, 'col', j, 'value', values(k)
                    k = k + 1
                end if
            end do
            row_ptr(i+1) = k
        end do

        ! Debugging prints
        print *, 'CSR row_ptr: ', row_ptr
        print *, 'CSR col_ind: ', col_ind
        print *, 'CSR values: ', values
    end subroutine convert_to_csr
end module sparse_matrix_module

program test_csr
    use sparse_matrix_module
    !use arpack
    implicit none
    integer, allocatable :: row_ptr(:), col_ind(:)
    complex*16, allocatable :: values(:), x(:), y(:)
    complex*16 :: HLoc(3, 3)
    integer :: nn, nev, ncv, lworkl, info, ido
    character(1) :: bmat
    character(2) :: which
    complex*16, allocatable :: resid(:), v(:,:), workd(:), workl(:), rwork(:), workev(:)
    integer, allocatable :: iparam(:), ipntr(:)
    complex*16, allocatable :: d(:), z(:,:)
    logical, allocatable :: select(:)
    complex*16 :: sigma

    ! Initialize matrix and vectors
    nn = 3
    HLoc = reshape([(1.0, 0.0), (0.0, 0.0), (2.0, 0.0), &
                    (0.0, 0.0), (3.0, 0.0), (0.0, 0.0), &
                    (4.0, 0.0), (0.0, 0.0), (5.0, 0.0)], shape=[3, 3])
    x = [(1.0, 0.0), (2.0, 0.0), (3.0, 0.0)]
    allocate(y(nn))

    ! Convert HLoc to CSR format
    call convert_to_csr(HLoc, nn, row_ptr, col_ind, values)

    ! ARPACK parameters
    nev = 2
    ncv = 2*nev
    lworkl = 3*ncv**2 + 5*ncv
    allocate(resid(nn), v(nn, ncv), workd(3*nn), workl(lworkl), rwork(ncv), d(nev), z(nn, nev), iparam(11), ipntr(14), select(ncv), workev(2*ncv))

    bmat = 'I'
    which = 'SM'  ! Smallest magnitude
    ido = 0
    iparam = 0
    iparam(1) = 1
    iparam(3) = 1000
    iparam(7) = 1
    sigma = (0.0, 0.0)

    ! ARPACK iteration
    do while (ido /= 99)
        call znaupd(ido, bmat, nn, which, nev, 1.0e-10, resid, ncv, v, nn, iparam, ipntr, workd, workl, lworkl, rwork, info)

        if (ido == -1 .or. ido == 1) then
            ! Perform sparse matrix-vector multiplication
            call sparse_matvec(nn, row_ptr, col_ind, values, &
                               workd(ipntr(1):ipntr(1)+nn-1), &
                               workd(ipntr(2):ipntr(2)+nn-1))
        else if (info /= 0) then
            print *, 'Error with znaupd, info = ', info
            stop
        end if
    end do

    if (info /= 0) then
        print *, 'Error with znaupd, info = ', info
        stop
    end if

    ! Extract eigenvalues and eigenvectors
    call zneupd(.false., 'A', select, d, z, nn, sigma, workev, bmat, nn, which, nev, 1.0e-10, resid, ncv, v, nn, iparam, ipntr, workd, workl, lworkl, rwork, info)
    if (info /= 0) then
        print *, 'Error with zneupd, info = ', info
        stop
    end if

    ! Print the eigenvalues
    print *, 'Eigenvalues: ', d

end program test_csr


