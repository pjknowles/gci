!>@brief The \c memory module supports the allocation, management and deallocation of Fortran pointers through several possible
!> implementations: Molpro's memory stack ('STACK'), or Fortran ALLOCATE ('HEAP')
!> Stack allocation should be used wherever possible, since it avoids memory fragmentation, and offers the possibility
!> to easily release on exit all the memory used in a particular code segment and its children.
!> Heap allocation can be used for data structures that need to persist.
!!
!! The principal routines are \c memory_allocate which creates a pointer to a new (possibly multidimensional) array,
!! and \c memory_release which destroys the previously allocated array,
!! and, in the case of stack implementation, resets the stack such that the array corresponding to the pointer given to it, and all
!! allocated subsequently, are destroyed.
!! Both these routines are overloaded interfaces to routines for 1-, 2-, 3- and 4-dimensional double-precision arrays;
!! \c memory_allocate_double is also provided, and is identical to \c memory_allocate.
!! \c memory_allocate_integer provides corresponding allocation of integer arrays, which can also be released
!! with \c memory_release.
!! Optionally, explicit lower bounds can be assigned to the pointer arrays.
!! There is a single generic \c memory_release that accepts any of the created double precision or integer pointers.
!!
!! \c memory_resize is provided for changing the size of an array. It works by allocating a new array, copying the
!! old to the new, and, in the case of HEAP implementation only, releasing the old, and returning the pointer to the
!! new in its place. \c memory_resize can also be used for STACK arrays but no attempt is made to free or reuse the original
!! storage unless the array is one-dimensional and either at the top of the stack or smaller than the original.
!!
!! Although the explicit dimension-dependent routines are documented here, they are declared private to the module,
!! and should be called via the PUBLIC GENERIC overloaded
!! forms \c memory_allocate, \c memory_allocate_integer,
!! \c memory_release, \c memory_resize
!! using the appropriate documented arguments.
!! Note that the pointers returned by \c memory_allocate and \c memory_allocate_integer, as well as those accepted
!! by \c memory_resize, are declared to have the \c CONTIGUOUS attribute. The calling routine should normally declare
!! the associated variable as \c CONTIGUOUS too; if this is not done with \c memory_resize, compilation will fail, typically
!! with the message that no matching subroutine definition can be found.
!!
!! The case of arrays of zero length is handled specially. \c
!! memory_allocate and \c memory_allocate_integer will simply return
!! a pointer to a zero-length array, and will not record this in any
!! tables. As a consequence, a subsequent \c memory_release on the
!! pointer has no effect at all, and in particular zero-length
!! allocations can not be used as a way of saving the state of the
!! stack. Instead, \c memory_save() is provided for this purpose, and
!! the integer that it delivers can be passed to a subsequent call to
!! \c memory_release, with the effect of deallocating all arrays
!! (stack and heap) created after the call to \c memory_save.
!!
!! The routines in this module are completely interoperable
!! with Molpro's historic stack interface, i.e. \c icorr and \c corlsr, which is
!! also implemented in the module, with unmodularised wrappers provided too.
!! However, using a consistent programming style for arrays is encouraged within a
!! module or subroutine.  To assist the communication with existing Molpro
!! code using the historic stack interface, functions \c memory_pointer_to_q
!! and \c memory_pointer_to_iq are provided to create pointer arrays pointing
!! to arrays allocated through \c icorr and \c corlsr. In addition,
!! the function \c memory_get_stack_position can be used to get the index of
!! a stack array allocated by this module on the stacks `q` or `iq`.
!!
!>Example of use:
!>\code
!> DOUBLE PRECISION, POINTER, DIMENSION(:)   :: array1
!> DOUBLE PRECISION, POINTER, DIMENSION(:,:) :: array2
!> array1 => memory_allocate(3)
!> array1=1; array1(2)=2; print *, array1 ! produces 1,2,1
!> array2 => memory_allocate(4,5) ! allocate a two-dimensional array
!> array2(2,3)=value
!> call memory_release(array1) ! resets the stack to where it was before array1 was allocated, and nullifies both pointers array1
!> and array2
!>\endcode
!>
!>Example of use of legacy routines:
!>\code
!> include "common/big"
!> ir1 = icorr(n1) ! allocate a double-precision array of length \c n1 starting at \c q(ir1)
!> ii2 = icorr(n2) ! allocate a integer array of length \c n2 starting at \c iq(ii2)
!> q(ir1:ir1-1+n1) = 0d0
!> iq(ii2:ii2-1+n2) = 0
!> call corlsr(ir1) ! resets the stack head to where it was before \c ir1 was allocated
!>                  ! with the effect that q(ir1:) and iq(ii2:) arrays are destroyed.
!>\endcode
MODULE memory
! dependence on Molpro:
! common/tapes for iout
! Error('message','place') error exit
 USE iso_c_binding, ONLY : c_loc, c_f_pointer, c_ptr, c_associated, c_size_t
 USE iso_fortran_env, ONLY : character_storage_size
 IMPLICIT NONE
 PRIVATE
 PUBLIC :: memory_initialize, memory_clean, memory_close
 PUBLIC :: memory_initialized
 PUBLIC :: memory_default_implementation
 PUBLIC :: memory_allocate, memory_allocate_double, memory_allocate_integer
 PUBLIC :: memory_allocate_character, memory_allocate_character_array!, memory_allocate_character8
 PUBLIC :: memory_release
 PUBLIC :: memory_reset_maximum_stack
 PUBLIC :: memory_resize
 PUBLIC :: memory_save
 PUBLIC :: memory_print_status, memory_remaining, memory_used
 PUBLIC :: memory_get_stack_position, memory_pointer_to_q, memory_pointer_to_iq, memory_register_stack_array
 PUBLIC :: memory_duplicate
 PUBLIC :: memory_top
 PUBLIC :: icorr, icori, corlsr, corlsi, icorrm, icorim, icorhw, icorhwres

!> \public Create a new double precision array on the stack, and return a 1-dimensional pointer to it.
 INTERFACE memory_allocate
  MODULE PROCEDURE memory_allocate_generic,&
       memory_allocate_double_1, memory_allocate_double_2, memory_allocate_double_3, memory_allocate_double_4,&
       memory_allocate_double_5
 END INTERFACE memory_allocate
!> \public Create a new double precision array on the stack, and return a 1-dimensional pointer to it.
 INTERFACE memory_allocate_double
  MODULE PROCEDURE memory_allocate_generic,&
       memory_allocate_double_1, memory_allocate_double_2, memory_allocate_double_3, memory_allocate_double_4,&
       memory_allocate_double_5
 END INTERFACE memory_allocate_double
!> \public Create a new integer array on the stack, and return a 1-dimensional pointer to it.
 INTERFACE memory_allocate_integer
  MODULE PROCEDURE memory_allocate_integer_generic,&
       memory_allocate_integer_1, memory_allocate_integer_2, memory_allocate_integer_3, memory_allocate_integer_4
 END INTERFACE memory_allocate_integer

!> \public Allocate a new array and fill it with the contents of an existing array
 INTERFACE memory_duplicate
  MODULE PROCEDURE memory_duplicate_double_1
  MODULE PROCEDURE memory_duplicate_integer_1
  MODULE PROCEDURE memory_duplicate_char_array
  MODULE PROCEDURE memory_duplicate_character
 END INTERFACE memory_duplicate

!> \public Resize a double precision array on the stack, and return a 1-dimensional pointer to it.
 INTERFACE memory_resize
  MODULE PROCEDURE &
       memory_resize_double_1, memory_resize_double_2, memory_resize_double_3, memory_resize_double_4
  MODULE PROCEDURE &
       memory_resize_integer_1, memory_resize_integer_2, memory_resize_integer_3, memory_resize_integer_4
  MODULE PROCEDURE &
       memory_resize_character, memory_resize_character_array
 END INTERFACE memory_resize

!> \public Release the memory allocated to a given pointer.
!! In the case of \c STACK implementation, also reset the top of the stack to be just below the given pointer,
!! thereby releasing memory allocated to other pointers allocated later than the given one.
 INTERFACE memory_release
  MODULE PROCEDURE memory_release_double_1, memory_release_double_2, memory_release_double_3, memory_release_double_4,&
       memory_release_double_5
  MODULE PROCEDURE memory_release_integer_1, memory_release_integer_2, memory_release_integer_3, memory_release_integer_4
  MODULE PROCEDURE memory_release_character, memory_release_character_array
  MODULE PROCEDURE memory_release_saved
 END INTERFACE memory_release

!> \public Get a pointer array pointing to the legacy `q` stack.
 INTERFACE memory_pointer_to_q
  MODULE PROCEDURE memory_pointer_to_q_1, memory_pointer_to_q_2, memory_pointer_to_q_3, memory_pointer_to_q_4
 END INTERFACE memory_pointer_to_q

!> \public Get a pointer array pointing to the legacy `iq` stack.
 INTERFACE memory_pointer_to_iq
  MODULE PROCEDURE memory_pointer_to_iq_1, memory_pointer_to_iq_2, memory_pointer_to_iq_3, memory_pointer_to_iq_4
 END INTERFACE memory_pointer_to_iq

!> \public Get the stack position allocated to a given pointer (for compatibility with subroutines using `q` or `iq`).
!! An error will be raised if the pointer is not allocated with the \c STACK implementation by this module.
 INTERFACE memory_get_stack_position
  MODULE PROCEDURE memory_get_stack_position_d1, memory_get_stack_position_d2, memory_get_stack_position_d3,&
       memory_get_stack_position_d4
  MODULE PROCEDURE memory_get_stack_position_i1, memory_get_stack_position_i2, memory_get_stack_position_i3,&
       memory_get_stack_position_i4
 END INTERFACE memory_get_stack_position

 INTERFACE find_memory_entry
  MODULE PROCEDURE find_memory_entry1, find_memory_entry2, find_memory_entry3, find_memory_entry4
  MODULE PROCEDURE find_memory_entry_integer1, find_memory_entry_integer2, find_memory_entry_integer3, find_memory_entry_integer4
  MODULE PROCEDURE find_memory_entry_character, find_memory_entry_character_arr
 END INTERFACE find_memory_entry

 INTEGER, SAVE, PUBLIC :: memory_print=-1 !< How verbose to be in reporting memory use.
 INTEGER, PARAMETER :: maximum_array_rank=5

 INTEGER, PARAMETER :: implementation_stack=1, implementation_heap=2
 CHARACTER(len=5), DIMENSION(2), PARAMETER :: implementations=['STACK','HEAP ']
 INTEGER, SAVE :: default_implementation=implementation_stack

 INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: zero_sized_integer
 DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, TARGET :: zero_sized_double

 LOGICAL, PUBLIC :: memory_count_heap_size !< whether to allow stack to grow to full buffer size, or (buffer size - heap allocated)


 TYPE memory_entry
  DOUBLE PRECISION, POINTER, DIMENSION(:), CONTIGUOUS :: p
  CHARACTER(len=9) :: data_type !< \c 'double' or \c 'integer' or \c 'character'
  INTEGER :: data_length !< number of bytes occupied by one element of the array
  INTEGER :: rank
  INTEGER, DIMENSION(maximum_array_rank) :: lbound, ubound
  CHARACTER(len=128) :: description=' '
  integer :: implementation=implementation_stack
  INTEGER :: level ! only for stack
  INTEGER :: size
  INTEGER :: handle=0 !< an identifier
 END TYPE memory_entry

#ifdef MOLPRO
#define loc molpro_loc
 INTEGER, EXTERNAL :: loc
#endif
 LOGICAL :: initialised=.FALSE.
 INTEGER, SAVE :: heap_size
 INTEGER, SAVE :: maximum_heap_size
 DOUBLE PRECISION, DIMENSION(:), POINTER, SAVE :: memory_stack
 INTEGER, SAVE :: stack_size
 INTEGER, SAVE :: maximum_stack_size
 INTEGER(kind=c_size_t), SAVE :: legacy_stack_offset

 INTEGER, SAVE :: n_segments
 TYPE(memory_entry), DIMENSION(:), SAVE, TARGET, ALLOCATABLE :: segments

CONTAINS

!> Set or report default implementation
 FUNCTION memory_default_implementation(implementation)
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation
  CHARACTER(len=5) :: memory_default_implementation
  memory_default_implementation = implementations(default_implementation)
  IF (PRESENT(implementation)) THEN
   if (implementation.eq.'STACK' .or. implementation.eq.'stack') default_implementation = implementation_stack
   if (implementation.eq.'HEAP' .or. implementation.eq.'heap') default_implementation = implementation_heap
  END IF
 END FUNCTION memory_default_implementation

!> Report remaining stack memory
 FUNCTION memory_remaining(datatype)
  INTEGER :: memory_remaining
  CHARACTER(len=*), INTENT(in), OPTIONAL :: datatype !< \c 'double' (default) or \c 'integer'
  INTEGER :: datalengt, i
  CHARACTER(len=9) :: datatyp
  datatyp='double'; IF (PRESENT(datatype)) datatyp=datatype
  DO i=1,len_TRIM(datatyp)
   IF (IACHAR(datatyp(i:i)).LE.IACHAR('Z').AND.IACHAR(datatyp(i:i)).GE.IACHAR('A')) &
   datatyp(i:i)=CHAR(IACHAR(datatyp(i:i))+IACHAR('a')-IACHAR('A'))
  ENDDO
  IF (datatyp.EQ.'double') THEN
   datalengt=storage_size(0.0d0)/8
  ELSE IF (datatyp.EQ.'integer') THEN
   datalengt=storage_size(0)/8
  ELSE IF (datatyp.EQ.'character') THEN
   datalengt=1
  ELSE
   CALL Error('Unknown datatype "'//TRIM(datatyp)//'"','memory')
  END IF
  i=stack_size
  IF (memory_count_heap_size) i=stack_size+heap_size
  memory_remaining = (SIZE(memory_stack)-i)*storage_SIZE(0.0d0)/(8*datalengt)
 END FUNCTION memory_remaining


!> Report used memory
 FUNCTION memory_used(implementation, maximum)
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'stack' (default) or 'heap'
  LOGICAL, INTENT(in), optional :: maximum !< if true (default is false) report the maximum memory allocated to date
  INTEGER :: memory_used
  IF (.NOT.PRESENT(implementation) .OR. implementation.EQ.'stack' .OR. implementation.EQ.'STACK') THEN
  memory_used = stack_size
  IF (PRESENT(maximum)) THEN
   IF (maximum) memory_used = maximum_stack_size
  END IF
 ELSE IF (implementation.EQ.'heap' .OR. implementation.EQ.'HEAP') THEN
  memory_used = heap_size
  IF (PRESENT(maximum)) THEN
   IF (maximum) memory_used = maximum_heap_size
  END IF
 ELSE
  CALL Error('Unknown argument "'//TRIM(implementation)//'"','memory::memory_used')
 END IF
 END FUNCTION memory_used

!> Print the state of the memory
 SUBROUTINE memory_print_status (maximum_depth,title)
  INTEGER, INTENT(in), OPTIONAL :: maximum_depth
  CHARACTER(len=*), INTENT(in), OPTIONAL :: title
#ifdef MOLPRO
  INCLUDE "common/tapes"
#else
  INTEGER, PARAMETER :: iout=6
#endif
  TYPE(memory_entry) :: s
  TYPE(memory_entry) :: h
  INTEGER :: k,l,levelr
  INTEGER :: n_entries
  IF (.NOT. initialised) RETURN
  WRITE (iout,112)
112 FORMAT(40(' ='))
113 FORMAT(40(' -'))
  WRITE (iout,'('' Memory manager status, default implementation '',A)') TRIM(implementations(default_implementation))
  IF (PRESENT(title)) THEN
   WRITE (iout,'(1X, A)') TRIM(title)
  END IF
   DO k=n_segments,1,-1
   END DO
  n_entries=0; DO k=1,n_segments; IF (segments(k)%implementation.EQ.implementation_heap) n_entries=n_entries+1; END DO
  IF (heap_size.GT.0 .OR. default_implementation .EQ. implementation_heap) THEN
   WRITE (iout,'('' Heap status:'',I6,'' entries, memory used='',A,'' doubles'')') n_entries, TRIM(IntegerString(heap_size))
   IF (n_entries.GT.0) WRITE (iout,'(27X,''Type'',T47,''Address'',T58,''Size  Rank Bounds'')')
   IF (n_entries.GT.0) WRITE (iout,113)
   IF (PRESENT(maximum_depth)) n_entries=MIN(n_entries,maximum_depth)
   DO k=n_segments,1,-1
    if (segments(k)%implementation.ne.implementation_heap) cycle
    IF (n_entries.LE.0) EXIT
    n_entries=n_entries-1
    h=segments(k)
     IF (h%rank.GT.0) THEN
110   FORMAT(1X,A24,T27,1X,A,T37,I17,T54,I8,I6,15(:,2X,A,':',A))
      WRITE (iout,110) TRIM(h%description),h%data_type,loc(h%p(1)),h%size,h%rank&
           ,(TRIM(IntegerString(h%lbound(l))),TRIM(IntegerString(h%ubound(l))),l=1,h%rank)
     ELSE
      WRITE (iout,110) TRIM(h%description),h%data_type,h%size*storage_size(0d0)/storage_size(0),h%rank
     END IF
   END DO
    WRITE (iout,112)
  END IF
  n_entries=0; DO k=1,n_segments; IF (segments(k)%implementation.EQ.implementation_stack) n_entries=n_entries+1; END DO
  IF (n_entries.GT.0 .OR. default_implementation.EQ.implementation_stack) THEN
   levelr = stack_size+1
   WRITE (iout,'('' Stack status: Remaining memory='',A,'' doubles ('',A,'' currently used, '',A,'' maximum used)'')'&
        ) TRIM(IntegerString(memory_remaining())),&
        TRIM(IntegerString(memory_used('STACK',.FALSE.))), TRIM(IntegerString(memory_used('STACK',.TRUE.)))
   IF (n_entries.GT.0) WRITE (iout,'(27X,''Type'',T40,''Depth'',T57,''Address'',T70,''Size  Rank Bounds'')')
   if (n_entries.gt.0) WRITE (iout,113)
   IF (PRESENT(maximum_depth)) n_entries=MIN(n_entries,maximum_depth)
   DO k=n_segments,1,-1
    IF (segments(k)%implementation.NE.implementation_stack) CYCLE
    IF (n_entries.LE.0) EXIT
    n_entries = n_entries-1
    s=segments(k)
     IF (s%rank.GT.0) THEN
111   FORMAT(1X,A24,T26,1X,A,T35,I10,T45,I19,T66,I8,I6,15(:,2X,A,':',A))
      WRITE (iout,111) TRIM(s%description),s%data_type,levelr-s%level,loc(s%p(1)),s%size,s%rank&
           ,(TRIM(IntegerString(s%lbound(l))),TRIM(IntegerString(s%ubound(l))),l=1,s%rank)
     ELSE
      WRITE (iout,111) TRIM(s%description),s%data_type,levelr-s%level,s%size,s%rank
     END IF
   END DO
   WRITE (iout,112)
  END IF
  RETURN
 END SUBROUTINE memory_print_status
 FUNCTION IntegerString(n)
  CHARACTER(len=24) :: IntegerString, resul
  INTEGER, INTENT(in) :: n
  WRITE (resul,'(I24)') n
  IntegerString = ADJUSTL(resul)
 END FUNCTION IntegerString
 SUBROUTINE upcase(l)
  CHARACTER(*), INTENT(inout) :: l
  INTEGER, PARAMETER :: ia = ICHAR('a')
  INTEGER, PARAMETER :: iz = ICHAR('z')
  INTEGER, PARAMETER :: ishift = ICHAR('A')-ICHAR('a')
  INTEGER :: i, j
  DO i=1,LEN(l)
   j=ICHAR(l(i:i))
   IF (j.GE.ia.AND.j.LE.iz) l(i:i)=CHAR(j+ishift)
  END DO
 END SUBROUTINE upcase

!> Initialise the memory management system
 SUBROUTINE memory_initialize(stack_length)
  INTEGER, intent(in) :: stack_length !< the amount of memory to allocate to the stack
  !CALL memory_close ! maybe could instead implement copying of the old stack, and keeping the heap
  IF (ASSOCIATED(memory_stack)) DEALLOCATE(memory_stack)
  ALLOCATE (memory_stack(stack_length))
  stack_size = 0
  maximum_stack_size = 0


  IF (.NOT. ALLOCATED(segments)) THEN
   heap_size = 0
   maximum_heap_size = 0
   n_segments=0
   ALLOCATE(segments(1024))
  END IF
  IF (.NOT. ALLOCATED(zero_sized_integer)) ALLOCATE(zero_sized_integer(0))
  IF (.NOT. ALLOCATED(zero_sized_double)) ALLOCATE(zero_sized_double(0))
  initialised = .TRUE.
  memory_count_heap_size = .TRUE. ! default is to constrain total heap+stack not to exceed \c stack_length
 END SUBROUTINE memory_initialize

!> Reset the maximum stack used statistic to the currently-used stack
 SUBROUTINE memory_reset_maximum_stack(level)
  INTEGER, INTENT(in), OPTIONAL :: level !< if given, the stack position desired; default is current stack size
  IF (PRESENT(level)) THEN
   maximum_stack_size = level
  ELSE
   maximum_stack_size = stack_size
 END IF
 END SUBROUTINE memory_reset_maximum_stack

!> Enquire on whether memory management system is active
 FUNCTION memory_initialized()
  LOGICAL :: memory_initialized
  memory_initialized = initialised
 END FUNCTION memory_initialized

!> Close the memory management system
 SUBROUTINE memory_close(only_clear,heap,stack)
  LOGICAL, INTENT(in), OPTIONAL :: only_clear !< If set, only destroy all memory segments, and do not continue to destroy buffer
  LOGICAL, INTENT(in), optional :: heap !< If set (default is \c .true.) close down heap
  LOGICAL, INTENT(in), optional :: stack !< If set (default is \c .true.) close down stack
  LOGICAL :: hp,stk
  INTEGER :: i
  IF (.NOT. initialised) RETURN
  hp = .TRUE.; IF (PRESENT(heap)) hp = heap
  stk = .TRUE.; IF (PRESENT(stack)) stk = stack
  IF (heap_size.GT.0 .AND. hp) THEN
   CALL Warning('?? WARNING: some heap entries remain on closing down memory system; posible memory leak','memory')
   CALL memory_print_status
  END IF
  DO i=1,n_segments
   IF (.NOT. ASSOCIATED(segments(n_segments)%p)) CYCLE
   IF (segments(n_segments)%implementation.EQ.implementation_heap .AND. .NOT.hp) CYCLE
   IF (segments(n_segments)%implementation.EQ.implementation_stack .AND. .NOT.stk) CYCLE
   NULLIFY(segments(n_segments)%p)
  END DO
  DO WHILE (n_segments.GT.0 .AND. .NOT. ASSOCIATED(segments(n_segments)%p))
   n_segments = n_segments-1
  END DO
  IF (stk) stack_size=0
  IF (hp) heap_size=0
  IF (PRESENT(only_clear)) THEN
   IF (only_clear) RETURN
  END IF
  IF (stk .AND. ASSOCIATED(memory_stack)) DEALLOCATE(memory_stack)
  IF (n_segments.GT.0) RETURN
  DEALLOCATE(segments)
  DEALLOCATE(zero_sized_double)
  DEALLOCATE(zero_sized_integer)
  initialised = .FALSE.
END SUBROUTINE memory_close

SUBROUTINE extend_segments
! add another entry in the segments table, making the table bigger if necessary
 TYPE(memory_entry), DIMENSION(:),  ALLOCATABLE :: old_segments
 n_segments=n_segments+1
 IF (n_segments.GT.SIZE(segments)) THEN
  ALLOCATE(old_segments(1:SIZE(segments)))
  old_segments = segments
  DEALLOCATE(segments)
  ALLOCATE(segments(SIZE(old_segments)*2))
  segments(:SIZE(old_segments)) = old_segments
 END IF
END SUBROUTINE extend_segments


!> Clean up the maintained list of pointers, throwing away all those that are above the current top of the stack
SUBROUTINE memory_clean(stack_only)
 LOGICAL, INTENT(in), OPTIONAL :: stack_only
  INTEGER :: levelr
  INTEGER :: k,l
  levelr = stack_size+1
! remove dead entries from the heap and stack
! first look to see if the end of the segment table can be trashed
  DO k=n_segments,1,-1
   IF (segments(k)%implementation.EQ.implementation_stack .AND. segments(k)%level+segments(k)%size.LE.levelr) EXIT
   IF (segments(k)%implementation.EQ.implementation_heap .AND. ASSOCIATED(segments(k)%p)) EXIT
  END DO
  DO l=k+1,n_segments
   IF (ASSOCIATED(segments(l)%p))  NULLIFY(segments(l)%p)
  END DO
  n_segments=MIN(n_segments,k)
  IF (PRESENT(stack_only)) THEN
   IF (stack_only) RETURN
  END IF
! now clean up the rest if there are holes
  k=0
  DO WHILE (k.LT.n_segments)
   k=k+1
   IF (segments(k)%implementation.EQ.implementation_stack) THEN
    IF (segments(k)%size.GT.0 .AND. segments(k)%level+segments(k)%size .GT. levelr) THEN
     IF (ASSOCIATED(segments(k)%p))  NULLIFY(segments(k)%p)
     n_segments = n_segments-1
     DO l=k,n_segments
      segments(l) = segments(l+1)
     END DO
     k=k-1
    END IF
   ELSE IF (segments(k)%implementation.EQ.implementation_heap) THEN
    IF ( ASSOCIATED(segments(k)%p) ) CYCLE
    n_segments = n_segments-1
     heap_size = heap_size - segments(k)%size
    DO l=k,n_segments
     segments(l) = segments(l+1)
    END DO
    k=k-1
   END IF
  END DO
 END SUBROUTINE memory_clean

 SUBROUTINE stack_check_overflow(n, description)
  INTEGER, INTENT(in) :: n
  CHARACTER(len=*), INTENT(in), OPTIONAL :: description
  IF (memory_remaining().GE.n) RETURN
  CALL memory_print_status
  IF (PRESENT(description)) THEN
   CALL error('Insufficient memory to allocate array '//description//' of length '&
        & //trim(IntegerString(n))//' 8-byte words','memory')
  ELSE
   CALL error('Insufficient memory to allocate a new array of length '&
        & //trim(IntegerString(n))//' 8-byte words','memory')
  END IF
 END SUBROUTINE stack_check_overflow

!> \public Create a new double precision array on the stack, and return a 1-dimensional pointer to it.
!! The size of the array is specified by the dimensions of, in general, a multidimensional data structure,
!! including the possibility to specify notional lower bounds for indexes.
!! The pointer returned is always a 1-dimensional double precision; however, the first dimension will be
!! adjusted to the appropriate multiple of \c n(1) if \c datataype or \c datalength are specified,
!! such that the right number of 8-byte words are allocated.
!! The calling routine will normally then recast the pointer to the desired type.
 FUNCTION memory_allocate_generic(n,lb,description,implementation,datatype,datalength,maximum)
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: memory_allocate_generic
  CONTIGUOUS :: memory_allocate_generic
  INTEGER, INTENT(in), DIMENSION(:) :: n !< the lengths of each dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL, DIMENSION(:) :: lb !< the lower bounds in each dimension (default (/1,1,.../)).
!! the lower bounds have no effect on the size of the matrix calculated, but are remembered.
  CHARACTER(len=*), INTENT(in), OPTIONAL :: description !< string describing the matrix
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' (default) or 'HEAP'
  CHARACTER(len=*), INTENT(in), OPTIONAL :: datatype !< 'double' (default) or 'integer' or 'character'.
  INTEGER, INTENT(in), OPTIONAL :: datalength !< in the case of datatype=='character', the character string length.
  LOGICAL, INTENT(in), OPTIONAL :: maximum !< whether to make this allocation count
!! as a candidate for the maximum stack usage (default \c .TRUE.)

  LOGICAL :: maximum_
  INTEGER, DIMENSION(LBOUND(n,1):UBOUND(n,1)) :: lbb,nn
  INTEGER :: i, level, ntotal
  INTEGER :: implem
  CHARACTER(len=9) :: datatyp
  integer :: datalengt
  maximum_ = .TRUE.; IF (PRESENT(maximum)) maximum_=maximum
  lbb=(/(1,i=LBOUND(n,1),UBOUND(n,1))/); IF (PRESENT(lb)) lbb=lb
  IF (PRESENT(lb) .AND. (UBOUND(lbb,1).NE.UBOUND(n,1) .OR. LBOUND(lbb,1).NE.LBOUND(n,1)))&
       CALL Error('Wrong shape lower bound','memory::memory_allocate')
  nn=n
  IF (PRESENT(datatype)) THEN
   datatyp=datatype
   DO i=1,len_TRIM(datatyp)
    IF (IACHAR(datatyp(i:i)).LE.IACHAR('Z').AND.IACHAR(datatyp(i:i)).GE.IACHAR('A')) &
         datatyp(i:i)=CHAR(IACHAR(datatyp(i:i))+IACHAR('a')-IACHAR('A'))
   ENDDO
   datalengt=1; IF (PRESENT(datalength)) datalengt=datalength
   IF (datatyp.EQ.'double') THEN
    datalengt=storage_SIZE(0.0d0)/8
   ELSE IF (datatyp.EQ.'integer') THEN
    datalengt=storage_SIZE(0)/8
   ELSE IF (datatyp.EQ.'character') THEN
    CONTINUE
   ELSE
    CALL Error('Unknown datatype "'//TRIM(datatyp)//'"','memory')
   END IF
   nn(1)=(n(1)*datalengt-1)/(storage_SIZE(0.0d0)/8)+1
  ELSE
   datatyp='double'
   datalengt=8
  END IF
  !CALL memory_clean(.TRUE.)
  IF (SIZE(n).GT.maximum_array_rank) CALL Error('Attempt to allocate array of too high rank','memory::memory_allocate')
  ntotal = PRODUCT(nn)
  IF (ntotal == 0) THEN
! do not use the system if zero length, but simply return an appropriate pointer
! memory_release cannot be used successfully on the result
    memory_allocate_generic => zero_sized_double
    RETURN
  END IF
  implem=default_implementation
  IF (PRESENT(implementation)) THEN
   if (implementation.eq.'STACK' .or. implementation.eq.'stack') implem = implementation_stack
   if (implementation.eq.'HEAP' .or. implementation.eq.'heap') implem = implementation_heap
  END IF
  IF (implem.EQ.implementation_heap) THEN
   level=0
   ALLOCATE(memory_allocate_generic(1:ntotal))
  ELSE ! default: stack
   CALL stack_check_overflow(ntotal, description)
   level = stack_size+1
   stack_size = stack_size + ntotal
   IF (maximum_) maximum_stack_size = MAX(maximum_stack_size, stack_size)
   memory_allocate_generic=>memory_stack(level:stack_size)
  END IF
  CALL extend_segments
  segments(n_segments)%data_type=datatyp
  segments(n_segments)%data_length=datalengt
  segments(n_segments)%p=>memory_allocate_generic
  segments(n_segments)%rank=UBOUND(n,1)-LBOUND(n,1)+1
  segments(n_segments)%lbound=lbb
  segments(n_segments)%ubound=n+lbb-1
  segments(n_segments)%size=ntotal
  segments(n_segments)%description=' '
  IF (PRESENT(description)) segments(n_segments)%description=description
  segments(n_segments)%implementation=implem
  IF (implem.EQ.implementation_stack) segments(n_segments)%level = level
  IF (implem.EQ.implementation_heap) heap_size = heap_size + segments(n_segments)%size
  maximum_heap_size = MAX(maximum_heap_size, heap_size)
 END FUNCTION memory_allocate_generic

!> \public Create a new 1-dimensional character array on the stack, and return a pointer to it.
 FUNCTION memory_allocate_character_array(n,lb,description,implementation,maximum)
  CHARACTER(len=1), POINTER, DIMENSION(:) :: memory_allocate_character_array
  CONTIGUOUS :: memory_allocate_character_array
  INTEGER, INTENT(in) :: n !< the dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL :: lb !< the lower bounds (default 1).
!! the lower bounds have no effect on the size of the matrix calculated, but are remembered.
  CHARACTER(len=*), INTENT(in), OPTIONAL :: description !< string describing the matrix
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' (default) or 'HEAP'
  LOGICAL, INTENT(in), OPTIONAL :: maximum !< whether to make this allocation count
!! as a candidate for the maximum stack usage (default \c .TRUE.)

  DOUBLE PRECISION, POINTER, DIMENSION(:) :: db
  CHARACTER(len=1), POINTER, DIMENSION(:) :: cb
  INTEGER :: lbb
  TYPE(c_ptr) :: cptr
  lbb=1; if (present(lb)) lbb=lb

  db => memory_allocate_generic((/n/),(/lbb/),description,implementation,'character',maximum=maximum)
  cptr = c_loc(db(1))
  CALL c_f_pointer(cptr,cb,shape=[n])
  memory_allocate_character_array(lbb:lbb+n-1) => cb
 END FUNCTION memory_allocate_character_array

!> \public Create a new character scalar on the stack, and return a pointer to it.
 FUNCTION memory_allocate_character(n,description,implementation,maximum)
  CHARACTER(len=:), POINTER :: memory_allocate_character
  INTEGER, INTENT(in) :: n !< the length of the result
!! the lower bounds have no effect on the size of the matrix calculated, but are remembered.
  CHARACTER(len=*), INTENT(in), OPTIONAL :: description !< string describing the matrix
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' (default) or 'HEAP'
  LOGICAL, INTENT(in), OPTIONAL :: maximum !< whether to make this allocation count
!! as a candidate for the maximum stack usage (default \c .TRUE.)

  DOUBLE PRECISION, POINTER, DIMENSION(:) :: db
  CHARACTER(len=1), POINTER, DIMENSION(:) :: cb
  TYPE(c_ptr) :: cptr

  db => memory_allocate_generic((/n/),(/1/),description,implementation,'character',maximum=maximum)
  cptr = c_loc(db(1))
  CALL c_f_pointer(cptr,cb,shape=[n])
  !memory_allocate_character => cb(1)
  CALL get_scalar_pointer(n,cb,memory_allocate_character)
 END FUNCTION memory_allocate_character
 SUBROUTINE get_scalar_pointer(scalar_len, scalar, ptr)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR
    INTEGER, INTENT(IN) :: scalar_len
    CHARACTER(KIND=C_CHAR,LEN=scalar_len), INTENT(IN), TARGET :: scalar(1)
    CHARACTER(:,KIND=C_CHAR), INTENT(OUT), POINTER :: ptr
    ptr => scalar(1)
  END SUBROUTINE get_scalar_pointer




! !> \public Create a new 1-dimensional character*8 array on the stack, and return a pointer to it.
!  FUNCTION memory_allocate_character8(n,description,implementation)
!   CHARACTER(len=8), POINTER, DIMENSION(:) :: memory_allocate_character8
!   CONTIGUOUS :: memory_allocate_character8
!   INTEGER, INTENT(in) :: n !< the dimension of the matrix
!   CHARACTER(len=*), INTENT(in), OPTIONAL :: description !< string describing the matrix
!   CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' (default) or 'HEAP'

!   character(len=1), POINTER, dimension(:) :: c1
!   DOUBLE PRECISION, POINTER, dimension(:) :: db
!   db => memory_allocate_generic((/n*8/),(/1/),description,implementation,'character')
!   PRINT *,size(db)
!   CALL c_f_pointer(c_loc(db),c1,shape=[n*8])
!  END FUNCTION memory_allocate_character8


!> \public Create a new integer array on the stack, and return a 1-dimensional pointer to it.
!! The size of the array is specified by the dimensions of, in general, a multidimensional data structure,
!! including the possibility to specify notional lower bounds for indexes.
 FUNCTION memory_allocate_integer_generic(n,lb,description,implementation,maximum)
  INTEGER, POINTER, DIMENSION(:) :: memory_allocate_integer_generic
  CONTIGUOUS :: memory_allocate_integer_generic
  INTEGER, INTENT(in), DIMENSION(:) :: n !< the lengths of each dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL, DIMENSION(:) :: lb !< the lower bounds in each dimension (default (/1,1,.../).
!! the lower bounds have no effect on the size of the matrix calculated, but are remembered.
  CHARACTER(len=*), INTENT(in), OPTIONAL :: description !< string describing the matrix
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' (default) or 'HEAP'
  LOGICAL, INTENT(in), OPTIONAL :: maximum !< whether to make this allocation count
!! as a candidate for the maximum stack usage (default \c .TRUE.)

  INTEGER, POINTER, DIMENSION(:) :: ip
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: db,db2
  TYPE(c_ptr) :: cptr, cptr2

  db => memory_allocate_generic(n,lb,description,implementation,'integer',maximum=maximum)
  cptr = c_loc(db(1))
  CALL c_f_pointer(cptr,ip,shape=[PRODUCT(n)])
  cptr2 = c_loc(ip(1))
  call c_f_pointer(cptr2,db2,shape=[product(n)])
  memory_allocate_integer_generic => ip
 END FUNCTION memory_allocate_integer_generic

!> \public Create a new double precision array on the stack, and return a 1-dimensional pointer to it.
 FUNCTION memory_allocate_double_1(n1,lb,description,implementation,maximum)
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: memory_allocate_double_1
  CONTIGUOUS :: memory_allocate_double_1
  INTEGER, INTENT(in) :: n1 !< the length of the matrix
  INTEGER, INTENT(in), OPTIONAL,DIMENSION(:) :: lb!< the lower bound of the index (default (/1/).
  CHARACTER(len=*), INTENT(in), OPTIONAL :: description !< string describing the matrix
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' (default) or 'HEAP'
  LOGICAL, INTENT(in), OPTIONAL :: maximum !< whether to make this allocation count
!! as a candidate for the maximum stack usage (default \c .TRUE.)
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: x
  INTEGER, DIMENSION(1) :: lbb
  IF (n1.LE.0) THEN
   memory_allocate_double_1 => zero_sized_double
   RETURN
  END IF
  IF (PRESENT(lb)) THEN
    lbb=lb
  ELSE
    lbb=(/1/)
  END IF
  x => memory_allocate_generic((/n1/),lbb,description,implementation,maximum=maximum)
  memory_allocate_double_1(lbb(1):lbb(1)-1+n1) => x
 END FUNCTION memory_allocate_double_1

!> \public Create a new double-precision 2-dimensional array on the stack, and return a 2-dimensional pointer to it.
 FUNCTION memory_allocate_double_2(n1,n2,lb,description,implementation,maximum)
  DOUBLE PRECISION, POINTER, DIMENSION(:,:) :: memory_allocate_double_2
  CONTIGUOUS :: memory_allocate_double_2
  INTEGER, INTENT(in) :: n1 !< the length of the first dimension of the matrix
  INTEGER, INTENT(in) :: n2 !< the length of the second dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL,DIMENSION(2) :: lb !< the lower bounds in each dimension (default (/1,1/).
  CHARACTER(len=*), INTENT(in), OPTIONAL :: description !< string describing the matrix
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' (default) or 'HEAP'
  LOGICAL, INTENT(in), OPTIONAL :: maximum !< whether to make this allocation count
!! as a candidate for the maximum stack usage (default \c .TRUE.)
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: x
  INTEGER, DIMENSION(2) :: lbb
  IF (PRESENT(lb)) THEN
    lbb=lb
  ELSE
    lbb=(/1,1/)
  END IF
  x => memory_allocate_generic((/n1,n2/),lbb,description,implementation,maximum=maximum)
  memory_allocate_double_2(lbb(1):lbb(1)-1+n1,lbb(2):lbb(2)-1+n2) => x
 END FUNCTION memory_allocate_double_2

!> \public Create a new double-precision 3-dimensional array on the stack, and return a 3-dimensional pointer to it.
 FUNCTION memory_allocate_double_3(n1,n2,n3,lb,description,implementation,maximum)
  DOUBLE PRECISION, POINTER, DIMENSION(:,:,:) :: memory_allocate_double_3
  CONTIGUOUS :: memory_allocate_double_3
  INTEGER, INTENT(in) :: n1 !< the length of the first dimension of the matrix
  INTEGER, INTENT(in) :: n2 !< the length of the second dimension of the matrix
  INTEGER, INTENT(in) :: n3 !< the length of the third dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL,DIMENSION(3) :: lb !< the lower bounds in each dimension (default (/1,1,1/).
  CHARACTER(len=*), INTENT(in), OPTIONAL :: description !< string describing the matrix
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' (default) or 'HEAP'
  LOGICAL, INTENT(in), OPTIONAL :: maximum !< whether to make this allocation count
!! as a candidate for the maximum stack usage (default \c .TRUE.)
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: x
  INTEGER, DIMENSION(3) :: lbb
  IF (PRESENT(lb)) THEN
    lbb=lb
  ELSE
    lbb=(/1,1,1/)
  END IF
  x =>memory_allocate_generic((/n1,n2,n3/),lbb,description,implementation,maximum=maximum)
  memory_allocate_double_3(lbb(1):lbb(1)-1+n1,lbb(2):lbb(2)-1+n2,lbb(3):lbb(3)-1+n3) => x
 END FUNCTION memory_allocate_double_3

!> \public Create a new double-precision 4-dimensional array on the stack, and return a 4-dimensional pointer to it.
 FUNCTION memory_allocate_double_4(n1,n2,n3,n4,lb,description,implementation,maximum)
  DOUBLE PRECISION, POINTER, DIMENSION(:,:,:,:) :: memory_allocate_double_4
  CONTIGUOUS :: memory_allocate_double_4
  INTEGER, INTENT(in) :: n1 !< the length of the first dimension of the matrix
  INTEGER, INTENT(in) :: n2 !< the length of the second dimension of the matrix
  INTEGER, INTENT(in) :: n3 !< the length of the third dimension of the matrix
  INTEGER, INTENT(in) :: n4 !< the length of the fourth dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL,DIMENSION(4) :: lb !< the lower bounds in each dimension (default (/1,1,1,1/).
  CHARACTER(len=*), INTENT(in), OPTIONAL :: description !< string describing the matrix
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' (default) or 'HEAP'
  LOGICAL, INTENT(in), OPTIONAL :: maximum !< whether to make this allocation count
!! as a candidate for the maximum stack usage (default \c .TRUE.)
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: x
  INTEGER, DIMENSION(4) :: lbb
  IF (PRESENT(lb)) THEN
    lbb=lb
  ELSE
    lbb=(/1,1,1,1/)
  END IF
  x =>memory_allocate_generic((/n1,n2,n3,n4/),lbb,description,implementation,maximum=maximum)
  memory_allocate_double_4(lbb(1):lbb(1)-1+n1,lbb(2):lbb(2)-1+n2,lbb(3):lbb(3)-1+n3,lbb(4):lbb(4)-1+n4) => x
 END FUNCTION memory_allocate_double_4


!> \public Create a new double-precision 5-dimensional array on the stack, and return a 5-dimensional pointer to it.
 FUNCTION memory_allocate_double_5(n1,n2,n3,n4,n5,lb,description,implementation,maximum)
  DOUBLE PRECISION, POINTER, DIMENSION(:,:,:,:,:) :: memory_allocate_double_5
  CONTIGUOUS :: memory_allocate_double_5
  INTEGER, INTENT(in) :: n1 !< the length of the first dimension of the matrix
  INTEGER, INTENT(in) :: n2 !< the length of the second dimension of the matrix
  INTEGER, INTENT(in) :: n3 !< the length of the third dimension of the matrix
  INTEGER, INTENT(in) :: n4 !< the length of the fourth dimension of the matrix
  INTEGER, INTENT(in) :: n5 !< the length of the fourth dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL,DIMENSION(5) :: lb !< the lower bounds in each dimension (default (/1,1,1,1,1/).
  CHARACTER(len=*), INTENT(in), OPTIONAL :: description !< string describing the matrix
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' (default) or 'HEAP'
  LOGICAL, INTENT(in), OPTIONAL :: maximum !< whether to make this allocation count
!! as a candidate for the maximum stack usage (default \c .TRUE.)
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: x
  INTEGER, DIMENSION(5) :: lbb
  IF (PRESENT(lb)) THEN
    lbb=lb
  ELSE
    lbb=(/1,1,1,1,1/)
  END IF
  x =>memory_allocate_generic((/n1,n2,n3,n4,n5/),lbb,description,implementation,maximum=maximum)
  memory_allocate_double_5(lbb(1):lbb(1)-1+n1,lbb(2):lbb(2)-1+n2,lbb(3):lbb(3)-1+n3,lbb(4):lbb(4)-1+n4,lbb(5):lbb(5)-1+n5) => x
 END FUNCTION memory_allocate_double_5

!> \public Create a new integer array on the stack, and return a 1-dimensional pointer to it.
 FUNCTION memory_allocate_integer_1(n1,lb,description,implementation,maximum)
  INTEGER, POINTER, DIMENSION(:) :: memory_allocate_integer_1
  CONTIGUOUS :: memory_allocate_integer_1
  INTEGER, INTENT(in) :: n1 !< the length of the matrix
  INTEGER, INTENT(in), OPTIONAL,DIMENSION(1) :: lb !< the lower bound of the index (default (/1/).
  CHARACTER(len=*), INTENT(in), OPTIONAL :: description !< string describing the matrix
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' (default) or 'HEAP'
  LOGICAL, INTENT(in), OPTIONAL :: maximum !< whether to make this allocation count
!! as a candidate for the maximum stack usage (default \c .TRUE.)
  INTEGER, POINTER, DIMENSION(:) :: x
  INTEGER, DIMENSION(1) :: lbb
  IF (n1.LE.0) THEN
   memory_allocate_integer_1 => zero_sized_integer
   RETURN
  END IF
  IF (PRESENT(lb)) THEN
    lbb=lb
  ELSE
    lbb=(/1/)
  END IF
  x => memory_allocate_integer_generic((/n1/),lbb,description,implementation,maximum=maximum)
  memory_allocate_integer_1(lbb(1):lbb(1)-1+n1) => x
 END FUNCTION memory_allocate_integer_1

!> \public Create a new integer 2-dimensional array on the stack, and return a 2-dimensional pointer to it.
 FUNCTION memory_allocate_integer_2(n1,n2,lb,description,implementation,maximum)
  INTEGER, POINTER, DIMENSION(:,:) :: memory_allocate_integer_2
  CONTIGUOUS :: memory_allocate_integer_2
  INTEGER, INTENT(in) :: n1 !< the length of the first dimension of the matrix
  INTEGER, INTENT(in) ::  n2 !< the length of the second dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL,DIMENSION(2) :: lb !< the lower bounds in each dimension (default (/1,1/).
  CHARACTER(len=*), INTENT(in), OPTIONAL :: description !< string describing the matrix
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' (default) or 'HEAP'
  LOGICAL, INTENT(in), OPTIONAL :: maximum !< whether to make this allocation count
!! as a candidate for the maximum stack usage (default \c .TRUE.)
  INTEGER, POINTER, DIMENSION(:) :: x
  INTEGER, DIMENSION(2) :: lbb
  IF (PRESENT(lb)) THEN
    lbb=lb
  ELSE
    lbb=(/1,1/)
  END IF
  x => memory_allocate_integer_generic((/n1,n2/),lbb,description,implementation,maximum=maximum)
  memory_allocate_integer_2(lbb(1):lbb(1)-1+n1,lbb(2):lbb(2)-1+n2) => x
 END FUNCTION memory_allocate_integer_2

!> \public Create a new integer 3-dimensional array on the stack, and return a 3-dimensional pointer to it.
 FUNCTION memory_allocate_integer_3(n1,n2,n3,lb,description,implementation,maximum)
  INTEGER, POINTER, DIMENSION(:,:,:) :: memory_allocate_integer_3
  CONTIGUOUS :: memory_allocate_integer_3
  INTEGER, INTENT(in) :: n1 !< the length of the first dimension of the matrix
  INTEGER, INTENT(in) :: n2 !< the length of the second dimension of the matrix
  INTEGER, INTENT(in) :: n3 !< the length of the third dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL,DIMENSION(3) :: lb !< the lower bounds in each dimension (default (/1,1,1/).
  CHARACTER(len=*), INTENT(in), OPTIONAL :: description !< string describing the matrix
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' (default) or 'HEAP'
  LOGICAL, INTENT(in), OPTIONAL :: maximum !< whether to make this allocation count
!! as a candidate for the maximum stack usage (default \c .TRUE.)
  INTEGER, POINTER, DIMENSION(:) :: x
  INTEGER, DIMENSION(3) :: lbb
  IF (PRESENT(lb)) THEN
    lbb=lb
  ELSE
    lbb=(/1,1,1/)
  END IF
  x =>memory_allocate_integer_generic((/n1,n2,n3/),lbb,description,implementation,maximum=maximum)
  memory_allocate_integer_3(lbb(1):lbb(1)-1+n1,lbb(2):lbb(2)-1+n2,lbb(3):lbb(3)-1+n3) => x
 END FUNCTION memory_allocate_integer_3

!> \public Create a new integer 4-dimensional array on the stack, and return a 4-dimensional pointer to it.
 FUNCTION memory_allocate_integer_4(n1,n2,n3,n4,lb,description,implementation,maximum)
  INTEGER, POINTER, DIMENSION(:,:,:,:) :: memory_allocate_integer_4
  CONTIGUOUS :: memory_allocate_integer_4
  INTEGER, INTENT(in) :: n1 !< the length of the first dimension of the matrix
  INTEGER, INTENT(in) :: n2 !< the length of the second dimension of the matrix
  INTEGER, INTENT(in) :: n3 !< the length of the third dimension of the matrix
  INTEGER, INTENT(in) :: n4 !< the length of the fourth dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL,DIMENSION(4) :: lb !< the lower bounds in each dimension (default (/1,1,1,1/).
  CHARACTER(len=*), INTENT(in), OPTIONAL :: description !< string describing the matrix
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' (default) or 'HEAP'
  LOGICAL, INTENT(in), OPTIONAL :: maximum !< whether to make this allocation count
!! as a candidate for the maximum stack usage (default \c .TRUE.)
  INTEGER, POINTER, DIMENSION(:) :: x
  INTEGER, DIMENSION(4) :: lbb
  IF (PRESENT(lb)) THEN
    lbb=lb
  ELSE
    lbb=(/1,1,1,1/)
  END IF
  x =>memory_allocate_integer_generic((/n1,n2,n3,n4/),lbb,description,implementation,maximum=maximum)
  memory_allocate_integer_4(lbb(1):lbb(1)-1+n1,lbb(2):lbb(2)-1+n2,lbb(3):lbb(3)-1+n3,lbb(4):lbb(4)-1+n4) => x
 END FUNCTION memory_allocate_integer_4

!> \public
 SUBROUTINE memory_release_double_1(p)
  DOUBLE PRECISION, POINTER, DIMENSION(:), INTENT(in) :: p
  INTEGER :: k
  TYPE(c_ptr) :: cptr1, cptr2
  IF (SIZE(p) == 0) RETURN
  cptr1 = c_loc(p(LBOUND(p,1)))
  DO k=n_segments,1,-1
   cptr2 = c_loc(segments(k)%p(1))
   IF ( ASSOCIATED(p,segments(k)%p) .OR. c_ASSOCIATED(cptr1,cptr2) ) THEN
    IF (segments(k)%implementation.EQ.implementation_stack) THEN
     stack_size=segments(k)%level-1
     maximum_stack_size = MAX(maximum_stack_size, stack_size)
    ELSE IF (segments(k)%implementation.EQ.implementation_heap) THEN
     heap_size = heap_size - segments(k)%size
     DEALLOCATE(segments(k)%p)
     NULLIFY(segments(k)%p)
     segments(k)%size=0
    ELSE
     CALL Error('Unexpected case reached','memory::memory_release')
    END IF
    CALL memory_clean(segments(k)%implementation.eq.implementation_stack)
    RETURN
   END IF
  END DO
! perhaps at the top of the stack
  cptr2 = c_loc(memory_stack(stack_size+1))
  IF (c_ASSOCIATED(cptr1,cptr2)) RETURN
! perhaps somewhere in stack, and we want to allow cavalier resetting of the top of the stack
  IF (loc(p(LBOUND(p,1))).GE.loc(memory_stack(1)) .AND. &
       loc(p(LBOUND(p,1))).LT.loc(memory_stack(SIZE(memory_stack)))) THEN
   stack_size = (loc(p(LBOUND(p,1)))-loc(memory_stack(1)))/(storage_SIZE(0d0)/storage_SIZE(' '))
   maximum_stack_size = MAX(maximum_stack_size, stack_size)
   CALL memory_clean(.TRUE.)
   RETURN
  END IF
! nothing found, so crash
  CALL memory_print_status(title='memory_release_double_1 processing pointer at '//IntegerString(INT(loc(p(LBOUND(p,1))))))
  PRINT *,loc(p(LBOUND(p,1))), loc(memory_stack(stack_size+1))
  CALL Error('Cannot find allocated pointer to release','memory::memory_release')
 END SUBROUTINE memory_release_double_1

!> \public
 SUBROUTINE memory_release_double_2(p)
  DOUBLE PRECISION, POINTER, DIMENSION(:,:), INTENT(in) :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: pd
  TYPE(c_ptr) :: cptr
  IF (SIZE(P) == 0) RETURN
! pgf95 can't cope without defining intermediate cptr like this
  cptr = c_loc(p(LBOUND(p,1),LBOUND(p,2)))
  CALL c_f_pointer(cptr,pd,shape=[SIZE(p)])
  CALL memory_release(pd)
 END SUBROUTINE memory_release_double_2

!> \public
 SUBROUTINE memory_release_double_3(p)
  DOUBLE PRECISION, POINTER, DIMENSION(:,:,:), INTENT(in) :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: pd
  TYPE(c_ptr) :: cptr
  IF (SIZE(P) == 0) RETURN
  cptr = c_loc(p(LBOUND(p,1),LBOUND(p,2),LBOUND(p,3)))
  CALL c_f_pointer(cptr,pd,shape=[SIZE(p)])
  CALL memory_release(pd)
 END SUBROUTINE memory_release_double_3

!> \public
 SUBROUTINE memory_release_double_4(p)
  DOUBLE PRECISION, POINTER, DIMENSION(:,:,:,:), INTENT(in) :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: pd
  TYPE(c_ptr) :: cptr
  IF (SIZE(P) == 0) RETURN
  cptr = c_loc(p(LBOUND(p,1),LBOUND(p,2),LBOUND(p,3),LBOUND(p,4)))
  CALL c_f_pointer(cptr,pd,shape=[SIZE(p)])
  CALL memory_release(pd)
 END SUBROUTINE memory_release_double_4

!> \public
 SUBROUTINE memory_release_double_5(p)
  DOUBLE PRECISION, POINTER, DIMENSION(:,:,:,:,:), INTENT(in) :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: pd
  TYPE(c_ptr) :: cptr
  IF (SIZE(P) == 0) RETURN
  cptr = c_loc(p(LBOUND(p,1),LBOUND(p,2),LBOUND(p,3),LBOUND(p,4),LBOUND(p,5)))
  CALL c_f_pointer(cptr,pd,shape=[SIZE(p)])
  CALL memory_release(pd)
 END SUBROUTINE memory_release_double_5


!> \public
 SUBROUTINE memory_release_integer_1(p)
  INTEGER, POINTER, DIMENSION(:), intent(in) :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: pd
  TYPE(c_ptr) :: cptr
  IF (SIZE(P) == 0) RETURN
  cptr = c_loc(p(LBOUND(p,1)))
  CALL c_f_pointer(cptr,pd,shape=[(SIZE(p)-1)/(storage_size(0d0)/storage_size(0))+1])
  CALL memory_release(pd)
 END SUBROUTINE memory_release_integer_1

!> \public
 SUBROUTINE memory_release_integer_2(p)
  INTEGER, POINTER, DIMENSION(:,:), INTENT(in) :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: pd
  TYPE(c_ptr) :: cptr
  IF (SIZE(P) == 0) RETURN
  cptr = c_loc(p(LBOUND(p,1),LBOUND(p,2)))
  CALL c_f_pointer(cptr,pd,shape=[(SIZE(p)-1)/(storage_size(0d0)/storage_size(0))+1])
  CALL memory_release(pd)
 END SUBROUTINE memory_release_integer_2

!> \public
 SUBROUTINE memory_release_integer_3(p)
  INTEGER, POINTER, DIMENSION(:,:,:), INTENT(in) :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: pd
  TYPE(c_ptr) :: cptr
  IF (SIZE(P) == 0) RETURN
  cptr = c_loc(p(LBOUND(p,1),LBOUND(p,2),LBOUND(p,3)))
  CALL c_f_pointer(cptr,pd,shape=[(SIZE(p)-1)/(storage_size(0d0)/storage_size(0))+1])
  CALL memory_release(pd)
 END SUBROUTINE memory_release_integer_3

!> \public
 SUBROUTINE memory_release_integer_4(p)
  INTEGER, POINTER, DIMENSION(:,:,:,:), INTENT(in) :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: pd
  TYPE(c_ptr) :: cptr
  IF (SIZE(p) == 0) RETURN
  cptr = c_loc(p(LBOUND(p,1),LBOUND(p,2),LBOUND(p,3),LBOUND(p,4)))
  CALL c_f_pointer(cptr,pd,shape=[(SIZE(p)-1)/(storage_size(0d0)/storage_size(0))+1])
  CALL memory_release(pd)
 END SUBROUTINE memory_release_integer_4

!> \public
 SUBROUTINE memory_release_character_array(p)
  CHARACTER(len=1), POINTER, DIMENSION(:), INTENT(in) :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: pd
  TYPE(c_ptr) :: cptr
  IF (SIZE(p) == 0) RETURN
  cptr = c_loc(p(LBOUND(p,1)))
  CALL c_f_pointer(cptr,pd,shape=[(SIZE(p)-1)/(storage_size(0d0)/character_storage_size)+1])
  CALL memory_release(pd)
 END SUBROUTINE memory_release_character_array

!> \public
 SUBROUTINE memory_release_character(p)
  CHARACTER(len=:), POINTER, INTENT(in) :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: pd
  CHARACTER(len=1), POINTER :: ptr
  TYPE(c_ptr) :: cptr
  IF (LEN(p) == 0) RETURN
  ptr => p(1:1) ! to circumvent bug in gfortran 4.8
  cptr = c_loc(ptr)
  CALL c_f_pointer(cptr,pd,shape=[(LEN(p)-1)/(storage_SIZE(0d0)/character_storage_size)+1])
  CALL memory_release(pd)
 END SUBROUTINE memory_release_character

!> Remember the state of the stack, with a view to coming
!!back to it via a call to \c memory_release
 FUNCTION memory_save()
  INTEGER :: memory_save
  memory_save = stack_size
 END FUNCTION memory_save

!> \public Release all stack memory allocated since call to \ref memory_save.
!!Call via the generic interface \c memory_release
 SUBROUTINE memory_release_saved(handle)
  INTEGER, INTENT(in) :: handle !< returned by a previous \ref memory_save
  IF (handle.GT.stack_size) CALL Error('Invalid attempt to release stack to a point above the current stack size','memory')
  stack_size = handle
  maximum_stack_size = MAX(maximum_stack_size, stack_size)
  CALL memory_clean
 END SUBROUTINE memory_release_saved

!> \public State whether an existing allocated pointer is at the top of the stack.
 FUNCTION memory_top(ptr)
  LOGICAL :: memory_top
  DOUBLE PRECISION, POINTER, INTENT(in), DIMENSION(:) :: ptr !< pointer to an array obtained from memory_allocate
  CONTIGUOUS :: ptr
  CLASS(memory_entry), POINTER :: t
  t => find_memory_entry(ptr)
  memory_top = t%implementation.EQ.implementation_stack .AND.  t%level+SIZE(t%p).EQ.stack_size+1
 END FUNCTION memory_top

!> \public Take an existing allocated pointer, and reallocate it, retaining any existing data that is within the new bounds.
!! If the existing pointer is on the heap, the old array will be released.
!! If the existing pointer is at the top of the stack, it will be adjusted in size without any copying of data,
!! provided that the lower bound is not changed too.
!! Similarly, if the existing pointer is on the stack and it is being reduced in size without a change in lower bound,
!! the existing storage will be reused without copying, and without releasing the unused memory.
 SUBROUTINE memory_resize_double_1(ptr, n1,lb,description,implementation,maximum)
  DOUBLE PRECISION, POINTER, INTENT(inout), DIMENSION(:) :: ptr !< pointer to an array obtained from memory_allocate
  CONTIGUOUS :: ptr
  INTEGER, INTENT(in) :: n1 !< the new length of the first dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL, DIMENSION(:) :: lb !< the new lower bounds in each dimension (default existing).
  CHARACTER(len=*), INTENT(in), OPTIONAL, TARGET :: description !< string describing the matrix; defaults to existing
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' or 'HEAP'; defaults to existing
  LOGICAL, INTENT(in), OPTIONAL :: maximum !< whether to make this allocation count
!! as a candidate for the maximum stack usage (default \c .TRUE.)

  LOGICAL :: maximum_
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: newptr
  CONTIGUOUS :: newptr
  CLASS(memory_entry), POINTER :: oldtype, newtype

  INTEGER, DIMENSION(1) :: lbb
  CHARACTER(len=:), POINTER :: desc
  INTEGER :: impl
  INTEGER :: i,j,k,l
  oldtype => find_memory_entry(ptr)
  maximum_ = .TRUE.; IF (PRESENT(maximum)) maximum_=maximum
  lbb=oldtype%lbound(:1); IF (PRESENT(lb)) lbb=lb
  desc => oldtype%description; IF (PRESENT(description)) desc => description
  impl=oldtype%implementation
  IF (PRESENT(implementation)) THEN
   if (implementation.eq.'STACK' .or. implementation.eq.'stack') impl = implementation_stack
   if (implementation.eq.'HEAP' .or. implementation.eq.'heap') impl = implementation_heap
  END IF
! attempt to do in-place if possible
  IF (impl.EQ.implementation_stack .AND. impl.EQ.oldtype%implementation .AND. &
       lbb(1).EQ.oldtype%LBOUND(1) .AND. (memory_top(oldtype%p) .OR. SIZE(oldtype%p).GE.n1) &
       ) THEN
   IF (oldtype%level+SIZE(oldtype%p).EQ.stack_size+1) THEN ! simply move the top of the stack
    stack_size = oldtype%level-1+n1
    IF (maximum_) maximum_stack_size = MAX(maximum_stack_size, stack_size)
   END IF
   oldtype%size = n1
   oldtype%LBOUND(1)=lbb(1)
   oldtype%UBOUND(1)=n1
   NULLIFY(oldtype%p)
   oldtype%p => memory_stack(oldtype%level:oldtype%level-1+n1)
   ptr => oldtype%p
  ELSE
   newptr => memory_allocate_double(n1,lbb,desc,implementations(impl))
! find oldtype again in case segments() was reallocated
   oldtype => find_memory_entry(ptr)
   newtype => find_memory_entry(newptr)

   DO i=MAX(oldtype%LBOUND(1),newtype%LBOUND(1)),MIN(oldtype%UBOUND(1),newtype%UBOUND(1))
    newptr(i)=ptr(i)
   END DO

   IF (oldtype%implementation.EQ.implementation_heap) CALL memory_release(ptr)
   ptr => newptr
  END IF
 END SUBROUTINE memory_resize_double_1

!> \public Take an existing allocated pointer, and reallocate it, retaining any existing data that is within the new bounds.
!! If the existing pointer is on the heap, the old array will be released.
 SUBROUTINE memory_resize_double_2(ptr, n1, n2,lb,description,implementation)
  DOUBLE PRECISION, POINTER, INTENT(inout), DIMENSION(:,:) :: ptr !< pointer to an array obtained from memory_allocate
  CONTIGUOUS :: ptr
  INTEGER, INTENT(in) :: n1 !< the new length of the first dimension of the matrix
  INTEGER, INTENT(in) :: n2 !< the new length of the second dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL, DIMENSION(:) :: lb !< the new lower bounds in each dimension (default existing).
  CHARACTER(len=*), INTENT(in), OPTIONAL, TARGET :: description !< string describing the matrix; defaults to existing
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' or 'HEAP'; defaults to existing

  DOUBLE PRECISION, POINTER, DIMENSION(:,:) :: newptr
  CONTIGUOUS :: newptr
  CLASS(memory_entry), POINTER :: oldtype, newtype

  INTEGER, DIMENSION(2) :: lbb
  CHARACTER(len=:), POINTER :: desc
  CHARACTER(len=5) :: impl
  INTEGER :: i,j,k,l
  oldtype => find_memory_entry(ptr)
  lbb=oldtype%lbound(:2); IF (PRESENT(lb)) lbb=lb
  desc => oldtype%description; IF (PRESENT(description)) desc => description
  impl=implementations(oldtype%implementation); IF (PRESENT(implementation)) impl = implementation
  newptr => memory_allocate_double(n1,n2,lbb,desc,impl)
! find oldtype again in case segments() was reallocated
  oldtype => find_memory_entry(ptr)
  newtype => find_memory_entry(newptr)

  DO j=MAX(oldtype%LBOUND(2),newtype%LBOUND(2)),MIN(oldtype%UBOUND(2),newtype%UBOUND(2))
   DO i=MAX(oldtype%LBOUND(1),newtype%LBOUND(1)),MIN(oldtype%UBOUND(1),newtype%UBOUND(1))
    newptr(i,j)=ptr(i,j)
   END DO
  END DO

  IF (oldtype%implementation.EQ.implementation_heap) CALL memory_release(ptr)
  ptr => newptr
 END SUBROUTINE memory_resize_double_2

!> \public Take an existing allocated pointer, and reallocate it, retaining any existing data that is within the new bounds.
!! If the existing pointer is on the heap, the old array will be released.
 SUBROUTINE memory_resize_double_3(ptr, n1, n2, n3 ,lb,description,implementation)
  DOUBLE PRECISION, POINTER, INTENT(inout), DIMENSION(:,:,:) :: ptr !< pointer to an array obtained from memory_allocate
  CONTIGUOUS :: ptr
  INTEGER, INTENT(in) :: n1 !< the new length of the first dimension of the matrix
  INTEGER, INTENT(in) :: n2 !< the new length of the second dimension of the matrix
  INTEGER, INTENT(in) :: n3 !< the new length of the third dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL, DIMENSION(:) :: lb !< the new lower bounds in each dimension (default existing).
  CHARACTER(len=*), INTENT(in), OPTIONAL, TARGET :: description !< string describing the matrix; defaults to existing
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' or 'HEAP'; defaults to existing

  DOUBLE PRECISION, POINTER, DIMENSION(:,:,:) :: newptr
  CONTIGUOUS :: newptr
  CLASS(memory_entry), POINTER :: oldtype, newtype

  INTEGER, DIMENSION(3) :: lbb
  CHARACTER(len=:), POINTER :: desc
  CHARACTER(len=5) :: impl
  INTEGER :: i,j,k,l
  oldtype => find_memory_entry(ptr)
  lbb=oldtype%lbound(:3); IF (PRESENT(lb)) lbb=lb
  desc => oldtype%description; IF (PRESENT(description)) desc => description
  impl=implementations(oldtype%implementation); IF (PRESENT(implementation)) impl = implementation
  newptr => memory_allocate_double(n1,n2,n3,lbb,desc,impl)
! find oldtype again in case segments() was reallocated
  oldtype => find_memory_entry(ptr)
  newtype => find_memory_entry(newptr)

  DO k=MAX(oldtype%LBOUND(3),newtype%LBOUND(3)),MIN(oldtype%UBOUND(3),newtype%UBOUND(3))
   DO j=MAX(oldtype%LBOUND(2),newtype%LBOUND(2)),MIN(oldtype%UBOUND(2),newtype%UBOUND(2))
    DO i=MAX(oldtype%LBOUND(1),newtype%LBOUND(1)),MIN(oldtype%UBOUND(1),newtype%UBOUND(1))
     newptr(i,j,k)=ptr(i,j,k)
    END DO
   END DO
  END DO

  IF (oldtype%implementation.EQ.implementation_heap) CALL memory_release(ptr)
  ptr => newptr
 END SUBROUTINE memory_resize_double_3

!> \public Take an existing allocated pointer, and reallocate it, retaining any existing data that is within the new bounds.
!! If the existing pointer is on the heap, the old array will be released.
 SUBROUTINE memory_resize_double_4(ptr, n1, n2, n3, n4 ,lb,description,implementation)
  DOUBLE PRECISION, POINTER, INTENT(inout), DIMENSION(:,:,:,:) :: ptr !< pointer to an array obtained from memory_allocate
  CONTIGUOUS :: ptr
  INTEGER, INTENT(in) :: n1 !< the new length of the first dimension of the matrix
  INTEGER, INTENT(in) :: n2 !< the new length of the second dimension of the matrix
  INTEGER, INTENT(in) :: n3 !< the new length of the third dimension of the matrix
  INTEGER, INTENT(in) :: n4 !< the new length of the third dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL, DIMENSION(:) :: lb !< the new lower bounds in each dimension (default existing).
  CHARACTER(len=*), INTENT(in), OPTIONAL, TARGET :: description !< string describing the matrix; defaults to existing
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' or 'HEAP'; defaults to existing

  DOUBLE PRECISION, POINTER, DIMENSION(:,:,:,:) :: newptr
  CONTIGUOUS :: newptr
  CLASS(memory_entry), POINTER :: oldtype, newtype

  INTEGER, DIMENSION(4) :: lbb
  CHARACTER(len=:), POINTER :: desc
  CHARACTER(len=5) :: impl
  INTEGER :: i,j,k,l
  oldtype => find_memory_entry(ptr)
  lbb=oldtype%lbound(:4); IF (PRESENT(lb)) lbb=lb
  desc => oldtype%description; IF (PRESENT(description)) desc => description
  impl=implementations(oldtype%implementation); IF (PRESENT(implementation)) impl = implementation
  newptr => memory_allocate_double(n1,n2,n3,n4,lbb,desc,impl)
! find oldtype again in case segments() was reallocated
  oldtype => find_memory_entry(ptr)
  newtype => find_memory_entry(newptr)

  DO l=MAX(oldtype%LBOUND(4),newtype%LBOUND(4)),MIN(oldtype%UBOUND(4),newtype%UBOUND(4))
   DO k=MAX(oldtype%LBOUND(3),newtype%LBOUND(3)),MIN(oldtype%UBOUND(3),newtype%UBOUND(3))
    DO j=MAX(oldtype%LBOUND(2),newtype%LBOUND(2)),MIN(oldtype%UBOUND(2),newtype%UBOUND(2))
     DO i=MAX(oldtype%LBOUND(1),newtype%LBOUND(1)),MIN(oldtype%UBOUND(1),newtype%UBOUND(1))
      newptr(i,j,k,l)=ptr(i,j,k,l)
     END DO
    END DO
   END DO
  END DO

  IF (oldtype%implementation.EQ.implementation_heap) CALL memory_release(ptr)
  ptr => newptr
 END SUBROUTINE memory_resize_double_4


!> \public Take an existing allocated pointer, and reallocate it, retaining any existing data that is within the new bounds.
!! If the existing pointer is on the heap, the old array will be released.
!! If the existing pointer is at the top of the stack, it will be adjusted in size without any copying of data,
!! provided that the lower bound is not changed too.
!! Similarly, if the existing pointer is on the stack and it is being reduced in size without a change in lower bound,
!! the existing storage will be reused without copying, and without releasing the unused memory.
 SUBROUTINE memory_resize_integer_1(ptr, n1,lb,description,implementation,maximum)
  INTEGER, POINTER, INTENT(inout), DIMENSION(:) :: ptr !< pointer to an array obtained from memory_allocate
  CONTIGUOUS :: ptr
  INTEGER, INTENT(in) :: n1 !< the new length of the first dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL, DIMENSION(:) :: lb !< the new lower bounds in each dimension (default existing).
  CHARACTER(len=*), INTENT(in), OPTIONAL, TARGET :: description !< string describing the matrix; defaults to existing
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' or 'HEAP'; defaults to existing
  LOGICAL, INTENT(in), OPTIONAL :: maximum !< whether to make this allocation count
!! as a candidate for the maximum stack usage (default \c .TRUE.)

  LOGICAL :: maximum_

  INTEGER, POINTER, DIMENSION(:) :: newptr
  CONTIGUOUS :: newptr
  CLASS(memory_entry), POINTER :: oldtype, newtype

  INTEGER, DIMENSION(1) :: lbb
  CHARACTER(len=:), POINTER :: desc
  CHARACTER(len=5) :: impl
  INTEGER :: i,j,k,l,n1r
  TYPE(c_ptr) :: cptr
  oldtype => find_memory_entry(ptr)
  maximum_ = .TRUE.; IF (PRESENT(maximum)) maximum_=maximum
  lbb=oldtype%lbound(:1); IF (PRESENT(lb)) lbb=lb
  desc => oldtype%description; IF (PRESENT(description)) desc => description
  impl=implementations(oldtype%implementation); IF (PRESENT(implementation)) impl = implementation
  n1r=((n1-1)/(storage_size(0d0)/storage_size(0))+1)
! attempt to do in-place if possible
  IF (impl.EQ.'STACK' .AND. impl.EQ.implementations(oldtype%implementation) .AND. &
       lbb(1).EQ.oldtype%LBOUND(1) .AND. (oldtype%level+SIZE(oldtype%p).EQ.stack_size+1 &
       .OR. SIZE(oldtype%p).GE.n1r) &
       ) THEN
   IF (oldtype%level+SIZE(oldtype%p).EQ.stack_size+1) THEN ! simply move the top of the stack
    stack_size = oldtype%level-1+n1r
    IF (maximum_) maximum_stack_size = MAX(maximum_stack_size, stack_size)
   END IF
   oldtype%size = n1r
   oldtype%LBOUND(1)=lbb(1)
   oldtype%UBOUND(1)=n1
   NULLIFY(oldtype%p)
   oldtype%p => memory_stack(oldtype%level:oldtype%level-1+n1r)
   cptr = c_loc(oldtype%p(LBOUND(oldtype%p,1)))
   CALL c_f_pointer(cptr,ptr,[n1])
  ELSE
   newptr => memory_allocate_integer(n1,lbb,desc,impl)
! find oldtype again in case segments() was reallocated
   oldtype => find_memory_entry(ptr)
   newtype => find_memory_entry(newptr)

   DO i=MAX(oldtype%LBOUND(1),newtype%LBOUND(1)),MIN(oldtype%UBOUND(1),newtype%UBOUND(1))
    newptr(i)=ptr(i)
   END DO

   IF (oldtype%implementation.EQ.implementation_heap) CALL memory_release(ptr)
   ptr => newptr
  END IF
  END SUBROUTINE memory_resize_integer_1

!> \public Take an existing allocated pointer, and reallocate it, retaining any existing data that is within the new bounds.
!! If the existing pointer is on the heap, the old array will be released.
 SUBROUTINE memory_resize_integer_2(ptr, n1, n2,lb,description,implementation)
  INTEGER, POINTER, INTENT(inout), DIMENSION(:,:) :: ptr !< pointer to an array obtained from memory_allocate
  CONTIGUOUS :: ptr
  INTEGER, INTENT(in) :: n1 !< the new length of the first dimension of the matrix
  INTEGER, INTENT(in) :: n2 !< the new length of the second dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL, DIMENSION(:) :: lb !< the new lower bounds in each dimension (default existing).
  CHARACTER(len=*), INTENT(in), OPTIONAL, TARGET :: description !< string describing the matrix; defaults to existing
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' or 'HEAP'; defaults to existing

  INTEGER, POINTER, DIMENSION(:,:) :: newptr
  CONTIGUOUS :: newptr
  CLASS(memory_entry), POINTER :: oldtype, newtype

  INTEGER, DIMENSION(2) :: lbb
  CHARACTER(len=:), POINTER :: desc
  CHARACTER(len=5) :: impl
  INTEGER :: i,j,k,l
  oldtype => find_memory_entry(ptr)
  lbb=oldtype%lbound(:2); IF (PRESENT(lb)) lbb=lb
  desc => oldtype%description; IF (PRESENT(description)) desc => description
  impl=implementations(oldtype%implementation); IF (PRESENT(implementation)) impl = implementation
  newptr => memory_allocate_integer(n1,n2,lbb,desc,impl)
! find oldtype again in case segments() was reallocated
  oldtype => find_memory_entry(ptr)
  newtype => find_memory_entry(newptr)

  DO j=MAX(oldtype%LBOUND(2),newtype%LBOUND(2)),MIN(oldtype%UBOUND(2),newtype%UBOUND(2))
   DO i=MAX(oldtype%LBOUND(1),newtype%LBOUND(1)),MIN(oldtype%UBOUND(1),newtype%UBOUND(1))
    newptr(i,j)=ptr(i,j)
   END DO
  END DO

  IF (oldtype%implementation.EQ.implementation_heap) CALL memory_release(ptr)
  ptr => newptr
 END SUBROUTINE memory_resize_integer_2

!> \public Take an existing allocated pointer, and reallocate it, retaining any existing data that is within the new bounds.
!! If the existing pointer is on the heap, the old array will be released.
 SUBROUTINE memory_resize_integer_3(ptr, n1, n2, n3 ,lb,description,implementation)
  INTEGER, POINTER, INTENT(inout), DIMENSION(:,:,:) :: ptr !< pointer to an array obtained from memory_allocate
  CONTIGUOUS :: ptr
  INTEGER, INTENT(in) :: n1 !< the new length of the first dimension of the matrix
  INTEGER, INTENT(in) :: n2 !< the new length of the second dimension of the matrix
  INTEGER, INTENT(in) :: n3 !< the new length of the third dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL, DIMENSION(:) :: lb !< the new lower bounds in each dimension (default existing).
  CHARACTER(len=*), INTENT(in), OPTIONAL, TARGET :: description !< string describing the matrix; defaults to existing
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' or 'HEAP'; defaults to existing

  INTEGER, POINTER, DIMENSION(:,:,:) :: newptr
  CONTIGUOUS :: newptr
  CLASS(memory_entry), POINTER :: oldtype, newtype

  INTEGER, DIMENSION(3) :: lbb
  CHARACTER(len=:), POINTER :: desc
  CHARACTER(len=5) :: impl
  INTEGER :: i,j,k,l
  oldtype => find_memory_entry(ptr)
  lbb=oldtype%lbound(:3); IF (PRESENT(lb)) lbb=lb
  desc => oldtype%description; IF (PRESENT(description)) desc => description
  impl=implementations(oldtype%implementation); IF (PRESENT(implementation)) impl = implementation
  newptr => memory_allocate_integer(n1,n2,n3,lbb,desc,impl)
! find oldtype again in case segments() was reallocated
  oldtype => find_memory_entry(ptr)
  newtype => find_memory_entry(newptr)

  DO k=MAX(oldtype%LBOUND(3),newtype%LBOUND(3)),MIN(oldtype%UBOUND(3),newtype%UBOUND(3))
   DO j=MAX(oldtype%LBOUND(2),newtype%LBOUND(2)),MIN(oldtype%UBOUND(2),newtype%UBOUND(2))
    DO i=MAX(oldtype%LBOUND(1),newtype%LBOUND(1)),MIN(oldtype%UBOUND(1),newtype%UBOUND(1))
     newptr(i,j,k)=ptr(i,j,k)
    END DO
   END DO
  END DO

  IF (oldtype%implementation.EQ.implementation_heap) CALL memory_release(ptr)
  ptr => newptr
 END SUBROUTINE memory_resize_integer_3

!> \public Take an existing allocated pointer, and reallocate it, retaining any existing data that is within the new bounds.
!! If the existing pointer is on the heap, the old array will be released.
 SUBROUTINE memory_resize_integer_4(ptr, n1, n2, n3, n4 ,lb,description,implementation)
  INTEGER, POINTER, INTENT(inout), DIMENSION(:,:,:,:) :: ptr !< pointer to an array obtained from memory_allocate
  CONTIGUOUS :: ptr
  INTEGER, INTENT(in) :: n1 !< the new length of the first dimension of the matrix
  INTEGER, INTENT(in) :: n2 !< the new length of the second dimension of the matrix
  INTEGER, INTENT(in) :: n3 !< the new length of the third dimension of the matrix
  INTEGER, INTENT(in) :: n4 !< the new length of the third dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL, DIMENSION(:) :: lb !< the new lower bounds in each dimension (default existing).
  CHARACTER(len=*), INTENT(in), OPTIONAL, TARGET :: description !< string describing the matrix; defaults to existing
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' or 'HEAP'; defaults to existing

  INTEGER, POINTER, DIMENSION(:,:,:,:) :: newptr
  CONTIGUOUS :: newptr
  CLASS(memory_entry), POINTER :: oldtype, newtype

  INTEGER, DIMENSION(4) :: lbb
  CHARACTER(len=:), POINTER :: desc
  CHARACTER(len=5) :: impl
  INTEGER :: i,j,k,l
  oldtype => find_memory_entry(ptr)
  lbb=oldtype%lbound(:4); IF (PRESENT(lb)) lbb=lb
  desc => oldtype%description; IF (PRESENT(description)) desc => description
  impl=implementations(oldtype%implementation); IF (PRESENT(implementation)) impl = implementation
  newptr => memory_allocate_integer(n1,n2,n3,n4,lbb,desc,impl)
! find oldtype again in case segments() was reallocated
  oldtype => find_memory_entry(ptr)
  newtype => find_memory_entry(newptr)

  DO l=MAX(oldtype%LBOUND(4),newtype%LBOUND(4)),MIN(oldtype%UBOUND(4),newtype%UBOUND(4))
   DO k=MAX(oldtype%LBOUND(3),newtype%LBOUND(3)),MIN(oldtype%UBOUND(3),newtype%UBOUND(3))
    DO j=MAX(oldtype%LBOUND(2),newtype%LBOUND(2)),MIN(oldtype%UBOUND(2),newtype%UBOUND(2))
     DO i=MAX(oldtype%LBOUND(1),newtype%LBOUND(1)),MIN(oldtype%UBOUND(1),newtype%UBOUND(1))
      newptr(i,j,k,l)=ptr(i,j,k,l)
     END DO
    END DO
   END DO
  END DO

  IF (oldtype%implementation.EQ.implementation_heap) CALL memory_release(ptr)
  ptr => newptr
 END SUBROUTINE memory_resize_integer_4

!> \public Take an existing allocated pointer, and reallocate it, retaining any existing data that is within the new bounds.
!! If the existing pointer is on the heap, the old array will be released.
!! If the existing pointer is at the top of the stack, it will be adjusted in size without any copying of data,
!! provided that the lower bound is not changed too.
!! Similarly, if the existing pointer is on the stack and it is being reduced in size without a change in lower bound,
!! the existing storage will be reused without copying, and without releasing the unused memory.
 SUBROUTINE memory_resize_character_array(ptr, n1,lb,description,implementation,maximum)
  CHARACTER(LEN=1), POINTER, INTENT(inout), DIMENSION(:) :: ptr !< pointer to an array obtained from memory_allocate
  CONTIGUOUS :: ptr
  INTEGER, INTENT(in) :: n1 !< the new length of the first dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL :: lb !< the new lower bound (default existing).
  CHARACTER(len=*), INTENT(in), OPTIONAL, TARGET :: description !< string describing the matrix; defaults to existing
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' or 'HEAP'; defaults to existing
  LOGICAL, INTENT(in), OPTIONAL :: maximum !< whether to make this allocation count
!! as a candidate for the maximum stack usage (default \c .TRUE.)

  LOGICAL :: maximum_

  CHARACTER(LEN=1), POINTER, DIMENSION(:) :: newptr
  CONTIGUOUS :: newptr
  CLASS(memory_entry), POINTER :: oldtype, newtype

  INTEGER :: lbb
  CHARACTER(len=:), POINTER :: desc
  CHARACTER(len=5) :: impl
  INTEGER :: i,j,k,l,n1r
  TYPE(c_ptr) :: cptr
  oldtype => find_memory_entry(ptr)
  maximum_ = .TRUE.; IF (PRESENT(maximum)) maximum_=maximum
  lbb=oldtype%lbound(1); IF (PRESENT(lb)) lbb=lb
  desc => oldtype%description; IF (PRESENT(description)) desc => description
  impl=implementations(oldtype%implementation); IF (PRESENT(implementation)) impl = implementation
  n1r=((n1-1)/(storage_size(0d0)/character_storage_size)+1)
! attempt to do in-place if possible
  IF (impl.EQ.implementations(implementation_stack) .AND. impl.EQ.implementations(oldtype%implementation) .AND. &
       lbb.EQ.oldtype%LBOUND(1) .AND. (oldtype%level+SIZE(oldtype%p).EQ.stack_size+1 &
       .OR. SIZE(oldtype%p).GE.n1r) &
       ) THEN
   IF (oldtype%level+SIZE(oldtype%p).EQ.stack_size+1) THEN ! simply move the top of the stack
    stack_size = oldtype%level-1+n1r
    if (maximum_) maximum_stack_size = MAX(maximum_stack_size, stack_size)
   END IF
   oldtype%size = n1r
   oldtype%LBOUND(1)=lbb
   oldtype%UBOUND(1)=n1+lbb-1
   NULLIFY(oldtype%p)
   oldtype%p => memory_stack(oldtype%level:oldtype%level-1+n1r)
   cptr = c_loc(oldtype%p(LBOUND(oldtype%p,1)))
   CALL c_f_pointer(cptr,ptr,[n1])
  ELSE
   newptr => memory_allocate_character_array(n1,lbb,desc,impl)
! find oldtype again in case segments() was reallocated
   oldtype => find_memory_entry(ptr)
   newtype => find_memory_entry(newptr)

   DO i=MAX(oldtype%LBOUND(1),newtype%LBOUND(1)),MIN(oldtype%UBOUND(1),newtype%UBOUND(1))
    newptr(i)=ptr(i)
   END DO

   IF (oldtype%implementation.EQ.implementation_heap) CALL memory_release(ptr)
   ptr => newptr
  END IF
  END SUBROUTINE memory_resize_character_array

!> \public Take an existing allocated pointer, and reallocate it, retaining any existing data that is within the new bounds.
!! If the existing pointer is on the heap, the old array will be released.
!! If the existing pointer is at the top of the stack, it will be adjusted in size without any copying of data,
!! provided that the lower bound is not changed too.
!! Similarly, if the existing pointer is on the stack and it is being reduced in size without a change in lower bound,
!! the existing storage will be reused without copying, and without releasing the unused memory.
 SUBROUTINE memory_resize_character(ptr, n1,description,implementation)
  CHARACTER(LEN=:), POINTER, INTENT(inout) :: ptr !< pointer to an array obtained from memory_allocate
  INTEGER, INTENT(in) :: n1 !< the new length of the first dimension of the matrix
  CHARACTER(len=*), INTENT(in), OPTIONAL, TARGET :: description !< string describing the matrix; defaults to existing
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' or 'HEAP'; defaults to existing

  CHARACTER(LEN=:), POINTER :: newptr
  CHARACTER(LEN=1), POINTER, DIMENSION(:) :: ptr1
  CONTIGUOUS :: ptr1
  CLASS(memory_entry), POINTER :: oldtype, newtype

  CHARACTER(len=:), POINTER :: desc
  CHARACTER(len=5) :: impl
  INTEGER :: i,j,k,l,n1r
  TYPE(c_ptr) :: cptr
  oldtype => find_memory_entry(ptr)
  desc => oldtype%description; IF (PRESENT(description)) desc => description
  impl=implementations(oldtype%implementation); IF (PRESENT(implementation)) impl = implementation
  n1r=((n1-1)/(storage_size(0d0)/character_storage_size)+1)
! attempt to do in-place if possible
  IF (impl.EQ.implementations(implementation_stack) .AND. impl.EQ.implementations(oldtype%implementation) .AND. &
       (oldtype%level+SIZE(oldtype%p).EQ.stack_size+1 &
       .OR. SIZE(oldtype%p).GE.n1r) &
       ) THEN
   IF (oldtype%level+SIZE(oldtype%p).EQ.stack_size+1) THEN ! simply move the top of the stack
    stack_size = oldtype%level-1+n1r
    maximum_stack_size = MAX(maximum_stack_size, stack_size)
   END IF
   oldtype%size = n1r
   oldtype%LBOUND(1)=1
   oldtype%UBOUND(1)=n1
   NULLIFY(oldtype%p)
   oldtype%p => memory_stack(oldtype%level:oldtype%level-1+n1r)
   cptr = c_loc(oldtype%p(1))
   CALL c_f_pointer(cptr,ptr1,[n1])
   CALL get_scalar_pointer(n1,ptr1,ptr)
  ELSE
   newptr => memory_allocate_character(n1,desc,impl)
! find oldtype again in case segments() was reallocated
   oldtype => find_memory_entry(ptr)
   newtype => find_memory_entry(newptr)

   DO i=1,MIN(LEN(ptr),LEN(newptr))
    newptr(i:i)=ptr(i:i)
   END DO

   IF (oldtype%implementation.EQ.implementation_heap) CALL memory_release(ptr)
   ptr => newptr
  END IF
  END SUBROUTINE memory_resize_character

!> \public Allocate a new integer array and fill it with the contents of an existing array
 FUNCTION memory_duplicate_integer_1(source,lb,description,implementation)
  INTEGER, POINTER, DIMENSION(:) :: memory_duplicate_integer_1
  CONTIGUOUS :: memory_duplicate_integer_1
  INTEGER, DIMENSION(:), INTENT(in) :: source !< the source matrix which gives the shape and contents
  INTEGER, INTENT(in), OPTIONAL,DIMENSION(:) :: lb!< the lower bound of the index (default (/1/).
  CHARACTER(len=*), INTENT(in), OPTIONAL :: description !< string describing the matrix
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' (default) or 'HEAP'
  memory_duplicate_integer_1 => memory_allocate_integer(SIZE(source),lb,description,implementation)
  memory_duplicate_integer_1 = source
 END FUNCTION memory_duplicate_integer_1

!> \public Allocate a new double precision array and fill it with the contents of an existing array
 FUNCTION memory_duplicate_double_1(source,lb,description,implementation)
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: memory_duplicate_double_1
  CONTIGUOUS :: memory_duplicate_double_1
  DOUBLE PRECISION, dimension(:), INTENT(in) :: source !< the source matrix which gives the shape and contents
  INTEGER, INTENT(in), OPTIONAL,DIMENSION(:) :: lb!< the lower bound of the index (default (/1/).
  CHARACTER(len=*), INTENT(in), OPTIONAL :: description !< string describing the matrix
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' (default) or 'HEAP'
  memory_duplicate_double_1 => memory_allocate(SIZE(source),lb,description,implementation)
  memory_duplicate_double_1 = source
 END FUNCTION memory_duplicate_double_1

!> \public Allocate a new character array and fill it with the contents of an existing character array
 FUNCTION memory_duplicate_char_array(source,lb,description,implementation)
  CHARACTER(len=1), DIMENSION(:), POINTER :: memory_duplicate_char_array
  CONTIGUOUS :: memory_duplicate_char_array
  CHARACTER(len=*), DIMENSION(:), INTENT(in) :: source !< the source character array
  INTEGER, INTENT(in), OPTIONAL :: lb!< the lower bound of the index (default (/1/).
  CHARACTER(len=*), INTENT(in), OPTIONAL :: description !< string describing the matrix
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' (default) or 'HEAP'
  memory_duplicate_char_array => memory_allocate_character_array(SIZE(source),lb,description,implementation)
  memory_duplicate_char_array = source
 END FUNCTION memory_duplicate_char_array

!> \public Allocate a new character scalar and fill it with the contents of an existing character scalar
 FUNCTION memory_duplicate_character(source,description,implementation)
  CHARACTER(len=:), POINTER :: memory_duplicate_character
  CHARACTER(len=*), INTENT(in) :: source !< the source character string
  CHARACTER(len=*), INTENT(in), OPTIONAL :: description !< string describing the matrix
  CHARACTER(len=*), INTENT(in), OPTIONAL :: implementation !< 'STACK' (default) or 'HEAP'
  memory_duplicate_character => memory_allocate_character(LEN(source),description,implementation)
  memory_duplicate_character = source
 END FUNCTION memory_duplicate_character

!> \public Register the location of the base array (normally common/big/q(1)) used for the old-style stack
!! allocation calls (\ref icorr etc).
 SUBROUTINE memory_register_stack_array(q)
  DOUBLE PRECISION, DIMENSION(1), INTENT(in) :: q
  legacy_stack_offset = (loc(memory_stack(1))-loc(q(1)))/(storage_SIZE(0d0)/storage_size(' '))
  IF (legacy_stack_offset.GT.HUGE(0)) CALL Error(&
       'Traditional stack memory allocation routines cannot be used because offset is too large to store','memory')
  !PRINT *, 'loc(memory_stack(1)) ',loc(memory_stack(1))
  !PRINT *, 'loc(q(1)) ',loc(q(1))
  !PRINT *, 'legacy_stack_offset ',legacy_stack_offset
 END SUBROUTINE memory_register_stack_array

!> \public Create a double precision 1-dimensional array pointer pointing to an existing array
!! at a given stack position.
 FUNCTION memory_pointer_to_q_1(ip, n1, lb) result(ptr)
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: ptr
  CONTIGUOUS :: ptr
  INTEGER, INTENT(in) :: ip !< Position on the `q` stack at which the array begins.
  INTEGER, INTENT(in) :: n1 !< the length of the first dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL, DIMENSION(1) :: lb !< the lower bounds in each dimension (default (/1/).
  INTEGER, DIMENSION(1) :: lbb
  IF (PRESENT(lb)) THEN
   lbb = lb
  ELSE
   lbb = (/1/)
  END IF
  ptr(lbb(1):lbb(1)-1+n1) => memory_stack(ip-legacy_stack_offset:ip-legacy_stack_offset-1+n1)
 END FUNCTION memory_pointer_to_q_1

!> \public Create a double precision 2-dimensional array pointer pointing to an existing array
!! at a given stack position.
 FUNCTION memory_pointer_to_q_2(ip, n1, n2, lb) result(ptr)
  DOUBLE PRECISION, POINTER, DIMENSION(:, :) :: ptr
  CONTIGUOUS :: ptr
  INTEGER, INTENT(in) :: ip !< Position on the `q` stack at which the array begins.
  INTEGER, INTENT(in) :: n1 !< the length of the first dimension of the matrix
  INTEGER, INTENT(in) :: n2 !< the length of the second dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL, DIMENSION(2) :: lb !< the lower bounds in each dimension (default (/1,1/).
  INTEGER, DIMENSION(2) :: lbb
  IF (PRESENT(lb)) THEN
    lbb=lb
  ELSE
    lbb=(/1,1/)
  END IF
  ptr(lbb(1):lbb(1)-1+n1,lbb(2):lbb(2)-1+n2) => memory_stack(ip-legacy_stack_offset:ip-legacy_stack_offset-1+n1*n2)
 END FUNCTION memory_pointer_to_q_2

!> \public Create a double precision 3-dimensional array pointer pointing to an existing array
!! at a given stack position.
 FUNCTION memory_pointer_to_q_3(ip, n1, n2, n3, lb) result(ptr)
  DOUBLE PRECISION, POINTER, DIMENSION(:, :, :) :: ptr
  CONTIGUOUS :: ptr
  INTEGER, INTENT(in) :: ip !< Position on the `q` stack at which the array begins.
  INTEGER, INTENT(in) :: n1 !< the length of the first dimension of the matrix
  INTEGER, INTENT(in) :: n2 !< the length of the second dimension of the matrix
  INTEGER, INTENT(in) :: n3 !< the length of the third dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL, DIMENSION(3) :: lb !< the lower bounds in each dimension (default (/1,1,1/).
  INTEGER, DIMENSION(3) :: lbb
  IF (PRESENT(lb)) THEN
    lbb=lb
  ELSE
    lbb=(/1,1,1/)
  END IF
  ptr(lbb(1):lbb(1)-1+n1,lbb(2):lbb(2)-1+n2,lbb(3):lbb(3)-1+n3) => &
       memory_stack(ip-legacy_stack_offset:ip-legacy_stack_offset-1+n1*n2*n3)
 END FUNCTION memory_pointer_to_q_3

!> \public Create a double precision 4-dimensional array pointer pointing to an existing array
!! at a given stack position.
 FUNCTION memory_pointer_to_q_4(ip, n1, n2, n3, n4, lb) result(ptr)
  DOUBLE PRECISION, POINTER, DIMENSION(:, :, :, :) :: ptr
  CONTIGUOUS :: ptr
  INTEGER, INTENT(in) :: ip !< Position on the `q` stack at which the array begins.
  INTEGER, INTENT(in) :: n1 !< the length of the first dimension of the matrix
  INTEGER, INTENT(in) :: n2 !< the length of the second dimension of the matrix
  INTEGER, INTENT(in) :: n3 !< the length of the third dimension of the matrix
  INTEGER, INTENT(in) :: n4 !< the length of the fourth dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL, DIMENSION(4) :: lb !< the lower bounds in each dimension (default (/1,1,1,1/).
  INTEGER, DIMENSION(4) :: lbb
  IF (PRESENT(lb)) THEN
    lbb=lb
  ELSE
    lbb=(/1,1,1,1/)
  END IF
  ptr(lbb(1):lbb(1)-1+n1,lbb(2):lbb(2)-1+n2,lbb(3):lbb(3)-1+n3,lbb(4):lbb(4)-1+n4) => &
       memory_stack(ip-legacy_stack_offset:ip-legacy_stack_offset-1+n1*n2*n3*n4)
 END FUNCTION memory_pointer_to_q_4

!> \public Create an integer 1-dimensional array pointer pointing to an existing array
!! at a given stack position.
 FUNCTION memory_pointer_to_iq_1(ip, n1, lb) RESULT(ptr)
  INTEGER, POINTER, DIMENSION(:) :: ptr
  CONTIGUOUS :: ptr
  INTEGER, INTENT(in) :: ip !< Position on the `iq` stack at which the array begins.
  INTEGER, INTENT(in) :: n1 !< the length of the first dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL, DIMENSION(1) :: lb !< the lower bounds in each dimension (default (/1/).
  INTEGER, DIMENSION(1) :: lbb
  INTEGER, POINTER, DIMENSION(:) :: x
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: db
  TYPE(c_ptr) :: cptr
  IF (PRESENT(lb)) THEN
   lbb = lb
  ELSE
   lbb = (/1/)
  END IF
  db => memory_pointer_to_q((ip-1)/(storage_size(0d0)/storage_size(0))+1, 1)
  cptr = c_loc(db(1))
  CALL c_f_pointer(cptr,x,SHAPE=[n1])
  ptr(lbb(1):lbb(1)-1+n1) => x
 END FUNCTION memory_pointer_to_iq_1

!> \public Create an integer 2-dimensional array pointer pointing to an existing array
!! at a given stack position.
 FUNCTION memory_pointer_to_iq_2(ip, n1, n2, lb) RESULT(ptr)
  INTEGER, POINTER, DIMENSION(:, :) :: ptr
  CONTIGUOUS :: ptr
  INTEGER, INTENT(in) :: ip !< Position on the `iq` stack at which the array begins.
  INTEGER, INTENT(in) :: n1 !< the length of the first dimension of the matrix
  INTEGER, INTENT(in) :: n2 !< the length of the second dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL, DIMENSION(2) :: lb !< the lower bounds in each dimension (default (/1,1/).
  INTEGER, DIMENSION(2) :: lbb
  INTEGER, POINTER, DIMENSION(:) :: x
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: db
  TYPE(c_ptr) :: cptr
  IF (PRESENT(lb)) THEN
   lbb = lb
  ELSE
   lbb = (/1,1/)
  END IF
  db => memory_pointer_to_q((ip-1)/(storage_size(0d0)/storage_size(0))+1, 1)
  cptr = c_loc(db(1))
  CALL c_f_pointer(cptr,x,SHAPE=[n1*n2])
  ptr(lbb(1):lbb(1)-1+n1,lbb(2):lbb(2)-1+n2) => x
 END FUNCTION memory_pointer_to_iq_2

!> \public Create an integer 3-dimensional array pointer pointing to an existing array
!! at a given stack position.
 FUNCTION memory_pointer_to_iq_3(ip, n1, n2, n3, lb) RESULT(ptr)
  INTEGER, POINTER, DIMENSION(:, :, :) :: ptr
  CONTIGUOUS :: ptr
  INTEGER, INTENT(in) :: ip !< Position on the `iq` stack at which the array begins.
  INTEGER, INTENT(in) :: n1 !< the length of the first dimension of the matrix
  INTEGER, INTENT(in) :: n2 !< the length of the second dimension of the matrix
  INTEGER, INTENT(in) :: n3 !< the length of the third dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL, DIMENSION(3) :: lb !< the lower bounds in each dimension (default (/1,1,1/).
  INTEGER, DIMENSION(3) :: lbb
  INTEGER, POINTER, DIMENSION(:) :: x
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: db
  TYPE(c_ptr) :: cptr
  IF (PRESENT(lb)) THEN
   lbb = lb
  ELSE
   lbb = (/1,1,1/)
  END IF
  db => memory_pointer_to_q((ip-1)/(storage_size(0d0)/storage_size(0))+1, 1)
  cptr = c_loc(db(1))
  CALL c_f_pointer(cptr,x,SHAPE=[n1*n2*n3])
  ptr(lbb(1):lbb(1)-1+n1,lbb(2):lbb(2)-1+n2,lbb(3):lbb(3)-1+n3) => x
 END FUNCTION memory_pointer_to_iq_3

!> \public Create an integer 4-dimensional array pointer pointing to an existing array
!! at a given stack position.
 FUNCTION memory_pointer_to_iq_4(ip, n1, n2, n3, n4, lb) RESULT(ptr)
  INTEGER, POINTER, DIMENSION(:, :, :, :) :: ptr
  CONTIGUOUS :: ptr
  INTEGER, INTENT(in) :: ip !< Position on the `iq` stack at which the array begins.
  INTEGER, INTENT(in) :: n1 !< the length of the first dimension of the matrix
  INTEGER, INTENT(in) :: n2 !< the length of the second dimension of the matrix
  INTEGER, INTENT(in) :: n3 !< the length of the third dimension of the matrix
  INTEGER, INTENT(in) :: n4 !< the length of the fourth dimension of the matrix
  INTEGER, INTENT(in), OPTIONAL, DIMENSION(4) :: lb !< the lower bounds in each dimension (default (/1,1,1,1/).
  INTEGER, DIMENSION(4) :: lbb
  INTEGER, POINTER, DIMENSION(:) :: x
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: db
  TYPE(c_ptr) :: cptr
  IF (PRESENT(lb)) THEN
   lbb = lb
  ELSE
   lbb = (/1,1,1,1/)
  END IF
  db => memory_pointer_to_q((ip-1)/(storage_size(0d0)/storage_size(0))+1, 1)
  cptr = c_loc(db(1))
  CALL c_f_pointer(cptr,x,SHAPE=[n1*n2*n3*n4])
  ptr(lbb(1):lbb(1)-1+n1,lbb(2):lbb(2)-1+n2,lbb(3):lbb(3)-1+n3,lbb(4):lbb(4)-1+n4) => x
 END FUNCTION memory_pointer_to_iq_4


!> \public
 FUNCTION memory_get_stack_position_d1(p) result(level)
  INTEGER :: level
  DOUBLE PRECISION, POINTER, DIMENSION(:), INTENT(in) :: p
  INTEGER :: k
  TYPE(c_ptr) :: cptr1, cptr2
  IF (SIZE(p) == 0) THEN
   level = stack_size+1
   RETURN
  END IF
  cptr1 = c_LOC(p(LBOUND(p,1)))
  DO k=n_segments,1,-1
   cptr2 = c_loc(segments(k)%p(1))
   IF ( ASSOCIATED(p,segments(k)%p) .OR. c_ASSOCIATED(cptr1,cptr2) ) THEN
    IF (segments(k)%implementation.EQ.implementation_stack) THEN
     level = segments(k)%level + legacy_stack_offset
    ELSE IF (segments(k)%implementation.EQ.implementation_heap) THEN
     CALL Error('Attempting to get the stack position of a heap memory pointer', 'memory::memory_get_stack_position')
    ELSE
     CALL Error('Unexpected case reached','memory::memory_release')
    END IF
    RETURN
   END IF
  END DO
  CALL memory_print_status
  CALL Error('Cannot find allocated pointer to release','memory::memory_get_stack_position')
 END FUNCTION memory_get_stack_position_d1

!> \public
 FUNCTION memory_get_stack_position_d2(p) result(level)
  INTEGER :: level
  DOUBLE PRECISION, POINTER, DIMENSION(:, :), INTENT(in) :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: pd
  TYPE(c_ptr) :: cptr
  IF (SIZE(P) == 0) THEN
   level = stack_size+1
   RETURN
  END IF
  cptr = c_loc(p(LBOUND(p,1),LBOUND(p,2)))
  CALL c_f_pointer(cptr,pd,shape=[SIZE(p)])
  level = memory_get_stack_position(pd)
 END FUNCTION memory_get_stack_position_d2

!> \public
 FUNCTION memory_get_stack_position_d3(p) result(level)
  INTEGER :: level
  DOUBLE PRECISION, POINTER, DIMENSION(:, :, :), INTENT(in) :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: pd
  TYPE(c_ptr) :: cptr
  IF (SIZE(P) == 0) THEN
   level = stack_size+1
   RETURN
  END IF
  cptr = c_loc(p(LBOUND(p,1),LBOUND(p,2),LBOUND(p,3)))
  CALL c_f_pointer(cptr,pd,shape=[SIZE(p)])
  level = memory_get_stack_position(pd)
 END FUNCTION memory_get_stack_position_d3

!> \public
 FUNCTION memory_get_stack_position_d4(p) result(level)
  INTEGER :: level
  DOUBLE PRECISION, POINTER, DIMENSION(:, :, :, :), INTENT(in) :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: pd
  TYPE(c_ptr) :: cptr
  IF (SIZE(P) == 0) THEN
   level = stack_size+1
   RETURN
  END IF
  cptr = c_loc(p(LBOUND(p,1),LBOUND(p,2),LBOUND(p,3),LBOUND(p,4)))
  CALL c_f_pointer(cptr,pd,shape=[SIZE(p)])
  level = memory_get_stack_position(pd)
 END FUNCTION memory_get_stack_position_d4

!> \public
 FUNCTION memory_get_stack_position_i1(p) result(level)
  INTEGER :: level
  INTEGER, POINTER, DIMENSION(:), INTENT(in) :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: pd
  TYPE(c_ptr) :: cptr
  IF (SIZE(P) == 0) THEN
   level = stack_size+1
   RETURN
  END IF
  cptr = c_loc(p(LBOUND(p,1)))
  CALL c_f_pointer(cptr,pd,shape=[(SIZE(p)-1)/(storage_size(0d0)/storage_size(0))+1])
  level = (memory_get_stack_position(pd) - 1) * (storage_size(0d0)/storage_size(0)) + 1
 END FUNCTION memory_get_stack_position_i1

!> \public
 FUNCTION memory_get_stack_position_i2(p) result(level)
  INTEGER :: level
  INTEGER, POINTER, DIMENSION(:, :), INTENT(in) :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: pd
  TYPE(c_ptr) :: cptr
  IF (SIZE(P) == 0) THEN
   level = stack_size+1
   RETURN
  END IF
  cptr = c_loc(p(LBOUND(p,1),LBOUND(p,2)))
  CALL c_f_pointer(cptr,pd,shape=[(SIZE(p)-1)/(storage_size(0d0)/storage_size(0))+1])
  level = (memory_get_stack_position(pd) - 1) * (storage_size(0d0)/storage_size(0)) + 1
 END FUNCTION memory_get_stack_position_i2

!> \public
 FUNCTION memory_get_stack_position_i3(p) result(level)
  INTEGER :: level
  INTEGER, POINTER, DIMENSION(:, :, :), INTENT(in) :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: pd
  TYPE(c_ptr) :: cptr
  IF (SIZE(P) == 0) THEN
   level = stack_size+1
   RETURN
  END IF
  cptr = c_loc(p(LBOUND(p,1),LBOUND(p,2),LBOUND(p,3)))
  CALL c_f_pointer(cptr,pd,shape=[(SIZE(p)-1)/(storage_size(0d0)/storage_size(0))+1])
  level = (memory_get_stack_position(pd) - 1) * (storage_size(0d0)/storage_size(0)) + 1
 END FUNCTION memory_get_stack_position_i3

!> \public
 FUNCTION memory_get_stack_position_i4(p) result(level)
  INTEGER :: level
  INTEGER, POINTER, DIMENSION(:, :, :, :), INTENT(in) :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: pd
  TYPE(c_ptr) :: cptr
  IF (SIZE(P) == 0) THEN
   level = stack_size+1
   RETURN
  END IF
  cptr = c_loc(p(LBOUND(p,1),LBOUND(p,2),LBOUND(p,3),LBOUND(p,4)))
  CALL c_f_pointer(cptr,pd,shape=[(SIZE(p)-1)/(storage_size(0d0)/storage_size(0))+1])
  level = (memory_get_stack_position(pd) - 1) * (storage_size(0d0)/storage_size(0)) + 1
 END FUNCTION memory_get_stack_position_i4


 FUNCTION find_memory_entry1(p)
  CLASS(memory_entry), POINTER :: find_memory_entry1
  DOUBLE PRECISION, POINTER, DIMENSION(:), INTENT(in) :: p
  CONTIGUOUS :: p
  INTEGER :: k
  TYPE(c_ptr) :: cptr1, cptr2
#ifdef MOLPRO
  INCLUDE "common/tapes"
#else
  INTEGER, PARAMETER :: iout=6
#endif
  cptr1 = c_LOC(p(LBOUND(p,1)))
  DO k=n_segments,1,-1
   find_memory_entry1 => segments(k)
   cptr2 = c_loc(segments(k)%p(1))
   IF (ASSOCIATED(p,find_memory_entry1%p) .OR. c_ASSOCIATED(cptr1,cptr2)) RETURN
  END DO
  CALL memory_print_status(title='Looking for pointer at address '//IntegerString(INT(loc(p(LBOUND(p,1))))))
  CALL Error ('Cannot locate pointer in heap and stack tables','memory::find_memory_entry')
 END FUNCTION find_memory_entry1

 FUNCTION find_memory_entry2(p)
  CLASS(memory_entry), POINTER :: find_memory_entry2
  DOUBLE PRECISION, POINTER, DIMENSION(:,:), INTENT(in) :: p
  CONTIGUOUS :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: p1
  CONTIGUOUS :: p1
  p1(1:SIZE(p)) => p(:,:)
  find_memory_entry2 => find_memory_entry1(p1)
 END FUNCTION find_memory_entry2

 FUNCTION find_memory_entry3(p)
  CLASS(memory_entry), POINTER :: find_memory_entry3
  DOUBLE PRECISION, POINTER, DIMENSION(:,:,:), INTENT(in) :: p
  CONTIGUOUS :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: p1
  CONTIGUOUS :: p1
  p1(1:SIZE(p)) => p(:,:,:)
  find_memory_entry3 => find_memory_entry1(p1)
 END FUNCTION find_memory_entry3

 FUNCTION find_memory_entry4(p)
  CLASS(memory_entry), POINTER :: find_memory_entry4
  DOUBLE PRECISION, POINTER, DIMENSION(:,:,:,:), INTENT(in) :: p
  CONTIGUOUS :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: p1
  CONTIGUOUS :: p1
  p1(1:SIZE(p)) => p(:,:,:,:)
  find_memory_entry4 => find_memory_entry1(p1)
 END FUNCTION find_memory_entry4


 FUNCTION find_memory_entry_integer1(p)
  CLASS(memory_entry), POINTER :: find_memory_entry_integer1
  INTEGER, POINTER, DIMENSION(:), INTENT(in) :: p
  CONTIGUOUS :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: pp
  CONTIGUOUS :: pp
  INTEGER :: l
  TYPE(c_ptr) :: cptr
  l = (SIZE(p)-1)/(storage_size(0d0)/storage_size(0)) + 1
  cptr = c_loc(p(1))
  CALL c_f_pointer(cptr,pp,shape=[l])
  find_memory_entry_integer1 => find_memory_entry1(pp)
 END FUNCTION find_memory_entry_integer1

 FUNCTION find_memory_entry_integer2(p)
  CLASS(memory_entry), POINTER :: find_memory_entry_integer2
  INTEGER, POINTER, DIMENSION(:,:), INTENT(in) :: p
  CONTIGUOUS :: p
  INTEGER, POINTER, DIMENSION(:) :: p1
  CONTIGUOUS :: p1
  p1(1:SIZE(p)) => p(:,:)
  find_memory_entry_integer2 => find_memory_entry_integer1(p1)
 END FUNCTION find_memory_entry_integer2

 FUNCTION find_memory_entry_integer3(p)
  CLASS(memory_entry), POINTER :: find_memory_entry_integer3
  INTEGER, POINTER, DIMENSION(:,:,:), INTENT(in) :: p
  CONTIGUOUS :: p
  INTEGER, POINTER, DIMENSION(:) :: p1
  CONTIGUOUS :: p1
  p1(1:SIZE(p)) => p(:,:,:)
  find_memory_entry_integer3 => find_memory_entry_integer1(p1)
 END FUNCTION find_memory_entry_integer3

 FUNCTION find_memory_entry_integer4(p)
  CLASS(memory_entry), POINTER :: find_memory_entry_integer4
  INTEGER, POINTER, DIMENSION(:,:,:,:), INTENT(in) :: p
  CONTIGUOUS :: p
  INTEGER, POINTER, DIMENSION(:) :: p1
  CONTIGUOUS :: p1
  p1(1:SIZE(p)) => p(:,:,:,:)
  find_memory_entry_integer4 => find_memory_entry_integer1(p1)
 END FUNCTION find_memory_entry_integer4

 FUNCTION find_memory_entry_character_arr(p)
  CLASS(memory_entry), POINTER :: find_memory_entry_character_arr
  CHARACTER(LEN=1), POINTER, DIMENSION(:), INTENT(in) :: p
  CONTIGUOUS :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: pp
  CONTIGUOUS :: pp
  INTEGER :: l
  TYPE(c_ptr) :: cptr
  l = (SIZE(p)-1)/(storage_size(0d0)/character_storage_size) + 1
  cptr = c_loc(p(LBOUND(p,1)))
  CALL c_f_pointer(cptr,pp,shape=[l])
  find_memory_entry_character_arr => find_memory_entry1(pp)
 END FUNCTION find_memory_entry_character_arr

 FUNCTION find_memory_entry_character(p)
  CLASS(memory_entry), POINTER :: find_memory_entry_character
  CHARACTER(LEN=:), POINTER,  INTENT(in) :: p
  DOUBLE PRECISION, POINTER, DIMENSION(:) :: pp
  CONTIGUOUS :: pp
  INTEGER :: l
  TYPE(c_ptr) :: cptr
  l = (LEN(p)-1)/(storage_SIZE(0d0)/character_storage_size) + 1
  cptr = c_loc(p(1:1))
  CALL c_f_pointer(cptr,pp,shape=[l])
  find_memory_entry_character => find_memory_entry1(pp)
 END FUNCTION find_memory_entry_character

! Legacy Fortran-77 routines
FUNCTION icorr(n)
 IMPLICIT NONE
 INTEGER :: icorr
 INTEGER, INTENT(in) :: n
 DOUBLE PRECISION, POINTER, DIMENSION(:) :: p
 icorr = stack_size + legacy_stack_offset + 1
 IF (n.NE.0) p => memory_allocate(ABS(n), implementation='STACK', maximum=n.GT.0)
END FUNCTION icorr

FUNCTION icori(n)
 IMPLICIT NONE
 INTEGER :: icori
 INTEGER, INTENT(in) :: n
 INTEGER, POINTER, DIMENSION(:) :: p
 icori = (stack_size + legacy_stack_offset)*(storage_SIZE(0d0)/storage_SIZE(0)) + 1
 IF (n.NE.0) p => memory_allocate_integer(ABS(n), implementation='STACK', maximum=n.GT.0)
END FUNCTION icori

SUBROUTINE corlsr(i)
 USE iso_c_binding, ONLY : c_ptr, c_loc, c_f_pointer
 INTEGER, INTENT(in) :: i
 DOUBLE PRECISION, POINTER, DIMENSION(:) :: p
 TYPE(c_ptr) :: cptr
 cptr = c_loc(memory_stack(i-legacy_stack_offset))
 CALL c_f_pointer(cptr,p,[1])
 CALL memory_release(p)
END SUBROUTINE corlsr

SUBROUTINE corlsi(i)
 USE iso_c_binding, ONLY : c_ptr, c_loc, c_f_pointer
 INTEGER, INTENT(in) :: i
 INTEGER, POINTER, DIMENSION(:) :: p
 TYPE(c_ptr) :: cptr
 cptr = c_LOC(memory_stack(1+(i-1)/(storage_SIZE(0d0)/storage_SIZE(0))-legacy_stack_offset))
 CALL c_f_pointer(cptr,p,[1])
 CALL memory_release(p)
END SUBROUTINE corlsi

FUNCTION icorrm()
 INTEGER :: icorrm
 icorrm = memory_remaining('double')
END FUNCTION icorrm

FUNCTION icorim()
 INTEGER :: icorim
 icorim = memory_remaining('integer')
END FUNCTION icorim

FUNCTION icorhw()
 INTEGER :: icorhw
 icorhw = memory_used('stack',.TRUE.)
END FUNCTION icorhw

SUBROUTINE icorhwres()
 CALL memory_reset_maximum_stack
END SUBROUTINE icorhwres


END MODULE memory

#ifndef MOLPRO
! outside module to avoid checking arguments
! provided in case not in compiler
FUNCTION loc(x)
 INTEGER :: loc
 DOUBLE PRECISION, INTENT(in) :: x
 loc = 0
END FUNCTION loc
#endif

! limited and dirty interface for C
! Simple C binding of byte allocator on stack
FUNCTION memory_allocate(n) BIND(C,name="memory_allocate")
 USE memory, ONLY : memory_allocate_character_array
 USE iso_c_binding, ONLY : c_loc, c_ptr, c_size_t
 INTEGER(kind=c_size_t), VALUE, INTENT(in) :: n
 TYPE(c_ptr) :: memory_allocate
 CHARACTER(len=1), DIMENSION(:), POINTER :: fp
 fp => memory_allocate_character_array(INT(n),implementation='STACK')
 memory_allocate = c_LOC(fp(1))
END FUNCTION memory_allocate

SUBROUTINE memory_release(p) BIND(C,name="memory_release")
 USE memory, ONLY : memory_releaseF => memory_release
 USE iso_c_binding, ONLY : c_loc, c_ptr, c_size_t, c_f_pointer
 TYPE(c_ptr), INTENT(in), VALUE :: p
 CHARACTER(len=1), POINTER, DIMENSION(:) :: fp
 CALL c_f_POINTER(p,fp,(/1/))
 CALL memory_releaseF(fp)
END SUBROUTINE memory_release

FUNCTION memory_save() BIND(C,name="memory_save")
 USE memory, ONLY : memory_saveF => memory_save
 USE iso_c_binding, ONLY : c_size_t
 INTEGER(kind=c_size_t) :: memory_save
 memory_save = INT(memory_saveF(),KIND(memory_save_c))
END FUNCTION memory_save

SUBROUTINE memory_release_saved(p) BIND(C,name="memory_release_saved")
 USE memory, ONLY : memory_release
 USE iso_c_binding, ONLY : c_size_t
 INTEGER(kind=c_size_t), INTENT(in), VALUE :: p
 CALL memory_release(INT(p))
END SUBROUTINE memory_release_saved

FUNCTION memory_used(maximum) BIND(C,name='memory_used')
 USE memory, ONLY : memory_usedF => memory_used
 USE iso_c_binding, ONLY : c_size_t, c_int, c_char
 INTEGER(c_int), INTENT(in), VALUE :: maximum !< if true, report the maximum memory allocated to date
 INTEGER(c_size_t) :: memory_used
 memory_used = (storage_SIZE(0d0)/storage_SIZE(' '))*memory_usedF('STACK',maximum=INT(maximum).NE.0)
END FUNCTION memory_used

FUNCTION memory_remaining() BIND(C,name='memory_remaining')
 USE memory, ONLY : memory_remainingF => memory_remaining
 USE iso_c_binding, ONLY : c_size_t, c_int, c_char
 INTEGER(c_size_t) :: memory_remaining
 memory_remaining = (storage_SIZE(0d0)/storage_SIZE(' '))*memory_remainingF('double')
END FUNCTION memory_remaining

SUBROUTINE memory_reset_maximum_stack(level) BIND(C, name='memory_reset_maximum_stack')
 USE iso_c_binding, ONLY : c_int64_t
 USE memory, ONLY : memory_reset_maximum_stackF => memory_reset_maximum_stack
 INTEGER(kind=c_int64_t), INTENT(in), VALUE :: level
 IF (level.LT.0) THEN
  CALL memory_reset_maximum_stackF
 ELSE
  CALL memory_reset_maximum_stackF(storage_SIZE(' ')*INT(level)/storage_SIZE(0d0))
 END IF
END SUBROUTINE memory_reset_maximum_stack

! unmodularised Fortran-77 wrappers
FUNCTION icorr(n)
 USE memory, only : m_icorr => icorr
 icorr = m_icorr(n)
END FUNCTION icorr

FUNCTION icori(n)
 USE memory, only : m_icori => icori
 icori = m_icori(n)
END FUNCTION icori

SUBROUTINE corlsr(i)
 USE memory, only : m_corlsr => corlsr
 call m_corlsr(i)
END SUBROUTINE corlsr

SUBROUTINE corlsi(i)
 USE memory, only : m_corlsi => corlsi
 call m_corlsi(i)
END SUBROUTINE corlsi

FUNCTION icorrm()
 USE memory, only : m_icorrm => icorrm
 icorrm = m_icorrm()
END FUNCTION icorrm

FUNCTION icorim()
 USE memory, only : m_icorim => icorim
 icorim = m_icorim()
END FUNCTION icorim

FUNCTION icorhw()
 USE memory, only : m_icorhw => icorhw
 icorhw = m_icorhw()
END FUNCTION icorhw

SUBROUTINE icorhwres()
 USE memory, only : m_icorhwres => icorhwres
 CALL m_icorhwres
END SUBROUTINE icorhwres

! outside module to avoid false positives from private module elements
SUBROUTINE memory_module_test(printlevel)
 USE memory
 IMPLICIT NONE
 INTEGER, INTENT(in) :: printlevel
 INTEGER :: i,j,k,l,ibase,iimpl,i2p
 INTEGER :: base,base2,memsave
 DOUBLE PRECISION, DIMENSION(:), POINTER :: zerotest,pd,pd2
 CONTIGUOUS :: zerotest,pd,pd2
 DOUBLE PRECISION, DIMENSION(:,:,:,:), POINTER :: r4
 CONTIGUOUS :: r4
 INTEGER, DIMENSION(:), POINTER :: pi
 INTEGER, DIMENSION(:,:), POINTER :: i2
 INTEGER, DIMENSION(:,:,:,:), POINTER :: i4
 CONTIGUOUS :: i4
 CHARACTER(5), DIMENSION(2) :: impls=(/'STACK','HEAP '/)
 CHARACTER(5) :: impl, impl_default
 INTEGER, DIMENSION(4) :: bounds=[4,3,2,1]
 CHARACTER(1), DIMENSION(:), POINTER :: cp
 CONTIGUOUS :: cp
 CHARACTER(len=:), pointer :: c
#ifdef MOLPRO
 INTEGER, EXTERNAL :: molpro_loc
 INCLUDE "common/big"
#else
 DOUBLE PRECISION, DIMENSION(1) :: q
 COMMON /big/ q
 INTEGER iq(1)
 EQUIVALENCE (q(1),iq(1))
#endif
 !WRITE (6,*) 'Test memory module, print level =',printlevel; CALL flush6
 !CALL memory_clean
 !PRINT *, 'memory_used() ',memory_used()
 !PRINT *, 'memory_remaining() ',memory_remaining()
 IF (printlevel.GT.1) CALL memory_print_status
 impl_default = memory_default_implementation()
 DO iimpl=1,2
  impl=memory_default_implementation(impls(iimpl))
  WRITE (6,*) ' Test implementation ',TRIM(impls(iimpl))
  ibase = memory_remaining()
  base = memory_save()
  memsave=memory_used()
  base2 = memory_save()
  zerotest => memory_allocate(0)
  IF (printlevel.GT.1) WRITE (6,*) 'zero-length pointer created, size=',SIZE(zerotest)
  IF (printlevel.GT.1) CALL memory_print_status(title='after creating one null')
  zerotest => memory_allocate(0)
  IF (printlevel.GT.1) CALL memory_print_status(title='after creating two nulls')
  CALL memory_release(zerotest)
  CALL memory_release(base2)
  IF (printlevel.GT.1) CALL memory_print_status(title='after releasing')
  IF (memory_used() .NE.memsave) CALL Error('Unexpected failure in memory_save/memory_release','memory')

! simple tests
  i2=>memory_allocate_integer(1000,1000)
  IF (printlevel.GT.9) CALL memory_print_status(title='after allocating two-dimensional integer array')
  !WRITE (6,*) ASSOCIATED(i2)
  call memory_release(i2)
  IF (printlevel.GT.9) CALL memory_print_status(title='after releasing two-dimensional integer array')
  pi=>memory_allocate_integer(1000,description='pi')
  if (printlevel.gt.9) call memory_print_status(title='after create pi')
  pd=>memory_allocate_double(1000,description='pd')
  if (printlevel.gt.9) call memory_print_status(title='after create pd')
  call memory_release(pd)
  if (printlevel.gt.9) call memory_print_status(title='after release pd')
  call memory_release(pi)
  if (printlevel.gt.9) call memory_print_status(title='after release pi')

  IF (printlevel.GT.1) CALL memory_print_status(title='after releasing one of the nulls')
  i4 => memory_allocate_integer(bounds(1),bounds(2),bounds(3),bounds(4),description='4-d integer array')
  IF (printlevel.GT.1) CALL memory_print_status(title='after creating i4')
  r4 => memory_allocate(bounds(1),bounds(2),bounds(3),bounds(4),description='4-d double array')
  IF (printlevel.GT.1) CALL memory_print_status(title='after creating d4')
  !WRITE (6,'(A,4(1X,i1,'':'',i1))') 'i4 bounds ',(LBOUND(i4,i),UBOUND(i4,i),i=1,4)
  !WRITE (6,'(A,4(1X,i1,'':'',i1))') 'r4 bounds ',(LBOUND(r4,i),UBOUND(r4,i),i=1,4)
  !call flush6
  DO l=LBOUND(i4,4),UBOUND(i4,4)
   DO k=LBOUND(i4,3),UBOUND(i4,3)
    DO j=LBOUND(i4,2),UBOUND(i4,2)
     DO i=LBOUND(i4,1),UBOUND(i4,1)
      i4(i,j,k,l) = val(i,j,k,l)
      r4(i,j,k,l) = val(i,j,k,l)
     END DO
    END DO
   END DO
  END DO
  IF (printlevel.GT.9) CALL memory_print_status(title='before resizing i4')
  CALL memory_resize(i4,5,4,3,2,lb=(/0,0,0,0/))
  IF (printlevel.GT.9) CALL memory_print_status(title='after resizing i4')
  CALL memory_resize(r4,5,4,3,3,lb=(/2,2,2,-1/))
  IF (printlevel.GT.9) CALL memory_print_status(title='after resizing r4')
  !WRITE (6,'(A,4(1X,i1,'':'',i1))') 'i4 bounds ',(LBOUND(i4,i),UBOUND(i4,i),i=1,4)
  !WRITE (6,'(A,4(1X,i1,'':'',i1))') 'r4 bounds ',(LBOUND(r4,i),UBOUND(r4,i),i=1,4)
  DO l=MAX(1,LBOUND(r4,4)),MIN(bounds(4),UBOUND(r4,4))
   DO k=MAX(1,LBOUND(r4,3)),MIN(bounds(3),UBOUND(r4,3))
    DO j=MAX(1,LBOUND(r4,2)),MIN(bounds(2),UBOUND(r4,2))
     DO i=MAX(1,LBOUND(r4,1)),MIN(bounds(1),UBOUND(r4,1))
      if (printlevel.gt.10) WRITE (6,*) 'i,j,k,l,r4,val ',i,j,k,l,r4(i,j,k,l),val(i,j,k,l)
      IF (r4(i,j,k,l) .NE. val(i,j,k,l)) CALL Error('Wrong double values after resize','memory_module_test')
     END DO
    END DO
   END DO
   END DO
  DO l=MAX(1,LBOUND(i4,4)),MIN(bounds(4),UBOUND(i4,4))
   DO k=MAX(1,LBOUND(i4,3)),MIN(bounds(3),UBOUND(i4,3))
    DO j=MAX(1,LBOUND(i4,2)),MIN(bounds(2),UBOUND(i4,2))
     DO i=MAX(1,LBOUND(i4,1)),MIN(bounds(1),UBOUND(i4,1))
      if (printlevel.gt.10) WRITE (6,*) 'i,j,k,l,i4,val ',i,j,k,l,i4(i,j,k,l),val(i,j,k,l)
      IF (i4(i,j,k,l) .NE. val(i,j,k,l)) CALL Error('Wrong integer values after resize','memory_module_test')
     END DO
    END DO
   END DO
  END DO
  IF (printlevel.GT.1) CALL memory_print_status(title='before releasing r4')
  CALL memory_release(r4)
  IF (printlevel.GT.1) CALL memory_print_status(title='after releasing r4')
  CALL memory_release(i4)
  IF (printlevel.GT.1) CALL memory_print_status(title='after releasing i4')

  cp => memory_allocate_character_array(9)
  IF (printlevel.GT.1) CALL memory_print_status(title='after allocating cp')
  cp=(/'A','b','3','d','e','f','g','h','i'/)
  IF (SIZE(cp).NE.9) CALL Error('Incorrect allocated character length: ','memory')
  !PRINT *,'cp1 ',LBOUND(cp,1),UBOUND(cp,1),cp,'9=',cp(9)
  CALL memory_resize(cp,24,lb=1)
  IF (printlevel.GT.1) CALL memory_print_status(title='after resizing cp')
  !PRINT *,'cp2 ',LBOUND(cp,1),UBOUND(cp,1),cp,'9=',cp(9)
  IF (cp(9).NE.'i') CALL Error('Something wrong with allocate_character','memory')
  DO i=LBOUND(cp,1),UBOUND(cp,1); cp(i)=CHAR(IACHAR('A')+i-1) ; ENDDO
   !PRINT *,'cp3 ',LBOUND(cp,1),UBOUND(cp,1),cp
  IF (printlevel.GT.2) CALL memory_print_status(title='before releasing cp')
  call memory_release(cp)
  IF (printlevel.GT.1) CALL memory_print_status(title='after releasing cp')

  c => memory_allocate_character(9,description='character(:)')
  IF (printlevel.GT.1) CALL memory_print_status(title='after allocating c')
  IF (LEN(c).NE.9) CALL Error('Incorrect allocated character length','memory')
  c='Ab3defghi'
  CALL memory_resize(c,24)
  IF (LEN(c).NE.24) CALL Error('Incorrect resized character length','memory')
  IF (printlevel.GT.1) CALL memory_print_status(title='after resizing c')
  IF (c(9:9).NE.'i') CALL Error('Something wrong with resize_character','memory')
  DO i=1,len(c); c(i:i)=CHAR(IACHAR('A')+i-1) ; ENDDO
   !PRINT *,'c3 ',c
  IF (printlevel.GT.2) CALL memory_print_status(title='before releasing c')
  call memory_release(c)
  IF (printlevel.GT.1) CALL memory_print_status(title='after releasing c')

  pd => memory_allocate(40,implementation='HEAP') ! Portland does some copying that goes wrong
!in the CALL to memory_duplicate IF a stack array is used as source
  pd=1d0
  pd2 => memory_duplicate(pd)
  IF (printlevel.GT.1) CALL memory_print_status(title='after duplicate')
  IF (.NOT. ALL(pd.EQ.pd2)) CALL Error('Faulty memory_duplicate','memory')
  call memory_release(pd2)
  call memory_release(pd)


!  PRINT *,'base=',base
  CALL memory_release(base)
  IF (printlevel.GT.1) CALL memory_print_status(title='after releasing base')
  if (printlevel.GT.1) WRITE (6,*) 'ibase, memory_remaining() ',ibase, memory_remaining()
  IF (ibase.NE.memory_remaining()) CALL Error ('Memory not released properly','memory_module_test')
 END DO

 ibase = icori(0)
 i2p = icori(100)
 call corlsi(ibase)

 WRITE (6,*) ' Test implementation of legacy stack routines'
 i2p = icori(100)
 CALL corlsi(i2p+50)
 IF (icori(0)-i2p.NE.50) CALL Error('stack adjustment failure','memory')
 CALL corlsi(i2p)

 WRITE (6,*) ' Test interoperability with the q/iq stack'
 i2p=icori(100**2)
 k=0
 DO j=1,100
  DO i=0,99
   iq(i2p+k)=100*i+j
   k=k+1
  END DO
 END DO
 i2 => memory_pointer_to_iq(i2p,100,100,[0,1])
 DO j=1,100
  DO i=0,99
   IF (i2(i,j) /= 100*i+j) CALL error('Incorrect array returned from memory_pointer_to_iq', 'memory')
  END DO
 END DO
 CALL corlsi(i2p)
 i2 => memory_allocate_integer(100,100,[0,1],implementation="STACK")
 DO j=1,100
  DO i=0,99
   i2(i,j)=100*i+j
  END DO
 END DO
 i2p = memory_get_stack_position(i2)
 k=0
 DO j=1,100
  DO i=0,99
   if (iq(i2p+k) /= 100*i+j) CALL error('Incorrect position returned from memory_get_stack_position', 'memory')
   k=k+1
  END DO
 END DO
 CALL memory_release(i2)


 impl=memory_default_implementation(impl_default)
 !WRITE (6,*) 'End of test memory module'
CONTAINS
!> \private
 FUNCTION val(i,j,k,l)
  INTEGER :: val
  INTEGER, INTENT(in) :: i,j,k,l
  val=i+j*10+k*100+l*1000
 END FUNCTION val
END SUBROUTINE memory_module_test

#ifndef MOLPRO
#ifndef NOMAIN
PROGRAM main
 USE memory
 CHARACTER(len=:), POINTER :: string
 DOUBLE PRECISION, DIMENSION(1) :: q
 COMMON /big/ q
 CALL memory_initialize(1001000)
 CALL memory_register_stack_array(q)
 CALL memory_module_test(0)
 string => memory_duplicate('abcde',implementation='HEAP')
 call memory_print_status
 CALL memory_close
END PROGRAM main
SUBROUTINE Error(msg,place)
 CHARACTER(*), INTENT(in) :: msg, place
 PRINT *,'Error: ',msg,':',place
 STOP
END SUBROUTINE Error
SUBROUTINE Warning(msg,place)
 CHARACTER(*), INTENT(in) :: msg, place
 PRINT *,'Warning: ',msg,':',place
 STOP
END SUBROUTINE Warning
#endif
#endif
