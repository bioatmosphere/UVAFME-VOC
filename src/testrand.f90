program testrand
     real randnum
     real urand

     do i=1,100000
          randnum=urand()
          print *, randnum
     enddo

end program testrand


function urand(lb,ub,seed)
     real                                 :: urand
     real, optional                       :: lb,ub
     real, optional                       :: seed

     real                                 :: rnd
     real                                 :: lower,upper
     integer                              :: isize, idate(8)
     integer, dimension(:), allocatable   :: iseed
     logical                              :: first=.true.

     ! Returns a uniformly-distributed random number in the range lb to ub.
     ! Currently sets the seed only on the first call.

     if ( first ) then
          call random_seed(size=isize)
          if (.not. allocated(iseed)) allocate(iseed(isize))

          if ( .not. present(seed) ) then
               call date_and_time(values=idate)
               call random_seed(get=iseed)
               iseed = iseed * (idate(8))      ! idate(8) contains millisecond
          else
               iseed=seed
          endif

          call random_seed(put=iseed)

          first=.false.
     endif

     if ( present(lb) ) then
          lower=lb
     else
          lower=0.0
     endif

     if ( present(ub) ) then
          upper = ub
     else
          upper = 1.0
     endif
     call random_number(rnd)

     urand = lower+(upper-lower)*rnd

     return
end function urand



function climrand(lb,ub)
     real                             :: urand
     real, optional                   :: lb,ub

     ! Used when it is desired to have a separate RNG for climate data
     ! Keep It Simple Stupid by George Marsaglia
     ! The  KISS (Keep It Simple Stupid) random number generator. Combines:
     ! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
     ! (2) A 3-shift shift-register generator, period 2^32-1,
     ! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
     !  Overall period>2^123;  Default seeds x,y,z,w.
     !
     integer KISS
     integer x,y,z,w,c,t
     ! Inline function
     m(y,k)=ieor(y,ishft(y,k))

     if ( present(lb) ) then
          lower=lb
     else
          lower=0.0
     endif

     if ( present(ub) ) then
          upper = ub
     else
          upper = 1.0
     endif

     x=x+545925293
     y=m(m(m(y,13),-17),5)
     t=z+w+c; z=w; c=ishft(t,-31); w=iand(t,2147483647)
     KISS=x+y+w

     climrand=lower+(upper-lower)*real(KISS)*4.656613e-10

     return
end function climrand
