program vm_trace

   implicit none

   real :: x1,x2		!  horizontal boundaries of model
   real :: z1,z2		!  vertical boundaries of model
   integer :: nx,nz		!  number of grid points in x and z direction
   integer :: nv		!  number of grid points in model (nx*nz)
   integer :: nm		!  number of grid points in all velocity slices (nr+1)*nv
   integer :: nr		!  number of reflecting boundaries
   integer :: nb		!  number of boundary points (nx*nr)
   real :: dx,dz		!  grid spacings
   real, allocatable :: zrf(:)	!  reflector depths
   integer, allocatable :: ird(:)   !  ID of reflector

   real, allocatable :: vel(:)	!  grids of velocities for all layers
   real, allocatable :: vir(:)	!  single grid of velocities

   real :: cmax			! max slowness (is this needed as a global variable?)
   real :: tinf			! infinite traveltime
   real :: star			! max distance between nodes in forward star

   real :: xp1,zp1		! lower left corner graph grid
   real :: xp2,zp2		! upper right corner graph grid


   integer :: nxp,nzp       ! graph grid dimensions
   real :: drp			! graph grid node spacing
   real :: drr          ! graph refl node spacing
   integer :: iseed		! seed for graph rgid randomization


   integer :: nn		! total number of graph nodes
   integer :: np		! number of nodes without source and refl.
   integer :: nxr		! number of interface nodes
   integer :: nps		! number of nodes for shots

   integer :: ins		! instrument ID
   real :: xin,zin		! location of instrument
   integer, allocatable :: isk(:)	! shot number with pick

   integer, allocatable :: fas(:)	! phase ID
   real, allocatable :: ttp(:)		! picked traveltime
   real, allocatable :: etp(:)		! pick uncertainty
   integer, allocatable :: ish(:)	! shot number with coord
   integer, allocatable :: ind(:)   ! shot file index for pick
   integer, allocatable :: igr(:)   ! node number for pick
   real, allocatable :: xsh(:)		! x-coord. shot
   real, allocatable :: zsh(:)		! x-coord. shot

   real, allocatable :: xgr(:)
   real, allocatable :: zgr(:)

   integer :: mrf		! index loop over model layer boundaries
   integer :: irf		! index loop over layers
   integer :: jrf		! index loop over phases (reflections, refractions)
   integer :: nsh		! number of shots
   integer :: npx		! number of picks

   integer :: nfr,nfl		! number of refractions and reflections

   integer :: ngr			! number of active grid nodes
   integer :: nhp			! size of heap
   integer, allocatable :: hp(:)	! pointer from grid to heap
   integer, allocatable :: qp(:)	! pointer from heap to grid
   integer, allocatable :: pd(:)	! downgoing graph
   integer, allocatable :: pu(:)	! upgoing graph
   real, allocatable :: tt(:)		! traveltimes

   integer, allocatable :: len(:)	! the length of a ray path
   real, allocatable :: tca(:)		! calc tt for ray
   real, allocatable :: xry(:),zry(:)	! ray path storage
   integer :: nrt			! total number ray points

   real, parameter :: fpr = 0.33	! spacing of ray sample points in drp
   real, parameter :: frr = 0.1		! ratio of drr and drp in graph method
   real, parameter :: fsr = 3.5		! ratio of star and drp in graph method

   real, allocatable :: sxr(:),szr(:)	! resampled ray

   real, allocatable :: grx(:),grz(:)	! gradient of ray
   real, allocatable :: hrx(:),hrz(:)	! conjugate gradient of ray
   real, allocatable :: frx(:),frz(:)	! copy of gradient of ray

   integer :: nry	! allocated space for one ray
   integer :: npa	! allocated space for rays



!		******************************************



   call read_vm				! read binary velocity model

   call nwk_prm				! read graph parameters

   call read_picks			! get all picks with specified phases


   call alloc_rays			! reserve space for ray paths

   do mrf = 0,nr			! loop over model layers

      call count_picks			! count refraction and reflection picks

      jrf = 0



      if (nps > 0) then         ! some shot points sampled layer mrf

         call graph_size		! calculate size of the graph

         call graph_nodes		! allocate and assign the graph nodes in x,z

     if (mrf==0) then
        print*,'Shortest path method for first arrivals '
     else
            print*,'Shortest path method to layer ',mrf
     end if

     call init_source(mrf)		! initiate the graph with a point source

         call fill_vel(mrf)		! select velocities for this calculation

         call loop_star(pd)		! propagate traveltimes down

     if (mrf==0) then		! mrf==0 represents first-arriving phase

        call build_rays		! direct arrivals use just one graph

         else

        do jrf = -mrf,mrf,2*mrf	! reflections / refractions in pairs

           if ((jrf<0 .and. nfl>0) .or. (jrf>0 .and. nfr>0)) then

                  if (jrf < 0) then	! for reflections we must active the layer
             irf = mrf      ! beneath the boundary mrf
          else
             irf = mrf+1
          end if

          print*,'Shortest path method back to surface for phase ',jrf

                  call fill_vel(irf)	! choose seismic velocities

                  call init_boundary(irf) ! initialize graph for upward calculation

                  call loop_star(pu)	! propagate traveltimes upward to shots

          print*

          call build_rays	! extract ray paths and traveltimes from graph

           end if

        end do

         end if

         call close_graph(mrf)	! free up memory for this graph

      end if        ! There are picks for this layer

   end do		! loop over layers

   call write_rays	! write all phases to disk

contains


   subroutine read_vm

      character(60) :: vmfile
      integer :: j

      integer :: m,ix,iz,i,ii
      real :: zi,v1,v2,v3,v4
      real :: x

      read(5,'(a)') vmfile

      open(unit=13,file=vmfile,access='stream')
      read(13) nx,nz,nr
      read(13) x1,x2
      read(13) z1,z2

      nv = nx*nz
      nm = nv*(nr+1)
      nb = nx*nr
      dx = (x2-x1)/(nx-1)
      dz = (z2-z1)/(nz-1)

      allocate(zrf(nb),ird(nb),vel(nm))

      read(13) (zrf(j),j=1,nb)
      read(13) (ird(j),j=1,nb)
      read(13) (vel(j),j=1,nm)
      close(13)

      allocate(vir(nv))		! fill in later

   end subroutine read_vm





   subroutine nwk_prm

      integer, parameter :: larg = 987654321
      integer :: isec

      read(5,*) cmax
      tinf = cmax*sqrt((x2-x1)*(x2-x1)+(z2-z1)*(z2-z1))
      read(5,*) drp     ! average grid graph node spacing

      drr = frr*drp     ! drr should be about 10x smaller than drp
      star = fsr*drp	! star should be about 3x larger than drp

      read(5,*) isec	! seed for random numbers
      iseed = larg - isec

      call got_rand	! test rand function.

   end subroutine nwk_prm


   subroutine my_rand(val)

     ! rand must be a user or compiler
     ! supplied random number generator
     ! that takes an integer as seed, and
     ! returns a value equally distributed
     ! between 0 and 1.


      real :: val
      real :: dice
      integer, parameter :: mul = 672182
      integer, parameter :: larg = 987654321
      integer :: irf,icnt
      real :: rand

      dice = rand(iseed)
      if (dice <= 0.0 .or.dice >= 1.0) then
         if (iseed == 0) iseed = larg
         irf = abs(iseed)
     icnt = 1
     do while (icnt < irf)
        icnt = 10*icnt
     end do
     dice = float(irf)/float(icnt)
      dice = (dice - 0.1)/0.9
      end if

      val = dice - 0.5

      iseed = iseed + nint(mul*val)
      if (iseed > larg) then
         iseed = iseed - larg
      elseif (iseed < -larg) then
        iseed = iseed + larg
      end if

   end subroutine my_rand






   subroutine got_rand

!	Test the rand function before using it to choose graph nodes.
!	Fortran function rand should give a uniformly distributed random
!	number between 0 and 1. If not, alert and stop the program.

      integer, parameter :: nbin = 10
      integer, parameter :: ndice = 5000
      integer :: ibn(nbin)
      real :: dbin
      real :: val,dice
      real :: vq1
      integer :: id,ib
      integer :: nmin,nmax
      real :: rand

      dbin = 1.0/nbin
      ibn = 0

      do id=1,ndice

     call my_rand(val)
     dice = 0.5 + val

     ib = int(dice/dbin)+1
     if (ib<1) ib = 1
     if (ib>nbin) ib = nbin
     ibn(ib) = ibn(ib) + 1

      end do

      nmax=0
      nmin=ndice
      do ib=1,nbin
         id = ibn(ib)
     if (id > nmax) nmax = id
     if (id < nmin) nmin = id
      end do

      if (2*nmin < nmax) then
         print*,'Function rand is a bad random number generator.'
     print*,'iseed ',iseed
     print*,'Largest 10% bin cnt = ',nmax,nmax/float(ndice)
     print*,'Smallest 10% bin cnt = ',nmin,nmin/float(ndice)
     print*,'Provide a better function rand.'
     print*

      end if

   end subroutine got_rand




   subroutine read_picks

      character(60) :: pick_file,shot_file
      real :: range
      integer :: instr
      integer :: j,k,m,i,ji,ii,jj,imin,jmin
      integer :: fi,fj,fia,fja,fim,ikm
      integer :: jsw,iswap
      real :: swap
      integer :: nshot
      integer :: npk
      integer :: is,id
      real :: xs,zs
      integer :: itop,ibot		! top and bottom raytracing
      integer :: ip1,ip2,iq1,iq2
      integer :: ik,fk
      real:: rk,tk,ek
      integer :: ipol,jpol

      integer, allocatable :: irn(:),ifx(:)
      integer, allocatable :: ift(:),ist(:)
      real, allocatable :: tpt(:),ept(:)

      integer :: fai,faj,famin

      read(5,*) itop,ibot		! shall and deep phase considered

      if (itop < 0) then		! in this calculation
         ip1 = itop
     iq1 = -itop
      else
         ip1 = -itop-1
     iq1 = itop
      end if
      if (ibot < 0) then
         ip2 = ibot
     iq2 = -ibot-1
      else
         ip2 = -ibot
     iq2 = ibot
      endif

      read(5,*) ins			! instrument # and location
      read(5,*) xin,zin

      if (xin < x1 .or. xin > x2 .or. zin < z1 .or. zin > z2) then
         print*,'instrument out of model bounds ',ins,xin,zin
     stop
      end if

      read(5,'(a)') pick_file		! read picks first to count them
      read(5,*) npk




      read(5,'(a)') shot_file	! read shot numbers and locations
      read(5,*) nshot


      open(22,file=shot_file)	! first scan file to count shots
      nsh = 0
      do k=1,nshot
         read(22,*) is,xs,zs
     if (xs >= x1 .and. xs <= x2 .and. zs >= z1 .and. zs <= z2) nsh=nsh+1
      end do
      close(22)
      open(22,file=shot_file)



      allocate(ish(nsh),xsh(nsh),zsh(nsh))  ! read shots again to store
      m = 0
      do k=1,nshot
         read(22,*) is,xs,zs
     if (xs >= x1 .and. xs <= x2 .and. zs >= z1 .and. zs <= z2) then

        m=m+1
        ish(m) = is
        xsh(m) = xs
        zsh(m) = zs

     end if
      end do
      close(22)




      open(21,file=pick_file)
      npx = 0
      do m=1,npk
         read(21,*) instr,ik,fk,rk,tk,ek

     id = -1
     do k=1,nsh
        if (ish(k) == ik) id = 1
     end do

     if (ins == instr .and. id > 0) then
        if (fk==0.or.(ip2<=fk.and.fk<=ip1).or.(iq1<=fk.and.fk<=iq2)) npx=npx+1
     end if
      end do
      close(21)



      open(21,file=pick_file)   ! after memory allocation, read picks to store

      allocate(ist(npx),ift(npx),tpt(npx),ept(npx),irn(npx),ifx(npx))
      j=0
      do m=1,npk
         read(21,*) instr,ik,fk,rk,tk,ek


     id = -1
     do k=1,nsh
        if (ish(k) == ik) id = 1
     end do

     if (ins == instr .and. id > 0) then
        if (fk==0 .or.(ip2<=fk.and.fk<=ip1).or.(iq1<=fk.and.fk<=iq2)) then
           j=j+1
           ist(j) = ik
           ift(j) = fk
           tpt(j) = tk
           ept(j) = ek
           irn(j) = j
        end if
     end if
      end do
      close(21)

      do i=1,npx
         ifx(i) = 1			! raytrace order?
      end do

      do j = 1,npx          ! loop over picks for this instrument (instr) and valid shots
                                !      sort by layer, phase, and shot number
         jmin = 999999999		! layer
     fim  = 999999999		! phase
     ikm  = 999999999       ! shot

     famin = 9999999

     jpol = 1

         do i =  1,npx
        fi = ift(i)
        fai =2*abs(fi)
        if (fi > 0) fai = fai + 1

        if (ifx(i) == 1) then

           if (fai < famin) then
              famin = fai
              imin = i
              fim = fi
              ikm = ist(i)
           elseif (famin == fai .and. ist(i) < ikm) then
              famin = fai
              imin = i
              fim = fi
              ikm = ist(i)
           end if
        end if
     end do

     irn(j) = imin
     ifx(imin) = 0

      end do				!	loop over j

      allocate(fas(npx),isk(npx),ttp(npx),etp(npx))

      do j=1,npx

     i = irn(j)
     isk(j) = ist(i)
     fas(j) = ift(i)
     ttp(j) = tpt(i)
     etp(j) = ept(i)

    ! write(6,'(3i5,2f9.3)') j,isk(j),fas(j),ttp(j),etp(j)

      end do



      deallocate(ist,ift,irn,tpt,ept)


   end subroutine read_picks






   subroutine alloc_rays

      real :: xm,zm
      xm = x2-x1
      zm = z2-z1

      nry =  nint((xm+2*zm)/(fpr*drp))
      npa = npx*nry
      allocate(xry(npa),zry(npa))
      allocate(len(npx),tca(npx))
      allocate(ind(nsh))   ! This array holds a temporary node # for each shot
      allocate(igr(npx))   ! This array holds the graph node # for each pick

      allocate(sxr(nry),szr(nry))	! resampled ray

      allocate(frx(nry),frz(nry))	! copy ray traveltime gradient
      allocate(grx(nry),grz(nry))	! ray traveltime gradient
      allocate(hrx(nry),hrz(nry))	! ray traveltime conjugate gradient

      nrt = 0		! init ray point counter

   end subroutine alloc_rays





   subroutine count_picks

!		Count picks of arrivals that engage model boundary mrf
      integer :: j,m
      integer :: jf,id

      nfr = 0
      nfl = 0
      do j=1,npx
         jf = fas(j)
     if (jf == -mrf) nfl = nfl+1	! number of reflections for graph
     if (jf ==  mrf) nfr = nfr+1	! number of refractions for graph
      end do

      nps = 0
      do m = 1,nsh
         id = -1	! initialize
     ind(m) = -1
     do j = 1,npx
        if (isk(j) == ish(m)) then
           jf = fas(j)
               if (abs(jf) == mrf) id = 1
        end if
     end do
     if (id == 1) then
        nps=nps+1
        ind(m) = nps	! ind numbers the shot locations
     end if         ! in the graph from 1 to nps
      end do

   end subroutine count_picks



   subroutine graph_size
!			Determine the graph size that will be used
!			in the following calculation.

      real :: xmin,xmax
      real :: zmin,zmax

      integer :: im
      integer :: ix
      integer :: j,m,id
      real :: x
      real :: zm


      xmin = x2		! initialize
      xmax = x1
      xp1 = xmin
      xp2 = xmax

      if (xin - star < xmin) xmin = xin - star
      if (xin + star > xmax) xmax = xin + star

      do m=1,nsh
         id = ind(m)
     if (id > 0) then
        if (xsh(m)-star < xmin) xmin = xsh(m)-star
        if (xsh(m)+star > xmax) xmax = xsh(m)+star
     end if
      end do

      zmin = z1
      if (mrf == 0 .or. mrf == nr) then
     zmax = z2
      else
     zmax = zin+star
     do ix=1,nx
        x = x1 + (ix-1)*dx
        if (x >= xmin .and. x <= xmax) then
           im = (mrf  )*nx + ix
           zm = zrf(im) + drp
           if (zm > zmax) zmax = zm
        end if
     end do
      end if

      if (xmin<x1) xmin = x1
      if (xmax>x2) xmax = x2
      if (zmin<z1) zmin = z1
      if (zmax>z2) zmax = z2

      nzp = nint((zmax-zmin)/drp)+1
      nxp = nint((xmax-xmin)/drp)+1

      nxr = nint((xmax-xmin)/drr) + 1
      xp1 = xmin
      xp2 = xmin + (nxr-1)*drr

      np = nxp*nzp		! number of grid nodes

      zp1 = zmin
      zp2 = zmin + (nzp-1)*drp

      if (mrf == 0) then
         nn = np
      else			! for other phases (mrf not 0),
         nn = np + nxr		! add refl/refr boundary nodes
      end if

      nn = nn + nps + 1		! source is separate graph node
                    ! each shot point (nps) is a graph node

   end subroutine graph_size



   subroutine graph_nodes

      integer :: j,m
      integer :: k
      integer :: ix,iz
      integer :: jx
      integer :: ik,jf
      integer :: id
      integer :: knot
      real :: xi,zi,xh,zh
      real :: val


      allocate(xgr(nn),zgr(nn))

       k = 0
       do ix = 1,nxp            ! loop over grid
      xi = xp1 + (ix-1)*drp
      if (xi < x1) xi = x1
      if (xi > x2) xi = x2

      do iz = 1,nzp
        zi = zp1 + (iz-1)*drp
        if (zi < z1) zi = z1
        if (zi > z2) zi = z2
            ! randomize the location of grid graph nodes (xh,zh)
        k=k+1	! increment graph number
        call my_rand(val)

        xh =  xi + 0.5*val*drp
        if (xh < xp1) xh = xp1
        if (xh > xp2) xh = xp2
        xgr(k) = xh
        call my_rand(val)

        zh =  zi + 0.5*val*drp
        if (zh < zp1) zh = zp1
        if (zh > zp2) zh = zp2
        zgr(k) = zh
     end do
       end do

      if (mrf /= 0) then	! model boundary is densely
     do ix=1,nxr		! sampled (drr small!) with graph nodes
            k = k+1
            xi = xp1 + (ix-1)*drr
        call get_zrp(xi,zi,mrf)
        xgr(k) = xi
        zgr(k) = zi
     end do
      end if

      knot = k		! "knot" is set at the current graph node #
      do j=1,npx    ! loop over picks
     igr(j) = -1
     jf = fas(j)
     if (abs(jf) == mrf) then	! if phase fits the graph
            ik = isk(j)
        do m = 1,nsh		! loop over shots
           if (ish(m) == ik) then   ! if shot matches
              id = ind(m)		! get graph # for this shot
              if (id < 0) then
             print*,'problem in graphnodes with node assignment'
             print*,m,ik,jf,id
             stop
          else
                 k = knot + id	! graph node number
             igr(j) = k		! pointer between pick and graph #
          end if
           end if
        end do
     end if
      end do

      do m = 1,nsh
         id = ind(m)
     if (id > 0) then
        k = knot + id
        xgr(k) = xsh(m)
        zgr(k) = zsh(m)
         end if
      end do
      k = knot + nps	! "nps" is the number of graph nodes added from shots

      k=k+1		! source location
      xgr(k) = xin
      zgr(k) = zin

   end subroutine graph_nodes




   subroutine init_source(irf)

      integer :: irf

      integer :: i,j,ir
      integer :: k
      integer :: ix,iz
      real :: xi,zi
      real :: zr

      integer :: krf

      allocate(qp(nn))

      allocate(tt(nn))

      allocate(pd(nn))

      allocate(pu(nn))



      do i=1,nn		! initialize
         qp(i) = -1
     pd(i) = i
     tt(i) = tinf
      end do

      allocate(hp(nn))

      if (irf == 0) then
         do i=1,nn-1
            hp(i) = -1
         end do
         ngr = nn-1
      else

     krf = irf
         xi = xgr(nn)
         zi = zgr(nn)
         do ir= irf,nr
            call get_zrp(xi,zr,ir)
        if (zr < zi) then
           krf = ir+1
        end if
         end do

     ngr = 0
     do i = 1,np		! grid nodes
        xi = xgr(i)
        zi = zgr(i)
        if (krf <= nr) then
           call get_zrp(xi,zr,krf)
        else
           zr = z2
        end if
        if (zr >= zi) then
           hp(i) = -1
           ngr = ngr + 1
        else
           hp(i) = -3
        end if
     end do

     do i = 1,nxr		! boundary nodes
            hp(np+i) = -1
        ngr=ngr+1
     end do

     do j = 1,nps		! shot nodes
        i = np+nxr+j
        xi = xgr(i)
        zi = zgr(i)
        call get_zrp(xi,zr,irf)
        if (zr >= zi) then
           hp(i) = -1
           ngr = ngr + 1
        else
           hp(i) = -3
        end if
     end do

      end if

      ! source location
      ngr=ngr+1
      i = nn		! Initialize the graph with the source
      nhp = 1		! heap size 1: The source
      qp(nhp) = i
      pd(i) = 0
      tt(i) = 0.0
      hp(i) = nhp
      k = nhp
      call sort_heap(k)



   end subroutine init_source



   subroutine slavg(xa,za,xb,zb,slog)

!		Rapid (i.e., less accurate) calculation of
!		average slowness along a line from "a" to "b".
!		Sampling of "nearby" grid points with equal
!		weights. Hence, the larger the grid spacing,
!		the larger the deviation from a true line integral

      real :: xa,za,xb,zb
      real :: slog

      integer :: ixa,iza,ixb,izb
      integer :: ixab,izab
      integer :: kx,kz,mx,mz
      integer :: igx,igz
      integer :: nsum
      real :: sums
      integer :: k
      integer :: ie,ik
      real :: vik,vmin
      real :: slok
      integer :: icut

      icut = 0
      vmin = 1.0/cmax

      ixa = nint((xa-x1)/dx)+1
      IF (ixa<1) ixa=1
      IF (ixa> nx) ixa = nx
      iza = nint((za-z1)/dz)+1
      IF (iza<1) iza=1
      IF (iza> nz) iza = nz

      ixb = nint((xb-x1)/dx)+1
      IF (ixb<1) ixb=1
      IF (ixb> nx) ixb = nx
      izb = nint((zb-z1)/dz)+1
      IF (izb<1) izb=1
      IF (izb> nz) izb = nz

      ixab = ixb-ixa
      izab = izb-iza

      mx = ixab
      mz = izab

      igx = 1
      igz = 1

      if (ixab < 0) then
         mx = -mx
         igx = -1
      endif

      if (izab < 0) then
         mz = -mz
         igz = -1
      endif

      sums = 0.0
      nsum = 0
      if (mx > mz) then

         ie = 2*mz-mx
       kz = iza
         do kx = ixa,ixa+mx

        ik = (kx-1)*nz + kz
        if (ie > 0) then
               kz=kz+igz
               ie = ie + 2*(mz-mx)
            else
               ie = ie + 2*mz
            end if

        vik = vir(ik)
        if (vik >= vmin) then
           slok = 1.0/vik
           sums = sums + slok
           nsum = nsum + 1
        else
           icut = 1
        end if

     end do
      else

     ie = 2*mx-mz
     kx = ixa
         do kz=iza,iza+mz
        ik = (kx-1)*nz + kz
        if (ie > 0) then
               kx=kx+igx
               ie = ie + 2*(mx-mz)
            else
               ie = ie + 2*mx
            end if

        vik = vir(ik)
        if (vik >= vmin) then
           slok = 1.0/vik
           sums = sums + slok
           nsum = nsum + 1
        else
           icut = 1
        end if

     end do

      end if

      if (icut == 1) then
         slog = cmax
      else
         slog = sums/nsum
      end if

   end subroutine slavg






   SUBROUTINE sort_heap(k)
!		This classic (I wish the programming style was classy) algorithm
!		sorts the heap of active nodes (size nhp) in a binary tree. The
!		top of the tree is always the node with the smallest traveltime (tt).

      INTEGER :: k

      INTEGER :: m
      INTEGER :: inw

!	*************************************

1     IF (k > nhp/2) GOTO 3

      m = 2*k
      IF (m == nhp) goto 2

      IF (tt(qp(m+1)) < tt(qp(m))) m=m+1
2     IF (tt(qp(k)) <= tt(qp(m))) GOTO 3

      hp(qp(m)) = k
      hp(qp(k)) = m
      inw = qp(m)
      qp(m) = qp(k)
      qp(k) = inw
      k = m
      GOTO 1

3     IF (k == 1) GOTO 4
      m = k/2
      IF (tt(qp(m)) <= tt(qp(k))) GOTO 4
      hp(qp(m)) = k
      hp(qp(k)) = m
      inw = qp(m)
      qp(m) = qp(k)
      qp(k) = inw
      k = m
      GOTO 3

4     CONTINUE

   END SUBROUTINE sort_heap




   subroutine loop_star(p)

      integer :: p(:)

      integer :: in
      real :: ti,tij
      real :: tnew
      integer :: im,j
      integer :: k
      real :: xi,zi
      real :: xj,zj
      real :: xij,zij,rij
      real :: cij
      integer :: nprint
      integer :: iper

      integer :: ik,ig

      nprint = ngr/20	! show progress every 5% of calculation

      print*,'Percentage / traveltime '

      do in = 1,ngr

         if (nhp > 0) then

        im = qp(1)		! set aside grid node with smallest t
            hp(qp(nhp)) = 1
        hp(qp(1)) = -2	! Take out of graph
        ti = tt(im)
        qp(1) = qp(nhp)
        qp(nhp) = -5	! Take out of graph
        nhp = nhp - 1

        if (nhp > 0) then
           k=1					! sort the heap
           call sort_heap(k)
        end if

            xi = xgr(im)
            zi = zgr(im)


        do j = 1,nn				! loop over receiving nodes

           if (hp(j) > -2 .and. im /= j) then

          xj = xgr(j)
          xij = abs(xi-xj)
          if (xij < star) then

              zj = zgr(j)
          zij = abs(zi-zj)
          if (zij < star) then

          rij = sqrt(xij*xij+zij*zij)
          if (rij < star .and. rij >= drp) then

             call slavg(xi,zi,xj,zj,cij)

             tij = rij*cij

                 tnew = ti + tij

             if (tnew < tt(j)) then	! need an update

                    tt(j) = tnew
                    p(j) = im
                    k = hp(j)

!                    in case node wasn't in heap yet (h=-1)
!                    add it, so nhp increases:
                    if (k == -1) then        ! increase heap
                       nhp = nhp +1
                       k = nhp
                       hp(j) = nhp
                       qp(nhp)= j
                end if           ! increase heap

!               rearrange the heap:
            call sort_heap(k)
             end if      ! node not excluded
          end if
          end if
          end if
           end if
        end do  ! end loop over receiving nodes

        hp(im) = -2

        IF (mod(in,nprint) == 0 )  then
           iper = nint(in*100.0/(ngr))
            write (6,'(i6,f10.3)') iper,ti
          ! write (6,'(i6,i8,f11.3)') iper,nhp,ti
        END IF

        qp(ngr+1-in) = im
     end if
      end do	! loop over in (fw stars)

      print*

   end subroutine loop_star



   subroutine init_boundary(irf)

      integer :: irf
      integer :: i,j,k
      real :: xi,zi
      real :: zr

      do i=1,nn
         qp(i) = -1
      end do

      do i=1,nn-1
         hp(i) = -1
     pu(i) = i
      end do

      do i = 1,np		! init grid nodes
         tt(i) = tinf
      end do

      do j = 1,nps		! init shot nodes
         i = np+nxr+j
         tt(i) = tinf
      end do

      ngr = nxr
      if (irf > 0 .and. irf <= nr) then

         do i = 1,np		! grid nodes
            xi = xgr(i)
            zi = zgr(i)
            call get_zrp(xi,zr,irf)
            if (zr >= zi) then
           hp(i) = -1
           ngr = ngr + 1
            else
           hp(i) = -3
            end if
         end do



     do j = 1,nps		! shot point nodes
        i = np+nxr+j
        xi = xgr(i)
            zi = zgr(i)
            call get_zrp(xi,zr,irf)
            if (zr >= zi) then
           hp(i) = -1
           ngr = ngr + 1
            else
           hp(i) = -3
             end if

     end do

      else
         ngr = nn-1
      end if

      hp(nn) = -3		! source node is not active
      pu(nn) = i
      tt(nn) = tinf

      nhp = 0
      do i = 1,nxr
         nhp=nhp+1
     qp(nhp) = np+i
     hp(np+i) = nhp
     pu(np+i) = 0
     k = nhp
     call sort_heap(k)

      end do

   end subroutine init_boundary




   subroutine build_rays

      integer :: j,m,k
      integer :: kj
      integer :: pp
      integer :: nj
      integer :: kfix
      integer :: kref

      integer, parameter :: nam = 5

      integer :: iam
      real :: to1,to2,to3,to4
      real :: xst,rdist

      real :: xj1,zj1,xj2,zj2



      kref = abs(jrf)

      print*,'building rays :'
      print*,'Instr /    Shot / Phase / Pts / Dist /  Tcalc /  Tpick /  Tbend /  dT '
      do j = 1,npx          ! loop over picks
         k = igr(j)

     if (k > 0 .and. jrf == fas(j)) then	! if connected to graph

        tca(j) = tt(k)
        nj = 1
        xry(nrt+nj) = xgr(k)	! first ray point is the shot location
        zry(nrt+nj) = zgr(k)

        xst = xgr(k)
        rdist = abs(xin-xst)

        if (hp(k) /= -2) then
           len(j) = 1
           print*,'shot oob ',j,isk(j),k,hp(k),xgr(k),zgr(k),tt(k)
        else

           if (jrf==0) then
          kj = pd(k)
           else
          kj = pu(k)
           end if

           do while (kj /= 0)
          nj=nj+1
          k=kj
          xry(nrt+nj) = xgr(k)
          zry(nrt+nj) = zgr(k)
          if (jrf == 0) then
                 kj = pd(k)
          else
                 kj = pu(k)
          end if
           end do

           kfix = nj

           if (jrf == 0) then
          len(j) = nj
          kfix = 0
           else
          kj = pd(k)
          do while (kj /= 0)
                 nj=nj+1
             k=kj
                 xry(nrt+nj) = xgr(k)
                 zry(nrt+nj) = zgr(k)
             kj = pd(k)
          end do
          len(j) = nj
           end if

           to1 = tca(j)

           do iam = 1,nam        ! a few loops
               call remove_kinks(j,kfix)
               call sample_ray(j,kfix)
              if (iam == 1) then
                 call calc_tt(j)
                 to2 = tca(j)
              end if
               call smooth_ray(j,kfix)
               call bend(j,kfix,kref)
               call remove_knots(j,kfix)

               if (len(j) > nry) then
                 print*,'ray path exceeds array size in buildrays '
                 print*,'ray group ',mrf,jrf
                 print*,'size ',nry,len(j)
                 stop
              endif

           end do

           xj1 = xry(nrt+1)
           xj2 = xry(nrt+len(j))
           zj1 = zry(nrt+1)
           zj2 = zry(nrt+len(j))


           call calc_tt(j)
           nrt=nrt+ len(j)

           if (nrt > npa) then
              print*,'rays out of space in buildrays '
              print*,'ray group ',mrf,jrf
              print*,'sizes ',nrt,npa
              stop
           endif


           to3 = tca(j) - to2
           to4 = ttp(j) - tca(j)



         if (mod(j,5)==1) &
         write(6,'(i5,i11,2i5,5f9.3)') ins,isk(j),fas(j),len(j),rdist,tca(j),ttp(j),to3,to4


        !   if (mod(j,10)==1) &
        !   write(6,'(i5,i11,2i5,7f9.3)') ins,isk(j),fas(j),len(j),rdist,tca(j),ttp(j),xj1,zj1,xj2,zj2

        end if
     end if
      end do

      print*

   end subroutine build_rays




   subroutine smooth_ray(j,kfix)

      integer :: j,kfix

      integer, parameter :: nsm = 3	! # sweeps
      real, parameter :: acc = 0.80	! how much of neighbors
      real :: xh,zh
      integer :: k,n


      do n = 1,nsm

         do k=1,len(j)
            sxr(k) = xry(nrt+k)
            szr(k) = zry(nrt+k)
         end do

     do k=2,len(j)-1
        if (k /= kfix) then


           xh = 0.5*(sxr(k-1)+sxr(k+1))
           zh = 0.5*(szr(k-1)+szr(k+1))

           xry(nrt+k) = (1.0-acc)*sxr(k) + acc*xh
           zry(nrt+k) = (1.0-acc)*szr(k) + acc*zh

        end if
     end do


      end do

   end subroutine smooth_ray



   subroutine remove_kinks(j,kfix)

      integer :: j
      integer :: kfix
      integer :: nkill
      integer :: n,nj
      integer :: ns
      real :: xv,zv
      real :: xw,zw
      real :: xi,zi
      real :: hx,hz
      real :: hr
      real :: fr
      integer :: nt,it
      integer :: nfix

      real :: xa,xb,xc
      real :: za,zb,zc
      real :: xab,xbc
      real :: zab,zbc
      real :: rab,rbc
      real :: cosb

      real, parameter :: angl = 150.0	! sharpest angle (180 = straight only).
      real :: cosc
      cosc = cos(3.141592*angl/180.0)

      nj = len(j)

      ns = 1
      nkill = 0
      sxr(ns) = xry(nrt+1)
      szr(ns) = zry(nrt+1)
      do n = 2,nj-1
         if (n /= kfix) then

        xa =  xry(nrt+n-1)
        xb =  xry(nrt+n)
        xc =  xry(nrt+n+1)

        za =  zry(nrt+n-1)
        zb =  zry(nrt+n)
        zc =  zry(nrt+n+1)

        xab = xa-xb
        xbc = xc-xb
        zab = za-zb
        zbc = zc-zb
        rab = sqrt(xab*xab + zab*zab)
        rbc = sqrt(xbc*xbc + zbc*zbc)

        cosb = (xab*xbc+zab*zbc)/(rab*rbc)
        if (cosb < cosc) then

           ns = ns+1
           sxr(ns) = xry(nrt+n)
               szr(ns) = zry(nrt+n)
        end if
     else
        ns = ns+1
        nfix = ns
        sxr(ns) = xry(nrt+n)
            szr(ns) = zry(nrt+n)
     end if

      end do

      ns = ns+1
      sxr(ns) = xry(nrt+nj)
      szr(ns) = zry(nrt+nj)
      len(j) = ns

      do n=1,len(j)
         xry(nrt+n) =  sxr(n)
         zry(nrt+n) =  szr(n)
      end do
      kfix=nfix

   end  subroutine remove_kinks




    subroutine remove_knots(j,kfix)

      integer :: j
      integer :: kfix
      integer :: nkill
      integer :: n,nj
      integer :: ns
      real :: xv,zv
      real :: xw,zw
      real :: xi,zi
      real :: hx,hz
      real :: hr
      real :: fr
      integer :: nt,it
      integer :: nfix

      real :: xa,xb
      real :: za,zb
      real :: xab
      real :: zab
      real :: rab
      real :: cosb

      real :: dmin,dcur

      dmin = 0.6*drp*fpr

      dcur = 0.0

      nj = len(j)

      ns = 1
      nkill = 0
      sxr(ns) = xry(nrt+1)
      szr(ns) = zry(nrt+1)
      do n = 2,nj-1
         if (n /= kfix) then

        xa =  xry(nrt+n-1)
        xb =  xry(nrt+n)
        za =  zry(nrt+n-1)
        zb =  zry(nrt+n)
        xab = xa-xb
        zab = za-zb
        rab = sqrt(xab*xab + zab*zab)
        dcur = dcur + rab

        if (dcur > dmin) then

           ns = ns+1
           sxr(ns) = xry(nrt+n)
               szr(ns) = zry(nrt+n)
           dcur = 0.0
        end if
     else
        ns = ns+1
        nfix = ns
        sxr(ns) = xry(nrt+n)
            szr(ns) = zry(nrt+n)
        dcur = 0.0
     end if

      end do

      ns = ns+1
      sxr(ns) = xry(nrt+nj)
      szr(ns) = zry(nrt+nj)
      len(j) = ns

      do n=1,len(j)
         xry(nrt+n) =  sxr(n)
         zry(nrt+n) =  szr(n)
      end do
      kfix = nfix

   end  subroutine remove_knots


   subroutine bend(j,kfix,kref)

      integer :: j,kfix
      integer :: kref

      integer, parameter :: icut = 5
      real, parameter :: gcrit = 0.07
      integer, parameter :: nw = 41
      integer, parameter :: nsweep = 17
      integer, parameter :: nmin = 4
      real, parameter :: wa = -0.9
      real, parameter :: wb = 2.4
      integer :: n
      real :: to2

      real :: wm,dw
      real :: gm
      real :: tm
      integer :: m
      real :: gmin,tmin
      real :: xadd,zadd
      real :: xnew,znew
      integer :: k
      integer :: icnt
      real :: tup
      real :: zr


      dw = (wb-wa)/(nw-1)

      frx = 0.0
      frz = 0.0
      grx = 0.0
      grz = 0.0
      hrx = 0.0
      hrz = 0.0

      tup = tca(j)
      icnt = 0

      do n=1,nsweep

         if (n>1) then
        do k=1,len(j)
           frx(k) = grx(k)
           frz(k) = grz(k)
        end do
     end if

     call ray_grad(j,kfix)

     call conj_grad(n,j)

     gmin = 0.0
     tmin = tup

     do m = 1,nw
            wm = wa + (m-1)*dw
        gm = wm*wm
        if (wm < 0) gm = -gm

        call probe_tt(j,kref,kfix,gm,tm)
        if (tm < tmin) then
           gmin = gm
           tmin = tm
        end if
     end do

     if (abs(gmin) <= gcrit) then
        icnt = icnt + 1
     else
        icnt = 0
     end if

     do k=1,len(j)

        xadd = gmin*hrx(k)
            zadd = gmin*hrz(k)
            xnew = xry(nrt+k)+xadd
            znew = zry(nrt+k)+zadd
        if (xnew < x1) xnew=x1
        if (xnew > x2) xnew=x2
        if (znew < z1) znew=z1
        if (znew > z2) znew=z2

        if (k == kfix .and. kref /=0) then
           call get_zrp(xnew,zr,kref)
           znew = zr
            end if
        xry(nrt+k) = xnew
        zry(nrt+k) = znew

     end do
     call calc_tt(j)
     tup = tca(j)

     if (icnt >= icut .and. n > nmin) exit

      end do

   end subroutine bend



   subroutine ray_grad(j,kfix)

      integer :: j
      integer :: kfix

      integer :: nj
      integer :: k
      real :: xa,xb,xc
      real :: za,zb,zc
      real :: xab,xbc,xac
      real :: zab,zbc,zac

      real :: hxac,hzac
      real :: rab,rbc
      real :: rac
      real :: rt
      real :: rdx,rdz
      real :: gx,gz
      real :: va,vb,vc
      real :: ua,ub,uc
      real :: uab,ubc

      real :: gtx,gtz
      real :: glx,glz
      real :: gpx,gpz

      real :: prg



      nj = len(j)

      do k=1,nj-1

         if (k>1.and. k<nj .and. k /= kfix) then

        xa = xry(nrt+k-1)
        xb = xry(nrt+k  )
        xc = xry(nrt+k+1)

        za = zry(nrt+k-1)
        zb = zry(nrt+k  )
        zc = zry(nrt+k+1)

        xab = xb - xa
        zab = zb - za

        xbc = xb - xc
        zbc = zb - zc

        xac = xc - xa
        zac = zc - za

        rab = sqrt(xab*xab + zab*zab)
        rbc = sqrt(xbc*xbc + zbc*zbc)
        rac = sqrt(xac*xac + zac*zac)

        hxac = xac/rac
        hzac = zac/rac

        rt = rab + rbc

        call grdu(xb,zb,gx,gz)

        call get_vel(xa,za,va)
        ua = 1.0/va
        call get_vel(xb,zb,vb)
        ub = 1.0/vb
        call get_vel(xc,zc,vc)
        uc = 1.0/vc

        uab = 0.5*(ua + ub)
        ubc = 0.5*(ub + uc)

        gtx = 0.5*rt*gx + xab*uab/rab + xbc*ubc/rbc
        gtz = 0.5*rt*gz + zab*uab/rab + zbc*ubc/rbc

        prg = gtx*hxac + gtz*hzac

        glx = prg*hxac
        glz = prg*hzac

        gpx = gtx - glx
        gpz = gtz - glz


        grx(k) = gpx
            grz(k) = gpz

     else
        grx(k) = 0.0
        grz(k) = 0.0

     end if

      end do

      if (kfix > 1 .and. kfix < nj) then

         grx(kfix) = 0.5*(grx(kfix-1) + grx(kfix+1))
         grz(kfix) = 0.5*(grz(kfix-1) + grz(kfix+1))
      end if

   end subroutine ray_grad



   subroutine conj_grad(n,j)

      integer :: n,j

      integer :: k
      real :: f1,f2
      real :: fac
      real :: px,pz

      if (n==1) then

         do k = 1,len(j)
            hrx(k) = -grx(k)
        hrz(k) = -grz(k)
         end do
      else

     f2 = 0.0
     do k=1,len(j)
       px = frx(k)
       pz = frz(k)
       f2 = f2 + px*px+pz*pz
     end do

     f1 = 0.0
     do k=1,len(j)
       px = grx(k)
       pz = grz(k)
       f1 = f1 + px*px+pz*pz
     end do

     fac = f1/f2

     do k=1,len(j)
        hrx(k) = fac*hrx(k) - grx(k)
        hrz(k) = fac*hrz(k) - grz(k)
     end do

      end if

   end subroutine conj_grad



   subroutine grdu(x,z,gx,gz)

!		approximate local velocity gradient
      real :: x,z,gx,gz

      integer :: ix,iz
      integer :: ixa,ixb
      integer :: iza,izb
      integer :: ia,ib
      real :: dxab,dzab
      real :: duab
      real :: va,vb
      real :: ua,ub

      ix = nint((x-x1)/dx)+1
      if (ix < 1) ix = 1
      if (ix > nx) ix=nx
      iz = nint((z-z1)/dz)+1
      if (iz < 1) iz = 1
      if (iz > nz) iz=nz

      ixa = ix-1
      if (ixa < 1) ixa = 1
      ixb = ix+1
      if (ixb  > nx) ixb = nx

      ia = (ixa-1)*nz+iz
      ib = (ixb-1)*nz+iz
      dxab = (ixb-ixa)*dx
      va = vir(ia)
      vb = vir(ib)

      ua = 1.0/va
      ub = 1.0/vb
      duab = ub - ua
      gx = duab/dxab

      iza = iz-1
      if (iza < 1) iza = 1
      izb = iz+1
      if (izb  > nz) izb = nz

      ia = (ix-1)*nz+iza
      ib = (ix-1)*nz+izb
      dzab = (izb-iza)*dz
      va = vir(ia)
      vb = vir(ib)

      ua = 1.0/va
      ub = 1.0/vb
      duab = ub - ua
      gz = duab/dzab

   end subroutine grdu




   subroutine calc_tt(j)
        !   recalculate traveltime with a line
        !   integral of model slowness along the ray path

      integer :: j	! ray #

      integer :: m	! ray point #
      real :: xa,za,xb,zb,xh,zh
      real :: dxab,dzab,drab
      real :: tray
      real :: dsum
      real :: slab
      real :: va,vb,vh

      tray = 0.0
      do m=2,len(j)

         if (m==2) then

        xa = xry(nrt+m-1)
            za = zry(nrt+m-1)
        call get_vel(xa,za,va)
        if (va < 1.0/cmax) va = 1.0/cmax
     else
        xa = xb
        za = zb
        va = vb
     end if

     xb = xry(nrt+m)
         zb = zry(nrt+m)
     call get_vel(xb,zb,vb)
     if (vb < 1.0/cmax) vb = 1.0/cmax

     xh = 0.5*(xa+xb)
     zh = 0.5*(za+zb)
     call get_vel(xh,zh,vh)
     if (vh < 1.0/cmax) vh = 1.0/cmax

     slab = 4.0/(va+2*vh+vb)	! average slowness

     dxab = xb-xa
     dzab = zb-za
     drab = sqrt(dxab*dxab + dzab*dzab)

     tray = tray + drab*slab

      end do

      tca(j) = tray

   end subroutine calc_tt


   subroutine probe_tt(j,kref,kfix,wm,tm)
        !   recalculate traveltime with a line
        !   integral of model slowness along the ray path

      integer :: j	! ray #
      integer :: kref
      integer :: kfix

      real :: wm	!
      real :: tm

      integer :: m	! ray point #
      real :: xa,za,xb,zb,xh,zh
      real :: dxab,dzab,drab
      real :: tray
      real :: dsum
      real :: slab
      real :: va,vb,vh
      real :: zr

    !  dsum = 0.0
      tray = 0.0
      do m=2,len(j)

         if (m==2) then

        xa = xry(nrt+m-1) + wm*hrx(m-1)
            za = zry(nrt+m-1) + wm*hrz(m-1)

        if (xa < x1) xa=x1
        if (xa > x2) xa=x2
        if (za < z1) za=z1
        if (za > z2) za=z2


        if (m-1 == kfix .and. kref /=0) then
           call get_zrp(xa,zr,kref)
           za = zr
        end if

        call get_vel(xa,za,va)
        if (va < 1.0/cmax) va = 1.0/cmax
     else
        xa = xb
        za = zb
        va = vb
     end if

     xb = xry(nrt+m) + wm*hrx(m)
         zb = zry(nrt+m) + wm*hrz(m)

     if (xb < x1) xb=x1
     if (xb > x2) xb=x2
     if (zb < z1) zb=z1
     if (zb > z2) zb=z2


     if (m == kfix .and. kref /=0) then
        call get_zrp(xb,zr,kref)
        zb = zr
     end if

     call get_vel(xb,zb,vb)
     if (vb < 1.0/cmax) vb = 1.0/cmax

     xh = 0.5*(xa+xb)
     zh = 0.5*(za+zb)
     call get_vel(xh,zh,vh)
     if (vh < 1.0/cmax) vh = 1.0/cmax

     slab = 4.0/(va+2*vh+vb)	! average slowness

     dxab = xb-xa
     dzab = zb-za
     drab = sqrt(dxab*dxab + dzab*dzab)

     tray = tray + drab*slab

      end do

      tm = tray

   end subroutine probe_tt





   subroutine sample_ray(j,kfix)
      integer :: j
      integer :: kfix
      integer :: nkill
      integer :: n,nj
      integer :: ns
      real :: xv,zv
      real :: xw,zw
      real :: xi,zi
      real :: hx,hz
      real :: hr
      real :: fr
      integer :: nt,it
      integer :: nfix

      real :: drs

      drs = 1.1*drp*fpr


      nfix = 1
      nj = len(j)
      ns = 0

      do n = 1,nj

         xw = xry(nrt+n)
     zw = zry(nrt+n)

         if (n > 1) then
        hx = xw-xv
        hz = zw-zv
        hr = sqrt(hx*hx+hz*hz)

        nt = nint(hr/drs) - 1

        if (nt > 0) then
           do it = 1,nt
              fr = float(it)/float(nt+1)
              xi = xv + fr*hx
          zi = zv + fr*hz
          ns = ns + 1
          sxr(ns) = xi
          szr(ns) = zi
         ! write(6,'(i5,2f9.3)') ns,xi,zi
           end do
        end if
     end if

     ns=ns+1

     if (n == kfix) nfix = ns
     sxr(ns) = xw
         szr(ns) = zw


     xv=xw
     zv=zw

      end do

      kfix = nfix
      nj = ns
      do n=1,nj
         xry(nrt+n) = sxr(n)
     zry(nrt+n) = szr(n)
      end do

      len(j) = nj



   end subroutine sample_ray



   subroutine get_zrp(xi,zr,irf)

      real :: xi,zr
      integer :: irf

      integer :: ix
      integer :: kx
      integer :: mx
      integer :: m
      integer :: jx

      real :: xj,xm
      real :: ex(0:1)
      real, parameter :: zero = 0.0
      real, parameter :: one = 1.0


      zr = zero

      jx = int((xi-x1)/dx)+1
      jx = max(1,jx)
      jx = min(nx-1,jx)
      xj = x1+(jx-1)*dx
      xm = xi-xj
      ex(1) = xm/dx
      ex(0) = one - ex(1)
      do kx=0,1
         mx=jx+kx
     m = (irf-1)*nx + mx
     zr = zr + ex(kx)*zrf(m)
      end do

   end subroutine get_zrp



   subroutine get_vel(xi,zi,vr)

      real :: xi,zi
      real :: vr

      integer :: ix,iz
      integer :: kx,kz
      integer :: mx,mz
      integer :: m
      integer :: jx,jz

      real :: xj,xm
      real :: zj,zm
      real :: ex(0:1), ez(0:1)
      real, parameter :: zero = 0.0
      real, parameter :: tiny = 0.001
      real, parameter :: one = 1.0

      real :: velm
      real :: w,wsum
      real :: vsum

      vr = zero

      jx = int((xi-x1)/dx)+1
      jx = max(1,jx)
      jx = min(nx-1,jx)
      xj = x1+(jx-1)*dx
      xm = xi-xj
      ex(1) = xm/dx
      ex(0) = one - ex(1)


      jz = int((zi-z1)/dz)+1
      jz = max(1,jz)
      jz = min(nz-1,jz)
      zj = z1+(jz-1)*dz
      zm = zi-zj
      ez(1) = zm/dz
      ez(0) = one - ez(1)

      wsum = zero
      vsum = zero
      do kx=0,1
         mx=jx+kx
     do kz=0,1
            mz=jz+kz
        m = (mx-1)*nz+mz
        w = ex(kx)*ez(kz)
        velm = vir(m)
        if (velm >= 1.0/cmax) then
          wsum = wsum + w
          vsum = vsum + w*velm
        end if
     end do
      end do

      if (wsum > tiny) then
         vr = vsum/wsum
      else
         vr = 1.0/cmax
      end if

   end subroutine get_vel



   subroutine fill_vel(irf)

      integer :: irf

      integer :: krf
      integer :: i,m

      integer :: ix,iz
      integer :: ir,ic
      real :: za,zb

      real :: zi,xi,zr
      integer :: ire

      krf = irf
      if (krf == 0) then
         krf = nr + 1
      else
         xi = xgr(nn)
         zi = zgr(nn)
         do ir= irf,nr
            call get_zrp(xi,zr,ir)
        if (zr < zi) then
           krf = ir+1
        end if
         end do
      end if

      do ix=1,nx
         do ir = 1,krf

        if (ir == 1) then
           za = z1
           zb = zrf(ic)
        elseif (ir <= nr) then
           ic = (ir-2)*nx+ix
           za = zrf(ic)
           zb = zrf(ic+nx)
        else
           ic = (nr-1)*nx+ix
           za = zrf(ic)
           zb = z2
        end if

        zb = z2

        do iz = 1,nz
           i = (ix-1)*nz+iz
           m = (ir-1)*nx*nz + i
           zi = z1 + (iz-1)*dz
           if (zi >= za .and. zi <= zb) vir(i) = vel(m)
            end do
     end do
      end do


   end subroutine fill_vel



   subroutine output_graph

      character(50) :: graph_file
      integer :: i
      real :: x,z,ti

      read(5,'(a)') graph_file
      open(unit=14,file=graph_file)

      do i  = 1,nn
     x = xgr(i)
     z = zgr(i)
     ti = tt(i)
     if (ti > 0.99*tinf) ti = -1.00
     write(14,'(3f11.4)') x,z,ti
      end do
      close(14)

   end subroutine output_graph




   subroutine close_graph(mrf)

      integer :: mrf

      deallocate(xgr,zgr)
      deallocate(tt)
      deallocate(hp)
      deallocate(qp)
      deallocate(pd)
      deallocate(pu)


   end subroutine close_graph



   subroutine write_rays

      character(80) :: rayfile
      integer :: inew
      integer :: inp
      integer :: num,inum
      integer ::  nwd,ipos
      integer :: nh,nq
      integer :: j,k,m
      real :: xi,zi,zr
      real :: xa,za,xb,zb,xj1,xj2,zj1,zj2

      integer :: ncnt

      num = 0
      k = 0

      read(5,*) inew



      read(5,'(a)') rayfile
      open(unit=17,file=rayfile,access='stream')

      ipos = 5
      if (inew == 0) then

     read(17) num

     nwd = 1
     ipos = 1 + 4*nwd

     do inum=1,num

        read(17,pos=ipos) inp,nh,nq

        nwd = nwd + 3 + 6*nh + 2*nq
        ipos = 1 + 4*nwd

     end do

      end if

      num = num+1

      ncnt=0

      do j=1,npx

         xj1 = xry(ncnt+1)
     xj2 = xry(ncnt+len(j))
     zj1 = zry(ncnt+1)
     zj2 = zry(ncnt+len(j))
       !  write(6,'(3i5,2i13,6f9.3)') ins,j,fas(j),isk(j),len(j),ttp(j),tca(j),xj1,zj1,xj2,zj2
     ncnt=ncnt+len(j)
      end do


      write(17,pos=1) num
      write(17,pos=ipos) ins,npx,nrt
      write(17) (isk(j),j=1,npx)
      write(17) (fas(j),j=1,npx)
      write(17) (ttp(j),j=1,npx)
      write(17) (etp(j),j=1,npx)
      write(17) (tca(j),j=1,npx)
      write(17) (len(j),j=1,npx)
      write(17) (xry(k),k=1,nrt)
      write(17) (zry(k),k=1,nrt)
      close(17)

   end subroutine write_rays


end program vm_trace
