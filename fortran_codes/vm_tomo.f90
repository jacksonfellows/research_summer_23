program vm_tomo

   implicit none

   character(20) :: matfile     ! scratch file for frechet matrix

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

   real :: cmax			! max slowness (is this needed as a global variable?)

   real, allocatable :: vin(:)	! velocity model on inversion grid
   real, allocatable :: rin(:)	! reflector boundaries on inversion spacing
   integer, allocatable :: iin(:)!  ID of reflector, resampled
   integer, allocatable :: ijn(:)!  ID of reflector, active nodes only


   real :: vscal
   real :: zscal

   real :: vpmin

   real :: reg0			! damping strength
   real :: reg1			! flattening strength
   real :: reg2			! smoothing strength
   real :: crf			! relative scaling of regularization boundary depths versus velocities

   integer :: npx			! number of picks
   integer, allocatable :: inj(:)	! instrument with pick
   integer, allocatable :: isk(:)	! shot number with pick
   integer, allocatable :: fas(:)	! phase ID
   real, allocatable :: ttp(:)		! picked traveltime
   real, allocatable :: etp(:)		! pick uncertainty
   integer, allocatable :: len(:)	! the length of a ray path
   real, allocatable :: tca(:)		! calc tt for ray
   real, allocatable :: xry(:),zry(:)	! ray path storage

   integer, allocatable :: ind(:)	! pointer from solution vector to model
   integer, allocatable :: ivc(:)   ! pointer from model to solution vector

   integer :: nxc,nzc
   real :: dxc
   real, allocatable :: dzc(:)

   integer :: itop,ibot
   integer :: na1,na2           ! selected bounds on reflected phases
   integer :: nb1,nb2           ! selected bounds on refracted phases
   integer :: nc1,nc2           ! range of velocity layers in inversion
   integer :: nd1,nd2           ! range of layer boundaries in inversion

   integer :: np		! number of model parameters
   integer :: nsol		! size of solution vector
   integer :: nfz		! number of nonzero elements in frechet matrix
   integer :: nnz		! number of nonzero elements in entire matrix
   integer :: neq		! number of equations
   integer :: ndat		! number of data used in inversion

   integer, allocatable :: nrow(:)
   integer, allocatable :: nro(:)

   integer, allocatable :: inz(:)
   real, allocatable :: anz(:)
   integer, allocatable :: imz(:)
   real, allocatable :: amz(:)


   logical, allocatable :: ynh(:)	! flag covered nodes
   logical, allocatable :: yng(:)	! flag active nodes
   real, allocatable :: dws(:)	! derivative weight sum

   real, allocatable :: dat(:)	! rhs vector with data residuals
   real, allocatable :: rhs(:)	! right hand side vector for inversion
   real, allocatable :: rhc(:)	! copy of right hand side vector
   real, allocatable :: rms(:)	! rhs misfit

   real, allocatable :: dm(:)
   real, allocatable :: vvec(:)
   real, allocatable :: wvec(:)


   real :: xreg,zreg		! extent of model space beyond ray coverage
   real :: asr			! scaling horizontal and vertical derivatives in smoothing

   real :: ch2n			! target chi-squared for inversion result
   logical :: echo		! diagnostics output switch
   real :: tofa			! first guess L-curve parameter
   integer :: ilm       ! max number of iterations in lsqr
   real :: rlm          ! cut-off tolerance for lsqr
   real :: dax			! lsqr built-in damping: set to zero

   real :: vars			! variance before inversion
   real :: varn			! variance after inversion
   real :: vard			! variance reduction
   real :: pen

   real :: rega,regb
   real :: ch2a,ch2b
   integer :: icur




   call read_vm

   call get_rays

   call inv_prm

   call sample_model

   call frechet

   call out_dws

   call model_coverage

   call out_mask

   call size_matrix

   call regularization

   call init_inverse

   call solve_eqn(rega,ch2a)

   call bracket

   call secant

   call new_sol

   call update_model

   call write_vm

contains


   subroutine get_rays

      character(80) :: rayfile
      integer :: nrt
      integer :: j,k,jt,kt,m

      integer :: num,inum
      integer :: ipos,nwd
      integer :: nh,nq,inp

      read(5,'(a)') rayfile
      open(unit=18,file=rayfile,access='stream')
      read(18) num

      nrt = 0
      npx = 0
      nwd = 1
      ipos = 1+4*nwd
      do inum = 1,num
         read(18,pos=ipos) inp,nh,nq
     npx = npx + nh
     nrt = nrt + nq
     nwd = nwd + 3 + 6*nh + 2*nq
     ipos = 1+4*nwd
      end do

      allocate(isk(npx),inj(npx),fas(npx),ttp(npx),etp(npx),tca(npx),len(npx))
      allocate(xry(nrt),zry(nrt))

      ipos = 1
      read(18,pos=ipos) num
      jt = 0
      kt = 0
      do inum = 1,num
         read(18) inp,nh,nq
     do j = 1,nh
        inj(jt+j) = inp
     end do
     do j = 1,nh
        read(18) isk(jt+j)
     end do
     do j = 1,nh
        read(18) fas(jt+j)
     end do
     do j = 1,nh
        read(18) ttp(jt+j)
     end do
     do j = 1,nh
        read(18) etp(jt+j)
     end do
         do j = 1,nh
        read(18) tca(jt+j)
     end do
     do j = 1,nh
        read(18) len(jt+j)
     end do
         do k = 1,nq
        read(18) xry(kt+k)
     end do
         do k = 1,nq
        read(18) zry(kt+k)
     end do

     jt=jt+nh
     kt=kt+nq

      end do
      close(18)

   end subroutine get_rays




   subroutine read_vm

      character(60) :: vmfile
      integer :: j

      integer :: m,ix,iz,i,ii
      real :: zi

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

   end subroutine read_vm




   subroutine inv_prm

      real :: f,z
      integer :: iz

      read(5,*) cmax

      read(5,*) nxc,nzc
      dxc=(x2-x1)/(nxc-1)


      allocate(dzc(nzc-1)) ! The grid spacing can be made to increase linearly with depth
      read(5,*) dzc(1)          ! vertical grid spacing near surface


      f = 2*(z2 - z1 - (nzc-1)*dzc(1))/((nzc-1)*(nzc-2))

      z = z1+dzc(1)
      do iz=2,nzc-1
         dzc(iz) = dzc(iz-1)+f
     z=z+dzc(iz)
    ! write(6,'(i5,2f9.3)') iz,z,dzc(iz)
      end do

      read(5,*) itop,ibot		! shallowest and deepest phase used in inversion



      if (itop < 0) then
         na1 = itop
     nb1 = -itop
     nc1 = -itop
     nd1 = -itop
      else
         na1 = -itop-1
     nb1 = itop
     nc1 = itop+1
     nd1 = itop+1
      end if
      if (ibot < 0) then
         na2 = ibot
     nb2 = -ibot-1
     nc2 = -ibot
     nd2 = -ibot
      else
         na2 = -ibot
     nb2 = ibot
     nc2 = ibot+1
     nd2 = ibot
      end if


      np = (nr+1)*nxc*nzc + nr*nxc


      read(5,*) vscal,zscal		! dimension scaling velocity and depth in inversion

      read(5,*) xreg,zreg		! horizontal and vert. reach of regularization
                        ! beyond ray coverage
      read(5,*) asr			! strength of horizontal vs vertical derivatives in regularization

      read(5,*) reg0,reg1,reg2
      read(5,*) crf			! relative strength regularization of model boundaries

      read(5,*) ch2n

      echo = .FALSE.			! no diagnostics to standard output
      ilm = 120				! max number of lsqr iterations
      rlm = 1.0e-6              ! tolerance level for breaking off lsqr iterations
      dax = 0.0

      read(5,*) vpmin			! lowest velocity that can be masked in velocity plot

   end subroutine inv_prm





   subroutine sample_model

      integer :: ir,ix,iz
      integer :: mx,mz
      integer :: kx,kz,ki
      integer :: mc,m,k
      integer :: i,j,jr
      real :: xi,zi
      real :: xm,zm
      real :: ex(0:1),ez(0:1)
      real :: samp
      !  real :: x,z

      allocate(vin(nxc*nzc*(nr+1)))
      allocate(rin(nxc*nr))

      allocate(iin(nxc*nr),ijn(nxc*nr))

      do ir = 1,nr+1
         do mx=1,nxc
        xm = x1 + (mx-1)*dxc
        ix = int((xm-x1)/dx)+1
        xi = x1 + (ix-1)*dx
        ex(1) = (xm-xi)/dx
        ex(0) = 1.0 - ex(1)

        zm = z1
        do mz = 1,nzc
           iz = int((zm-z1)/dz)+1
           zi = z1 + (iz-1)*dz
           ez(1) = (zm-zi)/dz
           ez(0) = 1.0 - ez(1)


           m = (mx-1)*nzc+mz
           mc = (ir-1)*nxc*nzc+m

           samp = 0.0

           do kx=0,1
              do kz=0,1
             k = (ix+kx-1)*nz+iz+kz
             i = (ir-1)*nx*nz+k
             samp = samp + ex(kx)*ez(kz)*vel(i)
          end do
           end do
           vin(mc) = samp

           if (mz < nzc) zm = zm + dzc(mz)
        end do
         end do
      end do



      do ir=1,nr
         do mx=1,nxc
        xm = x1 + (mx-1)*dxc
        ix = int((xm-x1)/dx)+1
        xi = x1 + (ix-1)*dx

        ex(1) = (xm-xi)/dx
        ex(0) = 1.0 - ex(1)

        mc = (ir-1)*nxc+mx

        samp = 0.0

        do kx=0,1
           k = ix+kx
           i = (ir-1)*nx+k
           samp = samp + ex(kx)*zrf(i)
           if (ex(kx) >= 0.5) iin(mc) = ird(i)
        end do
        rin(mc) = samp

     end do
      end do

      do mx=1,nxc
         do ir=1,nr

        mc = (ir-1)*nxc+mx
        i = iin(mc)

        if (ir < nd1 .or. ir > nd2) then

           ijn(mc) = -1
           if (i >= nd1 .and. i <= nd2) then
              m = (i-1)*nxc+mx
          ijn(m) = -1
           end if

        else

           if (i < nd1 .or. i > nd2) then
              ijn(mc) = -1
           else
              ijn(mc) = iin(mc)
           end if

            end if

         end do
      end do


   end subroutine sample_model




   subroutine frechet

      real, parameter :: tiny = 0.003

      integer :: i,ir,k
      integer :: j
      integer :: m
      real :: xa,xb,xc,xh,xd
      real :: za,zb,zc,zh,zd
      real :: dd
      real :: dxab,dxac,dzab,dzac
      real :: dxbc,dzbc,dbc
      real :: wgt
      real :: slope,argc
      real :: dab,dac

      integer :: kx,kz

      integer :: mm,irc,mcx

      logical :: yn

      real :: zp1,zp2
      integer :: irh
      integer :: ira

      integer :: ipr,idp

      integer :: mx,mz
      real :: xm,zm
      real :: ex(0:1),ez(0:1)

      real :: vh
      real :: zr,zrb,zrc,zra
      real :: dzrb,dzrc,dzra
      real :: dslo

      integer :: ip
      integer :: im,in
      real :: elm,val
      real :: cosb
      real :: cost,cos2t
      real :: v1,v2
      real :: du

      real :: scal
      integer :: iflag
      integer :: iv


      logical, allocatable :: ynr(:)
      real, allocatable :: row(:)

      real :: dtsum,dtj,dtj2,rmst
      integer :: ntsum

      matfile = "scratch.bin"
      open(unit=7,file=matfile,access='stream')

      allocate(nrow(npx),dat(npx))
      allocate(row(np),inz(np),anz(np),ynr(np),ynh(np),dws(np))


      do ip=1,np
         ynh(ip) = .false.
         dws(ip) = 0.0
      end do

      print*
      print*,'Calculating frechet matrix with ', npx,' data. '

      m = 0
      nfz = 0
      dtsum = 0.0
      ntsum = 0
      do j = 1,npx

         nrow(j) = 0
     dat(j) = 0.0

         do ip=1,np
        row(ip) = 0.0
        ynr(ip) = .false.
     end do

     if (fas(j) == 0) then
        yn = .true.
     else if (fas(j) < 0) then
        if (fas(j) <= na1 .and. fas(j) >= na2) then
           yn = .true.
        else
           yn = .false.
        end if
     else
        if (fas(j) >= nb1 .and. fas(j) <= nb2) then
           yn = .true.
        else
           yn = .false.
        end if
     end if

         if (yn) then
     if (len(j) > 2) then
        dtj = ttp(j) - tca(j)
        dtj2 = dtj*dtj
        dat(j) = dtj/etp(j)

        do k=1,len(j)-1



           xa = xry(m+k)
           xb = xry(m+k+1)
           za = zry(m+k)
           zb = zry(m+k+1)
           xd = xb-xa
           zd = zb-za
           dd = sqrt(xd*xd+zd*zd)

           xh = 0.5*(xa+xb)
           zh = 0.5*(za+zb)

           irh = 0
           do ir = nc1,nc2
              if (ir == 1) then
             zp1 = z1
          else
                 call get_zrp(xh,zp1,ir-1)
          end if

          if (ir == nr+1) then
             zp2 = z2
          else
                 call get_zrp(xh,zp2,ir)
          end if

          if (zh >= zp1 .and. zh <= zp2) then
             irh = ir
             exit
          end if
           end do


           if (irh /= 0) then

              call get_veli(xh,zh,irh,vh)
          scal = -1.0*vscal/(vh*vh*etp(j))

          mx = int((xh-x1)/dxc) + 1
          xm = x1 + (mx-1)*dxc
          ex(1) = (xh-xm)/dxc
          ex(0) = 1.0 - ex(1)

          zm = z1
          mz = 1
          do while (zm+dzc(mz) < zh .and. mz < nzc-1)
             zm = zm + dzc(mz)
                 mz=mz+1
          end do
          ez(1) = (zh-zm)/dzc(mz)
          ez(0) = 1.0 - ez(1)

          do kx = 0,1			! velocities
                 do kz = 0,1
            wgt = ex(kx)*ez(kz)
            im = (mx+kx-1)*nzc + mz+kz
            ip = (irh-1)*nxc*nzc + im
            val = wgt*dd*scal
            ynr(ip) = .true.
            ynh(ip) = .true.
            row(ip) = row(ip) + val
            dws(ip) = dws(ip) + abs(val)
             end do
          end do

           end if


           if (k > 1) then
!							look or reflection points
          xc = xry(m+k-1)
              zc = zry(m+k-1)

          ira = 0
          do ir = nd1,nd2

             mcx = int((xa-x1)/dx) + 1
             mm = (ir-1)*nx+mcx
             call get_zrp(xa,zr,ir)

             if (abs(zr-za) < tiny*dz .and. ird(mm) == ir) then
                ira = ir
            exit
                 end if
          end do

          if (ira /= 0) then
             call get_zrp(xb,zrb,ira)
             call get_zrp(xc,zrc,ira)
          else
             zrb = 0.0
             zrc = 0.0
          end if

          if (ira /= 0 .and. zc < zrc .and. zb < zrb) then

             dxab = xb-xa
             dzab = zb-za
             dab = sqrt(dxab*dxab+dzab*dzab)

             dxac = xc-xa
             dzac = zc-za
             dac = sqrt(dxac*dxac+dzac*dzac)

             cos2t = dxab*dxac + dzab*dzac
                 cos2t = cos2t/(dab*dac)

             if (abs(cos2t) < 0.9999) then

                cost = sqrt((1.0+cos2t)/2.0)

                dxbc = xb-xc
                dzbc = zrb-zrc
            if (abs(dxbc) > tiny*dx) then
                   slope = abs(dzbc)/dxbc
                   argc = atan(slope)
                           cosb = cos(argc)
             else
                cosb = 1.0
             end if

            mcx = int((xa-x1)/dx) + 1

            do irc = ira,1,-1	! find velocity right above node
               mm = (irc-1)*nx+mcx
               if (ird(mm) == ira) call get_veli(xa,za,irc,v1)
            end do

                if (cosb*cost > 0.1) then
                   mx = int((xa-x1)/dxc) + 1
                   xm = x1 + (mx-1)*dxc
                   ex(1) = (xa-xm)/dxc
                   ex(0) = 1.0 - ex(1)
                   elm = 2*cosb*cost/v1

                   do kx = 0,1
                      im = mx+kx
                  ipr =  (ira-1)*nxc + im
                      ip = (nr+1)*nxc*nzc + ipr

    !  if this node, or any other dependent node is not to be modified,
        ! it will not move in the inversion.

                  iflag = 1
                  do irc = 1,nr	! find velocity right above node
                     mm = (irc-1)*nx+mcx
                     if (ird(mm) == ira) then
                    if (irc < nd1 .or. irc > nd2) iflag = 0
                 end if
                  end do

                  idp = ijn(ipr)
                  if (idp == ira .and. iflag==1) then
                            ! if the boundary node is independent
                     val = ex(kx)*elm*zscal/etp(j)

                         ynr(ip) = .true.
                         ynh(ip) = .true.
                         row(ip) = row(ip) + val
                         dws(ip) = dws(ip) + abs(val)

                  end if

                   end do

                end if          ! cosb, cost > 0.1
             end if			! cos2t < 1
          end if			! ira /= 0



          do ir = nd1,nd2		! now look for pierce points

             call get_zrp(xb,zrb,ir)
             call get_zrp(xc,zrc,ir)
             dzrb = zrb - zb
             dzrc = zrc - zc

             mcx = int((xa-x1)/dx) + 1
             mm = (ir-1)*nx+mcx

             if (dzrb*dzrc < 0 .and. ird(mm) == ir) then

                call get_zrp(xa,zra,ir)
            dzra = zra - za
            if (abs(dzra) < abs(dzrb) .and. abs(dzra) < abs(dzrc)) then

               mcx = int((xa-x1)/dx) + 1
               do irc = ir,1,-1	! find velocity right above node
                  mm = (irc-1)*nx+mcx
                  if (ird(mm) == ir) call get_veli(xa,za,irc,v1)
               end do

               call get_veli(xa,za,ir+1,v2)

                       du = 1.0/v1 - 1.0/v2

               dxbc = xb-xc
               dzbc = abs(zb-zc)
               dbc = sqrt(dxbc*dxbc+dzbc*dzbc)
               cost = dzbc/dbc

               dxbc = xb-xc
                   dzbc = zrb-zrc

               if (abs(dxbc) > tiny*dx) then
                      slope = abs(dzbc)/dxbc
                      argc = atan(slope)
                              cosb = cos(argc)
               else
                  cosb = 1.0
               end if

               if (cost*cosb > 0.1) then

                  elm = cosb*cost*du
                      mx = int((xa-x1)/dxc) + 1
                      xm = x1 + (mx-1)*dxc
                      ex(1) = (xa-xm)/dxc
                      ex(0) = 1.0 - ex(1)
                      do kx = 0,1
                     im = mx+kx
                 ipr = (ir-1)*nxc + im
                     ip = (nr+1)*nxc*nzc + ipr
                 idp = ijn(ipr)

                 iflag = 1
                      do irc = 1,nr	! find velocity right above node
                         mm = (irc-1)*nx+mcx
                         if (ird(mm) == ir) then
                        if (irc < nd1 .or. irc > nd2) iflag = 0
                    end if
                      end do

                 if (idp == ir .and. iflag ==1) then

                        val = ex(kx)*elm*zscal/etp(j)


                        ynr(ip) = .true.
                        ynh(ip) = .true.
                        row(ip) = row(ip) + val
                        dws(ip) = dws(ip) + abs(val)
                 end if
                      end do

               end if
            end if
             end if
          end do	! loop over boundaries



           end if		! k > 1

        end do		! loop over ray points
        end if              ! if len(j) > 2

        iv = 0
        do ip=1,np
           if (ynr(ip)) then
          iv=iv+1
          inz(iv) = ip
          anz(iv) = row(ip)
           end if
        end do
        nrow(j) = iv
            nfz = nfz + iv

       if (mod(j,50) ==0 ) &
           !if (abs(etp(j)) < 0.01 .or. abs(dtj) > 1) &
          write(6,'(i6,i7,i10,i5,i8,4f9.3)')  j,inj(j),isk(j),fas(j),nrow(j),ttp(j),tca(j),etp(j),dtsum

            if (nrow(j) > 0) then
               dtsum = dtsum + dtj2
           ntsum = ntsum + 1
            endif

        write(7) (inz(i),i=1,iv)
        write(7) (anz(i),i=1,iv)


         end if

     m = m + len(j)

      end do			! loop over picks

      close(7)
      deallocate(ynr,row)



      rmst = sqrt(dtsum/ntsum)
      print*,'Root-mean-square misfit (sec) before inversion ',rmst
      print*,'Of ',npx,' picks, ',ntsum,'are used .'



   end subroutine frechet




   subroutine out_dws

      character(50) :: dwsv_file,dwsz_file
      integer :: ix,iz
      real :: x,z
      real :: za,zb,zr
      integer :: i,ip
      integer :: ir,ira
      real :: val

      read(5,'(a)') dwsv_file
      open(11,file=dwsv_file)
      do ix = 1,nxc
         x = x1 + (ix-1)*dxc
     z = z1
         do iz = 1,nzc

        i = (ix-1)*nzc+iz
        ira = 0
        do ir=1,nr+1
           if (ir == 1) then
              za = z1
           else
              call get_zrp(x,za,ir-1)
           endif

           if (ir == nr+1) then
              zb = z2
           else
              call get_zrp(x,zb,ir)
           end if

           if (z >= za .and. z <= zb) then
              ira = ir
          exit
           end if

        end do

        if (ira /= 0) then
           ip = (ira-1)*nxc*nzc + i
           if (vin(ip) < vpmin .or.dws(ip) > 0) then
              val = 0.0
           else
              val = -0.5
           end if
           write(11,'(2f10.4,f11.5,f9.3)') x,z,dws(ip),val
        end if

        if (iz < nzc) z=z+dzc(iz)
     end do
      end do
      close(11)


      read(5,'(a)') dwsz_file
      open(12,file=dwsz_file)

      do ir=1,nr
         do ix=1,nxc
        x= x1 + (ix-1)*dxc
        i = (ir-1)*nxc+ix
        ip = (nr+1)*nxc*nzc+i
        call get_zrp(x,zr,ir)
        write(12,'(2f10.4,f11.5)') x,zr,dws(ip)
     end do
       !  write(12,'(a1)') ">"
      end do
      close(12)

   end subroutine out_dws




   subroutine model_coverage

      integer :: ir
      integer :: ix,iz,i
      integer :: jx,jz,j
      integer :: im,jm,ip
      real :: xi,zi,xj,zj
      real :: xij,zij,rij
      integer :: isol


      allocate(yng(np))
      do ip=1,np
         yng(ip) = .false.
      end do

      do ir=nc1,nc2	! velocity layers

     do ix = 1,nxc
        xi = x1+(ix-1)*dxc
        zi = z1

        do iz = 1,nzc
           i = (ix-1)*nzc+iz
           im = (ir-1)*nxc*nzc + i

           do jx = 1,nxc
              xj = x1+(jx-1)*dxc
          xij = (xj-xi)/xreg

              zj = z1
              do jz = 1,nzc

             j = (jx-1)*nzc+jz
             jm = (ir-1)*nxc*nzc + j

             if (ynh(jm)) then
                zij = (zj-zi)/zreg
                rij = sqrt(xij*xij+zij*zij)
                if (rij <= 1.0) yng(im) = .true.
                 end if

             if (jz < nzc) zj=zj+dzc(jz)
          end do	! loop over jz
           end do		! loop over jx

           if (iz < nzc) zi=zi+dzc(iz)
        end do		! loop over iz
     end do			! loop over ix

      end do			! loop over ir


      do ir = nd1,nd2	! reflecting boundaries

         do ix = 1,nxc
        xi = x1+(ix-1)*dxc
        i = (ir-1)*nxc+ix
        im = (nr+1)*nxc*nzc+i
        if (ijn(i) == ir) then	! only independent boundary nodes

           do jx = 1,nxc

              j = (ir-1)*nxc+jx
              jm = (nr+1)*nxc*nzc+j

              if (ynh(jm)) then
                 xj = x1+(jx-1)*dxc
                 xij = abs(xj-xi)/xreg
             if (xij <= 1.0) then
                yng(im) = .true.
             end if
              end if

           end do		! loop over jx
        end if

     end do			! loop over ix

      end do			! loop over ir

      deallocate(ynh)

      nsol = 0
      do ip = 1,np
         if (yng(ip)) nsol=nsol+1
      end do

      allocate(ind(nsol),ivc(np))
      isol = 0
      do ip = 1,np
         if (yng(ip)) then
        isol=isol+1
        ind(isol) = ip
        ivc(ip) = isol
     else
        ivc(ip) = 0
     end if
      end do

   end subroutine model_coverage


   subroutine out_mask

      character(50) :: mask_file
      integer :: ir
      integer :: i,im
      integer :: ix,iz
      real :: xi,zi

      read(5,'(a)') mask_file
      open(17,file=mask_file)

      do ir=nc1,nc2	! velocity layers

     do ix = 1,nxc
        xi = x1+(ix-1)*dxc
        zi = z1

        do iz = 1,nzc
           i = (ix-1)*nzc+iz
           im = (ir-1)*nxc*nzc + i
           if (yng(im)) write(17,'(2f9.3)') xi,zi

           if (iz < nzc) zi=zi+dzc(iz)
        end do		! loop over iz
     end do			! loop over ix

      end do

      close(17)

   end subroutine out_mask


   subroutine size_matrix


      integer :: ix,iz,ir
      integer :: i,ia,ib,ic
      integer :: ip,ipa,ipb,ipc
      integer :: j,je
      integer :: nn
      integer :: m

      integer :: iia,iib,iic
      integer :: ica,icb,icc
      integer :: ina,inb,inc
      integer :: flag,flag1,flag2,flag3

      integer :: inm


      je = 0
      do j = 1,npx
         if (nrow(j) > 0) je = je + 1
      end do
      ndat = je

      nnz = nfz
      neq = ndat

      do ir = nc1,nc2	! velocity layers

     do ix=1,nxc			! damping velocities
        do iz = 1,nzc
           i = (ix-1)*nzc + iz
           ip = (ir-1)*nxc*nzc+i
           if (yng(ip)) then
              nnz=nnz+1
          neq=neq+1
           end if
        end do
     end do

     do ix=1,nxc			! flattening velocities	z
        do iz = 1,nzc-1
           ia = (ix-1)*nzc + iz
           ib = (ix-1)*nzc + iz+1
           ipa = (ir-1)*nxc*nzc+ia
           ipb = (ir-1)*nxc*nzc+ib
           if (yng(ipa) .and. yng(ipb)) then
              nnz=nnz+2
          neq=neq+1
           end if
        end do
     end do

     do ix=1,nxc-1			! flattening velocities	x
        do iz = 1,nzc
           ia = (ix-1)*nzc + iz
           ib = (ix  )*nzc + iz
           ipa = (ir-1)*nxc*nzc+ia
           ipb = (ir-1)*nxc*nzc+ib
           if (yng(ipa) .and. yng(ipb)) then
              nnz=nnz+2
          neq=neq+1
           end if
        end do
     end do

     do ix=1,nxc			! smoothing velocities z
        do iz = 2,nzc-1
           ia = (ix-1)*nzc + iz-1
           ib = (ix-1)*nzc + iz
           ic = (ix-1)*nzc + iz+1
           ipa = (ir-1)*nxc*nzc+ia
           ipb = (ir-1)*nxc*nzc+ib
           ipc = (ir-1)*nxc*nzc+ic
           if (yng(ipa) .and. yng(ipb) .and. yng(ipc)) then
              nnz=nnz+3
          neq=neq+1
           end if
        end do
     end do

     do ix=2,nxc-1			! smoothing velocities x
        do iz = 1,nzc
           ia = (ix-2)*nzc + iz
           ib = (ix-1)*nzc + iz
           ic = (ix  )*nzc + iz
           ipa = (ir-1)*nxc*nzc+ia
           ipb = (ir-1)*nxc*nzc+ib
           ipc = (ir-1)*nxc*nzc+ic
           if (yng(ipa) .and. yng(ipb) .and. yng(ipc)) then
              nnz=nnz+3
          neq=neq+1
           end if
        end do
     end do

      end do

      do ir = nd1,nd2	! model boundaries

     do ix=1,nxc			! damping boundaries
        i = (ir-1)*nxc+ix
        ip = (nr+1)*nxc*nzc + i
        if (yng(ip)) then
           nnz=nnz+1
           neq=neq+1
        end if
     end do

     do ix=1,nxc-1			! flattening boundaries	x
        ia = (ir-1)*nxc+ix
        ib = (ir-1)*nxc+ix+1

        iia = iin(ia)
        iib = iin(ib)

        ica = (iia-1)*nxc+ix
        ina = (nr+1)*nxc*nzc + ica

        icb = (iib-1)*nxc+ix+1
        inb = (nr+1)*nxc*nzc + icb

        ipa = (nr+1)*nxc*nzc + ia
        ipb = (nr+1)*nxc*nzc + ib


        flag = 0
        if (yng(ipa) .or. yng(ipb)) flag = 1

        flag1 = 0
        if (yng(ipa) .or. ijn(ia) > 0) flag1 = 1

        flag2 = 0
        if (yng(ipb) .or. ijn(ib) > 0) flag2 = 1


        if (flag == 1 .and. flag1 == 1 .and. flag2 == 1) then
           neq=neq+1
           if (yng(ina))  nnz=nnz+1
           if (yng(inb))  nnz=nnz+1
        end if

     end do

     do ix=2,nxc-1			! smoothing boundaries	x
        ia = (ir-1)*nxc+ix-1
        ib = (ir-1)*nxc+ix
        ic = (ir-1)*nxc+ix+1

        iia = iin(ia)
        iib = iin(ib)
        iic = iin(ic)

        ipa = (nr+1)*nxc*nzc + ia
        ipb = (nr+1)*nxc*nzc + ib
        ipc = (nr+1)*nxc*nzc + ic

        ica = (iia-1)*nxc+ix
        ina = (nr+1)*nxc*nzc + ica

        icb = (iib-1)*nxc+ix+1
        inb = (nr+1)*nxc*nzc + icb

        icc = (iib-1)*nxc+ix+1
        inc = (nr+1)*nxc*nzc + icb



        flag = 0
        if (yng(ipa) .or. yng(ipb) .or. yng(ipc)) flag = 1

        flag1 = 0
        if (yng(ipa) .or. ijn(ia) == -1) flag1 = 1

        flag2 = 0
        if (yng(ipb) .or. ijn(ib) == -1) flag2 = 1

        flag3 = 0
        if (yng(ipc) .or. ijn(ic) == -1) flag3 = 1

        if (flag == 1 .and. flag1 == 1 .and. flag2 ==1 .and. flag3 ==1) then
           neq=neq+1
           if (yng(ina))  nnz=nnz+1
           if (yng(inb))  nnz=nnz+1
           if (yng(inc))  nnz=nnz+1
        end if
     end do

      end do

      allocate(imz(nnz),amz(nnz),nro(neq),rhs(neq))

      je = 0
      do j=1,npx
         if (nrow(j) > 0) then
        je = je + 1
            nro(je) = nrow(j)
        rhs(je) = dat(j)
     end if
      end do

      open(unit=7,file=matfile,access='stream')
      m = 0
      je = 0
      do j=1,ndat
         nn = nro(j)

         read(7) (inz(i),i=1,nn)
     read(7) (anz(i),i=1,nn)

     do i = 1,nn
        inm = inz(i)
        if (inm >= 1 .and. inm <= np .and. m+i >= 1 .and. m+i <= nnz) then
           imz(m+i) = ivc(inz(i))
           amz(m+i) = anz(i)
        else
           print*,'bug in size matrix ',inm,np,m+1,nnz
           stop
        end if
     end do

     m = m + nn
      end do
      close(7)



      deallocate(nrow,dat)
      deallocate(inz,anz)

      print*
      print*,'Inversion has ',nsol,' active model parameters, and ',neq,' equations.'
      print*



   end subroutine size_matrix




   subroutine regularization

      integer :: ix,iz,ir
      integer :: i,ia,ib,ic
      integer :: ip,ipa,ipb,ipc
      integer :: j
      integer :: nn
      integer :: m
      integer :: knz
      integer :: jeq
      real :: dzci
      real :: vol
      real :: dvel,dref
      real :: gfac

      integer :: iia,iib,iic
      integer :: ica,icb,icc
      integer :: ina,inb,inc

      integer :: flag1,flag2,flag3,flag

      print*
      print*,'Constructing the regularization matrix with ',neq-ndat ,' equations. '


      jeq = ndat
      knz = nfz

      do ir = nc1,nc2	! velocity layers

     do ix=1,nxc			! damping velocities
        do iz = 1,nzc
           if (iz < nzc) then
              dzci = dzc(iz)
           else
              dzci = dzc(iz-1)
           end if
           vol = dxc*dzci

           i = (ix-1)*nzc + iz
           ip = (ir-1)*nxc*nzc+i
           if (yng(ip)) then
              knz=knz+1
          imz(knz) = ivc(ip)
          amz(knz) = reg0*vol
          jeq=jeq+1
          rhs(jeq) = 0.0
          nro(jeq) = 1

           end if
        end do
     end do

     do ix=1,nxc			! flattening velocities	z
        do iz = 1,nzc-1
           if (iz < nzc) then
              dzci = dzc(iz)
           else
              dzci = dzc(iz-1)
           end if
           vol = dxc*dzci
           gfac = reg1*vol/dzci

           ia = (ix-1)*nzc + iz
           ib = (ix-1)*nzc + iz+1
           ipa = (ir-1)*nxc*nzc+ia
           ipb = (ir-1)*nxc*nzc+ib
           dvel = (vin(ipb)-vin(ipa))/vscal

           if (yng(ipa) .and. yng(ipb)) then
              knz=knz+1
          imz(knz) = ivc(ipa)
          amz(knz) = gfac
          knz=knz+1
          imz(knz) = ivc(ipb)
          amz(knz) = -gfac
          jeq=jeq+1
          rhs(jeq) = gfac*dvel
          nro(jeq) = 2
           end if
        end do
     end do

     do ix=1,nxc-1			! flattening velocities	x
        do iz = 1,nzc
           if (iz < nzc) then
              dzci = dzc(iz)
           else
              dzci = dzc(iz-1)
           end if
           vol = dxc*dzci
           gfac = reg1*vol*asr/dxc

           ia = (ix-1)*nzc + iz
           ib = (ix  )*nzc + iz
           ipa = (ir-1)*nxc*nzc+ia
           ipb = (ir-1)*nxc*nzc+ib
           dvel = (vin(ipb)-vin(ipa))/vscal

           if (yng(ipa) .and. yng(ipb)) then
              knz=knz+1
          imz(knz) = ivc(ipa)
          amz(knz) = gfac
          knz=knz+1
          imz(knz) = ivc(ipb)
          amz(knz) = -gfac
          jeq=jeq+1
          rhs(jeq) = gfac*dvel
          nro(jeq) = 2
           end if
        end do
     end do

     do ix=1,nxc			! smoothing velocities z
        do iz = 2,nzc-1
           if (iz < nzc) then
              dzci = dzc(iz)
           else
              dzci = dzc(iz-1)
           end if
           vol = dxc*dzci
           gfac = reg2*vol/(dzci*dzci)

           ia = (ix-1)*nzc + iz-1
           ib = (ix-1)*nzc + iz
           ic = (ix-1)*nzc + iz+1
           ipa = (ir-1)*nxc*nzc+ia
           ipb = (ir-1)*nxc*nzc+ib
           ipc = (ir-1)*nxc*nzc+ic
           dvel = (vin(ipa)+vin(ipc)-2*vin(ipb))/vscal

           if (yng(ipa) .and. yng(ipb) .and. yng(ipc)) then
              knz=knz+1
          imz(knz) = ivc(ipa)
          amz(knz) = -gfac
          knz=knz+1
          imz(knz) = ivc(ipb)
          amz(knz) = 2*gfac
          knz=knz+1
          imz(knz) = ivc(ipc)
          amz(knz) = -gfac
          jeq=jeq+1
          rhs(jeq) = gfac*dvel
          nro(jeq) = 3

           end if
        end do
     end do

     do ix=2,nxc-1			! smoothing velocities x
        do iz = 1,nzc
           if (iz < nzc) then
              dzci = dzc(iz)
           else
              dzci = dzc(iz-1)
           end if
           vol = dxc*dzci
           gfac = reg2*vol*asr*asr/(dxc*dxc)

           ia = (ix-2)*nzc + iz
           ib = (ix-1)*nzc + iz
           ic = (ix  )*nzc + iz
           ipa = (ir-1)*nxc*nzc+ia
           ipb = (ir-1)*nxc*nzc+ib
           ipc = (ir-1)*nxc*nzc+ic
           dvel = (vin(ipa)+vin(ipc)-2*vin(ipb))/vscal

           if (yng(ipa) .and. yng(ipb) .and. yng(ipc)) then
              knz=knz+1
          imz(knz) = ivc(ipa)
          amz(knz) = -gfac
          knz=knz+1
          imz(knz) = ivc(ipb)
          amz(knz) = 2*gfac
          knz=knz+1
          imz(knz) = ivc(ipc)
          amz(knz) = -gfac
          jeq=jeq+1
          rhs(jeq) = gfac*dvel
          nro(jeq) = 3
           end if
        end do
     end do

      end do


      do ir = nd1,nd2	! model boundaries

         do ix=1,nxc			! damping boundary depth perturbations

        i = (ir-1)*nxc+ix
        ip = (nr+1)*nxc*nzc + i
        if (yng(ip)) then
           knz=knz+1
           imz(knz) = ivc(ip)
           amz(knz) = reg0*crf*dxc
           jeq=jeq+1
           rhs(jeq) = 0.0
           nro(jeq) = 1
        end if
     end do

     do ix=1,nxc-1			! flattening boundaries	x

        ia = (ir-1)*nxc+ix
        ib = (ir-1)*nxc+ix+1

        iia = iin(ia)
        iib = iin(ib)

        ica = (iia-1)*nxc+ix
        ina = (nr+1)*nxc*nzc + ica

        icb = (iib-1)*nxc+ix+1
        inb = (nr+1)*nxc*nzc + icb

        ipa = (nr+1)*nxc*nzc + ia
        ipb = (nr+1)*nxc*nzc + ib

        dref = (rin(icb)-rin(ica))/zscal
        gfac = reg1*crf*dxc*asr/dxc

        flag = 0
        if (yng(ipa) .or. yng(ipb)) flag = 1

        flag1 = 0
        if (yng(ipa) .or. ijn(ia) == -1) flag1 = 1

        flag2 = 0
        if (yng(ipb) .or. ijn(ib) == -1) flag2 = 1


        if (flag == 1 .and. flag1 ==1 .and. flag2 ==1) then

           jeq=jeq+1
           rhs(jeq) = gfac*dref
           nro(jeq) = 0

           if (yng(ina)) then
              knz=knz+1
          nro(jeq) = nro(jeq) + 1
          imz(knz) = ivc(ina)
          amz(knz) = gfac
           else
              rhs(jeq) = rhs(jeq) - gfac
           end if

           if (yng(inb)) then
              knz=knz+1
              nro(jeq) = nro(jeq) + 1
              imz(knz) = ivc(inb)
          amz(knz) = -gfac
           else
          rhs(jeq) = rhs(jeq) + gfac
           end if

        end if
     end do

     do ix=2,nxc-1			! smoothing boundaries	x

        ia = (ir-1)*nxc+ix-1
        ib = (ir-1)*nxc+ix
        ic = (ir-1)*nxc+ix+1

        iia = iin(ia)
        iib = iin(ib)
        iic = iin(ic)

        ica = (iia-1)*nxc+ix-1
        icb = (iia-1)*nxc+ix
        icc = (iib-1)*nxc+ix+1

        ipa = (nr+1)*nxc*nzc + ia
        ipb = (nr+1)*nxc*nzc + ib
        ipc = (nr+1)*nxc*nzc + ic

        ina = (nr+1)*nxc*nzc + ica
        inb = (nr+1)*nxc*nzc + icb
        inc = (nr+1)*nxc*nzc + icc


        flag = 0
        if (yng(ipa) .or. yng(ipb) .or. yng(ipc)) flag = 1

        flag1 = 0
        if (yng(ipa) .or. ijn(ia) == -1) flag1 = 1

        flag2 = 0
        if (yng(ipb) .or. ijn(ib) == -1) flag2 = 1

        flag3 = 0
        if (yng(ipc) .or. ijn(ic) == -1) flag3 = 1

        gfac = reg2*crf*dxc*asr*asr/(dxc*dxc)
        dref = (rin(ica)+rin(icc)-2*rin(icb))/zscal

        if (flag == 1 .and. flag1 == 1 .and. flag2 == 1 .and. flag3 ==1) then

           jeq=jeq+1
           rhs(jeq) = gfac*dref
           nro(jeq) = 0

           if (yng(ina)) then
              knz=knz+1
          nro(jeq) = nro(jeq) + 1
          imz(knz) = ivc(ina)
          amz(knz) = -gfac
           else
              rhs(jeq) = rhs(jeq) + gfac
           end if

           if (yng(inb)) then
              knz=knz+1
              nro(jeq) = nro(jeq) + 1
              imz(knz) = ivc(inb)
          amz(knz) = 2*gfac
           else
          rhs(jeq) = rhs(jeq) - 2*gfac
           end if

           if (yng(inc)) then
              knz=knz+1
              nro(jeq) = nro(jeq) + 1
              imz(knz) = ivc(inc)
          amz(knz) = -gfac
           else
          rhs(jeq) = rhs(jeq) + gfac
           end if

        end if
     end do

      end do

   end subroutine regularization





   subroutine init_inverse

      integer :: i,m
      integer :: j
      real :: chi2

      allocate(rhc(neq),anz(nnz),rms(neq))
      allocate(wvec(nsol),vvec(nsol))
      allocate(dm(nsol))


      chi2 = 0.0
      do j=1,ndat
         chi2 = chi2 + rhs(j)*rhs(j)
      end do
      chi2 = chi2/ndat

      print*
      print*,'starting with chi2 = ',chi2
      print*

      do m=1,nnz
         anz(m) = amz(m)
      end do

      do j = 1,neq
     rhc(j) = rhs(j)
      end do

      do i=1,nsol
         dm(i) = 0.0
     wvec(i) = 0.0
     vvec(i) = 0.0
      end do

      rega = 1.0		! set the first Tikonov parameter

   end subroutine init_inverse




   subroutine solve_eqn(reg,ch2)

      real :: reg,ch2

      call init_eqn(reg)

      call lsqr(neq,nsol,nnz,ilm,imz,nro,dax,rlm,dm,amz,rhs,vvec,wvec,echo)

      call misfit
      vard = (vars-varn)/vars
      ch2 = varn/ndat

     ! print*,'nominal ch2 after inversion ',ch2,'. Variance reduction ',vard


   end subroutine solve_eqn



   subroutine init_eqn(reg)

      real :: reg
      integer :: i,m,j,k

      m = 0
      vars = 0.0
      do j=1,neq

         if (j <= ndat) then
        rhs(j) = rhc(j)

        vars=vars+rhs(j)*rhs(j)
        do k=1,nro(j)
           amz(m+k) = anz(m+k)
        end do

     else
        do k=1,nro(j)
           amz(m+k) = anz(m+k)*reg
        end do
        rhs(j) = rhc(j)*reg
         end if
     m=m+nro(j)
      end do

      do i=1,nsol
         dm(i) = 0.0
      end do

   end subroutine init_eqn




   subroutine misfit

      real :: sums
      integer :: j

      integer :: naj,jj,lnz,i

      do j=1,neq
         rms(j) = -rhc(j)
      end do

      call avu(neq,nsol,nnz,imz,nro,anz,rms,dm)

      varn = 0.0
      do j=1,ndat
         varn = varn + rms(j)*rms(j)
      end do

     ! print*,'varn ',varn

      pen = 0.0
      do j = ndat+1,neq
         pen = pen + rms(j)*rms(j)
      end do

     ! print*,'pen ',pen



   end subroutine misfit




   subroutine bracket

      real, parameter :: sfac = 4.0
      integer, parameter :: nlim = 10
      integer :: nit
      real :: fac
      real :: swap

      print*,'First guess of lambda ',rega,' gives a chi2 of ', ch2a


      if (ch2a > ch2n) then
         regb = rega/sfac
      else
         regb = rega*sfac
      end if
      fac = 1.0
      ch2b = ch2a

      nit = 0
      do while (fac > 0.0 .and. nit < nlim)

         nit=nit+1
     call solve_eqn(regb,ch2b)

         icur = 1		! current model: 1=b, 0=a
         fac = (ch2n-ch2b)*(ch2n-ch2a)
         if (fac > 0 .and. nit < nlim) then
            rega = regb
            ch2a = ch2b
            icur = 0
            if (ch2a > ch2n) then
               regb = rega/sfac
            else
               regb=rega*sfac
            end if
         end if
      end do

!     Arrange bracket such that ch2a<ch2n<ch2b

      if (ch2a > ch2b) then
         swap=ch2b
         ch2b=ch2a
         ch2a=swap
         swap = regb
         regb = rega
         rega = swap
         icur = 0
      end if

      print*,'finished bracketing '
      Print*,'Lambda between ',rega,' and ',regb,' gives a chi2 between ',ch2a,' and ',ch2b

   end subroutine bracket





   subroutine secant

      integer, parameter :: nlim = 9
      real :: reg
      real :: ch2
      real :: cloa,clob
      real :: clox
      real :: sfac,cfac
      real, parameter :: close = 0.01      ! tolerance for fitting chisqn
      real, parameter :: tiny = 0.001
      integer :: nit

      cloa = abs(ch2n-ch2a)/ch2n
      clob = abs(ch2n-ch2b)/ch2n
      if (cloa < clob .and. cloa < close) THEN
         if (icur == 1) then
        call solve_eqn(rega,ch2a)
        icur=0
         end if
      else if (clob < close) then
         if (icur == 0) then
        call solve_eqn(regb,ch2b)
            icur=1
         end if
      end if

      print*,'Use the secant method to tune regularization parameter lambda :'

      if (ch2n < ch2a .or. ch2n > ch2b) then
         print*,'ch2n not bracketed in secant'
     write(6,'(3f11.4)') ch2a,ch2n,ch2b
     stop
      end if

      nit = 0
      do while (clob > close .and. cloa > close .and. nit < nlim)
         nit = nit+1
         sfac = abs(regb-rega)/(rega+regb)
         if (sfac > tiny) then
            reg = rega + (ch2n-ch2a)*(regb-rega)/(ch2b-ch2a)
            if (reg > regb .OR. reg < rega) then
               print*,'out of bracket in secant '
               reg = 0.5*(rega+regb)
            end if
            if (ch2b < ch2a) then
               print*,'trade-off has wrong sign in secant '
               print*,rega,ch2a
               print*,regb,ch2b
               stop
            end if

        call solve_eqn(reg,ch2)

            ! write(6,'(3f11.4)') rega,reg,regb

            clox = abs(ch2n-ch2)/ch2n

            if (clox > close) then

               if (ch2 > ch2a .and. ch2 < ch2b) then

              if (ch2n > ch2) then
                 rega = reg
                 ch2a = ch2
                 icur = 0
                 cloa = clox
              else
                 regb = reg
                 ch2b = ch2
                 icur = 1
                 clob = clox
              end if
          print*,'narrowing bracket, chi2 between ',ch2a,' and ',ch2b
               else
          print*,'problem in secant'
          print*,'iteration ',nit
          write(6,'(3f11.4)') rega,reg,regb
          write(6,'(3f11.4)') ch2a,ch2,ch2b
              cfac = (ch2-ch2b)*(ch2-ch2a)
          if (cfac < 0) then
             print*,'chi squared is bracketed'
          else
             print*,'chi squared is not bracketed'
          end if
          stop
               end if
            else

               print*,'data fit within tolerance for lambda = ',reg,' chi2 = ',ch2
               rega = reg
               regb = reg
               ch2a = ch2
               ch2b = ch2
               nit = nlim
            end if
         else

            nit = nlim
            print*,' Quitting on very steep chi squared curve in secant'

        reg = 0.5*(rega+regb)
            CALL solve_eqn(reg,ch2)
        rega = reg
            regb = reg
            ch2a = ch2
            ch2b = ch2
            nit = nlim
        icur = 2		! possibly a bad solution

     end if

      end do

   end subroutine secant





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



   subroutine get_veli(xi,zi,irf,vr)

      integer :: irf
      real :: xi,zi
      real :: vr

      integer :: ix,iz
      integer :: kx,kz
      integer :: mx,mz
      integer :: i,m
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
        i = (irf-1)*nx*nz+m
        w = ex(kx)*ez(kz)
        velm = vel(i)
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

   end subroutine get_veli





   SUBROUTINE lsqr(ndata,nnode,nsparse,itnlim0,ja,na,damp,redlim,x,ra,u,v,w,echo)

!		*************** DECLARATIONS ***********************

!		   ********* formal arguments **************

      INTEGER :: ndata      ! # equations
      INTEGER :: nnode      ! # unknowns
      INTEGER :: nsparse	! # non-zero elements
      INTEGER :: ja(*)      ! column index element ra(lnz), dim(nsparse)
      INTEGER :: na(*)      ! # of non-zero elements in row, dim(ndata)
      INTEGER :: itnlim0	! max # iterations performed
      REAL :: ra(*)		! non-zero matrix elements, dim(nsparse)
      REAL :: x(*)		! solution vector, dim(nnode)
      REAL :: damp		! damping
      REAL :: redlim		! reduction limit (scaled to 1.?)
      REAL :: u(*)		! rhs vector, dim(ndata)
      REAL :: v(nnode),w(nnode)

!           ************ local variables *************

      REAL, PARAMETER :: zero = 0.0
      REAL, PARAMETER :: one = 1.0
      INTEGER :: nout,istop,nhalf
      REAL :: otest1,xfit,xsfit

      INTEGER :: itnlim

      REAL :: atol,btol,conlim,anorm,acond,rnorm,arnorm,xnorm
!
      INTEGER :: i,itn,mod,nstop



      REAL :: alfa,bbnorm,beta,bnorm,cs,cs1,cs2,ctol,dampsq,ddnorm
      REAL :: delta,gamma,gambar,phi,phibar,psi,res1,res2,rho
      REAL :: rhobar,rhbar1,rhbar2,rhs,rtol,sn,sn1,sn2,t,tau,test1
      REAL :: test2,test3,theta,t1,t2,t3,xxnorm,z,zbar

      LOGICAL :: ECHO

!	******************* START SUBROUTINE ***********************


      atol=zero
      conlim = zero

      itnlim = itnlim0
      nout = 6
      nhalf = nnode/2+1

      IF (echo) WRITE(nout,20) ndata,nnode,damp,atol,conlim,btol,itnlim

      ctol=zero
      if (conlim > 6.e-7) ctol=one/conlim
      dampsq=damp*damp
      anorm=zero
      acond=zero
      arnorm=zero
      bbnorm=zero
      btol=zero
      conlim=zero
      ddnorm=zero
      res2=zero
      rnorm=zero
      xnorm=zero
      xxnorm=zero
      cs2=-one
      sn2=zero
      z=zero
      itn=0
      istop=0
      nstop=0
!
      DO i=1,nnode
         v(i) = zero
         x(i) = zero
      END DO

!             normalize u
      CALL normu(ndata,u,beta)

      CALL atuv(ndata,nnode,nsparse,ja,na,ra,u,v)

      CALL norml(nnode,v,alfa)

      CALL xscopy(1,1,nnode,v,w)

      rhobar = alfa
      phibar = beta
      bnorm = beta
      rnorm = beta
      arnorm = alfa*beta
!
      IF (arnorm > zero) THEN

         IF (echo) WRITE(nout,30)
         test1 = one
         otest1 = one
         test2 = alfa/beta
         xfit = one
         xsfit  =one

         IF(echo) WRITE(nout,50) &
             itn,xnorm,rnorm,test1,test2,anorm,acond
         IF(echo) WRITE(nout,60)

      END IF

      DO WHILE (istop == 0)

!            timing of lsqr?
!        dtime=time-tmold

         itn=itn+1
!
!            scale u
         CALL scale(ndata,u,-alfa)


         CALL avu(ndata,nnode,nsparse,ja,na,ra,u,v)


!         normalize u
         CALL normu(ndata,u,beta)
         bbnorm=bbnorm+alfa*alfa+beta*beta+dampsq
         CALL xsscal(1,nnode,-beta,v)

         CALL atuv(ndata,nnode,nsparse,ja,na,ra,u,v)


         CALL norml(nnode,v,alfa)
!
         rhbar2 = rhobar*rhobar+dampsq
         rhbar1 = SQRT(rhbar2)
         cs1 = rhobar/rhbar1
         sn1 = damp/rhbar1
         psi = sn1*phibar
         phibar = cs1*phibar
!
         rho = SQRT(rhbar2+beta*beta)
         cs = rhbar1/rho
         sn = beta/rho
         theta = sn*alfa
         rhobar = -cs*alfa
         phi = cs*phibar
         phibar = sn*phibar
         tau = sn*phi
!
         t1 = phi/rho
         t2 = -theta/rho
         t3 = one/rho
!
         DO i = 1,nnode
            t=w(i)
            x(i)=t1*t+x(i)
            w(i)=t2*t+v(i)
            t=(t3*t)
            t=t*t
            ddnorm=t+ddnorm
         END DO
!
!
         delta = sn2*rho
         gambar = -cs2*rho
         rhs = phi-delta*z
         zbar = rhs/gambar
         xnorm = SQRT(xxnorm+zbar*zbar)
         gamma = SQRT(gambar*gambar+theta*theta)
         cs2 = gambar/gamma
         sn2 = theta/gamma
         z = rhs/gamma
         xxnorm = xxnorm+z*z
!
         anorm = SQRT(bbnorm)
         acond = anorm*SQRT(ddnorm)
         res1 = phibar*phibar
         res2 = res2+psi*psi
         rnorm = SQRT(res1+res2)
         arnorm = alfa*ABs(tau)
!
         test1 = rnorm/bnorm
         test2 = arnorm/(anorm*rnorm)
         test3 = one/acond
         t1 = test1/(one+anorm*xnorm/bnorm)
         rtol = btol+atol*anorm*xnorm/bnorm
!
         t3 = one+test3
         t2 = one+test2
         t1 = one+t1

         IF (itn >= itnlim) istop = 7
         IF (t3 <= one) istop = 6
         IF (t2 <= one) istop = 5
         IF (t1 <= one) istop = 4
!
!        IF (test3 <= ctol) istop = 3
         IF (test2 <= atol) istop = 2
         IF (test1 <= rtol) istop = 1
!
         IF (itn <= 10 .OR. itn >= itnlim-10 &
                .OR. MOD(itn,10) == 0) THEN

            IF(echo) write(nout,50) &
                  itn,xnorm,rnorm,test1,test2,anorm,acond
         END IF

         IF (mod(itn,10) == 0) THEN
            if (echo) write(nout,60)
         END IF

         IF (echo .and. otest1-test1 < redlim) THEN
            write(6,*) 'lsqr: stopped otest1,test1, redlim', &
                     otest1,test1,redlim
            istop = 8
         END IF

         otest1=test1

      END DO
!
      t = one
      if (ndata > nnode) t=one*(ndata-nnode)
      if (dampsq > zero) t=one*ndata
      t = rnorm/SQRT(t)

      if (echo) write(nout,70) itn,istop

      itnlim = itn


!
      RETURN
!
20    format(//, &
      ' lsqr -- least squares solution of a*x=b',//, &
      10x,'the matrix has',i8,' rows and',i8,' columns',/,10x, &
      'the damping parameter is', e10.2,/,10x,'atol=',e10.2,'  conlim=',e10.2,'  btol=',e10.2,'  itnlim=',i10)
30    format(///,' itn',4x,'xnorm     residue   res(norm.)',3x,'incomp    anorm     acond ',/)
50    format(1x,i4,2e11.3,e14.6,e11.3,3e10.3)
60    format(1x)
70    format(/,' no. of iterations=',i6,8x,'istop=',i3)
666   format(///' lsqr stopped because the next iteration would hit the time limit'///)

   END SUBROUTINE lsqr







   SUBROUTINE avu(neq,nnode,nsparse,ja,na,ra,u,v)
!
!     calculates the product :
!		c(neq)=a(neq,nnode)*v(nnode)+u(neq)
!     matrix a is in 'lsqr-format'


      INTEGER :: nnode,neq,nsparse
      INTEGER :: ja(nsparse),na(neq)
      REAL :: ra(nsparse),v(nnode),u(neq)

      INTEGER :: jj
      REAL, PARAMETER :: zero=0.0
      INTEGER :: i,j,naj,lnz
      REAL :: sums

!	*****************************************************

      lnz=0
      DO j=1,neq
         sums=zero
         naj=na(j)
         DO jj=1,naj
            lnz=lnz+1
            i=ja(lnz)
            sums=sums+ra(lnz)*v(i)
         END DO
         u(j)=u(j)+sums
      END DO

   END SUBROUTINE avu



   SUBROUTINE atuv(ndata,nnode,nsparse,ja,na,ra,u,v)
!
!   calculates the product :
!             v(nnode)=a(ndata,nnode)^t*u(ndata)+v(nnode)

      INTEGER :: ndata,nnode,nsparse
      INTEGER :: ja(nsparse),na(ndata)
      REAL :: ra(nsparse),v(nnode),u(ndata)

      INTEGER :: j,lnz,naj,jj,i

!	***********************************************

      lnz=0
      DO j=1,ndata
         naj=na(j)
         DO jj=1,naj
            lnz=lnz+1
            i=ja(lnz)
        v(i)=v(i)+ra(lnz)*u(j)

     END DO
      END DO

   END SUBROUTINE atuv



   SUBROUTINE scale(ndata,u,fac)
      INTEGER :: ndata
      REAL :: u(*)
      REAL :: fac
      INTEGER :: i

      DO  i=1,ndata
         u(i)=u(i)*fac
      END DO
   END SUBROUTINE scale




   SUBROUTINE xscopy(i1,i2,n,a,b)

      INTEGER :: i1,i2
      INTEGER :: n
      REAL :: a(*)
      REAL :: b(*)

      INTEGER :: i

      DO i = i1,n
         b(i+i2-1)=a(i)
      END DO

   END SUBROUTINE xscopy





   FUNCTION xsdot(i1,i2,n,a,b) RESULT (sum)

      REAL, PARAMETER :: zero=0.
      INTEGER :: i1,i2
      INTEGER :: n
      REAL :: a(*)
      REAL :: b(*)
      INTEGER :: i

      REAL :: sum

      sum=zero
      DO  i=i1,n
         sum=sum+a(i)*b(i+i2-1)
      END DO

   END FUNCTION xsdot





   FUNCTION xsnrm2(i1,n,a) result (sum)

!            2-norm of a(i1->n)
      REAL, PARAMETER :: zero=0.0
      INTEGER :: i1
      INTEGER :: n
      REAL :: a(*)
      REAL :: sum
      INTEGER :: i

      sum = zero
      DO i = i1,n
     sum=sum+a(i)*a(i)
      END DO
      sum = SQRT(sum)

   END FUNCTION xsnrm2


   SUBROUTINE xsscal(i1,n,fac,a)
!             scale a(i1->n) with fac

      INTEGER :: i1
      INTEGER :: n
      REAL :: a(*)
      REAL :: fac
      INTEGER :: i

      DO i=i1,n
         a(i)=a(i)*fac
      END DO

   END SUBROUTINE xsscal



   SUBROUTINE normu(ndata,u,beta)
!             scales vector u (virtual on tape9) to length 1

      INTEGER :: ndata
      REAL :: u(*)
      REAL :: beta
      REAL, parameter :: one = 1.0

      beta = xsdot(1,1,ndata,u,u)
      beta = sqrt(beta)
      CALL scale(ndata,u,one/beta)

   END SUBROUTINE normu



   SUBROUTINE norml(n,x,fac)

      REAL, PARAMETER :: zero=0.
      REAL, PARAMETER :: one=1.
      INTEGER :: n
      REAL :: x(n)
      REAL :: fac

      fac = xsnrm2(1,n,x)
      IF (fac > zero) CALL xsscal(1,n,one/fac,x)

   END SUBROUTINE norml




   subroutine new_sol


      integer :: ix,iz
      integer :: jx,jz
      real :: xj,zj
      integer :: ixa,ixb
      real :: xa,xb
      real :: za,zb
      real :: xi,zi
      real :: rij,xij,zij
      real :: rmin
      integer :: ir,ik,ii,in
      integer :: i,mi,im
      integer :: j,mj,jm,ja,jb
      real :: dvel,dzrf

      do ir = nc1,nc2	! velocity layers with an update
         do ix = 1,nxc
        do iz = 1,nzc
           i = (ix-1)*nzc+iz
           mi = (ir-1)*nxc*nzc + i
           im = ivc(mi)
           if (im /= 0) then
              dvel = vscal*dm(im)
              vin(mi) = vin(mi) + dvel
           end if
        end do
     end do

      end do

      do ir = nc1,nc2	! fill in velocity coverage
         do ix = 1,nxc
        xi = x1 + (ix-1)*dxc
        zi = z1
        do iz = 1,nzc
           i = (ix-1)*nzc+iz
           mi = (ir-1)*nxc*nzc + i
           im = ivc(mi)

           if (im == 0) then
              rmin = sqrt((x2-x1)*(x2-x1)+asr*asr*(z2-z1)*(z2-z1))
              do jx = 1,nxc
             xj = x1 + (jx-1)*dxc
             xij = xi-xj
             zj = z1
             do jz = 1,nzc
                zij = zi-zj
            rij = sqrt(xij*xij + asr*asr*zij*zij)

            j = (jx-1)*nzc+jz
                    mj = (ir-1)*nxc*nzc + j
                    jm = ivc(mj)
            if (jm /= 0 .and. rij < rmin) then
               rmin = rij
               vin(mi) = vin(mj)
            end if
            if (jz < nzc) zj = zj + dzc(jz)
             end do
          end do
           end if

           if (iz < nzc) zi = zi + dzc(iz)
        end do			! loop over iz

     end do				! loop over ix

      end do

      ! link dependent reflector nodes to the inversion result:

      do ir= nd1,nd2
         do ix=1,nxc
        i = (ir-1)*nxc+ix
        ik = (nr+1)*nxc*nzc+i
        ii = iin(i)
        if (ii /= ir) then
           in = (ii-1)*nxc+ix
           mi = (nr+1)*nxc*nzc+in
           im = ivc(mi)
           if (im /= 0) then
              ivc(ik) = im
           end if
        end if
     end do

      end do


      do ir = nd1,nd2	! update reflector depths
         do ix = 1,nxc
        xi = x1 + (ix-1)*dxc
        i = (ir-1)*nxc+ix
        mi = (nr+1)*nxc*nzc+i
        im = ivc(mi)
        if (im /= 0) then

           dzrf = zscal*dm(im)

           rin(i) = rin(i) + dzrf

        end if
         end do
      end do





      do ir = nd1,nd2	! fill in reflector coverage

     do ix = 1,nxc
        xi = x1 + (ix-1)*dxc
        i = (ir-1)*nxc+ix
        mi = (nr+1)*nxc*nzc+i
        im = ivc(mi)
        if (im == 0) then
           xa = xi
           ixa = ix
           ja = im
           do while (ixa > 1 .and. ja==0)
              ixa = ixa -1
          xa = xa - dxc
              j = (ir-1)*nxc+ixa
              mj = (nr+1)*nxc*nzc+j
              ja = ivc(mj)
           end do
           if (ja /= 0) za = rin(j)

           xb = xi
           ixb = ix
           jb = im
           do while (ixb < nxc .and. jb==0)
              ixb = ixb +1
          xb = xb + dxc
              j = (ir-1)*nxc+ixb
              mj = (nr+1)*nxc*nzc+j
              jb = ivc(mj)
           end do
           if (jb /= 0) zb = rin(j)

           if (ja /= 0 .and. jb /= 0 .and. ixa /= ixb) then
              zi = ((xi-xa)*zb+(xb-xi)*za)/(xb-xa)
              rin(i) = zi
           else if (ja == 0 .and. jb /= 0) then
              rin(i) = zb
           else if (ja /= 0 .and. jb == 0) then
              rin(i) = za
           end if

        end if
         end do
      end do

   end subroutine new_sol




   subroutine update_model

      integer :: ix,iz
      real :: xi,zi
      integer :: jx,jz
      integer :: kx,kz
      integer :: mx,mz
      integer :: m
      real :: xj,zj
      integer :: i,j,ii,ij
      integer :: im,jm
      integer :: ir,jr
      real :: dzn
      real :: ex(0:1),ez(0:1)
      real :: wgt
      real :: samp

      do ir = nc1,nc2

     do ix=1,nx
        xi = x1 + (ix-1)*dx
        jx = int((xi-x1)/dxc)+1
        if (jx < 1) jx = 1
        if (jx > nxc-1) jx = nxc-1
        xj = x1+(jx-1)*dxc
        ex(1) = (xi-xj)/dxc
        ex(0) = 1.0 - ex(1)

        do iz=1,nz
           zi = z1 + (iz-1)*dz
           i = (ix-1)*nz+iz
           im = (ir-1)*nx*nz+i

           zj = z1
           jz = 1
           dzn = dzc(jz)
           do while (zj+dzn < zi)
              zj = zj + dzn
              jz=jz+1
          if (jz < nzc) dzn = dzc(jz)
           end do
           if (jz > nzc-1) jz = nzc-1
           ez(1) = (zi-zj)/dzn
           ez(0) = 1.0 - ez(1)

           samp = 0.0
           do kx = 0,1
              mx = jx+kx
              do kz = 0,1
             mz=jz+kz
             m = (mx-1)*nzc+mz
             jm = (ir-1)*nxc*nzc + m
             wgt = ex(kx)*ez(kz)
             samp = samp + wgt*vin(jm)
          end do
           end do
               vel(im) = samp
        end do		! loop over ix
         end do			! loop over iz

      end do			! loop over ir

      do ir = nd1,nd2
         do ix = 1,nx
        xi = x1 + (ix-1)*dx
        jx = int((xi-x1)/dxc)+1
        if (jx < 1) jx = 1
        if (jx > nxc-1) jx = nxc-1
        xj = x1+(jx-1)*dxc
        ex(1) = (xi-xj)/dxc
        ex(0) = 1.0 - ex(1)

        i = (ir-1)*nx+ix
        samp = 0.0
        do kx=0,1
           mx = jx+kx
           jm = (ir-1)*nxc+mx
           wgt = ex(kx)
           samp = samp + wgt*rin(jm)
            end do

        zrf(i) = samp

     end do
      end do

      do ix = 1,nx
         do ir=nr,2,-1
        ii = (ir-1)*nx+ix
        do jr= 1,ir-1
           ij = (jr-1)*nx+ix
           if (zrf(ii) < zrf(ij)) zrf(ii) =  zrf(ij)
        end do
     end do
      end do

   end subroutine update_model




   subroutine write_vm

      character(80) :: vmfile
      integer :: j

      read(5,'(a)') vmfile


     !

      open(unit=13,file=vmfile,access='stream')
      write(13) nx,nz,nr
      write(13) x1,x2
      write(13) z1,z2
      write(13) (zrf(j),j=1,nb)
      write(13) (ird(j),j=1,nb)
      write(13) (vel(j),j=1,nm)
      close(13)

   end subroutine write_vm



end program vm_tomo
