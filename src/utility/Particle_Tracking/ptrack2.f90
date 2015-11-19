!   Copyright 2014 College of William and Mary
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!										
!	SCHISM Particle tracking code for data format 5 (SZ). 				
!	Routines brought from SELFE/SCHISM:						
!	cpp, quicksearch, intersect2, signa, area_coord, levels			
!       Warning: indices in 2D arrays are not swapped (i.e., (np,nv)).	
!                 Also interpolation is along S-coord in tracking.		
!       No quads yet.
!										
!	Input: hgrid.gr3, vgrid.in, particle.bp, *hvel.64, *vert.63, *elev.61	
!	Input particle.bp:							
!	  (1) description;						
!	  (2) nscreen;								
!	  (3) ibf: forward (=1) or backward (=-1) tracking.			
!         (4) istiff: stiff (fixed distance frm f.s.; 1) or not (0).		
!	  (5) ics,slam0,sfea0: coordinate system and center for CPP projection
!             (same as in param.in);
!	  (6) h0,rnday,dtm,nspool,ihfskip,ndeltp: min. depth, # of days, 	
!	      time step used in the original run, nspool and ihfskip used 	
!             in the run (see param.in); # of sub-division in the tracking;
!             !!!WARNING: ihfskip must be equal to the actual # of steps 	
!                         contained in a file, in the case there is only 1 file.
!             Also: rnday may be smaller than the original run.
!	  (7) nparticle: # of particles;					
!	  (8) idp(i),st_p(i),xpar(i),ypar(i),zpar0(i): particle id, start time (s),
!		starting x,y, and z relative to the instant f.s. (<=0).		
!										
!	Output: particle.pth, fort.16 (screen), fort.11 (fatal errors).		
!										
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ifort -Bstatic -O2 -assume byterecl -o ptrack2 ptrack2.f90 ../../Utility/UtilLib/compute_     zcor.f90  ../../Utility/UtilLib/schism_geometry.f90

!   pgf90 -O2 -mcmodel=medium  -Bstatic -o ptrack2 ptrack2.f90  ../../Utility/UtilLib/compute_     zcor.f90  ../../Utility/UtilLib/schism_geometry.f90


!...  Data type consts
      module kind_par
        implicit none
        integer, parameter :: sng_kind1=4
        integer, parameter :: dbl_kind1=8
        real(kind=dbl_kind1), parameter :: small1=1.e-6 !small non-negative number; must be identical to that in global
      end module kind_par


!...  definition of variables
!...
!
!************************************************************************
!     			mnp < mne < mns					*
!************************************************************************
!
      module global
        implicit none
        integer, parameter :: sng_kind=4
        integer, parameter :: dbl_kind=8

!...  	Dimensioning parameters
!       integer, parameter :: mnp=130000
!       integer, parameter :: mne=260000
!       integer, parameter :: mns=400000
!       integer, parameter :: mnv=31
!       integer, parameter :: mnei=40 !neighbor
        integer, parameter :: nbyte=4
!       integer, parameter :: mnout=100 !max. # of output files
!      	integer, parameter :: mnope=6 !# of open bnd segements
!       	integer, parameter :: mnond=20000 !max. # of open-bnd nodes on each segment
!      	integer, parameter :: mnland=50 !# of land bnd segements
!      	integer, parameter :: mnlnd=10000 !max. # of land nodes on each segment
!      	integer, parameter :: mnoe=20000 !max. # of open-bnd elements on each segment
!      	integer, parameter :: mnosd=20000 !max # of open-bnd sides on each segment
!      	integer, parameter :: mnbfr=9 !# of forcing freqs.
!      	integer, parameter :: itmax=5000 !# of iteration for itpack solvers used for dimensioning
!      	integer, parameter :: nwksp=6*mne+4*itmax !available work space for itpack solvers
!      	integer, parameter :: mirec=1109000000 !997000000) !max. record # to prevent output ~> 4GB
      	real(kind=dbl_kind), parameter :: zero=1.e-5 !small postive number in lieu of 0; usually used to check areas 
        real(kind=dbl_kind), parameter :: small1=1.e-6 !small non-negative number; must be identical to that in kind_par
        real(kind=dbl_kind), parameter :: pi=3.1415926d0 

!...  	Important variables
      	integer :: np,ne,ns,nvrt,mnei,ibf,istiff,ivcor,kz,nsig
      	real(kind=dbl_kind) :: h0,rho0,dt
        real,save :: h_c,theta_b,theta_f,h_s !s_con1

!...    Output handles
        character(len=48) :: start_time,version
        character(len=12) :: ifile_char
        integer,save :: ihot,nrec,nspool,igmp,noutgm,ifile,noutput,ifort12(100)
!        character(len=48), dimension(mnout) :: outfile,variable_nm,variable_dim
!        integer :: ihot,nrec,nspool,igmp,noutgm,ifile,noutput,ifort12(100)
!        integer, dimension(mnout) :: ichan,irec,iof
!        real(kind=dbl_kind), dimension(mnout) :: vpos
        
        integer,save, allocatable :: nne(:),indel(:,:),idry(:),idry_e(:),idry_e0(:)
        integer,save, allocatable :: kbp(:),kbs(:),kbe(:),kbp00(:),isbnd(:)

        real(kind=dbl_kind),save, allocatable :: x(:),y(:),dp(:),eta1(:),eta2(:),eta3(:),hmod(:)
        real(kind=dbl_kind),save, allocatable :: area(:),xctr(:),yctr(:)
        real(kind=dbl_kind),save, allocatable :: snx(:),sny(:),distj(:),dps(:),dldxy(:,:,:)

        real(kind=dbl_kind),save, allocatable :: zpar0(:)
        !For interface with util routines
        real,save, allocatable :: ztot(:),sigma(:),xcj(:),ycj(:),sigma_lcl(:,:)

        integer,save, allocatable :: i34(:),elnode(:,:),ic3(:,:),elside(:,:),isdel(:,:),isidenode(:,:)
        integer,save, allocatable :: icum1(:,:),icum2(:,:,:)
        integer,save :: nxq(3,4,4),nodel(3)

        real(kind=dbl_kind),save, allocatable :: z(:,:)
        real(kind=dbl_kind),save, allocatable :: uu1(:,:),vv1(:,:),ww1(:,:),uu2(:,:),vv2(:,:),ww2(:,:)


!...  1D arrays
!       integer :: nne(mnp),nnp(mnp),idry(mnp),idry_e(mne),idry_e0(mne),iback(mnp)
!       integer :: kbp(mnp),kbs(mns),kbe(mne),kbp00(mnp)

!       real(kind=dbl_kind), dimension(mnp) :: x,y,dp,eta1,eta2,eta3,hmod(mnp)
!       real(kind=dbl_kind), dimension(mne) :: area,radiel,xctr,yctr
!       real(kind=dbl_kind), dimension(mns) :: snx,sny,distj,xcj,ycj,dps

!       real(kind=dbl_kind) :: sigma(mnv),cs(mnv),dcs(mnv),ztot(mnv)
!       real(kind=dbl_kind) :: zpar0(10*mnp)

!...  2D arrays
!       integer :: elnode(3,mne),ic3(3,mne),elside(3,mne)
!       integer :: nx(3,2),indel(mnei,mnp),isdel(2,mns), &
!          &isidenode(2,mns),iself(mnei,mnp)

!       real(kind=dbl_kind) :: ssign(3,mne),z(mnp,mnv),icum1(mnp,mnv),icum2(mnp,mnv,2)

!       real(kind=dbl_kind), dimension(mnp,mnv) :: uu1,vv1,ww1, uu2,vv2,ww2

      end module global

!...  Main program
      program ptrack
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)

      real(kind=sng_kind) :: floatout,floatout2
      real*8, allocatable :: xpar(:),ypar(:),zpar(:),st_p(:),upar(:),vpar(:),wpar(:)
      real*8, allocatable :: ztmp(:),ztmp2(:),dhfx(:),dhfy(:),dhfz(:),grdx(:),grdy(:), &
     &grdz(:),amas(:),wndx(:),wndy(:)
      integer, allocatable :: ielpar(:),levpar(:),iabnorm(:),ist(:),inbr(:)
      character(len=25), allocatable :: idp(:)
      real*8 :: vxl(4,2),vyl(4,2),vzl(4,2),vxn(4),vyn(4),vzn(4),arco(3), &
     &dx(10),dy(10),dz(10),val(4,2),vbl(4,2),vcl(4,2),vdl(4,2),van(4),vcn(4),vdn(4),vwx(4),vwy(4)
      real, allocatable :: zlcl(:)
      integer :: nodel2(3)

!     real(kind=dbl_kind), dimension(10*mnp) :: xpar,ypar,zpar,st_p,upar,vpar,wpar
!     integer, dimension(10*mnp) :: ielpar,levpar,iabnorm
!     character(len=25) :: idp(10*mnp)
!     real(kind=dbl_kind) :: vxl(4,2),vyl(4,2),vzl(4,2),vxn(4),vyn(4),vzn(4),staint(4)
!     dimension arco(3),ztmp(mnv),ztmp2(mnv)


      !Cyclical index
      do k=3,4 !elem. type
        do i=1,k  !local index
          do j=1,k-1 !offset
            nxq(j,i,k)=i+j
            if(nxq(j,i,k)>k) nxq(j,i,k)=nxq(j,i,k)-k
            if(nxq(j,i,k)<1.or.nxq(j,i,k)>k) then
              write(*,*)'nx wrong',i,j,k,nxq(j,i,k)
              stop
            endif
          enddo !j
        enddo !i
      enddo !k

      ifort12=0

!...  Read in particles
      open(95,file='particle.bp',status='old')
      read(95,*)
!      read(95,'(a48)') data_format
      read(95,*) nscreen
      read(95,*) ibf
      if(iabs(ibf)/=1) then
        write(*,*)'Wrong ibf',ibf
        stop
      endif
      read(95,*) istiff !1: fixed distance from F.S.
      if(istiff/=0.and.istiff/=1) then
        write(*,*)'Wrong istiff',istiff
        stop
      endif
      read(95,*) ics,slam0,sfea0
      slam0=slam0/180*pi
      sfea0=sfea0/180*pi
      read(95,*) h0,rnday,dtm,nspool,ihfskip,ndeltp !# of sub-divisions
      if(mod(ihfskip,nspool).ne.0) then
        write(*,*)'ihfskip must be a multiple of nspool'
        stop
      endif
      nrec=ihfskip/nspool !# of records (steps) per stack

      read(95,*) nparticle
      allocate(zpar0(nparticle),xpar(nparticle),ypar(nparticle),zpar(nparticle),zpar0(nparticle), &
     &st_p(nparticle),idp(nparticle),ielpar(nparticle),levpar(nparticle),upar(nparticle), &
     &vpar(nparticle),wpar(nparticle),iabnorm(nparticle),ist(nparticle),inbr(nparticle), &
     &dhfx(nparticle),dhfy(nparticle),dhfz(nparticle),grdx(nparticle),grdy(nparticle), &
     &grdz(nparticle),amas(nparticle),wndx(nparticle),wndy(nparticle),stat=istat)
      if(istat/=0) stop 'Failed to alloc (1)'

      levpar=-99
      iabnorm=0 !abnormal tracking exit flag

      dt=dtm*nspool !output time step
      st_m=(ibf+1)/2*rnday*86400 !initialize min. or max. starting time
      do i=1,nparticle
        if(ics.eq.1) then
!	  zpar0: relative to f.s.
          read(95,*)idp(i),st_p(i),xpar(i),ypar(i),zpar0(i)
        else
          read(95,*)idp(i),st_p(i),xparl,yparl,zpar0(i)
          xparl=xparl/180*pi
          yparl=yparl/180*pi
          call cpp(xpar(i),ypar(i),xparl,yparl,slam0,sfea0)
        endif
        if(st_p(i)<0.or.st_p(i)>rnday*86400) then
          write(11,*)'Starting time for particle ',i,' out of range:',st_p(i)
          stop
        endif
        if(zpar0(i)>0) then
          write(11,*)'Starting z-coord above f.s.',i
          stop
        endif
        if(ibf*st_p(i)<ibf*st_m) st_m=st_p(i)
      enddo !i
      close(95)

!...  Read in h- and v-grid and compute geometry
!...
      if(ics.eq.1) then
        open(14,file='hgrid.gr3',status='old')
      else
        open(14,file='hgrid.ll', status='old')
      endif
      read(14,*) 
      read(14,*) ne,np
      allocate(x(np),y(np),dp(np),hmod(np),i34(ne),elnode(4,ne),area(ne), &
     &nne(np),idry(np),idry_e(ne),idry_e0(ne),kbp(np),kbe(ne),kbp00(np), &
     &eta1(np),eta2(np),eta3(np),xctr(ne),yctr(ne),isbnd(np), &
     &dldxy(3,2,ne),stat=istat)
      if(istat/=0) stop 'Failed to alloc (2)'


      do i=1,np
        if(ics.eq.1) then
          read(14,*) j,x(i),y(i),dp(i)
        else if(ics.eq.2) then
          read(14,*) j,xlon,ylat,dp(i)
          ylat=ylat/180*pi
          xlon=xlon/180*pi
          call cpp(x(i),y(i),xlon,ylat,slam0,sfea0)
        endif
        hmod(i)=min(dp(i),h_s)
      enddo !i=1,np

      do i=1,ne
        read(14,*) j,i34(i),elnode(1:i34(i),i)

        n1=elnode(1,i)
        n2=elnode(2,i)
        n3=elnode(3,i)
        area(i)=signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3))

        !Derivative of shape function for triangles only
        !For quds, use nodes 1-3
        do j=1,3
          id1=j+1
          id2=j+2
          if(id1>3) id1=id1-3
          if(id2>3) id2=id2-3
          dldxy(j,1,i)=(y(elnode(id1,i))-y(elnode(id2,i)))/2/area(i)
          dldxy(j,2,i)=(x(elnode(id2,i))-x(elnode(id1,i)))/2/area(i)
        enddo !j

        if(i34(i)==4) then
          n4=elnode(4,i)
          area(i)=area(i)+signa(x(n1),x(n3),x(n4),y(n1),y(n3),y(n4))
        endif

        if(area(i)<=0) then
          write(11,*)'Negative area at',i
          stop
        endif
      enddo !i=1,ne

      close(14)
!...  End fort.14

!     vgrid
      open(19,file='vgrid.in',status='old')
      read(19,*)ivcor
      read(19,*)nvrt
      rewind(19)
      allocate(sigma_lcl(nvrt,np),z(np,nvrt),icum1(np,nvrt),icum2(np,nvrt,2), &
     &uu1(np,nvrt),vv1(np,nvrt),ww1(np,nvrt),uu2(np,nvrt),vv2(np,nvrt),ww2(np,nvrt), &
     &ztmp(nvrt),ztmp2(nvrt),zlcl(nvrt),sigma(nvrt),ztot(nvrt),stat=istat)
      if(istat/=0) stop 'Failed to alloc (3)'
      call get_vgrid('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)
      !kbp has been assigned for ivcor=1
!                                                                             *
!                                                                             *
!******************************************************************************
!                                                                             *
!     			Compute geometry 				      *
!                                                                             *
!******************************************************************************
!                                                                             *
!                                                                             *
!     We also need elem ball
      nne=0
      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
        enddo
      enddo
      mnei=maxval(nne)

      allocate(indel(mnei,np),stat=istat)
      if(istat/=0) stop 'Failed to alloc. indel'
      nne=0
      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
          if(nne(nd)>mnei) then
            write(11,*)'Too many neighbors',nd
            stop
          endif
          indel(nne(nd),nd)=i
        enddo
      enddo !i

      call compute_nside(np,ne,i34,elnode,ns)
      if(nscreen.eq.1) write(*,*) 'There are',ns,' sides in the grid...'
      write(16,*) 'There are',ns,' sides in the grid...'

!     Allocate side-related arrays
      allocate(ic3(4,ne),elside(4,ne),isdel(2,ns),isidenode(2,ns),xcj(ns),ycj(ns), &
     &snx(ns),sny(ns),distj(ns),kbs(ns),dps(ns),stat=istat)
      if(istat/=0) stop 'Failed to alloc (4)'
!     Then compute the rest of side related arrays with additional
!     inputs (xnd,ynd) (x,y coordinates of each node)
      call schism_geometry(np,ne,ns,real(x),real(y),i34,elnode,ic3,elside,isdel,isidenode,xcj,ycj)

      !Remaining side arrays
      do i=1,ns
        nd1=isidenode(1,i)
        nd2=isidenode(2,i)
        dps(i)=(dp(nd1)+dp(nd2))/2
        distj(i)=dsqrt((x(nd2)-x(nd1))**2+(y(nd2)-y(nd1))**2)
        if(distj(i)==0) then
          write(11,*)'Zero side',i
          stop
        endif

        thetan=datan2(x(nd1)-x(nd2),y(nd2)-y(nd1))
        snx(i)=dcos(thetan)
        sny(i)=dsin(thetan)
      enddo !i

!...  compute centers of each triangle 
!...
      do i=1,ne
        xctr(i)=0
        yctr(i)=0
        do j=1,i34(i)
          xctr(i)=xctr(i)+x(elnode(j,i))/i34(i)
          yctr(i)=yctr(i)+y(elnode(j,i))/i34(i)
        enddo !j
      enddo !i=1,ne

      if(nscreen.eq.1) write(*,*)'done computing geometry...'
      write(16,*)'done computing geometry...'

!...  Read in header 
      nt=rnday*86400/dt !total # of iterations
      ifile=1 !for st_m=0
      do i=1,nt/nrec+1
        if(st_m/dt>(i-1)*nrec.and.st_m/dt<=i*nrec) then
          ifile=i
          exit
        endif
      enddo
      if(ibf==1) then
        iths=(ifile-1)*nrec+1 !iteration # for the 1st time step output in ifile
      else !ibf=-1
        iths=ifile*nrec !iteration # for the last time step output in ifile
      endif
      write(ifile_char,'(i12)') ifile
      open(60,file=trim(ifile_char)//'_elev.61',access='direct',recl=nbyte)
      open(61,file=trim(ifile_char)//'_hvel.64',access='direct',recl=nbyte)
      open(62,file=trim(ifile_char)//'_vert.63',access='direct',recl=nbyte)

      irec00=5*48/nbyte+5+7+nvrt+2 !elev
!     Read initial bottom index
      do m=1,np
        read(60,rec=irec00+4)kbp00(m)
        irec00=irec00+4
      enddo !m=1,np
      do m=1,ne
        irec00=irec00+1
        do mm=1,i34(m)
          irec00=irec00+1
        enddo !mm
      enddo !m

!      irec00=irec00+4*ne

!     Initialize kbp for levels()
!     kbp=kbp00

!...  Compute record # offset for a node and level for 3D outputs
!...
      icount1=0
      icount2=0
      do i=1,np
        do k=max0(1,kbp00(i)),nvrt
          do m=1,2
            icount2=icount2+1
            icum2(i,k,m)=icount2
          enddo !m
          icount1=icount1+1
          icum1(i,k)=icount1
        enddo !k
      enddo !i=1,np

      if(ibf==1) then
        irec01=irec00 !end of the header for elev.61
        irec02=irec00 !hvel.64
        irec03=irec00 !vert.63
      else !ibf=-1
        irec01=irec00+(2+2*np)*nrec+1 !last record in ifile plus 1 (for backward reading) 
        irec02=irec00+(2+np+icum2(np,nvrt,2))*nrec+1
        irec03=irec00+(2+np+icum1(np,nvrt))*nrec+1
      endif
      irec1=irec01
      irec2=irec02
      irec3=irec03

!...  Compute initial elements for particle tracking
!...
      lp1: do i=1,nparticle
        do k=1,ne
          call pt_in_poly2(i34,x(elnode(1:i34(k),k)),y(elnode(1:i34(k),k)),xpar(i),ypar(i),inside)
          if(inside/=0) then
            ielpar(i)=k
            cycle lp1
          endif
        enddo !k=1,ne
        write(11,*)'Cannot find init. element for particle',i
        stop
      end do lp1 !i=1,nparticle

      open(95,file='particle.pth',status='replace')
      write(95,*)'Drogues'
      if(ibf==1) then
        write(95,*) nt-iths+1
      else
        write(95,*) iths
      endif

      if(nscreen.eq.1) write(*,*)'done initialization...'
      write(16,*)'done initialization...'


!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!									!
!	       Time iteration						!
!									!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

      if(ibf==1) then
        it2=nt
      else
        it2=1
      endif

!     Initialize for output before particle moving
      upar=0; vpar=0; wpar=0
      zpar=zpar0

      do it=iths,it2,ibf
!--------------------------------------------------------------------------
      time=it*dt
      
!...  Read in elevation and vel. info
      if((ibf==1.and.it>nrec*ifile).or.(ibf==-1.and.it<=nrec*(ifile-1))) then
        ifile=ifile+ibf
        write(ifile_char,'(i12)') ifile
        close(60)
        close(61)
        close(62)
        open(60,file=trim(ifile_char)//'_elev.61',access='direct',recl=nbyte)
        open(61,file=trim(ifile_char)//'_hvel.64',access='direct',recl=nbyte)
        open(62,file=trim(ifile_char)//'_vert.63',access='direct',recl=nbyte)
        irec1=irec01
        irec2=irec02
        irec3=irec03
      endif 

      if(ibf==1) then
        irec1=irec1+np+2
        irec2=irec2+np+2
        irec3=irec3+np+2

        do i=1,np
          read(60,rec=irec1+1) floatout
          irec1=irec1+1
          eta2(i)=floatout

          do k=max0(1,kbp00(i)),nvrt
            read(61,rec=irec2+1) floatout
            read(61,rec=irec2+2) floatout2
            irec2=irec2+2
            uu2(i,k)=floatout
            vv2(i,k)=floatout2
            read(62,rec=irec3+1) floatout
            irec3=irec3+1
            ww2(i,k)=floatout
          enddo !k
        enddo !i 
      else !ibf=-1
        do i=np,1,-1
          read(60,rec=irec1-1) floatout
          irec1=irec1-1
          eta2(i)=floatout

          do k=nvrt,max0(1,kbp00(i)),-1
            read(61,rec=irec2-2) floatout
            read(61,rec=irec2-1) floatout2
            irec2=irec2-2
            uu2(i,k)=floatout
            vv2(i,k)=floatout2
            read(62,rec=irec3-1) floatout
            irec3=irec3-1
            ww2(i,k)=floatout
          enddo !k
        enddo !i 

        if(irec1-np-2<=0.or.irec2-np-2<=0.or.irec3-np-2<=0) then
          write(*,*)'Negative record #:',irec1-np-2,irec2-np-2,irec3-np-2
          stop
        endif
        irec1=irec1-np-2
        irec2=irec2-np-2
        irec3=irec3-np-2
      endif !ibf

!...  Store info for first step
      if(it==iths) then
        uu1=uu2; vv1=vv2; ww1=ww2; eta1=eta2
      endif 

!...  Compute elevation eta3
      do i=1,np
        if(ibf==1) then
          eta3(i)=eta1(i)
        else
          eta3(i)=eta2(i)
        endif
      enddo !i

!...  Compute z-cor
      call levels

      if(nscreen.eq.1) write(*,*)'begin ptrack...'
      write(16,*)'begin ptrack...'

!...  Particle tracking
      write(95,*) time,nparticle
      do i=1,nparticle
        eta_p=0; dp_p=0 !for output before moving
        if((ibf==1.and.time<=st_p(i)).or.(ibf==-1.and.time>st_p(i)-dt)) go to 449

        pt=dt !tracking time step
!...    Initialize starting level 
        if(levpar(i)==-99) then !just started
          pt=(time-st_p(i))*ibf
          if(pt<=0) then
            write(*,*)'Tracking step negative:',pt
            stop
          endif
          iel=ielpar(i)
          if(idry_e(iel)==1) then !dry
            levpar(i)=-1
          else !wet
!           call area_coord(iel,xpar(i),ypar(i),arco)
            call pt_in_poly3(i34(iel),x(elnode(1:i34(iel),iel)),y(elnode(1:i34(iel),iel)), &
     &xpar(i),ypar(i),arco,nodel)

            do k=kbe(iel),nvrt
              ztmp2(k)=0
              do j=1,3
                nd=elnode(j,iel)
                ztmp2(k)=ztmp2(k)+z(nd,max(k,kbp(nd)))*arco(j)
              enddo !j
            enddo !k
            zpar(i)=max(zpar0(i)+ztmp2(nvrt),ztmp2(kbe(iel))) !zpar0<=0
            jlev=0
            do k=kbe(iel),nvrt-1
              if(zpar(i)>=ztmp2(k).and.zpar(i)<=ztmp2(k+1)) then
                jlev=k+1
                exit
              endif
            enddo !k
            if(jlev==0) then
              write(11,*)'Cannot find an init. level:',i,zpar(i),(ztmp2(k),k=kbe(iel),nvrt)
              stop
            endif
            levpar(i)=jlev
          
            upar(i)=0; vpar(i)=0; wpar(i)=0
            do j=1,3
              nd=elnode(j,iel)
              upar(i)=upar(i)+uu2(nd,jlev)*arco(j)
              vpar(i)=vpar(i)+vv2(nd,jlev)*arco(j)
              wpar(i)=wpar(i)+ww2(nd,jlev)*arco(j)
            enddo !j

          endif !wet
        endif !levpar=-99

!	Wetting/drying
        if(idry_e(ielpar(i))==1) then
          levpar(i)=-1
          go to 449
        endif

        nnel=ielpar(i)
        jlev=levpar(i)
!       Rewetted elements
        if(jlev==-1) then !element nnel wet
          jlev=nvrt
          zpar(i)=sum(eta3(elnode(1:i34(nnel),nnel)))/i34(nnel)
        endif
  
!	Tracking
        x0=xpar(i)
        y0=ypar(i)
        z0=zpar(i)
        nnel0=nnel
        jlev0=jlev
        dtb=pt/ndeltp
        do idt=1,ndeltp
          if(ibf==1) then
            trat=real(idt)/ndeltp
          else
            trat=real(idt-1)/ndeltp
          endif
          xt=x0+ibf*dtb*upar(i)
          yt=y0+ibf*dtb*vpar(i)
          zt=z0+ibf*dtb*wpar(i)
          call quicksearch(1,idt,i,nnel0,jlev0,dtb,x0,y0,z0,xt,yt,zt,nnel,jlev, &
     &nodel2,arco,zrat,nfl,eta_p,dp_p,ztmp,kbpl,ist(i),inbr(i),rnds,pbeach)

!	  nnel not dry
!	  Interpolate in time
          do j=1,i34(nnel)
            nd=elnode(j,nnel)
            do l=1,2
              lev=jlev+l-2
              vxl(j,l)=uu1(nd,lev)*(1-trat)+uu2(nd,lev)*trat
              vyl(j,l)=vv1(nd,lev)*(1-trat)+vv2(nd,lev)*trat
              vzl(j,l)=ww1(nd,lev)*(1-trat)+ww2(nd,lev)*trat
            enddo !l
          enddo !j

!	  Interpolate in vertical 
          do j=1,3
            vxn(j)=vxl(j,2)*(1-zrat)+vxl(j,1)*zrat
            vyn(j)=vyl(j,2)*(1-zrat)+vyl(j,1)*zrat
            vzn(j)=vzl(j,2)*(1-zrat)+vzl(j,1)*zrat
          enddo !j

!	  Interpolate in horizontal
          upar(i)=0; vpar(i)=0; wpar(i)=0
          do j=1,3
            upar(i)=upar(i)+vxn(j)*arco(j)
            vpar(i)=vpar(i)+vyn(j)*arco(j)
            wpar(i)=wpar(i)+vzn(j)*arco(j)
          enddo !j

          if(nfl==1) then
            iabnorm(i)=1
            go to 404
          endif

          x0=xt
          y0=yt
          z0=zt
          nnel0=nnel
          jlev0=jlev
        enddo !idt=1,ndeltp

404     xpar(i)=xt
        ypar(i)=yt
        zpar(i)=zt
        ielpar(i)=nnel
        levpar(i)=jlev

449     continue

        if(ics==2) then
           call cppinverse(xout, yout, xpar(i), ypar(i), slam0, sfea0)
           xout = xout * 180.0 / pi
           yout = yout * 180.0 / pi
        else
           xout = xpar(i)
           yout = ypar(i)
        endif

        write(95,*)i,real(xpar(i)),real(ypar(i)),real(zpar(i)-eta_p) !drogue format for xmvis6s; no extra lines after this
!        write(95,*)idp(i),xpar(i),ypar(i)
!        write(95,*)zpar(i),ielpar(i),levpar(i),eta_p,dp_p,iabnorm(i)
!        write(95,*)upar(i),vpar(i),wpar(i),(eta3(elnode(l,ielpar(i))),l=1,3)
!        write(95,'(e12.4)')zpar(i)-eta_p

!!        write(95,'(2e14.4)')time,ztmp2(nvrt)
!!        write(*,'(2e14.4)')time,zpar(i)-eta3(ielpar(i))
      enddo !i=1,nparticle

!...  Store info for next step
      uu1=uu2; vv1=vv2; ww1=ww2; eta1=eta2

      if(nscreen.eq.1) write(*,*)'Time=',time
      write(16,*)'Time=',time

!--------------------------------------------------------------------------
      enddo !it

      if(nscreen.eq.1) write(*,*)'Completed'
      write(16,*)'Completed'

      stop
      end



!******************************************************************************
!                                                                             *
!    Transform from lon,lat (rlambda,phi) coordinates into CPP coordinates.     *
!    Lon,Lat must be in radians.                                              *
!                                                                             *
!******************************************************************************

      subroutine cpp(x,y,rlambda,phi,rlambda0,phi0)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z), integer(i-n)

      r=6378206.4
      x=r*(rlambda-rlambda0)*cos(phi0)
      y=phi*r

      return
      end

      subroutine cppinverse(rlambda,phi,x,y,rlambda0,phi0)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z), integer(i-n)

      r=6378206.4
      rlambda=x / (r * cos(phi0)) + rlambda0
      phi=y/r

      return
      end


!
!********************************************************************************
!										*
!     Program to detect if two segments (1,2) and (3,4) have common pts   	*
!     Assumption: the 4 pts are distinctive.					*
!     The eqs. for the 2 lines are: X=X1+(X2-X1)*tt1 and X=X3+(X4-X3)*tt2.	*
!     Output: iflag: 0: no intersection or colinear; 1: exactly 1 intersection.	*
!     If iflag=1, (xin,yin) is the intersection.				*
!										*
!********************************************************************************
!

      subroutine intersect2(x1,x2,x3,x4,y1,y2,y3,y4,iflag,xin,yin,tt1,tt2)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z), integer(i-n)
      real(kind=dbl_kind1), parameter :: zero1=0.0 !small positive number or 0

      real(kind=dbl_kind1), intent(in) :: x1,x2,x3,x4,y1,y2,y3,y4
      integer, intent(out) :: iflag
      real(kind=dbl_kind1), intent(out) :: xin,yin,tt1,tt2

      tt1=-1000
      tt2=-1000
      iflag=0
      delta=(x2-x1)*(y3-y4)-(y2-y1)*(x3-x4)
      delta1=(x3-x1)*(y3-y4)-(y3-y1)*(x3-x4)
      delta2=(x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)

      if(delta.ne.0.0d0) then
        tt1=delta1/delta
        tt2=delta2/delta
        if(tt1.ge.-zero1.and.tt1.le.1+zero1.and.tt2.ge.-zero1.and.tt2.le.1+zero1) then
          iflag=1
          xin=x1+(x2-x1)*tt1
          yin=y1+(y2-y1)*tt1
        endif
      endif

      return
      end


       function signa(x1,x2,x3,y1,y2,y3)
       use kind_par
       implicit real(kind=dbl_kind1)(a-h,o-z),integer(i-n)
 
       signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2
       
       return
       end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                       !
!       Compute area coordinates of pt (xt,yt), which must be inside element nnel.      !
!       Impose bounds for area coordinates.                                             !
!                                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine area_coord(nnel,xt,yt,arco)
      use global, only : dbl_kind,elnode,area,x,y
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)
      integer, intent(in) :: nnel
      real(kind=dbl_kind), intent(in) :: xt,yt
      real(kind=dbl_kind), intent(out) :: arco(3)

      n1=elnode(1,nnel)
      n2=elnode(2,nnel)
      n3=elnode(3,nnel)
      arco(1)=signa(xt,x(n2),x(n3),yt,y(n2),y(n3))/area(nnel)
      arco(2)=signa(x(n1),xt,x(n3),y(n1),yt,y(n3))/area(nnel)
      arco(1)=dmax1(0.0d0,dmin1(1.0d0,arco(1)))
      arco(2)=dmax1(0.0d0,dmin1(1.0d0,arco(2)))
      if(arco(1)+arco(2)>1) then
        arco(3)=0
        arco(1)=dmin1(1.d0,dmax1(0.d0,arco(1)))
        arco(2)=1-arco(1)
      else
        arco(3)=1-arco(1)-arco(2)
      endif

      return
      end


!
!********************************************************************
!								    *
!	Routine to update z-coordinates and wetting and drying      *
!								    *
!********************************************************************
!
      subroutine levels
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)

      dimension idry_new(np) !,out2(12+mnv)
      real :: zlcl(nvrt)

!...  z-coor. for nodes
!...  
      do i=1,np
        if(dp(i)+eta3(i)<=h0) then !dry
          idry_new(i)=1 
          if(dp(i)>=h_s) then
            write(11,*)'Deep depth dry:',i
            stop
          endif
          kbp(i)=0
        else !wet
          idry_new(i)=0

          if(ivcor==2) then
            call zcor_SZ(real(dp(i)),real(eta3(i)),real(h0),h_s,h_c,theta_b, &
     &theta_f,kz,nvrt,ztot,sigma,zlcl,idry_tmp,kbpl)
            if(idry_tmp==1) then
              write(11,*)'Impossible dry (7):',i,idry_tmp,dp(i),eta1(i),eta2(i),eta3(i),kbpl
              stop
            endif
            !Cannot use kbp00 b/cos wet/dry
            kbp(i)=kbpl
            z(i,kbpl:nvrt)=zlcl(kbp(i):nvrt)
          else if(ivcor==1) then
            z(i,kbp(i):nvrt)=(eta3(i)+dp(i))*sigma_lcl(kbp(i):nvrt,i)+eta3(i)
          else
            write(11,*)'Unknown ivcor:',ivcor
            stop
          endif

          if(1==2) then
!--------------------------------------------------------------

!         S-levels
          do k=kz,nvrt
            kin=k-kz+1

            if(hmod(i)<=h_c) then
              if(ifort12(12)==0) then
                ifort12(12)=1
                write(12,*)'Initial depth too shallow for S:',i,hmod(i),h_c
              endif
             !iback(i)=1
              z(i,k)=sigma(kin)*(hmod(i)+eta3(i))+eta3(i)
            else if(eta3(i)<=-h_c-(hmod(i)-h_c)*theta_f/s_con1) then !hmod(i)>h_c>=0
              write(11,*)'Pls choose a larger h_c (1):',eta3(i),h_c
              stop
            else
!             z(i,k)=eta3(i)*(1+sigma(kin))+h_c*sigma(kin)+(hmod(i)-h_c)*cs(kin)
            endif
          enddo !k=kz,nvrt

!         z-levels
          if(dp(i)<=h_s) then
            kbp(i)=kz
          else !bottom index will never change
            if(kbp(i)>=kz.or.kbp(i)<1) then
              write(11,*)'Impossible 92:',kbp(i),kz,i
              stop
            endif
            z(i,kbp(i))=-dp(i)
            do k=kbp(i)+1,kz-1
              z(i,k)=ztot(k)
            enddo !k
          endif
!--------------------------------------------------------------
          endif !1==2


          do k=kbp(i)+1,nvrt
            if(z(i,k)-z(i,k-1)<=0) then
              write(11,*)'Inverted z-levels at:',i,k,z(i,k)-z(i,k-1),eta3(i),hmod(i)
              stop
            endif
          enddo !k
        endif !wet ot dry
      enddo !i=1,np

!...  Set wet/dry flags for element; element is "dry" if one of nodes is dry; conversely, 
!...  an element is wet if all nodes are wet (and all sides are wet as well)
!...  Weed out fake we nodes; a node is wet if and only if at least one surrounding element is wet
!...
!      idry_e0=idry_e !save
      idry=1 !dry unless wet
      kbe=0
      do i=1,ne
       idry_e(i)=maxval(idry_new(elnode(1:i34(i),i)))
        if(idry_e(i)==0) then
          idry(elnode(1:i34(i),i))=0
          kbe(i)=minval(kbp(elnode(1:i34(i),i)))
        endif
      enddo !i

      return
      end

!
!********************************************************************************
!	
!     Straightline search algorithm. Initially nnel0 is an element that 	
!     encompasses (x0,y0). iloc=0: do not nudge initial pt; iloc=1: nudge.	 
!     Input: iloc,idt,ipar,nnel0,x0,y0,z0,xt,yt,zt,jlev0, time, and uu2,vv2,ww2 for 	
!	abnormal cases;								
!     Output the updated end pt (xt,yt,zt) (if so), nnel1, jlev1, area          
!       coordinates, vertical ratio, a flag nfl, and local elevation and depth. 
!     nfl=1 if a bnd or dry element is hit and vel. there is small,		 
!	or death trap is reached.						
!										
!********************************************************************************
!
      subroutine quicksearch(iloc,idt,ipar,nnel0,jlev0,time,x0,y0,z0,xt,yt,zt,nnel1,jlev1, &
     &nodel2,arco,zrat,nfl,etal,dp_p,ztmp,kbpl,ist2,inbr2,rnds,pbeach)

      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)

      integer, intent(in) :: iloc,idt,ipar,nnel0,jlev0
      real(kind=dbl_kind), intent(in) :: time,x0,y0,z0,rnds,pbeach
      integer, intent(out) :: nnel1,jlev1,nodel2(3),nfl,kbpl,ist2,inbr2
      real(kind=dbl_kind), intent(inout) :: xt,yt,zt
      real(kind=dbl_kind), intent(out) :: arco(3),zrat,etal,dp_p,ztmp(nvrt)

      !Local
      real :: zlcl(nvrt)
      logical ::  ltmp1,ltmp2

      if(iloc>1) then
        write(11,*)'iloc > 1'
        stop
      endif
      if(idry_e(nnel0)==1) then
        write(11,*)'Starting element is dry'
        stop
      endif

      nfl=0
      trm=time !time remaining

!     Starting element nel
      nel=nnel0

!     An interior pt close to (x0,y0) to prevent underflow for iloc >=1.
      if(iloc==0) then
        xcg=x0
        ycg=y0
      else if(iloc==1) then
        xcg=(1-1.0d-4)*x0+1.0d-4*xctr(nel)
        ycg=(1-1.0d-4)*y0+1.0d-4*yctr(nel)
      endif

      call pt_in_poly2(i34(nel),x(elnode(1:i34(nel),nel)),y(elnode(1:i34(nel),nel)),xcg,ycg,inside)
!      aa=0
!      aa1=0
!      do i=1,3
!        n1=elnode(i,nel)
!        n2=elnode(nxq(1,i,i34(nel)),nel)
!        aa=aa+dabs(signa(x(n1),x(n2),x0,y(n1),y(n2),y0))
!        aa1=aa1+dabs(signa(x(n1),x(n2),xt,y(n1),y(n2),yt))
!      enddo !i
!      ae=dabs(aa-area(nel))/area(nel)
!      if(ae>small1) then

      if(inside==0) then
        write(11,*)'(x0,y0) not in nnel0 initially',nnel0,xcg,ycg
        stop
      endif

      call pt_in_poly2(i34(nel),x(elnode(1:i34(nel),nel)),y(elnode(1:i34(nel),nel)),xt,yt,inside)
!      ae=dabs(aa1-area(nel))/area(nel)
!      if(ae<small1) then
      if(inside/=0) then
        nnel1=nel
        go to 400
      endif

!     (xt,yt) not in nel, and thus (xcg,ycg) and (xt,yt) are distinctive
      pathl=dsqrt((xt-xcg)**2+(yt-ycg)**2)
      if(pathl==0) then
        write(11,*)'Zero path',x0,y0,xt,yt,xcg,ycg
        stop
      endif

!     Starting edge nel_j
      nel_j=0
      do j=1,i34(nel)
        jd1=elnode(nxq(1,j,i34(nel)),nel)
        jd2=elnode(nxq(2,j,i34(nel)),nel)
        call intersect2(xcg,xt,x(jd1),x(jd2),ycg,yt,y(jd1),y(jd2),iflag,xin,yin,tt1,tt2)
        if(iflag==1) then
          nel_j=j
          exit
        endif
      enddo !j=1,3
      if(nel_j==0) then
        write(11,*)'Found no intersecting edges I:',nel,xcg,ycg,xt,yt,ae
        stop
      endif

      zin=z0 !intialize
      it=0
      loop4: do
!----------------------------------------------------------------------------------------
      it=it+1
      if(it>1000) then
        if(ifort12(3)==0) then
          ifort12(3)=1
          write(12,*)'Death trap reached',idt,id0
        endif
        nfl=1
        xt=xin
        yt=yin
        zt=zin
        nnel1=nel
        exit loop4
      endif
      md1=elnode(nxq(1,nel_j,i34(nel)),nel)
      md2=elnode(nxq(2,nel_j,i34(nel)),nel)
      
!     Compute z position 
      dist=sqrt((xin-xt)**2+(yin-yt)**2)
!      if(dist/pathl.gt.1+1.0d-4) then
!        write(11,*)'Path overshot'
!        stop
!      endif
      tmp=min(1.d0,dist/pathl)
      zin=zt-tmp*(zt-zin)
      trm=trm*dist/pathl !time remaining
      
      pathl=sqrt((xin-xt)**2+(yin-yt)**2)
      if(pathl==0.or.trm==0) then
        write(11,*)'Target reached'
        stop
      endif

      lit=0 !flag
!     For horizontal exit and dry elements, compute tangential vel.,
!     update target (xt,yt,zt) and continue.
      if(ic3(nel_j,nel)==0.or.idry_e(max(1,ic3(nel_j,nel)))==1) then
        lit=1
        isd=elside(nel_j,nel)
        if(isidenode(1,isd)+isidenode(2,isd)/=md1+md2) then
          write(11,*)'Wrong side'
          stop
        endif

        !Set status flag
        if(ic3(nel_j,nel)==0) then
          ltmp1=isbnd(md1)>0.or.isbnd(md2)>0 !open bnd
          ltmp2=isbnd(md1)==-1.and.isbnd(md2)==-1.and.rnds>=pbeach !permanently stranded @ land bnd
          if(ltmp1.or.ltmp2) then
            if(ltmp1) then
              ist2=2
            else
              ist2=-1 !permanently stranded
            endif
            nfl=1
            xt=xin
            yt=yin
            zt=zin
            nnel1=nel
            exit loop4
          endif
        else !internal side with a dry elem.
          if(rnds>=pbeach) then !% exceeded
            ist2=-2  !stranded onshore
            inbr2=ic3(nel_j,nel)
            nfl=1
            xt=xin
            yt=yin
            zt=zin
            nnel1=nel
            exit loop4
          endif
        endif !ic3

!       Nudge intersect (xin,yin), and update starting pt
        eps=1.e-2
        xin=(1-eps)*xin+eps*xctr(nel)
        yin=(1-eps)*yin+eps*yctr(nel)
        xcg=xin
        ycg=yin

        vtan=-(uu2(md1,jlev0)+uu2(md2,jlev0))/2*sny(isd)+(vv2(md1,jlev0)+vv2(md2,jlev0))/2*snx(isd)
        vtan=max(vtan,1.d-2) !to prevent getting stuck
        xvel=-vtan*sny(isd)
        yvel=vtan*snx(isd)
        zvel=(ww2(md1,jlev0)+ww2(md2,jlev0))/2
        xt=xin-xvel*trm
        yt=yin-yvel*trm
        zt=zin-zvel*trm
        hvel=dsqrt(xvel**2+yvel**2)
        if(hvel<1.e-4) then
          write(11,*)'Impossible (5):',hvel
!          nfl=1
!          xt=xin
!          yt=yin
!          zt=zin
!          nnel1=nel
!          exit loop4
        endif
        pathl=hvel*trm
      endif !abnormal cases

!     Search for nel's neighbor with edge nel_j, or in abnormal cases, the same element
      if(lit==0) nel=ic3(nel_j,nel) !next front element

!      aa=0
!      do i=1,3
!        k1=elnode(i,nel)
!        k2=elnode(nxq(1,i,i34(nel)),nel)
!        aa=aa+dabs(signa(x(k1),x(k2),xt,y(k1),y(k2),yt))
!      enddo !i
!      ae=dabs(aa-area(nel))/area(nel)
!      if(ae<small1) then

      call pt_in_poly2(i34(nel),x(elnode(1:i34(nel),nel)),y(elnode(1:i34(nel),nel)),xt,yt,inside)
      if(inside/=0) then
        nnel1=nel
        exit loop4
      endif

!     Next intersecting edge
      do j=1,i34(nel)
         jd1=elnode(nxq(1,j,i34(nel)),nel)
         jd2=elnode(nxq(2,j,i34(nel)),nel)
!        For abnormal case, same side (border side) cannot be hit again
         if(jd1==md1.and.jd2==md2.or.jd2==md1.and.jd1==md2) cycle
         call intersect2(xcg,xt,x(jd1),x(jd2),ycg,yt,y(jd1),y(jd2),iflag,xin,yin,tt1,tt2)
         if(iflag==1) then
           nel_j=j !next front edge          
           cycle loop4
         endif
      enddo !j
      write(11,*)'Failed to find next edge I:',lit,xin,yin,xt,yt,nel,md1,md2,idt,id0,ae
      stop
!----------------------------------------------------------------------------------------
      end do loop4 

400   continue
!     No vertical exit from domain
      if(idry_e(nnel1)==1) then
        write(11,*)'Ending element is dry'
        stop
      endif

!     Compute area & sigma coord.
      !call area_coord(nnel1,xt,yt,arco)
      call pt_in_poly3(i34(nnel1),x(elnode(1:i34(nnel1),nnel1)),y(elnode(1:i34(nnel1),nnel1)), &
     &xt,yt,arco,nodel2)
      n1=elnode(nodel2(1),nnel1)
      n2=elnode(nodel2(2),nnel1)
      n3=elnode(nodel2(3),nnel1)
      etal=eta3(n1)*arco(1)+eta3(n2)*arco(2)+eta3(n3)*arco(3)
      dep=dp(n1)*arco(1)+dp(n2)*arco(2)+dp(n3)*arco(3)
      dp_p=dep
      if(etal+dep<h0) then
        write(11,*)'Weird wet element in quicksearch:',nnel1,eta3(n1),eta3(n2),eta3(n3)
        stop
      endif

      if(istiff==1) zt=etal+zpar0(ipar)

!     Compute z-levels
      if(ivcor==2) then
        call zcor_SZ(real(dep),real(etal),real(h0),h_s,h_c,theta_b, &
     &theta_f,kz,nvrt,ztot,sigma,zlcl,idry_tmp,kbpl)
        ztmp(kbpl:nvrt)=zlcl(kbpl:nvrt)
      else if(ivcor==1) then
        kbpl=nvrt+1 !local bottom index (maybe degerate)
        do j=1,3
          nd=elnode(nodel2(j),nnel1)
          if(kbp(nd)<kbpl) kbpl=kbp(nd)
        enddo !j
        ztmp(kbpl)=-dep !(etal+dep)*sigma_lcl(kbpl:nvrt,i)+etal
        ztmp(nvrt)=etal
        do k=kbpl+1,nvrt-1
          ztmp(k)=0
          do j=1,3
            nd=elnode(nodel2(j),nnel1)
            tmp=(eta3(nd)+dp(nd))*sigma_lcl(max(kbp(nd),k),nd)+eta3(nd)
            ztmp(k)=ztmp(k)+tmp*arco(j)
          enddo !j
        enddo !k
      else
        write(11,*)'Unknown ivcor:',ivcor
        stop
      endif

      if(1==2) then
!--------------------------------------------------------
      do k=kz,nvrt
        kin=k-kz+1
        hmod2=min(dep,h_s)
        if(hmod2<=h_c) then
          ztmp(k)=sigma(kin)*(hmod2+etal)+etal
        else if(etal<=-h_c-(dep-h_c)*theta_f/s_con1) then
          write(11,*)'Pls choose a larger h_c (2):',etal,h_c
          stop
        else
!          ztmp(k)=etal*(1+sigma(kin))+h_c*sigma(kin)+(hmod2-h_c)*cs(kin)
        endif

!       Following to prevent underflow
        if(k==kz) ztmp(k)=-hmod2
        if(k==nvrt) ztmp(k)=etal
      enddo !k

      if(dep<=h_s) then
        kbpl=kz
      else !z levels
!       Find bottom index
        kbpl=0
        do k=1,kz-1
          if(-dep>=ztot(k).and.-dep<ztot(k+1)) then
            kbpl=k
            exit
          endif
        enddo !k
        if(kbpl==0) then
          write(11,*)'Cannot find a bottom level at foot:',dep
          stop
        endif
        ztmp(kbpl)=-dep
        do k=kbpl+1,kz-1
          ztmp(k)=ztot(k)
        enddo !k
      endif
!--------------------------------------------------------
      endif !1==2

      do k=kbpl+1,nvrt
        if(ztmp(k)-ztmp(k-1)<=0) then
          write(11,*)'Inverted z-level in quicksearch:',nnel1,etal,dep,ztmp(k)-ztmp(k-1)
          stop
        endif
      enddo !k

      if(zt<=ztmp(kbpl)) then
        !Avoid getting stuck at bottom
        zt=ztmp(kbpl+1)
        zrat=0
        jlev1=kbpl+1
      else if(zt>=ztmp(nvrt)) then
        zt=ztmp(nvrt)
        zrat=0
        jlev1=nvrt
      else
        jlev1=0
        do k=kbpl,nvrt-1
          if(zt>=ztmp(k).and.zt<=ztmp(k+1)) then 
            jlev1=k+1
            exit
          endif
        enddo !k
        if(jlev1==0) then
          write(11,*)'Cannot find a vert. level:',zt,etal,dep
          write(11,*)(ztmp(k),k=kbpl,nvrt)
          stop
        endif
        zrat=(ztmp(jlev1)-zt)/(ztmp(jlev1)-ztmp(jlev1-1))
      endif

      if(zrat<0.or.zrat>1) then
        write(11,*)'Sigma coord. wrong (4):',jlev1,zrat
        stop
      endif

!      if(kbpl==kz) then !in pure S region
!        ss=(1-zrat)*sigma(jlev1-kz+1)+zrat*sigma(jlev1-kz)
!      else
!        ss=-99
!      endif

      end subroutine quicksearch

!======================================================================
      subroutine pt_in_poly2(i34,x,y,xp,yp,inside)
!     (Double-precision) Routine to perform point-in-polygon
!     (triangle/quads) test.
!     Inputs:
!            i34: 3 or 4 (type of elem)
!            x(i34),y(i34): coord. of polygon/elem. (counter-clockwise)
!            xp,yp: point to be tested
!     Outputs:
!            inside: 1, inside
      use global, only : small1
      implicit real*8(a-h,o-z)
      integer, intent(in) :: i34
      real*8, intent(in) :: x(i34),y(i34),xp,yp
      integer, intent(out) :: inside

      real*8 :: swild(i34)

      inside=0
      do j=1,i34
        j1=j+1
        if(j1>i34) j1=j1-i34
        swild(j)=signa(x(j),x(j1),xp,y(j),y(j1),yp)
      enddo !j
      ae=minval(swild(1:i34))
      if(ae>-small1) inside=1

      end subroutine pt_in_poly2
     
!======================================================================
      subroutine pt_in_poly3(i34,x,y,xp,yp,arco,nodel)
!     (Double-precision) Routine to perform point-in-polygon
!     (triangle/quads) test with assumption that it's inside, and calculate the area coord.
!     (for quad, split it into 2 triangles and return the 3 nodes and
!     area coord.)
!     Inputs:
!            i34: 3 or 4 (type of elem)
!            x(i34),y(i34): coord. of polygon/elem. (counter-clockwise)
!            xp,yp: point to be tested
!     Outputs:
!            arco(3), nodel(3) : area coord. and 3 local node indices (valid only if inside)
      implicit real*8(a-h,o-z)
      integer, intent(in) :: i34
      real*8, intent(in) :: x(i34),y(i34),xp,yp
      integer, intent(out) :: nodel(3)
      real*8, intent(out) :: arco(3)

      !Local
      integer :: list(3)
      real*8 :: ar(2),swild(2,3)

      !Areas of up to 2 triangles
      ar(1)=signa(x(1),x(2),x(3),y(1),y(2),y(3))
      ar(2)=0 !init
      if(i34==4) ar(2)=signa(x(1),x(3),x(4),y(1),y(3),y(4))
      if(ar(1)<=0.or.i34==4.and.ar(2)<=0) then
        print*, 'Negative area:',i34,ar,x,y
        stop
      endif

      aa=0
      ae_min=huge(1.0d0)
      do m=1,i34-2 !# of triangles
        if(m==1) then
          list(1:3)=(/1,2,3/) !local indices
        else !quads
          list(1:3)=(/1,3,4/)
        endif !m
        do j=1,3
          j1=j+1
          j2=j+2
          if(j1>3) j1=j1-3
          if(j2>3) j2=j2-3
          swild(m,j)=signa(x(list(j1)),x(list(j2)),xp,y(list(j1)),y(list(j2)),yp) !temporary storage
          aa=aa+abs(swild(m,j))
        enddo !j=1,3

        ae=abs(aa-ar(m))/ar(m)
        if(ae<=ae_min) then
          ae=ae_min
          nodel(1:3)=list(1:3)
          arco(1:3)=swild(m,1:3)/ar(m)
          arco(1)=max(0.d0,min(1.d0,arco(1)))
          arco(2)=max(0.d0,min(1.d0,arco(2)))
          if(arco(1)+arco(2)>1) then
            arco(3)=0
            arco(2)=1-arco(1)
          else
            arco(3)=1-arco(1)-arco(2)
          endif
        endif
      enddo !m

      end subroutine pt_in_poly3

!======================================================================
