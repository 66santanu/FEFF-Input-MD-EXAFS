        program feffInp 
        implicit none
        integer:: nframes,natom 
        integer:: j_ion,dimn,a,mina,natomtype
        integer:: i,j,k,m,icount,a1
        real:: bx,by,bz,minr,rj,center,distance,x,y,z
        real,dimension(3):: Xvec 
        real,allocatable,dimension(:,:):: xt
        real,allocatable,dimension(:):: rt
        integer,allocatable,dimension(:):: index_atom,inp_atom_index
        character(len=5), allocatable,dimension(:) :: ATOM,inp_atom_name

        ! Reading input parameters 
        read(*,*) natom,j_ion,natomtype
        allocate(inp_atom_index(natomtype),inp_atom_name(natomtype))
        do i=1,natomtype
           read(*,*) inp_atom_index(i),inp_atom_name(i)
           enddo 
        read(*, *) bx, by, bz
        ! *******************************************************
        allocate(ATOM(natom))
        dimn=3
        allocate(xt(natom,dimn))
        allocate(rt(natom))
        allocate(index_atom(natom))
        !********************************************************

        open(unit=20,file='Structure.xyz')
        read(20,*)
        read(20,*)
        
        ! Reading the positions
        do i=1,natom 
           read(20,*) ATOM(i),(xt(i,m),m=1,3)
        end do
        close(20)

        ! Distance calculations from the select ion  
        do i=1,natom 
           x=xt(i,1)-xt(j_ion,1)
           y=xt(i,2)-xt(j_ion,2)
           z=xt(i,3)-xt(j_ion,3)

           x=x-(bx*Anint(x/bx))
           y=y-(by*Anint(y/by))
           z=z-(bz*Anint(z/bz))

           rt(i)=sqrt(x**2+y**2+z**2)
           index_atom(i) = i
        end do 
 
        do i=1,natom-1
           minr = 1000.0;
           mina = -1;
           do a = i,natom
              if(rt(a)<minr) then
                 minr=rt(a)
                 mina=a
              endif
           enddo 
           rj=rt(i)
           rt(i)=minr;
           rt(mina)=rj;

           j=index_atom(i);
           index_atom(i)=index_atom(mina);
           index_atom(mina)=j;
         enddo

         open(unit=40,file="feff-conf.xyz")
         icount=0;
         center=0.0;
         write(40,*) "ATOM"
         write(40,"(3(F8.4,2X),2X,I2,2X,A5,2X,F8.4,2X,I5)") center,center,center,icount,ATOM(index_atom(1)),center,icount 
         do a=2,natom
            x=xt(index_atom(a),1)-xt(index_atom(1),1)
            y=xt(index_atom(a),2)-xt(index_atom(1),2)
            z=xt(index_atom(a),3)-xt(index_atom(1),3)

            x=x-(bx*Anint(x/bx))
            y=y-(by*Anint(y/by))
            z=z-(bz*Anint(z/bz))

            Xvec(1)=x
            Xvec(2)=y
            Xvec(3)=z

            distance=sqrt(x**2+y**2+z**2)

            do k=2,natomtype
             if(ATOM(index_atom(a))==inp_atom_name(k)) then
               a1=a-1     
               write(40,"(3(F8.4,2X),2X,I2,2X,A5,2X,F8.4,2X,I5)") (Xvec(i),i=1,3),inp_atom_index(k),ATOM(index_atom(a)),distance,a1
             endif        
            enddo   
         enddo
         write(40,*) "END"
         close(40)
                     
         deallocate(xt,rt,index_atom,inp_atom_index,inp_atom_name)
!CONTAINS

!!!!! Difference vector with pbc
!Function vecPbc(vector1, vector2, bx, by, bz)
!  implicit none
!  real, dimension(3) :: vector1, vector2
!  real,dimension(3):: vecPbc
!  real               :: rx, ry, rz, bx, by, bz
!  
!  rx = vector1(1) - vector2(1)
!  ry = vector1(2) - vector2(2)
!  rz = vector1(3) - vector2(3)
!  
!  vecPbc(1) = rx - bx * Anint(rx/bx)
!  vecPbc(2) = ry - by * Anint(ry/by)
!  vecPbc(3) = rz - bz * Anint(rz/bz)
!end function vecPbc

!!!!! distance with pbc
!Function distPbc(vector1, vector2, bx, by, bz)
!  implicit none
!  real, dimension(3) :: vector1, vector2
!  real               :: rx, ry, rz, distPbc, bx, by, bz
!
!  rx = vector1(1) - vector2(1)
!  ry = vector1(2) - vector2(2)
!  rz = vector1(3) - vector2(3)
!
!  rx = rx - bx * Anint(rx/bx)
!  ry = ry - by * Anint(ry/by)
!  rz = rz - bz * Anint(rz/bz)

!  distPbc = sqrt(rx**2+ry**2+rz**2)
!end function distPbc

end program feffInp 
