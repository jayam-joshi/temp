c  first code for one dimensional
       implicit real*8(a-h,o-z)
      character*8 name1(901)
      dimension col(9000),a0(9000),sumx(9000),sumy(9000),sum2(9000)
      dimension x(9000,2),y(9000,2),del(9000),sumx1(9000)
      dimension theta(9000,2),sumy2(9000),density(1019),densityi(1019),
     1 avesum2(6000000,20),p_acf(6000000,20)
     
cccccccccccccccccccccccccccc change the dimension of denity,densityi when system size is changed!
     
       integer ndim1(9000)
       real a,c,dt,vim,msd3,L,ii,jj,acf_density
       integer inew, iold, iseed
       character(100) :: filename_1, filename_2
 
       itst=0
c       
       nt=10 ! constant for the system
       L=16.0032d0 ! system size, change when changing packing fraction
       ii=16.0032d0
       jj=16.0032d0
       ax=0.0
       ay=0.0
       v=1.0d0
       n=40
       dt=0.001   ! time step
       D_R=0.1d0 !stregth of the gaussian white noise
       a00=0.2d0 ! particle radius
       ntime=5000000 ! total time steps
       nsnap=ntime/99
       am=-100.0! stregth of the repulsion force
       nens=3
       inew=2  ! for the particle
       iold=1 ! previous position of the particle
       iseed=-137481
        nrr=40
        
       atime1=100

       open(unit=1,file="fnamerhol")
       read(1,*)(name1(i),i=1,900)

       open(unit=38,file="msd_0.5_v1.dat")   ! change packing fraction here
       open(unit=11,file="dacf_0.5_v1.dat")
 

        pi=2.0*asin(1.0) ! how to define pi

        npro2=0
        npro=0
        npro1=0
        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do iens=1,nens
        npro2=npro2+1
        iseed = iseed + 2*iens
        ndim=0
        write(filename_1, '(A,I0,A)') 'msd_0.5_v1_ens', iens, '.dat' ! change packing fraction here
        open(unit=50+iens, file=trim(filename_1))
        
        write(filename_2, '(A,I0,A)') 'd_0.5_v1_ens', iens, '.dat'
        open(unit=100+iens, file=trim(filename_2))
 
       do i=1,nint(50*ii)
       do j=1,nint(50*jj)      

       if(ndim<1019) then    ! change no of particles here
        ndim=ndim+1
        x(ndim,iold)=i*0.02
        y(ndim,iold)=j*0.02
        theta(ndim,iold)=2.0*pi*ran2(iseed)
        endif
        
        enddo
        enddo
        
        packfrac=pi*(a00*a00*ndim)/(L*L)

           do k=1,ndim
           sumx1(k)=0.0d0
           sumy2(k)=0.0d0
           enddo

        do itime=1,ntime
           k1=0 
           
           do i=1,ndim
           density(i) = 0.0d0
           enddo

          avesum2(itime,iens)=0.0d0
          p_acf(itime,iens)=0.0d0
          
          do k=1,ndim
           ndim1(k)=0    ! particles in contact with the particle_ density
           sumfx=0.0d0
           sumfy=0.0d0
           do j=1,ndim
          
           if(j.ne.k) then
            disx=(x(j,iold)-x(k,iold))
            if(abs(disx).ge.0.5*L) disx=L-
     1         (FLOOR((disx/(0.5*L))))

            disy=(y(j,iold)-y(k,iold))
             if(abs(disy).ge.0.5*L) disy=L-
     1         (FLOOR((disy/(0.5*L))))

            dis=sqrt(disx**2+disy**2)
            disxunit=(disx)/(dis)
            disyunit=(disy)/(dis)

            sigma=2*a00
          if(dis.le.sigma) then
          ndim1(k)=ndim1(k)+1
          afx=((sigma-dis))*disxunit
          afy=((sigma-dis))*disyunit
           else
           afx=0.0d0
           afy=0.0d0
           ndim1(k)=ndim1(k)
           endif
           
           sumfx=sumfx+afx
           sumfy=sumfy+afy
           
           endif

          enddo ! here the particle jnekloop end here.

          
          a=ran2(iseed)
          c=ran2(iseed)
          b=sqrt(-2*alog(a))*cos(2*pi*c)
           


          x(k,inew)=x(k,iold)+(v*cos(theta(k,iold))
     1      +sumfx*am)*dt
          y(k,inew)=y(k,iold)+(v*sin(theta(k,iold))
     1       +sumfy*am)*dt

          dixx=(v*cos(theta(k,iold))+sumfx*am)*dt
     1       
          diyy=(v*sin(theta(k,iold))+sumfy*am)*dt
     1          
          theta(k,inew)=theta(k,iold)+
     1       b*(sqrt(2*D_R))*(sqrt(dt))
         
          sumx1(k)=sumx1(k)+dixx
          sumy2(k)=sumy2(k)+diyy
          sum2(k)=(sumx1(k))**2+(sumy2(k))**2

          if(theta(k,inew).ge.2*pi) theta(k,inew)=theta(k,inew)-2*pi
          if(theta(k,inew).lt.0.0) theta(k,inew)=theta(k,inew)+2*pi
          if(x(k,inew).le.0.0) x(k,inew)=x(k,inew)+L
          if(x(k,inew).gt.L) x(k,inew)=x(k,inew)-L ! boundary condition 
          if(y(k,inew).le.0.0) y(k,inew)=y(k,inew)+L
          if(y(k,inew).gt.L) y(k,inew)=y(k,inew)-L

                            
          avesum2(itime,iens)=avesum2(itime,iens)+sum2(k) 
                   
          density(k) = ndim1(k)/dfloat(6)           ! local density
          enddo
          
         
         avesum2(itime,iens)=avesum2(itime,iens)/dfloat(ndim)
         density = density - sum(density)/dfloat(ndim)   ! density fluctuation
         
         if(mod(itime,50).eq.0) then
         write(50+iens,*)itime*dt,avesum2(itime,iens)     ! Storing msd data for each ensemble
         endif
         
                  
         if(itime.ge.4900000) then
         if(itime.eq.4900000) then
          do k=1,ndim
           densityi(k) = density(k)
          enddo
         endif
         p_acf(itime,iens) = sum(densityi*density)/dfloat(ndim)  ! density auto-correlation
         if(mod(itime,50).eq.0) write(100+iens,*)(itime-4900000)*dt, p_acf(itime, iens)   !storing density-correlation with time for each ensemble
         endif
         
           ir=0
           if(iens.eq.1) then 
           do ip=1,(ntime-itst)/nsnap
           mtime=itst+ip*nsnap

            ir=ir+1
           if(itime.eq.mtime) then

            npro=npro+1
            
            open(unit=31,file=name1(npro))

            do k=1,ndim
              
            write(31,*)x(k,inew),y(k,inew),a00

            enddo

           endif
           
           enddo
           endif
        itemp=iold
        iold=inew
        inew=itemp
           
     
      
         enddo !time loop ends here
          
         enddo !ensemble loop end here
        
         msd3=0.0
         do itime=1,ntime
         acf_density=0.0d0
         do iens=1,nens
         msd3=msd3+avesum2(itime,iens)
         acf_density = acf_density + p_acf(itime,iens)
         enddo
         msd3=msd3/dfloat(nens)
         acf_density = acf_density/dfloat(nens)
         if(mod(itime,50).eq.0) write(38,*)itime*dt,msd3
         if(itime.ge.4900000) then
         if(mod(itime,50).eq.0) write(11,*)(itime-4900000)*dt,acf_density
         endif
         
         enddo

          stop
          end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       FUNCTION ran2(idum)
       INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
       REAL*8 ran2,AM,EPS,RNMX
       PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     * IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     * NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
       INTEGER idum2,j,k,iv(NTAB),iy
       SAVE iv,iy,idum2
       DATA idum2/123456789/, iv/NTAB*0/, iy/0/
       if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 11 j=NTAB+8,1,-1
           k=idum/IQ1
           idum=IA1*(idum-k*IQ1)-k*IR1
           if (idum.lt.0) idum=idum+IM1
           if (j.le.NTAB) iv(j)=idum
11       continue
         iy=iv(1)
       endif
       k=idum/IQ1
       idum=IA1*(idum-k*IQ1)-k*IR1
       if (idum.lt.0) idum=idum+IM1
       k=idum2/IQ2
       idum2=IA2*(idum2-k*IQ2)-k*IR2
       if (idum2.lt.0) idum2=idum2+IM2
       j=1+iy/NDIV
       iy=iv(j)-idum2
       iv(j)=idum
       if(iy.lt.1)iy=iy+IMM1
       ran2=min(AM*iy,RNMX)

       return
       END
