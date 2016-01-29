
*---------------------------
        program grain
*---------------------------

*     Integration numerique des equations de Hodgkin-Huxley
*     modele du grain - 5 variables - ca equation de bilan
*     modele reduit de Maex & De Schutter

        implicit none


        double precision y(5), f(5), yout(5)
        double precision g_na, g_k, g_ca, g_k_ca, g_l
        double precision v_na, v_k, v_l, v_ca
        double precision c_m, I_inj
        double precision FK, A, d,ff,k_ca
        double precision dt, t_max, t_print, t_jump
        double precision  m_inf, h_inf, n_inf,
     1		s_inf, kca_inf
        double precision tol, t, t_sol, dt_print, adtol
        double precision max1, max2, maxt1,maxt2

        integer  imax, nm
        integer n ,q, i, i_max,k, kk

        external fcn

        common/cond/ g_na, g_k, g_ca, g_k_ca, g_l
        common/pot_equilibre/v_na, v_k, v_l, v_ca
        common/k_cell/c_m, I_inj
        common/gating/m_inf, n_inf
        common/ca/FK, A, d, ff, k_ca
        common/steps/dt, t_max, t_print, dt_print, kk
        common/maxi/nm
        common/maxr/max1,max2,maxt1,maxt2


*        lecture donnees
*        open(3, file='gfreq.dat', status='unknown')
*        read(3,*)dt, t_max, t_print
*        close(3)



*        assignation des valeurs des constantes


*       conductances en nS

        g_na=172.0d0
        g_k=28.0d0*1.0
        g_ca=2.9d0*20.0d0
        g_k_ca=56.5d0*1.6
*        g_k_ca=56.5d0*1.0
*        g_l=0.0d0
        g_l=0.0971d0
*        g_l=0.5d0
*       sinon dv/dt toujours <0 si I_inj>0


*        potentiels d inversion en mV
        v_na=57.7d0
        v_k=-90.d0
        v_ca=80.d0
        v_l=-65.d0
*        v_l=120.d0


*       c_m en pF
        c_m=3.14d0

*       F en C/micromol
        FK=0.096500d0
*       surface cellule en um^2 (4*pi*r^2)
        A=314
*       epaisseur couche sous-membranaire en um
        d=0.084
*      rapport Ca libre/Ca lie
*        ff=0.1d0
        ff=0.01d0
*      extrusion calcium en ms^-1
        k_ca=10.0d0


        imax=30
*      I_inj=0.06535097d0
        I_inj=0.0d0
        do i=1, imax
            max1=100
            max2=100
            nm=0
            kk=1

*       parametres pour la routine d'integration
        t_max=10000.d0
		t_print=0.d0
		t_jump=6000.d0
        dt=0.005d0
        dt_print=0.005d0
        t=0.0d0
        adtol=0.01d0
        i_max=int(t_max/dt)
        q=1


        write(6,*)I_inj

*       conditions initiales
        y(1)=-60.d0
        y(2)=0.1d0
        y(3)=0.1d0
        y(4)=0.1d0
        y(5)=0.1d0

        yout(1)=-60.d0
        yout(2)=0.1d0
        yout(3)=0.1d0
        yout(4)=0.1d0
        yout(5)=0.1d0



        open(18)
        do n=1,i_max

           call adams(yout, t, dt, yout,fcn, adtol)
           call output(t, yout)
           t=t+dt
*		   if(t.ge.t_jump)then
* 		   I_inj=0.0d0
*		   endif
        end do
        close(18)
		I_inj=I_inj-1d0
        end do
        close(17)
      end


*-------------------------------------------------------
      subroutine adams(y, t, h, yout,fcn, adtol)
*-------------------------------------------------------
      implicit none


      double precision y(5), yt(5), yold(5),
     1     ynew(5),yout(5)
      double precision f(5), ftemp(5)
      double precision t, th, h, adtol
      integer i, j, k

*     prediction: Euler
      call fcn(t,y,f)

            do k=1,5
                yt(k)=y(k)+h*f(k)
                yold(k)=yt(k)
            end do

      th=t+h

*     correction: trapezes
      call fcn(th,yt,ftemp)

            do k=1,5
                ynew(k)=y(k)+h*(f(k)+ftemp(k))/2
            end do




            do k=1,5
               yout(k)=ynew(k)
            end do




      return
      end





*--------------------------------------------------------
      subroutine fcn(t,y,f)
*--------------------------------------------------------
      implicit none

        double precision y(5), f(5)
        double precision g_na, g_k, g_ca, g_k_ca, g_l
        double precision v_na, v_k, v_l, v_ca
        double precision c_m, I_inj
        double precision FK, A, d, ff, k_ca
        double precision dt, t_max, t_print, t, dt_print
        double precision I_na, I_k, I_ca, I_k_ca, I_l
	double precision  m_inf, h_inf, n_inf,
     1		s_inf, kca_inf
	double precision  tau_h,
     1		tau_s, tau_t, tau_kca
	double precision alpha_m, beta_m
	double precision alpha_h, beta_h
	double precision alpha_n, beta_n
	double precision alpha_s, beta_s
	double precision alpha_kca, beta_kca
        integer kk

        common/cond/g_na, g_k, g_ca, g_k_ca, g_l
  	    common/pot_equilibre/v_na, v_k, v_l, v_ca
	    common/k_cell/c_m, I_inj
	    common/gating/m_inf, n_inf
        common/ca/FK, A, d, ff, k_ca
        common/steps/dt, t_max, t_print, dt_print, kk

* gating Ina
	m_inf=alpha_m(y(1))/(alpha_m(y(1))+beta_m(y(1)))

	h_inf=alpha_h(y(1))/(alpha_h(y(1))+beta_h(y(1)))
	tau_h=1.0/(alpha_h(y(1))+beta_h(y(1)))
	if(tau_h.lt.0.045d0)then
		tau_h=0.045d0
	endif


* gating Ikdr
	n_inf=alpha_n(y(1))/(alpha_n(y(1))+beta_n(y(1)))


* gating Ica
	s_inf=alpha_s(y(1))/(alpha_s(y(1))+beta_s(y(1)))
	tau_s=1/(alpha_s(y(1))+beta_s(y(1)))

* gating Ikca
	kca_inf=alpha_kca(y(1),y(5))/(alpha_kca(y(1),y(5))
     1    +beta_kca(y(1),y(5)))
        tau_kca=1/(alpha_kca(y(1),y(5))
     1    +beta_kca(y(1),y(5)))


* currents

	I_na=g_na*m_inf**3*y(2)*(y(1)-v_na)
	I_k=g_k*n_inf**4*1*(y(1)-v_k)
	I_ca=g_ca*y(3)**2*1*(y(1)-v_ca)
      	I_k_ca=g_k_ca*y(4)*(y(1)-v_k)
	I_l=g_l*(y(1)-v_l)

* Equations

*     	dv/dt
	f(1)=-1/c_m*(I_na+I_k+I_ca+I_k_ca+I_l+I_inj)
*	dh/dt
	f(2)=(h_inf-y(2))/tau_h
*	ds/dt
	f(3)=(s_inf-y(3))/tau_s
*	dkca_act/dt
	f(4)=(kca_inf-y(4))/tau_kca
*	dcai/dt
	f(5)=ff*(-I_ca/(2*FK*A*d)-k_ca*y(5))

      return
      end

*--------------------------------------------------------
      subroutine output(t_sol, y)
*--------------------------------------------------------
      implicit none


        double precision y(5), f(5)
        double precision g_na, g_k, g_ca, g_k_ca, g_l
        double precision v_na, v_k, v_l, v_ca
        double precision c_m, I_inj
        double precision FK, A, d, ff, k_ca
        double precision dt, t_max, t_print, t_sol, dt_print
        double precision I_na, I_k, I_ca, I_k_ca, I_l
	    double precision  m_inf, h_inf, n_inf,
     1		s_inf, kca_inf
	    double precision  tau_h,
     1		tau_s, tau_t, tau_kca
	    double precision alpha_m, beta_m
	    double precision alpha_h, beta_h
	    double precision alpha_n, beta_n
	    double precision alpha_s, beta_s
	    double precision alpha_kca, beta_kca

	    double precision max1,max2,maxt1,maxt2
	    integer nm, kk



        common/cond/g_na, g_k, g_ca, g_k_ca, g_l
  	    common/pot_equilibre/v_na, v_k, v_l, v_ca
	    common/k_cell/c_m, I_inj
	    common/gating/m_inf, n_inf
        common/ca/FK, A, d, ff,k_ca
        common/steps/dt, t_max, t_print, dt_print, kk
	    common/maxi/nm
	    common/maxr/max1,max2,maxt1,maxt2

	I_na=g_na*m_inf**3*y(2)*(y(1)-v_na)
	I_k=g_k*n_inf**4*1*(y(1)-v_k)
	I_ca=g_ca*y(3)**2*1*(y(1)-v_ca)
	I_k_ca=g_k_ca*y(4)*(y(1)-v_k)
	I_l=g_l*(y(1)-v_l)

        if(t_sol.ge.t_print)then
        if(abs(t_sol-(kk*dt_print)).lt.0.001d0) then
*           	write(18,10)t_sol-t_print, y(1), y(5)
*     1		,I_na, I_k, I_ca, I_k_ca, I_l
        	call flush(18)
			kk=kk+1
* pour detecter amplitude max
           if((max2.gt.y(1)).and.(max2.gt.max1)
     1       .and.(max2.ge.-20.d0))then
	   if(nm.eq.0)then
              		maxt1=t_sol
			write(6,*)'maxt1',maxt1
			write(6,*)'max1',max1
			nm=1

	   elseif(nm.eq.1)then
		maxt2=t_sol
		write(6,*)'maxt2',maxt2
		write(6,*)'max2',max2

	        if(max2.ge.-20.d0)then
     			write(17,*)I_inj,1000.d0/(maxt2-maxt1),t_sol
			call flush(17)
			write(6,*)I_inj,1000.d0/(maxt2-maxt1)
			nm=2
		endif
	   endif
	 endif
           endif
           max1=max2
           max2=y(1)

        endif

*      t_sol=t_sol+0.01

 10   format(15(e15.6))
      return
      end



* gating Ina
*--------------------------------------------------------
	function alpha_m(v)
*--------------------------------------------------------
	implicit none
        double precision v, alpha_m
	alpha_m=7.5d0*dexp(0.081d0*(v+39.d0))
	return
	end
*--------------------------------------------------------
	function beta_m(v)
*--------------------------------------------------------
	implicit none
        double precision v, beta_m
	beta_m=7.5d0*dexp(-0.066d0*(v+39.d0))
	return
	end
*--------------------------------------------------------
	function alpha_h(v)
*--------------------------------------------------------
	implicit none
        double precision v, alpha_h
	alpha_h=0.6d0*dexp(-0.089d0*(v+50.d0))
	return
	end
*--------------------------------------------------------
	function beta_h(v)
*--------------------------------------------------------
	implicit none
        double precision v, beta_h
	beta_h=0.6d0*dexp(0.089d0*(v+50.d0))
	return
	end



* gating Ikdr
*--------------------------------------------------------
	function alpha_n(v)
*--------------------------------------------------------
	implicit none
        double precision v, alpha_n
	alpha_n=0.85d0*dexp(0.073d0*(v+38.d0))
	return
	end
*--------------------------------------------------------
	function beta_n(v)
*--------------------------------------------------------
	implicit none
        double precision v, beta_n
	beta_n=0.85d0*dexp(-0.018d0*(v+38.d0))
	return
	end



* gating Ica
*--------------------------------------------------------
	function alpha_s(v)
*--------------------------------------------------------
	implicit none
        double precision v, alpha_s
	alpha_s=8.0d0/(1+dexp(-0.072d0*(v-5.d0)))
	return
	end
*--------------------------------------------------------
	function beta_s(v)
*--------------------------------------------------------
	implicit none
        double precision v, beta_s
	beta_s=0.1d0*(v+8.9)/(dexp(0.2d0*(v+8.9d0))-1)
	return
	end


* gating Ikca
*--------------------------------------------------------
	function alpha_kca(v,ca)
*--------------------------------------------------------
	implicit none
        double precision v,ca, alpha_kca
	alpha_kca=12.5d0/(1+(1.5*0.1*dexp(-0.085d0*v)/ca))
	return
	end


*--------------------------------------------------------
	function beta_kca(v,ca)
*--------------------------------------------------------
	implicit none
        double precision v, ca,  beta_kca
	beta_kca=7.5d0/(1+(ca/(0.15*0.1*dexp(-0.077d0*v))))
	return
	end
