      SUBROUTINE lobo_dir(beta,gamma,m_res,
     &p0_res,px_res,py_res,pz_res,pt_res,p0_el_lab,px_el_lab,
     &py_el_lab,pz_el_lab,p0_po_lab,px_po_lab,py_po_lab,pz_po_lab)
                    
      IMPLICIT NONE
      
c**********************************************************************!

      real*8 pi,mass_electron,mass_muon
      real*8 phi_el,x,sin_theta_el,theta_el,phi_po,theta_po,beta_x_res
      real*8 p0_el,p0_po,p_el,p_po,px_el,py_el,pz_el,px_po,py_po,pz_po
      real*8 beta_y_res,beta_z_res,beta_res,gamma_res,p0_el_cm,px_el_cm
      real*8 py_el_cm,pz_el_cm,p0_po_cm,px_po_cm,py_po_cm,pz_po_cm
      real*8 beta_x_urqmd,beta_y_urqmd,beta_z_urqmd,beta_urqmd
      real*8 gamma_urqmd,p0_el_lab,px_el_lab,py_el_lab,pz_el_lab
      real*8 p0_po_lab,px_po_lab,py_po_lab,pz_po_lab,x_gstar,y_gstar
      real*8 z_gstar,t_gstar,x_nucleon,y_nucleon,z_nucleon,gamma,beta,s
      real*8 m_res,mass_nucleon,t_nucleon,x_el,y_el,z_el,x_po,y_po,z_po
      real*8 p0_res,px_res,py_res,pz_res,pt_res
    
c... Muon / Electron
      logical dimuon
      common /mu_el/ dimuon

      parameter(pi=3.141592654d0)
      parameter(mass_electron=0.000511d0)
      parameter(mass_muon=0.1057d0)

c... Random Function
      real*8 ranff
      external ranff

      integer seedinit
      common /randomf/ seedinit

c**********************************************************************!

      phi_el=2d0*pi*ranff(seedinit)
      x=ranff(seedinit)                                      
      sin_theta_el=(x-0.5)*2d0                        
      theta_el=asin(sin_theta_el)+pi/2d0

      phi_po=phi_el + pi                            
      if(phi_po.gt.2*pi) then                               
      phi_po=phi_po-2*pi
      end if

      if(theta_el.le.pi) then
      theta_po=pi-theta_el
      else
      theta_po=0
      end if

      x_el=sin(theta_el)*cos(phi_el)               
      y_el=sin(theta_el)*sin(phi_el)
      z_el=cos(theta_el)

      x_po=sin(theta_po)*cos(phi_po)               
      y_po=sin(theta_po)*sin(phi_po) 
      z_po=cos(theta_po)


      p0_el=m_res/2d0                               
      p0_po=m_res/2d0

c-------------------------------------------!
      if(.not.(dimuon)) then                !
c... DIELECTRONS                            !
      p_el=sqrt(p0_el**2-mass_electron**2)  !      
      p_po=sqrt(p0_po**2-mass_electron**2)  !
      endif                                 ! 
                                            !
      if(dimuon) then                       !
c... DIMUONS                                !
      p_el=sqrt(p0_el**2-mass_muon**2)      !  
      p_po=sqrt(p0_po**2-mass_muon**2)      !
      endif                                 !
c-------------------------------------------!

      px_el=x_el*p_el                              
      py_el=y_el*p_el 
      pz_el=z_el*p_el

      px_po=x_po*p_po                              
      py_po=y_po*p_po
      pz_po=z_po*p_po


      beta_x_res=-px_res/p0_res
      beta_y_res=-py_res/p0_res
      beta_z_res=-pz_res/p0_res

      beta_res=sqrt(beta_x_res**2+beta_y_res**2+beta_z_res**2)

      gamma_res=1d0/sqrt(1-beta_res**2)

      p0_el_cm=gamma_res*p0_el-beta_x_res*gamma_res*px_el-
     &beta_y_res*gamma_res*py_el-beta_z_res*gamma_res*pz_el
     
      px_el_cm=(-beta_x_res*gamma_res)*p0_el+
     &(1+(gamma_res-1)*((beta_x_res*beta_x_res)/
     &(beta_res*beta_res)))*px_el+
     &(gamma_res-1)*((beta_x_res*beta_y_res)/
     &(beta_res*beta_res))*py_el+
     &(gamma_res-1)*((beta_x_res*beta_z_res)/
     &(beta_res*beta_res))*pz_el

      py_el_cm=(-beta_y_res*gamma_res)*p0_el+
     &(gamma_res-1)*((beta_x_res*beta_y_res)/
     &(beta_res*beta_res))*px_el+
     &(1+(gamma_res-1)*((beta_y_res*beta_y_res)/
     &(beta_res*beta_res)))*py_el+
     &(gamma_res-1)*((beta_y_res*beta_z_res)/
     &(beta_res*beta_res))*pz_el 

      pz_el_cm=(-beta_z_res*gamma_res)*p0_el+
     &(gamma_res-1)*((beta_x_res*beta_z_res)/
     &(beta_res*beta_res))*px_el+
     &(gamma_res-1)*((beta_y_res*beta_z_res)/
     &(beta_res*beta_res))*py_el+
     &(1+(gamma_res-1)*((beta_z_res*beta_z_res)/
     &(beta_res*beta_res)))*pz_el 


      p0_po_cm=gamma_res*p0_po-beta_x_res*gamma_res*px_po- 
     &beta_y_res*gamma_res*py_po -beta_z_res*gamma_res*pz_po
      px_po_cm=(-beta_x_res*gamma_res)*p0_po+
     &(1+(gamma_res-1)*((beta_x_res*beta_x_res)/
     &(beta_res*beta_res)))*px_po+
     &(gamma_res-1)*((beta_x_res*beta_y_res)/
     &(beta_res*beta_res))*py_po+
     &(gamma_res-1)*((beta_x_res*beta_z_res)/
     &(beta_res*beta_res))*pz_po

      py_po_cm=(-beta_y_res*gamma_res)*p0_po+
     &(gamma_res-1)*((beta_x_res*beta_y_res)/
     &(beta_res*beta_res))*px_po+
     &(1+(gamma_res-1)*((beta_y_res*beta_y_res)/
     &(beta_res*beta_res)))*py_po+
     &(gamma_res-1)*((beta_y_res*beta_z_res)/
     &(beta_res*beta_res))*pz_po 

      pz_po_cm=(-beta_z_res*gamma_res)*p0_po+
     &(gamma_res-1)*((beta_x_res*beta_z_res)/
     &(beta_res*beta_res))*px_po+
     &(gamma_res-1)*((beta_y_res*beta_z_res)/
     &(beta_res*beta_res))*py_po+
     &(1+(gamma_res-1)*((beta_z_res*beta_z_res)/
     &(beta_res*beta_res)))*pz_po 

      beta_x_urqmd=0
      beta_y_urqmd=0
      beta_z_urqmd=-beta

      beta_urqmd=beta
      gamma_urqmd=gamma

      p0_el_lab=gamma_urqmd*p0_el_cm-beta_x_urqmd*gamma_urqmd*px_el_cm- 
     &beta_y_urqmd*gamma_urqmd*py_el_cm-beta_z_urqmd*gamma_urqmd*
     &pz_el_cm

      px_el_lab=(-beta_x_urqmd*gamma)*p0_el_cm + 
     &(1+(gamma_urqmd-1)*((beta_x_urqmd*beta_x_urqmd)/
     &(beta_urqmd*beta_urqmd)))*px_el_cm + 
     &(gamma_urqmd-1)*((beta_x_urqmd*beta_y_urqmd)/
     &(beta_urqmd*beta_urqmd))*py_el_cm + 
     &(gamma_urqmd-1)*((beta_x_urqmd*beta_z_urqmd)/
     &(beta_urqmd*beta_urqmd))*pz_el_cm

      py_el_lab=(-beta_y_urqmd*gamma)*p0_el_cm+ 
     &(gamma_urqmd-1)*((beta_x_urqmd*beta_y_urqmd)/
     &(beta_urqmd*beta_urqmd))*px_el_cm+ 
     &(1+(gamma_urqmd-1)*((beta_y_urqmd*beta_y_urqmd)/
     &(beta_urqmd*beta_urqmd)))*py_el_cm+
     &(gamma_urqmd-1)*((beta_y_urqmd*beta_z_urqmd)/
     &(beta_urqmd*beta_urqmd))*pz_el_cm 

      pz_el_lab=(-beta_z_urqmd*gamma)*p0_el_cm+
     &(gamma_urqmd-1)*((beta_x_urqmd*beta_z_urqmd)/
     &(beta_urqmd*beta_urqmd))*px_el_cm+
     &(gamma_urqmd-1)*((beta_y_urqmd*beta_z_urqmd)/
     &(beta_urqmd*beta_urqmd))*py_el_cm+
     &(1+(gamma_urqmd-1)*((beta_z_urqmd*beta_z_urqmd)/
     &(beta_urqmd*beta_urqmd)))*pz_el_cm 
  
      p0_po_lab=gamma_urqmd*p0_po_cm-beta_x_urqmd*gamma_urqmd*px_po_cm- 
     &beta_y_urqmd*gamma_urqmd*py_po_cm -beta_z_urqmd*gamma_urqmd*
     &pz_po_cm

      px_po_lab=(-beta_x_urqmd*gamma)*p0_po_cm+
     &(1+(gamma_urqmd-1)*((beta_x_urqmd*beta_x_urqmd)/
     &(beta_urqmd*beta_urqmd)))*px_po_cm+ 
     &(gamma_urqmd-1)*((beta_x_urqmd*beta_y_urqmd)/
     &(beta_urqmd*beta_urqmd))*py_po_cm+
     &(gamma_urqmd-1)*((beta_x_urqmd*beta_z_urqmd)/
     &(beta_urqmd*beta_urqmd))*pz_po_cm
     
      py_po_lab=(-beta_y_urqmd*gamma)*p0_po_cm+
     &(gamma_urqmd-1)*((beta_x_urqmd*beta_y_urqmd)/
     &(beta_urqmd*beta_urqmd))*px_po_cm+
     &(1+(gamma_urqmd-1)*((beta_y_urqmd*beta_y_urqmd)/
     &(beta_urqmd*beta_urqmd)))*py_po_cm+
     &(gamma_urqmd-1)*((beta_y_urqmd*beta_z_urqmd)/
     &(beta_urqmd*beta_urqmd))*pz_po_cm 

      pz_po_lab=(-beta_z_urqmd*gamma)*p0_po_cm+
     &(gamma_urqmd-1)*((beta_x_urqmd*beta_z_urqmd)/
     &(beta_urqmd*beta_urqmd))*px_po_cm+
     &(gamma_urqmd-1)*((beta_y_urqmd*beta_z_urqmd)/
     &(beta_urqmd*beta_urqmd))*py_po_cm+
     &(1+(gamma_urqmd-1)*((beta_z_urqmd*beta_z_urqmd)/
     &(beta_urqmd*beta_urqmd)))*pz_po_cm 

c**********************************************************************!

      end
