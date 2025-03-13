c234567**1*********2*********3*********4*********5*********6*********7**X
c.. SUBROUTINE PEL_PPO_DIST
c..
c..   This subroutine returns the momenta of the electron and positron
c..
c*****|****************************************************************|X

      SUBROUTINE pel_ppo_dist(beta,m_pair,
     &p0_pair,px_pair,py_pair,pz_pair,p0_el_lab,px_el_lab,
     &py_el_lab,pz_el_lab,p0_po_lab,px_po_lab,py_po_lab,pz_po_lab)         

      implicit none

      include 'defs.f'

c-----|----------------------------------------------------------------|X
                  
      real*8 phi_el,x,sin_theta_el,theta_el,phi_po,theta_po,beta_x_pair
      real*8 p0_el,p0_po,p_el,p_po,px_el,py_el,pz_el,px_po,py_po,pz_po
      real*8 beta_y_pair,beta_z_pair,beta_pair,gamma_pair,p0_el_cm,px_el_cm
      real*8 py_el_cm,pz_el_cm,p0_po_cm,px_po_cm,py_po_cm,pz_po_cm
      real*8 beta_x_urqmd,beta_y_urqmd,beta_z_urqmd,beta_urqmd
      real*8 gamma_urqmd,p0_el_lab,px_el_lab,py_el_lab,pz_el_lab
      real*8 p0_po_lab,px_po_lab,py_po_lab,pz_po_lab,x_gstar,y_gstar
      real*8 z_gstar,t_gstar,x_nucleon,y_nucleon,z_nucleon,gamma,beta,s
      real*8 m_pair,mass_nucleon,t_nucleon,x_el,y_el,z_el,x_po,y_po,z_po
      real*8 p0_pair,px_pair,py_pair,pz_pair,pt_pair

      real*8 ranff
c      external ranff

c-----|----------------------------------------------------------------|X
    
c.. Determine the mass of the lepton-type 
      if(dimuon) then
         m_lept=mass_muon
      else !dielectron
         m_lept=mass_electron
      endif

c      write(0,*)'lepton mass',m_lept ! Debug only 
c      write(0,*)'beta,m,p0l,pxl,pyl,pzl',beta,m_pair,p0_pair,px_pair,py_pair,
c     &pz_pair ! Debug only 

      phi_el=2d0*pi*ranff(seedinit)
      x=ranff(seedinit)     

c      write(0,*)'phi_x,x',phi_el,x !Debug only
                                 
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

      gamma=1.d0/sqrt(1.0d0-beta**2)

c      write(0,*)'gamma',gamma ! Debug only

      x_el=sin(theta_el)*cos(phi_el)               
      y_el=sin(theta_el)*sin(phi_el)
      z_el=cos(theta_el)

      x_po=sin(theta_po)*cos(phi_po)               
      y_po=sin(theta_po)*sin(phi_po) 
      z_po=cos(theta_po)


      p0_el=m_pair/2d0                               
      p0_po=m_pair/2d0

      p_el=sqrt(p0_el**2-m_lept**2)       
      p_po=sqrt(p0_po**2-m_lept**2)         

      px_el=x_el*p_el                              
      py_el=y_el*p_el 
      pz_el=z_el*p_el

      px_po=x_po*p_po                              
      py_po=y_po*p_po
      pz_po=z_po*p_po


      beta_x_pair=-px_pair/p0_pair
      beta_y_pair=-py_pair/p0_pair
      beta_z_pair=-pz_pair/p0_pair

c      write(0,*)'betaxpair,betaypair,betazpair',beta_x_pair,beta_y_pair,beta_z_pair !Debug only

      beta_pair=sqrt(beta_x_pair**2+beta_y_pair**2+beta_z_pair**2)

      gamma_pair=1.0d0/sqrt(1.0d0-beta_pair**2)

      p0_el_cm=gamma_pair*p0_el-beta_x_pair*gamma_pair*px_el-
     &beta_y_pair*gamma_pair*py_el-beta_z_pair*gamma_pair*pz_el

c      write(0,*)'betapair,gammapair,betaxpair,betaypair,betazpair',beta_pair,gamma_pair,
c     &beta_x_pair,beta_y_pair,beta_z_pair !Debug only

c      write(0,*)'p0_el,px_el,py_el,pz_el',p0_el,px_el,py_el,pz_el !Debug only
     
      px_el_cm=(-beta_x_pair*gamma_pair)*p0_el+
     &(1+(gamma_pair-1)*((beta_x_pair*beta_x_pair)/
     &(beta_pair*beta_pair)))*px_el+
     &(gamma_pair-1)*((beta_x_pair*beta_y_pair)/
     &(beta_pair*beta_pair))*py_el+
     &(gamma_pair-1)*((beta_x_pair*beta_z_pair)/
     &(beta_pair*beta_pair))*pz_el

      py_el_cm=(-beta_y_pair*gamma_pair)*p0_el+
     &(gamma_pair-1)*((beta_x_pair*beta_y_pair)/
     &(beta_pair*beta_pair))*px_el+
     &(1+(gamma_pair-1)*((beta_y_pair*beta_y_pair)/
     &(beta_pair*beta_pair)))*py_el+
     &(gamma_pair-1)*((beta_y_pair*beta_z_pair)/
     &(beta_pair*beta_pair))*pz_el 

      pz_el_cm=(-beta_z_pair*gamma_pair)*p0_el+
     &(gamma_pair-1)*((beta_x_pair*beta_z_pair)/
     &(beta_pair*beta_pair))*px_el+
     &(gamma_pair-1)*((beta_y_pair*beta_z_pair)/
     &(beta_pair*beta_pair))*py_el+
     &(1+(gamma_pair-1)*((beta_z_pair*beta_z_pair)/
     &(beta_pair*beta_pair)))*pz_el 


      p0_po_cm=gamma_pair*p0_po-beta_x_pair*gamma_pair*px_po- 
     &beta_y_pair*gamma_pair*py_po-beta_z_pair*gamma_pair*pz_po
     
      px_po_cm=(-beta_x_pair*gamma_pair)*p0_po+
     &(1+(gamma_pair-1)*((beta_x_pair*beta_x_pair)/
     &(beta_pair*beta_pair)))*px_po+
     &(gamma_pair-1)*((beta_x_pair*beta_y_pair)/
     &(beta_pair*beta_pair))*py_po+
     &(gamma_pair-1)*((beta_x_pair*beta_z_pair)/
     &(beta_pair*beta_pair))*pz_po

      py_po_cm=(-beta_y_pair*gamma_pair)*p0_po+
     &(gamma_pair-1)*((beta_x_pair*beta_y_pair)/
     &(beta_pair*beta_pair))*px_po+
     &(1+(gamma_pair-1)*((beta_y_pair*beta_y_pair)/
     &(beta_pair*beta_pair)))*py_po+
     &(gamma_pair-1)*((beta_y_pair*beta_z_pair)/
     &(beta_pair*beta_pair))*pz_po 

      pz_po_cm=(-beta_z_pair*gamma_pair)*p0_po+
     &(gamma_pair-1)*((beta_x_pair*beta_z_pair)/
     &(beta_pair*beta_pair))*px_po+
     &(gamma_pair-1)*((beta_y_pair*beta_z_pair)/
     &(beta_pair*beta_pair))*py_po+
     &(1+(gamma_pair-1)*((beta_z_pair*beta_z_pair)/
     &(beta_pair*beta_pair)))*pz_po 

c      write(0,*)'CM,electron:p0,px,py,pz',p0_el_cm,px_el_cm,
c     &py_el_cm,pz_el_cm
c      write(0,*)'CM,positron:p0,px,py,pz',p0_po_cm,px_po_cm,
c     &py_po_cm,pz_po_cm

c-----|----------------------------------------------------------------|X

      beta_x_urqmd=0.0d0
      beta_y_urqmd=0.0d0
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
 
c-----|----------------------------------------------------------------|X
      
      end

c*****|****************************************************************|X
