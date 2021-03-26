
submodule (phase:plastic) dislotwingnd
 
 real(pReal),                                    parameter,           private :: &
   kB = 1.38e-23_pReal                                                                                 !< Boltzmann constant in J/Kelvin

 type, private :: tParameters
   real(pReal) :: &
     mu                   = 1.0_pReal, &                                                            !< equivalent shear modulus
     nu                   = 1.0_pReal, &                                                            !< equivalent shear Poisson's ratio
     c1                   = 1.0_pReal, &                                                            !< material constant for passing stress
     c2                   = 1.0_pReal, &                                                            !< material constant for the formation of ssd dislo
     c3                   = 1.0_pReal                                                               !< material constant for the annihilation of ssd dislo
   real(pReal),                  dimension(:),     allocatable :: & 
     tau_0, &                                                                                       !< strength due to elements in solid solution
     b_sl, &                                                                                        !< absolute length of burgers vector [m] for each slip system
     Delta_F,&                                                                                      !< activation energy for glide [J] for each slip system
     Delta_V,&                                                                                      !< activation wolume for glide [b^3] for each slip system
     v0                                                                                             !< dislocation velocity prefactor [m/s] for each slip system
   real(pReal),                  dimension(:,:),   allocatable :: &
     M_sl, &                                                                                        !!!< slip direction m
     T_sl
   real(pReal),                  dimension(:,:,:), allocatable :: &
     P_sl
   integer :: & 
     sum_N_sl                                                                                       !< total number of active slip system
   integer(kind(undefined_ID)),  dimension(:),     allocatable :: &
     outputID                                                                                       !< ID of each post result output
 end type
 
 type, private :: tDislotwingndState
   real(pReal),                  dimension(:,:),   pointer :: &
     rho_mob, &
     rho_gnd_edge, &
     rho_gnd_screw, &
     gamma_sl
 end type tDislotwingndState 

 type, private :: tDislotwingndMicrostructure
   real(pReal),                  dimension(:,:),   allocatable :: &
      tau_pass
 end type tDislotwingndMicrostructure   

!--------------------------------------------------------------------------------------------------
! containers for parameters and state
 type(tParameters),                 allocatable, dimension(:), private :: param
 type(tDislotwingndState),          allocatable, dimension(:), private :: &
   dotState, &
   state
 type(tDislotwingndMicrostructure), allocatable, dimension(:), private :: microstructure

contains

!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function plastic_dislotwingnd_init() result(myPlasticity)

 logical, dimension(:), allocatable :: myPlasticity
 integer :: &
   ph, i, &
   Nmembers, &
   sizeState, sizeDotState, &
   startIndex, endIndex
 
 integer,     dimension(:), allocatable :: &
    N_sl
  real(pReal), allocatable, dimension(:) :: &
    rho_mob_0, &                                                                                    !< initial unipolar dislocation density per slip system
    rho_gnd_edge_0, &                                                                                 !< initial gnd_edge dislocation density ?special position?
    rho_gnd_screw_0
  character(len=pStringLen) :: &
    extmsg = ''
  class(tNode), pointer :: &
    phases, &
    phase, &
    mech, &
    pl


  myPlasticity = plastic_active('dislotwingnd')
  if(count(myPlasticity) == 0) return

  print'(/,a)', ' <<<+-  phase:mechanical:plastic:dislotwin init  -+>>>'
  print'(a,i0)', ' # phases: ',count(myPlasticity); flush(IO_STDOUT)

  print*, 'Ma and Roters, Acta Materialia 52(12):3603–3612, 2004'
  print*, 'https://doi.org/10.1016/j.actamat.2004.04.012'//IO_EOL

  print*, 'Roters et al., Computational Materials Science 39:91–95, 2007'
  print*, 'https://doi.org/10.1016/j.commatsci.2006.04.014'//IO_EOL

  print*, 'Wong et al., Acta Materialia 118:140–151, 2016'
  print*, 'https://doi.org/10.1016/j.actamat.2016.07.032'


  phases => config_material%get('phase')
  allocate(param(phases%length))
  allocate(state(phases%length))
  allocate(dotState(phases%length))
  allocate(dependentState(phases%length))


  do ph = 1, phases%length
    if(.not. myPlasticity(ph)) cycle

    associate(prm => param(ph), dot => dotState(ph), stt => state(ph), dst => dependentState(ph))

    phase => phases%get(ph)
    mech  => phase%get('mechanics')
    pl  => mech%get('plasticity')

#if defined (__GFORTRAN__)
    prm%output = output_asStrings(pl)
#else
    prm%output = pl%get_asStrings('output',defaultVal=emptyStringArray)
#endif

    ! This data is read in already in lattice
    prm%mu  = lattice_mu(ph)
    prm%nu  = lattice_nu(ph)
 
!--------------------------------------------------------------------------------------------------
! slip related parameters
   N_sl      = pl%get_asInts('N_sl',defaultVal=emptyIntArray)
   prm%sum_N_sl  = sum(abs(N_sl))
   slipActive: if (prm%sum_N_sl > 0) then
     prm%P_sl = lattice_SchmidMatrix_slip(N_sl,phase%get_asString('lattice'),&
                                          phase%get_asFloat('c/a',defaultVal=0.0_pReal))
     prm%M_sl = lattice_slip_direction(N_sl,phase%get_asString('lattice_structure'),&                           !!!slip direction m
                                       phase%get_asFloat('c/a',defaultVal=0.0_pReal))
     prm%T_sl = lattice_slip_transverse(N_sl,phase%get_asString('lattice_structure'),&                           !!!normal slip plane m
                                         phase%get_asFloat('c/a',defaultVal=0.0_pReal))
     
     rho_mob_0                = pl%get_asFloats('rho_mob_0',      requiredSize=size(N_sl))
     rho_gnd_edge_0           = pl%get_asFloats('rho_gnd_edge_0', requiredSize=size(N_sl), &
                                                      defaultVal=[(0.0_pReal, i=1,size(N_sl))])                       !first guessing
     rho_gnd_screw_0          = pl%get_asFloats('rho_gnd_screw_0',requiredSize=size(N_sl), &
                                                       defaultVal=[(0.0_pReal, i=1,size(N_sl))])                       !first guessing
     prm%v_0                  = pl%get_asFloats('v_0',            requiredSize=size(N_sl))
     prm%b_sl                 = pl%get_asFloats('b_sl',           requiredSize=size(N_sl))
     prm%Q_s                  = pl%get_asFloats('Q_s',            requiredSize=size(N_sl))
     prm%tau_0                = pl%get_asFloats('tau_0',          requiredSize=size(N_sl))
     
     prm%Delta_V              = pl%get_asFloat('activ_volume') * prm%b_sl**3.0_pReal
     prm%c1                   = pl%get_asFloat('c1_passsingstress')
     prm%c2                   = pl%get_asFloat('c2_rhoformation')
     prm%c3                   = pl%get_asFloat('c3_rhoathermal')   
     
     ! expand: family => system
     rho_mob_0        = math_expand(rho_mob_0,       N_sl)
     prm%tau_0        = math_expand(prm%tau_0,       N_sl)
     prm%v0           = math_expand(prm%v0,          N_sl)
     prm%b_sl         = math_expand(prm%b_sl,        N_sl)
     prm%Delta_F      = math_expand(prm%Delta_F,     N_sl)
     prm%Delta_V      = math_expand(prm%Delta_V,     N_sl)
     
     ! sanity checks
     if (any(rho_mob_0        <  0.0_pReal))         extmsg = trim(extmsg)//' rho_mob_0'
     if (any(prm%tau_0        <  0.0_pReal))         extmsg = trim(extmsg)//' tau_0'
     if (any(prm%v0           <  0.0_pReal))         extmsg = trim(extmsg)//' v0'
     if (any(prm%b_sl         <= 0.0_pReal))         extmsg = trim(extmsg)//' b_sl'
     if (any(prm%Delta_F      <= 0.0_pReal))         extmsg = trim(extmsg)//' Delta_F'
     if (any(prm%Delta_V      <= 0.0_pReal))         extmsg = trim(extmsg)//' Delta_V'

   else slipActive
     rho_mob_0 = emptyRealArray
     allocate(prm%b_sl,prm%Q_s,prm%v_0,source=emptyRealArray)
   endif slipActive

!--------------------------------------------------------------------------------------------------
! allocate state arrays
   Nmembers = count(material_phaseAt2 == ph)
   sizeDotState = size(['rho_mob      ', &
                        'rho_gnd_edge ', &
                        'rho_gnd_screw', &
                        'gamma_sl     ']) * prm%sum_N_sl                   !!!now only for slip
   sizeState = sizeDotState

   call phase_allocateState(plasticState(ph),Nmembers,sizeState,sizeDotState,0)                                     

!--------------------------------------------------------------------------------------------------
! locally defined state aliases and initialization of state0 and aTolState
   startIndex = 1
   endIndex   = prm%sum_N_sl
   stt%rho_mob=>plasticState(p)%state(startIndex:endIndex,:)
   stt%rho_mob= spread(prm%rho_mob_0,2,Nmembers)
   dot%rho_mob=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(ph)%atol(startIndex:endIndex) = pl%get_asFloat('atol_rho',defaultVal=1.0_pReal)
   if (any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pReal)) extmsg = trim(extmsg)//' atol_rho'
   
   startIndex = endIndex + 1
   endIndex   = endIndex + prm%sum_N_sl
   sst%rho_gnd_edge=>plasticState(p)%dotState(startIndex:endIndex,:)
   dot%rho_gnd_edge=>plasticState(p)%dotState(startIndex:endIndex,:)
   
   startIndex = endIndex + 1
   endIndex   = endIndex + prm%sum_N_sl
   stt%rho_gnd_screw=>plasticState(p)%state(startIndex:endIndex,:)
   dot%rho_gnd_screw=>plasticState(p)%dotState(startIndex:endIndex,:)
   
   startIndex = endIndex + 1
   endIndex   = endIndex + prm%sum_N_sl
   stt%gamma_sl=>plasticState(p)%state(startIndex:endIndex,:)
   dot%gamma_sl=>plasticState(p)%dotState(startIndex:endIndex,:)
   plasticState(ph)%atol(startIndex:endIndex) = 1.0e-2_pReal
   ! global alias
   plasticState(ph)%slipRate        => plasticState(ph)%dotState(startIndex:endIndex,:)

   allocate(dst%tau_pass              (prm%sum_N_sl,Nmembers),source=0.0_pReal)
   
   plasticState(ph)%state0 = plasticState(ph)%state                                                ! ToDo: this could be done centrally
   end associate
!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
   if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'(dislotwingnd)')

 enddo

end function plastic_dislotwingnd_init

!--------------------------------------------------------------------------------------------------
!> @brief calculates plastic velocity gradient and its tangent
!--------------------------------------------------------------------------------------------------
module subroutine dislotwingnd_LpAndItsTangent(Lp,dLp_dMp,Mp,T,ph,me)

 real(pReal), dimension(3,3),     intent(out) :: Lp
 real(pReal), dimension(3,3,3,3), intent(out) :: dLp_dMp
 real(pReal), dimension(3,3),     intent(in)  :: Mp
 integer,                         intent(in)  :: ph,me
 real(pReal),                     intent(in)  :: T
 
 integer :: i,k,l,m,n
 real(pReal) :: &
    ddot_gamma_dtau
 real(pReal), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma_sl,ddot_gamma_dtau_slip
 
 associate(prm => param(ph), stt => state(ph))

 
 Lp = 0.0_pReal
 dLp_dMp = 0.0_pReal  

 call kinetics_slip(Mp,T,ph,me,dot_gamma_sl,ddot_gamma_dtau_slip)
 slipContribution: do i = 1, prm%sum_N_sl
   Lp = Lp + dot_gamma_sl(i)*prm%P_sl(1:3,1:3,i)
   forall (k=1:3,l=1:3,m=1:3,n=1:3) &
     dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                      + ddot_gamma_dtau_slip(i) * prm%P_sl(k,l,i) * prm%P_sl(m,n,i)
 enddo slipContribution
 
 end associate
 
end subroutine dislotwingnd_LpAndItsTangent

!--------------------------------------------------------------------------------------------------
!> @brief calculates the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
module subroutine dislotwingnd_dotState(Mp,T,ph,me,el)

 real(pReal), dimension(3,3),  intent(in):: &
   Mp                                                                                               !< Mandel stress
 real(pReal),                  intent(in) :: &
   T                                                                                                !< temperature at integration point
 integer,                      intent(in) :: &
   ph, &
   me
 real(pReal) :: tau
 real(pReal), dimension(param(ph)%sum_N_sl) :: &
   dot_gamma_sl, &
   dot_rho_ssd_formation, &
   dot_rho_ssd_athermal_anni
 real(pReal), dimension(3,param(ph)%sum_N_sl) :: &
   grad_gamma_sl
 integer :: i
 
 associate(prm => param(ph),    stt => state(ph), &
           dot => dotstate(ph), dst => microstructure(ph))

 call kinetics_slip(Mp,T,ph,me,dot_gamma_sl)
 dot%gamma_sl(:,me) = abs(dot_gamma_sl) 

 SSD: do i = 1, prm%sum_N_sl
   tau = math_tensordot(Mp,prm%P_sl(1:3,1:3,i))
   dot_rho_ssd_formation(i)     = prm%c2*abs(dot_gamma_sl(i))/(prm%mu*prm%b_sl(i)**2.0_pReal)
   dot_rho_ssd_athermal_anni(i) = prm%c3*stt%rho_mob(i,me)*abs(dot_gamma_sl(i))
 enddo SSD

 dot%rho_mob(:,me) = dot_rho_ssd_formation - dot_rho_ssd_athermal_anni

 call getgradient_for_gamma_sl(ph, el)

 GND: do i = 1, prm%sum_N_sl
  dot%rho_gnd_edge (i,me)  = -math_inner(grad_gamma_sl(1:3,i),prm%M_sl(1:3,i))/prm%b_sl(i)
  dot%rho_gnd_screw(i,me)  =  math_inner(grad_gamma_sl(1:3,i),prm%T_sl(1:3,i))/prm%b_sl(i)         
 enddo GND
 
 end associate
 
end subroutine dislotwingnd_dotState

!--------------------------------------------------------------------------------------------------
!@brief evaluate gradient field of shear strain
!--------------------------------------------------------------------------------------------------
module subroutine getgradient_for_gamma_sl(ph, el)

 use spectral_utilities
 use discretization_grid
 
 integer,       intent(in) :: &
   ph,
   el
 integer :: &
   eli, ip, &
   k, j, i, &
   phi, me, &
   n
 real(pReal), dimension(param(ph)%sum_N_sl,grid(1),grid(2),grid3) :: &
    gamma_sl_forall
 real(pReal), dimension(3,param(ph)%sum_N_sl,discretization_Nelems) :: &   
    grad_gamma_sl_forall
  
 !$OMP PARALLEL DO PRIVATE(ph,me)
 eli = 0
 do k = 1, grid3; do j = 1, grid(2); do i = 1, grid(1)
   eli = eli + 1
     do ip = 1, size(material_phaseMemberAt,2)
       phi = material_phaseAt(1,el)
       me = material_phaseMemberAt(1,ip,eli)
        
       gamma_sl_forall(:,i,j,k)=plasticState(phi)%state(phi)%gamma_sl(:,me)

    enddo
 enddo;enddo;enddo   

 do n = 1, size(gamma_sl_forall,1)
   scalarField_real = 0.0_pReal
   scalarField_real(1:grid(1),1:grid(2),1:grid3) = gamma_sl_forall(n,1:grid(1),1:grid(2),1:grid3)         
   call utilities_FFTscalarForward
   call utilities_fourierScalarGradient                                                              !< calculate gradient of shear rate field
   call utilities_FFTvectorBackward
   
   eli = 0
    do k = 1, grid3;  do j = 1, grid(2);  do i = 1,grid(1)
       eli = eli + 1
       grad_gamma_sl_forall(1:3,n,eli) = vectorField_real(1:3,i,j,k)                 !save the gardient field into grad_gamma_sl
    enddo;enddo;enddo
 enddo
   grad_gamma_sl(1:3,:) =grad_gamma_sl_forall(1:3,:,el)
end subroutine getgradient_for_gamma_sl

!--------------------------------------------------------------------------------------------------
!> @brief calculates derived quantities from state
!--------------------------------------------------------------------------------------------------
module subroutine dislotwingnd_dependentState(T,ph,me)

 integer,       intent(in) :: &
   ph, &
   me
 real(pReal),   intent(in) :: &
   T

 integer :: &
   i
 associate(prm => param(ph),&
           stt => state(ph),&
           dst => microstructure(ph))  
           
 !* threshold stress for dislocation motion
 do i = 1 , prm%sum_N_sl
  rho_gnd(i)  = sqrt(sst%rho_gnd_edge(i,me)**2.0_pReal + sst%rho_gnd_screw(i,me)**2.0_pReal))
  dst%tau_pass(i,me) = prm%c1*prm%mu*prm%b_sl(i)*sqrt(stt%rho_mob(i,me)+rho_gnd(i))                      !!!now not consider the interaction
 enddo
 
 end associate

end subroutine dislotwingnd_dependentState

!--------------------------------------------------------------------------------------------------
!> @brief writes results to HDF5 output file
!--------------------------------------------------------------------------------------------------
module subroutine plastic_dislotwingnd_results(ph,group)

  integer,          intent(in) :: ph
  character(len=*), intent(in) :: group

  integer :: o

 associate(prm => param(ph), stt => state(ph), dst => dependentState(ph))
  outputsLoop: do o = 1,size(prm%output)
    select case(trim(prm%output(o)))

      case('rho_mob')
        if(prm%sum_N_sl>0) call results_writeDataset(group,stt%rho_mob,trim(prm%output(o)), &
                                                     'mobile dislocation density','1/m²')
      case('rho_gnd_edge')
        if(prm%sum_N_sl>0) call results_writeDataset(group,stt%rho_gnd_edge,trim(prm%output(o)), &
                                                     'edge_gnd dislocation density','1/m²')
      case('rho_gnd_screw')
        if(prm%sum_N_sl>0) call results_writeDataset(group,stt%rho_gnd_screw,trim(prm%output(o)), &
                                                     'screw_gnd dislocation density','1/m²')
      
      case('gamma_sl')
        if(prm%sum_N_sl>0) call results_writeDataset(group,stt%gamma_sl,trim(prm%output(o)), &
                                                     'plastic shear','1')
      case('tau_pass')
        if(prm%sum_N_sl>0) call results_writeDataset(group,dst%tau_pass,trim(prm%output(o)), &
                                                     'passing stress for slip','Pa')

    end select
  enddo outputsLoop
  end associate

end subroutine plastic_dislotwingnd_results

!--------------------------------------------------------------------------------------------------
!> @brief Shear rates on slip systems, their derivatives with respect to resolved stress and the
!  resolved stresss
!> @details Derivatives and resolved stress are calculated only optionally.
! NOTE: Against the common convention, the result (i.e. intent(out)) variables are the last to
! have the optional arguments at the end
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics_slip(Mp,T,ph,me, &
                              dot_gamma_sl,ddot_gamma_dtau_slip,tau_slip)

 real(pReal), dimension(3,3),  intent(in) :: &
   Mp                                                                                               !< Mandel stress
 real(pReal),                  intent(in) :: &
   T                                                                                                !< temperature
 integer,                      intent(in) :: &
   ph, &
   me
   
 real(pReal), dimension(param(ph)%sum_N_sl), intent(out) :: &
   dot_gamma_sl
 real(pReal), dimension(param(ph)%sum_N_sl), optional, intent(out) :: &
   ddot_gamma_dtau_slip, &
   tau_slip
 real(pReal), dimension(param(ph)%sum_N_sl) :: &
   ddot_gamma_dtau

 real(pReal), dimension(param(ph)%sum_N_sl) :: &
   tau, &
   tau_eff, &
   ActEnergyRatio, & 
   ActVolumeRatio, &
   sinh_tau, &
   dsinh_tau_dtau 
 
 integer :: i
 
 associate(prm => param(ph), stt => state(ph), dst => microstructure(ph))

 do i = 1, prm%sum_N_sl
   tau(i) = math_tensordot(Mp,prm%P_sl(1:3,1:3,i))
 enddo

 tau_eff = abs(tau)-dst%tau_pass(:,me)-prm%tau_0 
 
  do i = 1, prm%sum_N_sl
    if (tau_eff(i)>tol_math_check) then
      ActEnergyRatio(i)  = exp(-prm%Delta_F(i)/(kB*T))
      ActVolumeRatio(i)  = prm%Delta_V(i)/(kB*T)
      sinh_tau(i)        = sinh(ActVolumeRatio(i)*sign(tau_eff(i),tau(i))
      
      dot_gamma_sl(i)    = stt%rho_mob(i,me)*prm%b_sl(i)*prm%v0(i)*ActEnergyRatio(i)*sinh_tau(i)
      
      dsinh_tau_dtau(i)  = cosh(ActVolumeRatio(i)*tau_eff(i))*ActVolumeRatio(i)
      
      ddot_gamma_dtau(i) = stt%rho_mob(i,me)*prm%b_sl(i)*prm%v0(i)*ActEnergyRatio(i)*dsinh_tau_dtau(i)
    
     else
     
      dot_gamma_sl(i)    = 0.0_pReal
      ddot_gamma_dtau(i) = 0.0_pReal
      
     end if 
  enddo
 end associate
 
 if(present(ddot_gamma_dtau_slip)) ddot_gamma_dtau_slip = ddot_gamma_dtau
 if(present(tau_slip))             tau_slip             = tau
 
end subroutine kinetics_slip

end submodule dislotwingnd 
