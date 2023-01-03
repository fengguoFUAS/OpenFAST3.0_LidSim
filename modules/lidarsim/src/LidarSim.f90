    MODULE LidarSim

    USE LidarSim_Types
    USE LidarSim_Subs
    USE NWTC_Library
    !USE InflowWind
    !USE InflowWind_Subs
    USE InflowWind_Types
    !USE AeroDyn14_Types
    

    IMPLICIT NONE
    PRIVATE

    TYPE(ProgDesc), PARAMETER   ::  IfW_Ver = ProgDesc( 'LidarSim', 'v0.20', '12-December-2019' )
    PUBLIC                      ::  LidarSim_Init                                
    PUBLIC                      ::  LidarSim_CalcOutput    
    PUBLIC                      ::  LidarSim_End

    CONTAINS

    !#########################################################################################################################################################################

    SUBROUTINE LidarSim_Init(InitInp, y, p, Ifw_p,InitOutData, ErrStat, ErrMsg )
	
    IMPLICIT                                NONE
    CHARACTER(*),                           PARAMETER       ::  RoutineName="LidarSim_Init"
    
    TYPE(LidarSim_InitInputType),           INTENT(IN   )   ::  InitInp             ! Input data for initialization routine
    TYPE(LidarSim_OutputType),              INTENT(  OUT)   ::  y                   ! Output data for the lidar module
    TYPE(LidarSim_ParameterType),           INTENT(INOUT)   ::  p                   ! Parameter data for the lidar module
    TYPE(InflowWind_ParameterType),         INTENT(INOUT)   ::  Ifw_p               !< Inflow Parameters deliver to lidar
    TYPE(LidarSim_InitOutputType),          INTENT(  OUT)   ::  InitOutData         ! Data to initialize the outputs
    INTEGER(IntKi),                         INTENT(  OUT)   ::  ErrStat             !< Error status of the operation
    CHARACTER(*),                           INTENT(  OUT)   ::  ErrMsg              !< Error message if ErrStat /= ErrID_None

    !Local Variables
    TYPE(LidarSim_InputFile)                               	::  InputFileData      !< Structure to load the input file data into
    CHARACTER(1024)                                         ::  RootFileName
    CHARACTER(1024)                                         ::  EchoFileName
   

    ! Temporary variables for error handling
    INTEGER(IntKi)                                          ::  TmpErrStat          !< temporary error message
    CHARACTER(ErrMsgLen)                                    ::  TmpErrMsg           
    
    
    ! Initial error values
    ErrStat        =  0
    ErrMsg         =  ""  
    
    RootFileName  = InitInp%RootName
    IF (LEN_TRIM(RootFileName) == 0) CALL GetRoot( InitInp%InputInitFile, RootFileName )
    EchoFileName  = TRIM(RootFileName)//".ech"

    
    ! Reads the Config from the Input file and writes it into the LidarSim_InputFile data 
    CALL LidarSim_ReadInputFile(InitInp%InputInitFile,EchoFileName , InputFileData, TmpErrStat, TmpErrMsg)!EchoFileName
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    
	 ! needs to be further improved to show which line is wrong
    IF (ErrStat >= AbortErrLev) THEN
        ErrMsg = "error in reading the input file of lidar module. Check the format of the input file!"
        RETURN
    END IF
	
    
    !Transfering InputFileData to the p
    p%MeasurementMaxSteps   =   CEILING(REAL(NINT(InputFileData%t_measurement_interval*100000))/REAL(NINT(InitInp%DT*100000))) !NINT to remove float precision errors. Back to REAL, otherwise the divion ignores everything behind the decima point. Ceiling to round up to next integer
    p%LidarPosition_N(1)    =   InputFileData%LidarPositionX_N
    p%LidarPosition_N(2)    =   InputFileData%LidarPositionY_N
    p%LidarPosition_N(3)    =   InputFileData%LidarPositionZ_N
    p%URef                  =   InputFileData%URef
    p%MAXDLLChainOutputs    =   InputFileData%MAXDLLChainOutputs
    p%SpinnerMountedFlag    =   Inputfiledata%SpinnerMountedFlag
	p%NearestInterpFlag     =   Inputfiledata%NearestInterpFlag
     !p%GatesPerBeam          =   InputFileData%GatesPerBeam   see initialize the measuring points section, if it is Cartesian coordinate it should not have more than one gate
    
   
 
    !Creates the static rotationmatrix from the lidar system to the nacelle system
    CALL LidarSim_CreateRotationMatrix(InputFileData%RollAngle_N,InputFileData%PitchAngle_N,&
    InputFileData%YawAngle_N, p%LidarOrientation_N)

    !----- Calls different Subroutines to initialize the measuring points   
    IF(InputFileData%TrajectoryType == 0) THEN
        CALL LidarSim_InitMeasuringPoints_Cartesian(p, InputFileData, TmpErrStat, TmpErrMsg)   	! Calls Routine to initialize cartesian coordinate inputs
        CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         p%GatesPerBeam          =   1
    ELSEIF(InputFileData%TrajectoryType == 1)THEN
        CALL LidarSim_InitMeasuringPoints_Spherical(p, InputFileData, TmpErrStat, TmpErrMsg )  	! Calls Routine to initialize spherical coordinate inputs
        CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         p%GatesPerBeam          =   InputFileData%GatesPerBeam
    END IF
    
    !----- Calls different Subroutines to initialize the weighting points
    IF(InputFileData%WeightingType == 0) THEN                                                   ! Single point
        CALL AllocAry( p%WeightingDistance,1, 'p%WeightingDistance', TmpErrStat, TmpErrMsg )
        CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
        CALL AllocAry( p%Weighting,1, 'p%Weighting', TmpErrStat, TmpErrMsg )
        CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
        p%WeightingDistance(1) = 0
        p%Weighting(1) = 1
    ELSEIF(InputFileData%WeightingType == 1) THEN                                               ! Calls Routine to initialize the weighting with gaussian distribution
        CALL LidarSim_InitializeWeightingGauss(p, InputFileData, TmpErrStat, TmpErrMsg )       
        CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    ELSEIF(InputFileData%WeightingType == 2) THEN
        CALL LidarSim_InitializeWeightingManual(p, InputFileData, TmpErrStat, TmpErrMsg )      	! Calls Routine to initialize with manual weighting settings
        CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    ENDIF
    
        !----- Calls Subroutines to initialize the wind evolution
    IF(InputFileData%EvolutionFlag == .True.) THEN                                                   
        CALL LidarSim_InitializeWindEvolution(p, InputFileData, TmpErrStat, TmpErrMsg )      	! Calls Routine to initialize wind evolution
        CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
        P%WINDFILEFORMAT	= 		IFW_P%WINDTYPE
    ELSEIF (InputFileData%EvolutionFlag == .False.) THEN 
        ! deliver the flow field from inflow module, since the Taylor frozen theory applied
        CALL LidarSim_InitializeFrozenflow(p, Ifw_p,InputFileData, TmpErrStat, TmpErrMsg )    
        CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
        
    END IF
    
    IF(InputFileData%BladeBlockageFlag == .True.) THEN                                               
        CALL LidarSim_InitializeBladeBlockage(p,InputFileData, TmpErrStat, TmpErrMsg )    
        CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    ENDIF
    
   ! ----- Calls Subroutines to initialize the measurement availability (lidar data availability) 
    IF(InputFileData%AvailabilityFlag == .True.) THEN                                                   
        CALL LidarSim_InitializeAvailability(p, InputFileData, TmpErrStat, TmpErrMsg )      	! Calls Routine to initialize wind evolution
        CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    ENDIF

    CALL LidarSim_InitializeOutputs(y,p, InitOutData, InputFileData, TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    
    !initialize variables and outputs
    p%MeasurementCurrentStep = -1                                                               !< there was no measurement yet
    p%LastMeasuringPoint = 1                                                                    !< First measurement point
    p%NextBeamID = 0
    p%MeasurementTimeStep    = 0                                                                ! First time step of lidar measurement
    
    
    
    CALL LidarSim_DestroyInputFile(InputFileData, ErrStat, ErrMsg)                             	! Calls to destory the data from the inputfile. Important data has to be transfered to the parameter data before

    ErrStat = 		TMPERRSTAT
    ERRMSG	=       TMPERRMSG	                                                                                                                                                                                                                                                  
    
    END SUBROUTINE LidarSim_Init

    !#########################################################################################################################################################################

    SUBROUTINE LidarSim_CalcOutput (Time, y, p, u,ErrStat, ErrMsg )
    
    IMPLICIT                                    NONE
    CHARACTER(*),                               PARAMETER           ::  RoutineName="LidarSim_CalcOutput"
	   
    REAL(DbKi),                                 INTENT(IN   )       ::  Time                !< Current simulation time in seconds
    TYPE(LidarSim_ParameterType),              	INTENT(INOUT)       ::  p
    TYPE(LidarSim_OutputType),                 	INTENT(INOUT)       ::  y                   !< Outputs computed at Time (IN for mesh reasons and data allocation)
    TYPE(LidarSim_InputType),                  	INTENT(IN   )       ::  u                   !< Inputs from other Modules (e.g. ElastoDyn)
    
    !Data for CalcOutput of IfW_Subs and Blockage
    !TYPE(InflowWind_ParameterType),             INTENT(IN   )       ::  IfW_p                       !< Parameters
    !TYPE(AD14_InputType),                       INTENT(IN   )       ::  AD14_u                      ! Inputs from aerodynamic module
    !TYPE(AD14_ParameterType),                   INTENT(IN   )       ::  AD14_p                      ! Parameters from aerodynamic module
    !TYPE(InflowWind_ContinuousStateType),       INTENT(IN   )       ::  IfW_ContStates              !< Continuous states at Time
    !TYPE(InflowWind_DiscreteStateType),         INTENT(IN   )       ::  IfW_DiscStates              !< Discrete states at Time
    !TYPE(InflowWind_ConstraintStateType),       INTENT(IN   )       ::  IfW_ConstrStates            !< Constraint states at Time
    !TYPE(InflowWind_OtherStateType),            INTENT(IN   )       ::  IfW_OtherStates             !< Other/optimization states at Time
    !TYPE(InflowWind_MiscVarType),               INTENT(INOUT)       ::  IfW_m                       !< Misc variables for optimization (not copied in glue code)    
    INTEGER(IntKi),                             INTENT(  OUT)       ::  ErrStat                     !< Error status of the operation
    CHARACTER(*),                               INTENT(  OUT)       ::  ErrMsg                      !< Error message if ErrStat /= ErrID_None
    
    
    
    !Local Variables
    !TYPE(InflowWind_InputType)                                      ::  InputForCalculation         !Data Field needed for the calculation of the windspeed
    !TYPE(InflowWind_OutputType)                                     ::  OutputForCalculation        !datafield in which the calculated speed is stored
    REAL(ReKi)                                                      ::  UnitVector(3)               !Line of Sight Unit Vector
    REAL(ReKi)                                                      ::  MeasuringPosition_I(3)      !Transformed Measuring Position
    REAL(ReKi)                                                      ::  LidarPosition_I(3)          !Transformed Lidar Position
    REAL(ReKi)                                                      ::  Vlos                        !Line of sight speed
    REAL(ReKi)                                                      ::  BladeBlockageStatus         !Status of lidar beam blockage 1 for block, 0 for not
    INTEGER(IntKi)                                                  ::  LoopGatesPerBeam            !Counter to loop through all gate points of a line
    INTEGER(IntKi)                                                  ::  LoopCounter
    
    
    ! Temporary variables for error handling
    INTEGER(IntKi)                                                  ::  TmpErrStat        
    CHARACTER(ErrMsgLen)                                            ::  TmpErrMsg
    
    !Initialize error values
    ErrStat        =  0
    ErrMsg         =  ""
    BladeBlockageStatus = 0 ! it is always 0 if we do not consider
    
    
    IF(p%MeasurementCurrentStep>=p%MeasurementMaxSteps .OR. p%MeasurementCurrentStep == -1)THEN         !Check if there must be a new measurement     !(NINT(Time*1000)-NINT(p%t_last_measurement*1000)) >= NINT(p%t_measurement_interval*1000)
        p%MeasurementCurrentStep = 0
        
        LidarPosition_I = LidarSim_TransformLidarToInertial( u%NacelleMotion,p, (/0.0,0.0,0.0/) )                               !Calculation of the lidar positon ( 0 / 0 / 0 ) in the lidar coordinate system
        
        IF (p%SpinnerMountedFlag==.True.) THEN
            LidarPosition_I = LidarSim_TransformLidarToInertial( u%HubMotion,p, (/0.0,0.0,0.0/) )
        ELSE
            LidarPosition_I = LidarSim_TransformLidarToInertial( u%NacelleMotion,p, (/0.0,0.0,0.0/) )   
        END IF
    
        CALL LidarSim_CalculateIMU(p, y, u)                ! Caltulation of Inertial measurement unit
        
            IF (P%BladeBlockageFlag ==.True.) THEN    ! if consider blade blockage and then we check   
                
                IF (p%SpinnerMountedFlag==.True.) THEN
                    MeasuringPosition_I = LidarSim_TransformLidarToInertial(u%HubMotion,p,p%MeasuringPoints_L(:,p%LastMeasuringPoint))
                ELSE 
                    MeasuringPosition_I = LidarSim_TransformLidarToInertial(u%NacelleMotion,p,p%MeasuringPoints_L(:,p%LastMeasuringPoint)) ! Calculate the Measuringpoint coordinate in the initial system, here we only need the first measuring point, since all the measuring points are in the same ray (line)
                END IF
                
                IF (P%AeroynMode == 14) THEN
                 CALL LidarSim_CheckBladeBlockage_14(p,y,u,LidarPosition_I,MeasuringPosition_I,BladeBlockageStatus)  ! check blockage with AD14
                ELSEIF  (P%AeroynMode == 15) THEN
                 CALL LidarSim_CheckBladeBlockage_15(p,y,u,LidarPosition_I,MeasuringPosition_I,BladeBlockageStatus)  ! check blockage with AD15
                 !BladeBlockageStatus = 10
                END IF
            ELSE
                !BladeBlockageStatus = 0 ! it is always 0 if we do not consider
            END IF
        
            DO LoopGatesPerBeam = 0,p%GatesPerBeam-1
            
                !Transform measuring and lidar position to the inertial system
                IF (p%SpinnerMountedFlag==.True.) THEN
                    MeasuringPosition_I = LidarSim_TransformLidarToInertial(u%HubMotion,p,p%MeasuringPoints_L(:,p%LastMeasuringPoint+LoopGatesPerBeam))
                ELSE 
                    MeasuringPosition_I = LidarSim_TransformLidarToInertial(u%NacelleMotion,p,p%MeasuringPoints_L(:,p%LastMeasuringPoint+LoopGatesPerBeam)) ! Calculate the Measuringpoint coordinate in the initial system
                END IF
                
    
                !Line of Sight
                UnitVector    =   MeasuringPosition_I - LidarPosition_I             !Calculation of the Line of Sight Vector          
                UnitVector    =   UnitVector/NORM2(UnitVector)                      !=>Magnitude = 1
    
                 !Calculation of the wind speed at the calculated position or set to -99 if the measurement is not available
                IF (P%AVAILABILITYFLAG ==.True.) THEN ! if consider data availability
                    
                    
                    !NINT(1)
                    IF ((P%AVAIDATA(p%MeasurementTimeStep+1,p%LastMeasuringPoint+LoopGatesPerBeam)==1) .AND. (TIME<=P%AVAITIMEEND) .AND. (BladeBlockageStatus == 0))THEN
                         CALL LidarSim_CalculateVlos( p, UnitVector, Vlos, MeasuringPosition_I, LidarPosition_I, Time,TmpErrStat, TmpErrMsg) !Calculation of the line of sight wind speed
                         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)    
                         CALL LidarSim_SetOutputs(y,p,Vlos,BladeBlockageStatus,UnitVector,LoopGatesPerBeam,Time)    !Set all outputs to the output variable

                    ELSE
                        Vlos = -99   
                         CALL LidarSim_SetOutputs(y,p,Vlos,BladeBlockageStatus,UnitVector,LoopGatesPerBeam,Time)    !Set all outputs to the output variable
                         y%AllOutputs( 25 + p%GatesPerBeam + LoopGatesPerBeam ) = -99   ! set the LOS to a identifiable wrong value
                    END IF
    
                ELSE
                    IF (BladeBlockageStatus == 0) THEN
                        CALL LidarSim_CalculateVlos( p, UnitVector, Vlos, MeasuringPosition_I, LidarPosition_I, Time,TmpErrStat, TmpErrMsg) !Calculation of the line of sight wind speed
                         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)    
                         CALL LidarSim_SetOutputs(y,p,Vlos,BladeBlockageStatus,UnitVector,LoopGatesPerBeam,Time)    !Set all outputs to the output variable
                    ELSE
                         Vlos = -99   
                         CALL LidarSim_SetOutputs(y,p,Vlos,BladeBlockageStatus,UnitVector,LoopGatesPerBeam,Time)    !Set all outputs to the output variable
                         y%AllOutputs( 25 + p%GatesPerBeam + LoopGatesPerBeam ) = -99   ! set the LOS to a identifiable wrong value, remove the impact by nacelle motion
                    END IF
               
                END IF
    
            ENDDO                
        
        !Choosing which measuring point has to be calculated
        IF(p%LastMeasuringPoint+p%GatesPerBeam > SIZE(p%MeasuringPoints_L,2))THEN                      
            p%LastMeasuringPoint = 1        ! already reached the last point before ? => start over from the beginning
            p%NextBeamID = 0
        ELSE
            p%LastMeasuringPoint = p%LastMeasuringPoint + p%GatesPerBeam
            p%NextBeamID = p%NextBeamID + 1
        END IF
       p%MeasurementTimeStep = p%MeasurementTimeStep + 1 
         
    ELSE                                            !Set NewData signals to zero
        IF(ANY(p%ValidOutputs == 22)) THEN
            DO LoopCounter = 1,SIZE(p%ValidOutputs)
                IF(p%ValidOutputs(LoopCounter) == 22) THEN
                    y%WriteOutput(LoopCounter) = 0
                END IF
            END DO
        END IF
        y%SwapOutputs(1) = 0
    ENDIF
    p%MeasurementCurrentStep = p%MeasurementCurrentStep + 1
    
    END SUBROUTINE LidarSim_CalcOutput

    !#########################################################################################################################################################################

    SUBROUTINE LidarSim_End( y, p, u, ErrStat, ErrMsg)
	
    IMPLICIT                                 NONE    
    CHARACTER(*),              				   PARAMETER       :: RoutineName="LidarSim_End"
	
    TYPE(LidarSim_InputType),              	INTENT(INOUT)   ::  u           !< Input data for initialization
    TYPE(LidarSim_ParameterType),          	INTENT(INOUT)   ::  p           !< Parameters
    TYPE(LidarSim_OutputType),             	INTENT(INOUT)   ::  y           !< Output data
    
    ! Error Handling
    INTEGER( IntKi ),                        INTENT(  OUT)   :: ErrStat      !< error status
    CHARACTER(*),                            INTENT(  OUT)   :: ErrMsg       !< error message
    
    
    ErrStat = ErrID_None
    ErrMsg = ""
    
    CALL LidarSim_DestroyOutput(y, ErrStat, ErrMsg)
    CALL LidarSim_DestroyParam(p, ErrStat, ErrMsg)
	
    END SUBROUTINE LidarSim_End
    
    !#########################################################################################################################################################################
    
    END MODULE LidarSim