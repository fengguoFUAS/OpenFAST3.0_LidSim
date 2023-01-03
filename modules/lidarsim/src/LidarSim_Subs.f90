    MODULE LidarSim_Subs

    USE LidarSim_Types
    USE NWTC_Library
    !USE InflowWind_Subs
    USE InflowWind_Types   !needed to deliver the flow field from InflowWind module during initialization, only in the frozen case 
    !USE IfW_BladedFFWind
    !USE AeroDyn14_Types
    

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC  ::  LidarSim_ReadInputFile
    PUBLIC  ::  LidarSim_InitMeasuringPoints_Cartesian
    PUBLIC  ::  LidarSim_TransformLidarToInertial
    PUBLIC  ::  LidarSim_InitMeasuringPoints_Spherical
    PUBLIC  ::  LidarSim_CreateRotationMatrix
    PUBLIC  ::  LidarSim_InitializeWeightingGauss
    PUBLIC  ::  LidarSim_InitializeWeightingManual
    PUBLIC  ::  LidarSim_InitializeBladeBlockage
    PUBLIC  ::  LidarSim_CalculateVlos
    PUBLIC  ::  LidarSim_InitializeOutputs
    PUBLIC  ::  LidarSim_SetOutputs
    PUBLIC  ::  LidarSim_CalculateIMU
    PUBLIC  ::  LidarSim_InitializeWindEvolution
    PUBLIC  ::  LidarSim_CheckBladeBlockage_14
    PUBLIC  ::  LidarSim_CheckBladeBlockage_15
    PUBLIC  ::  LidarSim_InitializeAvailability
    PUBLIC  ::  LidarSim_InitializeFrozenflow
    PUBLIC  ::  LidarSim_CalculateUVW
    PUBLIC  ::  LidSim_TS_Bladed_FFWind_CalcOutput
    
    CONTAINS
    
	
!#########################################################################################################################################################################
    
    SUBROUTINE LidarSim_ReadInputFile(InputInitFile, EchoFileName, InputFileData, ErrStat, ErrMsg)
	
    IMPLICIT                                NONE
    CHARACTER(*),                           PARAMETER       ::  RoutineName="LidarSim_ReadInputFile"
    
    CHARACTER(1024),                        INTENT(IN   )   ::  InputInitFile       !< Name of the Input File
    CHARACTER(*),                           INTENT(IN   )   ::  EchoFileName        !< name of the echo file 
    TYPE(LidarSim_InputFile),               INTENT(INOUT)   ::  InputFileData
    INTEGER(IntKi),                         INTENT(  OUT)   ::  ErrStat             !< Error status of the operation
    CHARACTER(*),                           INTENT(  OUT)   ::  ErrMsg              !< Error message if ErrStat /= ErrID_None
    
    ! Local variables
    REAL(ReKi)                                              ::  RotationAngle       !< Variable to temporary store the rotation angle (roll, pitch or yaw)
    REAL(ReKi)                                              ::  Rotations(3,3,3)    !< DCMs for roll, pitch and yaw
    INTEGER(IntKi)                                          ::  UnitInput           !< Unit number for the input file
    INTEGER(IntKi)                                          ::  TemporaryFileUnit   !< Unit number for the same input file. opened twice to check for out commented variables
    INTEGER(IntKi)                                          ::  UnitEcho            !< The local unit number for this module's echo file
    LOGICAL                                                 ::  CommentLine         !< True if line is commented out
    CHARACTER(1024)                                         ::  ReadLine            !< Temporary variable to Read a Line in
    INTEGER(IntKi)                                          ::  CounterNumberOfPoints_Cartesian !< Loop counter for the cartesian coordinates
    INTEGER(IntKi)                                          ::  CounterNumberOfPoints_Spherical !< Loop counter for the spherical coordinates    
    INTEGER(IntKi)                                          ::  CounterNumberManualWeighting    !< Loop counter for the manual weighting points
    INTEGER(IntKi)                                          ::  TmpErrStat
    CHARACTER(ErrMsgLen)                                    ::  TmpErrMsg           !< temporary error message
    INTEGER(IntKi)                                          ::  ErrStatIO           !< temporary error for read commands
    
    ! Initialization 
    ErrStat        =  0
    ErrMsg         =  ""
    UnitEcho       = -1
    
    ! Allocate OutList space
    CALL AllocAry( InputFileData%OutList, 18, "InflowWind Input File's OutList", TmpErrStat, TmpErrMsg ) !Max additional output parameters = 18
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
 
    ! Allocate evolution filename space
    !CALL AllocAry(InputFileData%EvolutionFilenameRoot, 1,"Rootname of the up-stream full-field wind file to use (.wnd, .sum)", TmpErrStat, TmpErrMsg ) 
    !CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    !IF (ErrStat >= AbortErrLev) THEN
    !    CALL Cleanup()
    !    RETURN
    !END IF

	
    !-------------------------------------------------------------------------------------------------
    ! Open the file
    !-------------------------------------------------------------------------------------------------

    CALL GetNewUnit( UnitInput, TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL OpenFInpFile( UnitInput, TRIM(InputInitFile), TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    CALL GetNewUnit( TemporaryFileUnit, TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL OpenFInpFile( TemporaryFileUnit, TRIM(InputInitFile), TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
	
    !-------------------------------------------------------------------------------------------------
    ! File header
    !-------------------------------------------------------------------------------------------------
	
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadCom( UnitInput, InputInitFile, 'Lidar version', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
  
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg)
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    CALL ReadCom( UnitInput, InputInitFile, 'description', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadCom( UnitInput, InputInitFile, '-------------', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF

	
    !-------------------------------------------------------------------------------------------------
    ! General settings
    !-------------------------------------------------------------------------------------------------    
	
    ! Echo on/off
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%Echo, 'Echo', ' Echo input data to <RootName>.ech (flag)', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    IF ( InputFileData%Echo ) THEN
        CALL OpenEcho ( UnitEcho, TRIM(EchoFileName), TmpErrStat, TmpErrMsg )
        CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
        IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp()
            RETURN
        END IF
    
        REWIND(UnitInput)
        REWIND(TemporaryFileUnit)
    
        CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
        CALL ReadCom( UnitInput, InputInitFile, 'Lidar version', TmpErrStat, TmpErrMsg )
        CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    
        CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
        CALL ReadCom( UnitInput, InputInitFile, 'description', TmpErrStat, TmpErrMsg )
        CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    
        CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
        CALL ReadCom( UnitInput, InputInitFile, '-------------', TmpErrStat, TmpErrMsg )
        CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    
        CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
        CALL ReadVar ( UnitInput, InputInitFile, InputFileData%Echo, 'Echo', 'Echo input data to <RootName>.ech (flag)', TmpErrStat, TmpErrMsg)
        CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    
    END IF
    
    ! MAXDLLChainOutputs
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%MAXDLLChainOutputs, 'MAXDLLChainOutputs', 'Number of entries in the avrSWAP reserved for the DLL chain', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
	
    !-------------------------------------------------------------------------------------------------
    ! Lidar Configuration
    !-------------------------------------------------------------------------------------------------
	
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadCom( UnitInput, InputInitFile, 'Lidar input file separator line', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF    
    
    ! Trajectory Type
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName ) 
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%TrajectoryType, 'TrajectoryType', 'Switch : {0 = cartesian coordinates; 1 = spherical coordinates}', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    ! Weighting Type
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%WeightingType, 'WeightingType', 'Switch : {0 = single point; 1 = gaussian distribution}', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
        
    ! Position of the lidar system
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%LidarPositionX_N, 'LidarPositionX_N', 'Position of the lidar system (X coordinate) [m]', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%LidarPositionY_N, 'LidarPositionY_N', 'Position of the lidar system (Y coordinate) [m]', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%LidarPositionZ_N, 'LidarPositionZ_N', 'Position of the lidar system (Z coordinate) [m]', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    ! Rotation of the lidar system    
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%RollAngle_N, 'Roll_N', 'Roll angle between the Nacelle and the lidar coordinate system', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    InputFileData%RollAngle_N = InputFileData%RollAngle_N * (Pi_D/180)    
    
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%PitchAngle_N, 'Pitch_N', 'Pitch angle between the Nacelle and the lidar coordinate system', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    InputFileData%PitchAngle_N = InputFileData%PitchAngle_N * (Pi_D/180)

    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%YawAngle_N, 'Yaw_N', 'Yaw Pitch angle between the Nacelle and the lidar coordinate system', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    InputFileData%YawAngle_N = InputFileData%YawAngle_N * (Pi_D/180)
    
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%URef, 'URef', 'Mean u-component wind speed at the reference height', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%GatesPerBeam, 'GatesPerBeam', 'Amount of gates per point', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
	
    !-------------------------------------------------------------------------------------------------
    ! Measurement settings
    !-------------------------------------------------------------------------------------------------
	
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%t_measurement_interval, 't_measurement_interval', 'Time between each measurement [s]', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
       
	
	 !-------------------------------------------------------------------------------------------------
    ! Cartesian coordinates
    !-------------------------------------------------------------------------------------------------
	
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadCom( UnitInput, InputInitFile, 'cartesian coordinates', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    ! Number of cartesian points
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%NumberOfPoints_Cartesian, 'NumberOfPoints_Cartesian', 'Amount of Points [-]', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    ! Table header
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadCom( UnitInput, InputInitFile, 'Table Header', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    CALL AllocAry( InputFileData%X_Cartesian_L, InputFileData%NumberOfPoints_Cartesian, 'X Coordinate', TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( InputFileData%Y_Cartesian_L, InputFileData%NumberOfPoints_Cartesian, 'Y Coordinate', TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( InputFileData%Z_Cartesian_L, InputFileData%NumberOfPoints_Cartesian, 'Z Coordinate', TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    
    ! Loop through table
    DO CounterNumberOfPoints_Cartesian = 1,InputFileData%NumberOfPoints_Cartesian
        CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
        CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
        IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
        END IF
        READ(UnitInput,*,IOSTAT=ErrStatIO) InputFileData%X_Cartesian_L(CounterNumberOfPoints_Cartesian),&
        InputFileData%Y_Cartesian_L(CounterNumberOfPoints_Cartesian),InputFileData%Z_Cartesian_L(CounterNumberOfPoints_Cartesian)
        IF( ErrStatIO > 0 ) THEN
            CALL SetErrStat(ErrID_Fatal,'Error reading cartesian coordinates',ErrStat,ErrMsg,RoutineName)
            CALL Cleanup()
            RETURN
        END IF
    END DO
    
	 
    !-------------------------------------------------------------------------------------------------
    ! Spherical coordinates
    !-------------------------------------------------------------------------------------------------
	
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadCom( UnitInput, InputInitFile, 'spherical coordinates', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    ! Number of spherical points
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%NumberOfPoints_Spherical, 'NumberOfPoints_Spherical', 'Amount of Points [-]', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
        
    ! Table header
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadCom( UnitInput, InputInitFile, 'Table Header', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    CALL AllocAry( InputFileData%Azimuth, InputFileData%NumberOfPoints_Spherical, 'Azimuth', TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( InputFileData%Elevation, InputFileData%NumberOfPoints_Spherical, 'Elevation', TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( InputFileData%Range, InputFileData%NumberOfPoints_Spherical,InputFileData%GatesPerBeam , 'Range', TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    
    ! Loop through spherical data table
    DO CounterNumberOfPoints_Spherical = 1, InputFileData%NumberOfPoints_Spherical
        ! Reading the number, azimuth and elevation in the corresponding variables.
		  ! The left over parameters (1..n) are read into the range variable
        ! This allows multiple range measurements with the same azimuth / elevation
        CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
        CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
        IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
        END IF
        READ(UnitInput,*,IOSTAT=ErrStatIO) InputFileData%Azimuth(CounterNumberOfPoints_Spherical),&
        InputFileData%Elevation(CounterNumberOfPoints_Spherical),InputFileData%Range(CounterNumberOfPoints_Spherical,:) 
        IF( ErrStatIO > 0 ) THEN
            CALL SetErrStat(ErrID_Fatal,'Error reading spherical coordinates',ErrStat,ErrMsg,RoutineName)
            CALL Cleanup()
            RETURN
        END IF
    END DO
    InputFileData%Azimuth   = InputFileData%Azimuth * (Pi_D/180)
    InputFileData%Elevation = InputFileData%Elevation * (Pi_D/180)
	
    
    !-------------------------------------------------------------------------------------------------
    ! Gaussian distribution
    !-------------------------------------------------------------------------------------------------    
	
    ! Description
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadCom( UnitInput, InputInitFile, 'Weighting function Gauss', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF

    ! FWHM
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%FWHM, 'FWHM', 'Width of half maximum', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    ! Number of points to evaluate
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%PointsToEvaluate, 'PointsToEvaluate', 'points evaluated to "integrate" (odd number so there is a point in the peak', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    
    !-------------------------------------------------------------------------------------------------
    ! Manual distribution
    !-------------------------------------------------------------------------------------------------    
	
    ! Description
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadCom( UnitInput, InputInitFile, 'Weighting function Manual', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    ! Number of manual weighting points
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%ManualWeightingPoints, 'ManualWeightingPoints', 'Number manual weightingpoints', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    ! Table header
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadCom( UnitInput, InputInitFile, 'Table header', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    CALL AllocAry( InputFileData%ManualWeightingDistance, InputFileData%ManualWeightingPoints, 'ManualWeightingDistance', TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( InputFileData%ManualWeighting, InputFileData%ManualWeightingPoints, 'ManualWeighting', TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    
    ! Loop through manual weighting data table
    DO CounterNumberManualWeighting = 1, InputFileData%ManualWeightingPoints
        CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
        CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
        IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
        END IF
        READ(UnitInput,*,IOSTAT=ErrStatIO) InputFileData%ManualWeightingDistance(CounterNumberManualWeighting), InputFileData%ManualWeighting(CounterNumberManualWeighting)
        IF( ErrStatIO > 0 ) THEN
            CALL SetErrStat(ErrID_Fatal,'Error reading manual weighting',ErrStat,ErrMsg,RoutineName)
            CALL Cleanup()
            RETURN
        END IF
    END DO
    
    
    ! Wind Evolution
    
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadCom( UnitInput, InputInitFile, 'Wind Evolution', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    ! flag whether we consider evolution
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%EvolutionFlag, 'EvolutionFlag', 'Flag whether evolution considered', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    ! Root file name of the evolving turbulence data
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    InputFileData%EvolutionFilenameRoot = ""
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%EvolutionFilenameRoot,&
    'EvolutionFilenameRoot', 'Rootname of the up-stream full-field wind file to use (.wnd, .sum)', TmpErrStat, TmpErrMsg,UnitEcho )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    
    ! Availability
    
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadCom( UnitInput, InputInitFile, 'Measurement availability', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    ! flag whether we consider availiability
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%AvailabilityFlag, 'AvailabilityFlag', 'Flag whether availability considered', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    
    ! Root file name of the measurement availability data
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    InputFileData%AvailabilityFilenameRoot = ""
    CALL ReadVar ( UnitInput, InputInitFile, InputFileData%AvailabilityFilenameRoot,&
    'AvailabilityFilenameRoot', 'Rootname of the lidar availability time series data', TmpErrStat, TmpErrMsg,UnitEcho )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    
    ! blade blockage
    call lidarsim_skipcomments(temporaryfileunit, unitinput, tmperrstat, tmperrmsg, unitecho)
    call seterrstat( tmperrstat, tmperrmsg, errstat, errmsg, routinename )
    if (errstat >= aborterrlev) then
        call cleanup()
        return
    end if
    call readcom( unitinput, inputinitfile, 'blade blockage', tmperrstat, tmperrmsg )
    call seterrstat( tmperrstat, tmperrmsg, errstat, errmsg, routinename )
    if (errstat >= aborterrlev) then
        call cleanup()
        return
    end if
    
    ! flag whether we consider blade blockage
    call lidarsim_skipcomments(temporaryfileunit, unitinput, tmperrstat, tmperrmsg, unitecho)
    call seterrstat( tmperrstat, tmperrmsg, errstat, errmsg, routinename )
    if (errstat >= aborterrlev) then
        call cleanup()
        return
    end if
    call readvar ( unitinput, inputinitfile, inputfiledata%bladeblockageflag, 'bladeblockageflag', 'flag whether blade blockage considered', tmperrstat, tmperrmsg )
    call seterrstat( tmperrstat, tmperrmsg, errstat, errmsg, routinename )
    if (errstat >= aborterrlev) then
        call cleanup()
        return
    end if
    
    ! flag whether we consider spinner-mounted lidar
    call lidarsim_skipcomments(temporaryfileunit, unitinput, tmperrstat, tmperrmsg, unitecho)
    call seterrstat( tmperrstat, tmperrmsg, errstat, errmsg, routinename )
    if (errstat >= aborterrlev) then
        call cleanup()
        return
    end if
    call readcom( unitinput, inputinitfile, 'Spinner Mounted', tmperrstat, tmperrmsg )
    call seterrstat( tmperrstat, tmperrmsg, errstat, errmsg, routinename )
    if (errstat >= aborterrlev) then
        call cleanup()
        return
    end if
    
    call lidarsim_skipcomments(temporaryfileunit, unitinput, tmperrstat, tmperrmsg, unitecho)
    call seterrstat( tmperrstat, tmperrmsg, errstat, errmsg, routinename )
    if (errstat >= aborterrlev) then
        call cleanup()
        return
    end if
    call readvar ( unitinput, inputinitfile, inputfiledata%SpinnerMountedFlag, 'SpinnerMountedFlag', 'flag whether spinner-mounted considered', tmperrstat, tmperrmsg )
    call seterrstat( tmperrstat, tmperrmsg, errstat, errmsg, routinename )
    if (errstat >= aborterrlev) then
        call cleanup()
        return
    end if
    
        ! flag whether nearest interpolation should be used
    call lidarsim_skipcomments(temporaryfileunit, unitinput, tmperrstat, tmperrmsg, unitecho)
    call seterrstat( tmperrstat, tmperrmsg, errstat, errmsg, routinename )
    if (errstat >= aborterrlev) then
        call cleanup()
        return
    end if
    call readcom( unitinput, inputinitfile, 'Nearest Interpolation', tmperrstat, tmperrmsg )
    call seterrstat( tmperrstat, tmperrmsg, errstat, errmsg, routinename )
    if (errstat >= aborterrlev) then
        call cleanup()
        return
    end if
    
    call lidarsim_skipcomments(temporaryfileunit, unitinput, tmperrstat, tmperrmsg, unitecho)
    call seterrstat( tmperrstat, tmperrmsg, errstat, errmsg, routinename )
    if (errstat >= aborterrlev) then
        call cleanup()
        return
    end if
    call readvar ( unitinput, inputinitfile, inputfiledata%NearestInterpFlag, 'NearestInterpFlag', 'Flag to use nearest interpolation. By default, linear interpolation is used.', tmperrstat, tmperrmsg )
    call seterrstat( tmperrstat, tmperrmsg, errstat, errmsg, routinename )
    if (errstat >= aborterrlev) then
        call cleanup()
        return
    end if
	
    ! OutList
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadCom( UnitInput, InputInitFile, 'Outlist headline', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    CALL ReadCom( UnitInput, InputInitFile, 'Outlist', TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    ! Read OutputList
    CALL LidarSim_ReadOutputList (TemporaryFileUnit, UnitInput, InputInitFile, InputFileData%OutList, InputFileData%NumOuts, 'OutList',"List of user-requested output channels", TmpErrStat, TmpErrMsg,-1 )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
    IF (ErrStat >= AbortErrLev) THEN
        CALL Cleanup()
        RETURN
    END IF
    
    CALL Cleanup()
    
    RETURN
	
    CONTAINS
	
    !-------------------------------------------------------------------------------------------------
    SUBROUTINE Cleanup()
	
    CLOSE ( UnitInput )
    CLOSE ( TemporaryFileUnit )
    IF ( InputFileData%Echo ) THEN
        CLOSE(UnitEcho)
    END IF
	
    END SUBROUTINE Cleanup
	!-------------------------------------------------------------------------------------------------
    
    END SUBROUTINE LidarSim_ReadInputFile
    
	
!#########################################################################################################################################################################
   
    SUBROUTINE LidarSim_ReadOutputList (TemporaryFileUnit, UnitInput, FileName, OutputArray, ReadNumberOutputs, VariableName, VariableDescribtion, ErrStat, ErrMsg, UnitEcho )
    
    IMPLICIT      NONE
    CHARACTER(*), PARAMETER           ::  RoutineName="LidarSim_ReadOutputList"
	
    INTEGER,      INTENT(  OUT)       :: ReadNumberOutputs                          !< Length of the array that was actually read.
    INTEGER,      INTENT(IN   )       :: TemporaryFileUnit                          !< Temporary unit for skipping comments
    INTEGER,      INTENT(IN   )       :: UnitInput                                  !< I/O unit for input file.
    INTEGER,      INTENT(IN   )       :: UnitEcho                                   !< I/O unit for echo file (if > 0).
    INTEGER,      INTENT(  OUT)       :: ErrStat                                    !< Error status
    CHARACTER(*), INTENT(  OUT)       :: ErrMsg                                     !< Error message
    CHARACTER(*), INTENT(  OUT)       :: OutputArray(:)                             !< Character array being read (calling routine dimensions it to max allowable size).
    CHARACTER(*), INTENT(IN   )       :: FileName                                   !< Name of the input file.
    CHARACTER(*), INTENT(IN   )       :: VariableDescribtion                        !< Text string describing the variable.
    CHARACTER(*), INTENT(IN   )       :: VariableName                               !< Text string containing the variable name.

    ! Local variables
    INTEGER                           :: MaxOutputs                                  ! Maximum length of the array being read
    INTEGER                           :: NumberWords                                 ! Number of words contained on a line
    INTEGER(IntKi)                    :: TmpErrStat                                  !< Temporary error status
    CHARACTER(ErrMsgLen)              :: TmpErrMsg                                   !< temporary error message
    CHARACTER(1000)                   :: OutLine                                     ! Character string read from file, containing output list
    CHARACTER(3)                      :: EndOfFile

    ! Initialization
    ErrStat = ErrID_None
    ErrMsg  = ''
    MaxOutputs  = SIZE(OutputArray)
    ReadNumberOutputs = 0
    OutputArray = ''

    ! Read in all of the lines containing output parameters and store them in OutputArray(:).
    ! The end of this list is specified with the line beginning with END.
    DO
    CALL LidarSim_SkipComments(TemporaryFileUnit, UnitInput, TmpErrStat, TmpErrMsg, UnitEcho)
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    IF (ErrStat >= AbortErrLev) THEN
        RETURN
    END IF
    CALL ReadVar ( UnitInput, FileName, OutLine, VariableName, VariableDescribtion, ErrStat, ErrMsg )
    IF ( ErrStat >= AbortErrLev ) RETURN

    EndOfFile = OutLine(1:3)            ! EndOfFile is the 1st 3 characters of OutLine
    CALL Conv2UC( EndOfFile )           ! Convert EndOfFile to upper case
    IF ( EndOfFile == 'END' )  EXIT     ! End of OutList has been reached; therefore, exit this DO

    NumberWords = CountWords( OutLine ) ! The number of words in OutLine.

    ReadNumberOutputs = ReadNumberOutputs + NumberWords  ! The total number of output channels read in so far.

    ! Check to see if the maximum # allowable in the array has been reached.
    IF ( ReadNumberOutputs > MaxOutputs )  THEN

    ErrStat = ErrID_Fatal
    ErrMsg = 'ReadOutputList:The maximum number of output channels allowed is '//TRIM( Int2LStr(MaxOutputs) )//'.'
    RETURN

    ELSE

    CALL GetWords ( OutLine, OutputArray((ReadNumberOutputs - NumberWords + 1):ReadNumberOutputs), NumberWords )

    END IF
    END DO

    RETURN
	
    END SUBROUTINE LidarSim_ReadOutputList
    
	
!#########################################################################################################################################################################    
    
    SUBROUTINE LidarSim_SkipComments(TemporaryFileUnit, UnitInput, ErrStat, ErrMsg, UnitEcho)
	
    IMPLICIT            NONE
    CHARACTER(*),       PARAMETER                   ::  RoutineName="LidarSim_SkipComments"

    INTEGER(IntKi),     INTENT(IN   )               ::  TemporaryFileUnit   !Unit number to look ahead for comments in the input file
    INTEGER(IntKi),     INTENT(IN   )               ::  UnitInput           !Unit number of the "normal" input file
    INTEGER(IntKi),     INTENT(  OUT)               ::  ErrStat             !< Error status of the operation
    CHARACTER(*),       INTENT(  OUT)               ::  ErrMsg              !< Error message if ErrStat /= ErrID_None
    INTEGER(IntKi),     INTENT(IN   ), OPTIONAL     ::  UnitEcho            !< The local unit number for this module's echo file
	
    ! Local variables
    INTEGER(IntKi)                                  ::  ErrStatIO           !Error Status of the read commands
    CHARACTER(1024)                                 ::  TemporaryRead       !string to read a line in
    LOGICAL                                         ::  Commented           ! true if line is commented, false otherwise
    
	! Initialization
    ErrStat = 0
    ErrMsg = ''
    
    READ(TemporaryFileUnit,*, IOSTAT = ErrStatIO) TemporaryRead
    IF(ErrStatIO > 0) THEN
        CALL SetErrStat(ErrID_Fatal,'Error checking for comments in the temporary input file',ErrStat,ErrMsg,RoutineName)
        return
    ENDIF
    
    IF ( PRESENT(UnitEcho) )  THEN
        IF ( UnitEcho > 0 ) &
        WRITE (UnitEcho,'(A)')  TRIM(TemporaryRead) !Writes read line to Echo file
    END IF
    
    IF(TemporaryRead(1:1) == '!') THEN
        Commented = .TRUE.
    ELSE
        Commented = .FALSE.
    ENDIF
    
    DO while( Commented == .TRUE. )
        READ(UnitInput,*,IOSTAT = ErrStatIO)    !Skip commented line in unit Input
        IF(ErrStatIO > 0) THEN
            CALL SetErrStat(ErrID_Fatal,'Error checking for comments in the original input file',ErrStat,ErrMsg,RoutineName)
            return
        ENDIF
        READ(TemporaryFileUnit,*, IOSTAT = ErrStatIO) TemporaryRead !Checks if next line is commented again
        IF(ErrStatIO > 0) THEN
            CALL SetErrStat(ErrID_Fatal,'Error checking for comments in the temporary input file',ErrStat,ErrMsg,RoutineName)
            return
        ENDIF
        
        IF ( PRESENT(UnitEcho) )  THEN
            IF ( UnitEcho > 0 ) &
            WRITE (UnitEcho,'(A)')  TRIM(TemporaryRead) !Writes read line to Echo file
        END IF
        
        IF(TemporaryRead(1:1) == '!') THEN
            Commented = .TRUE.
        ELSE
            Commented = .FALSE.
        ENDIF
    END DO
        
    END SUBROUTINE
    
	
!#########################################################################################################################################################################    
    
    SUBROUTINE LidarSim_CreateRotationMatrix(Roll_N, Pitch_N, Yaw_N, LidarOrientation_N)
	
    IMPLICIT         NONE
	 CHARACTER(*),    PARAMETER       ::  RoutineName="LidarSim_CreateRotationMatrix"
	
    REAL(ReKi), 	   INTENT(IN   )   ::  Roll_N                      !Roll Rotation
    REAL(ReKi), 	   INTENT(IN   )   ::  Pitch_N                     !Pitch Rotation
    REAL(ReKi), 	   INTENT(IN   )   ::  Yaw_N                       !Yaw Rotation
    REAL(ReKi), 	   INTENT(INOUT)   ::  LidarOrientation_N(3,3)     !Output Rotation matrix
    
    ! Local variables
    REAL(ReKi)                       ::  Rotations(3,3,3)            !Temporary Rotation matrices
    
    ! Roll Rotation
    Rotations(1,1,1) = 1
    Rotations(1,2,1) = 0
    Rotations(1,3,1) = 0
    
    Rotations(2,1,1) = 0
    Rotations(2,2,1) = COS(Roll_N)
    Rotations(2,3,1) = SIN(Roll_N)
    
    Rotations(3,1,1) = 0
    Rotations(3,2,1) = - SIN(Roll_N)
    Rotations(3,3,1) =   COS(Roll_N)
    
    ! Pitch Rotation
    Rotations(1,1,2) =   COS(Pitch_N)
    Rotations(1,2,2) = 0
    Rotations(1,3,2) = - SIN(Pitch_N)
    
    Rotations(2,1,2) = 0
    Rotations(2,2,2) = 1
    Rotations(2,3,2) = 0
    
    Rotations(3,1,2) = SIN(Pitch_N)
    Rotations(3,2,2) = 0
    Rotations(3,3,2) = COS(Pitch_N)
    
    ! Yaw Rotation
    Rotations(1,1,3) = COS(Yaw_N)
    Rotations(1,2,3) = SIN(Yaw_N)
    Rotations(1,3,3) = 0
    
    Rotations(2,1,3) = - SIN(Yaw_N)
    Rotations(2,2,3) =   COS(Yaw_N)
    Rotations(2,3,3) = 0
    
    Rotations(3,1,3) = 0
    Rotations(3,2,3) = 0
    Rotations(3,3,3) = 1
    
    ! Combining rotations
    LidarOrientation_N = MATMUL( MATMUL( Rotations(:,:,3),Rotations(:,:,2) ), Rotations(:,:,1) )
    
    END SUBROUTINE LidarSim_CreateRotationMatrix
    
	
!#########################################################################################################################################################################
    
    SUBROUTINE LidarSim_InitMeasuringPoints_Cartesian(p, InputFileData, ErrStat, ErrMsg)
	
    IMPLICIT                            NONE
    CHARACTER(*),                       PARAMETER       ::  RoutineName="LidarSim_InitMeasuringPoints_Cartesian"
	
    TYPE(LidarSim_ParameterType),       INTENT(INOUT)   ::  p                   !parameter data (destination of the InputFileData)
    TYPE(LidarSim_InputFile),           INTENT(IN   )   ::  InputFileData       !data read from the input file
    INTEGER(IntKi),                     INTENT(  OUT)   ::  ErrStat             !< Error status of the operation
    CHARACTER(*),                       INTENT(  OUT)   ::  ErrMsg              !< Error message if ErrStat /= ErrID_None
	
    ! Local variables
    INTEGER(IntKi)                                      ::  LoopCounter         !< counter to run through all cartesian coordinates data
    INTEGER(IntKi)                                      ::  TmpErrStat          !< Temporary error status
    CHARACTER(ErrMsgLen)                                ::  TmpErrMsg           !< temporary error message
    
	 ! Initialization
    ErrStat        =  0
    ErrMsg         =  ""
    
    CALL AllocAry( p%MeasuringPoints_L, 3, SIZE(InputFileData%X_Cartesian_L), 'MeasuringPoints_L', TmpErrStat, TmpErrMsg )     !Allocate the array size for n=CounterChannelNumbers measuringpositions
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( p%MeasuringPoints_Spherical_L, 3,  SIZE(InputFileData%X_Cartesian_L), 'MeasuringPoints_Spherical_L', TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    
    DO LoopCounter = 1,SIZE(InputFileData%X_Cartesian_L)
        p%MeasuringPoints_L(1,LoopCounter) = InputFileData%X_Cartesian_L(LoopCounter)
        p%MeasuringPoints_L(2,LoopCounter) = InputFileData%Y_Cartesian_L(LoopCounter)
        p%MeasuringPoints_L(3,LoopCounter) = InputFileData%Z_Cartesian_L(LoopCounter)
        p%MeasuringPoints_Spherical_L(:,LoopCounter) = LidarSim_Cartesian2Spherical(InputFileData%X_Cartesian_L(LoopCounter),InputFileData%Y_Cartesian_L(LoopCounter),InputFileData%Z_Cartesian_L(LoopCounter))
    END DO
	
    END SUBROUTINE LidarSim_InitMeasuringPoints_Cartesian
   
   
!#########################################################################################################################################################################
    
    SUBROUTINE LidarSim_InitMeasuringPoints_Spherical(p, InputFileData, ErrStat, ErrMsg)
	
    IMPLICIT                            NONE
    CHARACTER(*),                       PARAMETER       ::  RoutineName="LidarSim_InitMeasuringPoints_Spherical"

    TYPE(LidarSim_ParameterType),       INTENT(INOUT)   ::  p                       !parameter data (destination of the InputFileData)
    TYPE(LidarSim_InputFile),           INTENT(IN   )   ::  InputFileData           !data read from the input file
    INTEGER(IntKi),                     INTENT(  OUT)   ::  ErrStat                 !< Error status of the operation
    CHARACTER(*),                       INTENT(  OUT)   ::  ErrMsg                  !< Error message if ErrStat /= ErrID_None
    
    ! Local variables
    INTEGER(IntKi)                                      ::  OuterLoopCounter        !counter for looping through the coordinate data
    INTEGER(IntKi)                                      ::  InnerLoopCounter        !counter for looping through the multiple range gates
    INTEGER(IntKi)                                      ::  CounterChannelNumbers   !counts the amount of channels
    INTEGER(IntKi)                                      ::  TmpErrStat              !< Temporary error status
    CHARACTER(ErrMsgLen)                                ::  TmpErrMsg               !< temporary error message

	 ! Initialization
    ErrStat        =  0
    ErrMsg         =  ""
    
    CALL AllocAry( p%MeasuringPoints_L, 3,InputFileData%NumberOfPoints_Spherical*InputFileData%GatesPerBeam, 'MeasuringPoints_L', TmpErrStat, TmpErrMsg )     !Allocate the array size for n=CounterChannelNumbers measuringpositions
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( p%MeasuringPoints_Spherical_L, 3,InputFileData%NumberOfPoints_Spherical*InputFileData%GatesPerBeam, 'MeasuringPoints_Spherical_L', TmpErrStat, TmpErrMsg )  
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)

    OuterLoopCounter = 1
    InnerLoopCounter = 1
    CounterChannelNumbers = 1    
        
    DO WHILE( OuterLoopCounter <= InputFileData%NumberOfPoints_Spherical)
        DO WHILE( InnerLoopCounter <= InputFileData%GatesPerBeam)
                IF ( InputFileData%Range(OuterLoopCounter,InnerLoopCounter) /= 0 )THEN  !Range gates mustn't be 0. => Divide by 0 !
                    p%MeasuringPoints_L(:,CounterChannelNumbers) = &   !Transformation from the spherical to cartesian coordinates
                        LidarSim_Spherical2Cartesian(InputFileData%Azimuth(OuterLoopCounter),InputFileData%Elevation(OuterLoopCounter),InputFileData%Range(OuterLoopCounter,InnerLoopCounter))
                    p%MeasuringPoints_Spherical_L(:,CounterChannelNumbers) = &
                        (/ InputFileData%Range(OuterLoopCounter,InnerLoopCounter),InputFileData%Azimuth(OuterLoopCounter),InputFileData%Elevation(OuterLoopCounter)/)
                    CounterChannelNumbers = CounterChannelNumbers + 1
                ELSE
                    CALL SetErrStat(ErrID_Fatal,"Range gates must not be 0",ErrStat,ErrMsg,RoutineName)
                ENDIF
            InnerLoopCounter = InnerLoopCounter + 1
        END DO
        InnerLoopCounter = 1
        OuterLoopCounter = OuterLoopCounter + 1
    END DO
	
    END SUBROUTINE LidarSim_InitMeasuringPoints_Spherical

	
!#########################################################################################################################################################################
    
    FUNCTION LidarSim_Spherical2Cartesian(Azimuth, Elevation, Range)
	
    IMPLICIT        NONE
	 CHARACTER(*),   PARAMETER       ::  RoutineName="LidarSim_Spherical2Cartesian"
	
    REAL(ReKi),     INTENT(IN   )   ::  Azimuth                             !Azimuth angle
    REAL(ReKi),     INTENT(IN   )   ::  Elevation                           !Elevation angle
    REAL(ReKi),     INTENT(IN   )   ::  Range                               !range gate
    REAL(ReKi),     DIMENSION (3)   ::  LidarSim_Spherical2Cartesian        !Output : x,y,z coordinates
    
    LidarSim_Spherical2Cartesian(1)  =   Range*COS(Elevation)*COS(Azimuth)   !x
    LidarSim_Spherical2Cartesian(2)  =   Range*COS(Elevation)*SIN(Azimuth)   !y
    LidarSim_Spherical2Cartesian(3)  =   Range*SIN(Elevation)                !z
	
    END FUNCTION LidarSim_Spherical2Cartesian
    
	
!#########################################################   
    
    FUNCTION LidarSim_Cartesian2Spherical(X, Y, Z)
	
    IMPLICIT        NONE
	 CHARACTER(*),   PARAMETER       ::  RoutineName="LidarSim_Cartesian2Spherical"
	
    REAL(ReKi),     INTENT(IN   )   ::  X                             !Azimuth angle
    REAL(ReKi),     INTENT(IN   )   ::  Y                           !Elevation angle
    REAL(ReKi),     INTENT(IN   )   ::  Z                               !range gate
    REAL(ReKi),     DIMENSION (3)   ::  LidarSim_Cartesian2Spherical       !Output : x,y,z coordinates

    LidarSim_Cartesian2Spherical(1)  =   SQRT(X**2+Y**2+Z**2)   !range
    
    IF(LidarSim_Cartesian2Spherical(1) > 0 ) THEN
        LidarSim_Cartesian2Spherical(3)  =   ASIN(Z/LidarSim_Cartesian2Spherical(1))   !y
        
        LidarSim_Cartesian2Spherical(2) = ATAN2(Y,X)
    ELSE
        LidarSim_Cartesian2Spherical(2) = 0
        LidarSim_Cartesian2Spherical(3) = 0
    ENDIF
	
    END FUNCTION LidarSim_Cartesian2Spherical
    
	
!#############################################
    
    FUNCTION LidarSim_TransformLidarToInertial(AttachBodyMotion, p, MeasuringPoint_L)
	
    IMPLICIT        NONE
	 CHARACTER(*),   PARAMETER       	 ::  RoutineName="LidarSim_TransformLidarToInertial"
	
    REAL(ReKi)                          ::  LidarSim_TransformLidarToInertial(3)    !Output calculated transformation from the lidar coord. sys. to the inertial system
    TYPE(MeshType)                      ::  AttachBodyMotion                           !Data describing the motion of the nacelle coord. sys.
    TYPE(LidarSim_ParameterType)        ::  p                                       !Parameter data 
    REAL(ReKi)                          ::  MeasuringPoint_L(3)                     !point which needs to be transformed
    
    ! local variables
    REAL(ReKi)                          :: PositionNacelle_I(3)                     !local variable to save the current Nacelle position (in the inerital cord. sys.)
    
    PositionNacelle_I = AttachBodyMotion%Position(:,1) + AttachBodyMotion%TranslationDisp(:,1)
       
    LidarSim_TransformLidarToInertial = PositionNacelle_I +  MATMUL(TRANSPOSE(AttachBodyMotion%Orientation(:,:,1)),(p%LidarPosition_N + MATMUL( p%LidarOrientation_N,MeasuringPoint_L ) ) )
  
    END FUNCTION LidarSim_TransformLidarToInertial
    
	
!###########################################
    
    SUBROUTINE LidarSim_InitializeWeightingManual(p, InputFileData, ErrStat, ErrMsg)
    
    IMPLICIT                        NONE
    CHARACTER(*),                   PARAMETER       ::  RoutineName="LidarSim_InitializeWeightingManual"
    
    TYPE(LidarSim_ParameterType),   INTENT(INOUT)   ::  p                   ! parameter data to write results in
    TYPE(LidarSim_InputFile),       INTENT(IN   )   ::  InputFileData       ! Inputdata from the input file
    INTEGER(IntKi),                 INTENT(  OUT)   ::  ErrStat             !< Temporary error status
    CHARACTER(*),                   INTENT(  OUT)   ::  ErrMsg              !< temporary error message
    
    ! local variables
    INTEGER(IntKi)          ::  TmpErrStat
    CHARACTER(ErrMsgLen)    ::  TmpErrMsg           !< temporary error message
    
    TmpErrStat = 0
    TmpErrMsg = ''

    CALL AllocAry( p%WeightingDistance,SIZE(InputFileData%ManualWeightingDistance), 'p%WeightingDistance', TmpErrStat, TmpErrMsg )  !Allocating the needed size for the weighting distance vector
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( p%Weighting,SIZE(InputFileData%ManualWeighting), 'p%Weighting', TmpErrStat, TmpErrMsg )                          !Allocating the needed size for the weighting vector
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    p%WeightingDistance = InputFileData%ManualWeightingDistance                                                                     !writing the input distances in the parameter
    p%Weighting(:) = InputFileData%ManualWeighting(:) / SUM(InputFileData%ManualWeighting)                                          !writing the input weighting in the parameter and normalizing it (sum = 1)
    
    END SUBROUTINE LidarSim_InitializeWeightingManual
    
    
!##########################################  

    SUBROUTINE LidarSim_InitializeWeightingGauss(p, InputFileData, ErrStat, ErrMsg)
    
    IMPLICIT                        NONE
    CHARACTER(*),                   PARAMETER       ::  RoutineName="LidarSim_InitializeWeightingGauss"
    
    TYPE(LidarSim_ParameterType),   INTENT(INOUT)   ::  p                   ! parameter data to write results in
    TYPE(LidarSim_InputFile),       INTENT(IN   )   ::  InputFileData       ! Inputdata from the input file
    INTEGER(IntKi),                 INTENT(  OUT)   ::  ErrStat             !< Error status of the operation
    CHARACTER(*),                   INTENT(  OUT)   ::  ErrMsg              !< Error message if ErrStat /= ErrID_None
    
    ! local variables
    INTEGER(IntKi)          ::  Counter                                     !Loopcounter for every point to evaluate
    REAL(ReKi)              ::  Dist                                        !Distance between each evaluation point
    INTEGER(IntKi)          ::  TmpErrStat                                  !< Temporary error status
    CHARACTER(ErrMsgLen)    ::  TmpErrMsg                                   !< temporary error message
    
    ErrStat =   0
    ErrMsg  =   ''
    Dist = 2*InputFileData%FWHM/(InputFileData%PointsToEvaluate+1)          !Distance between each weighting point
    
    CALL AllocAry( p%WeightingDistance,InputFileData%PointsToEvaluate, 'p%WeightingDistance', TmpErrStat, TmpErrMsg )   !Allocating the needed size for the weighting distance vector
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( p%Weighting,InputFileData%PointsToEvaluate, 'p%Weighting', TmpErrStat, TmpErrMsg )                   !Allocating the needed size for the weighting vector
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    
    DO Counter=1,InputFileData%PointsToEvaluate
        p%WeightingDistance(Counter) = (Counter) * Dist + (-InputFileData%FWHM)             !Creating the distance vector
        p%Weighting(Counter) = ( (2*SQRT(LOG(2.0)))/(InputFileData%FWHM*SQRT(Pi_D)) ) *&      !Calculation of the gaussian distribution
        EXP( -((p%WeightingDistance(Counter)**2)*4*LOG(2.0))/(InputFileData%FWHM**2) ) !&
    END DO
    p%Weighting = p%Weighting / SUM( p%Weighting )

    END SUBROUTINE LidarSim_InitializeWeightingGauss
    
! #################################################
    
    SUBROUTINE LidarSim_InitializeBladeBlockage(p, InputFileData, ErrStat, ErrMsg)
    
    IMPLICIT                        NONE
    CHARACTER(*),                   PARAMETER       ::  RoutineName="LidarSim_InitializeBladeBlockage"
    
    TYPE(LidarSim_ParameterType),   INTENT(INOUT)   ::  p                   ! parameter data to write results in
    TYPE(LidarSim_InputFile),       INTENT(IN   )   ::  InputFileData       ! Inputdata from the input file
    INTEGER(IntKi),                 INTENT(  OUT)   ::  ErrStat             !< Error status of the operation
    CHARACTER(*),                   INTENT(  OUT)   ::  ErrMsg              !< Error message if ErrStat /= ErrID_None
    
    ! local variables
    INTEGER(IntKi)          ::  TmpErrStat                                  !< Temporary error status
    CHARACTER(ErrMsgLen)    ::  TmpErrMsg                                   !< temporary error message
    
    ErrStat =   0
    ErrMsg  =   ''

     P%BladeBlockageFlag	= INPUTFILEDATA%BladeBlockageFlag	! store the flag in parameter
        CALL AllocAry( P%BladeELMCL, P%NELM,'Blade chord length for lidar sim.', TmpErrStat, TmpErrMsg )
        CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
    END SUBROUTINE LidarSim_InitializeBladeBlockage
    
!##########################################
   
    ! this routine read the binary wind data of upstream positions, the code is mostly similarly to ifw_BladedFFWind module.
    ! fot the initial implementing of wind evolution only bladed style file is supported when we generate the 4D turbulence
    ! we also output the data as the bladed style
    SUBROUTINE LidarSim_InitializeWindEvolution(p, InputFileData, ErrStat, ErrMsg)
    
    IMPLICIT                        NONE
    CHARACTER(*),                   PARAMETER       ::  RoutineName="LidarSim_InitializeWindEvolution"
    
    TYPE(LidarSim_ParameterType),   INTENT(INOUT)   ::  p           ! parameter data to write results in
    TYPE(LidarSim_InputFile),       INTENT(IN   )   ::  InputFileData       ! Inputdata from the input file
    INTEGER(IntKi),                 INTENT(  OUT)   ::  ErrStat             !< Error status of the operation
    CHARACTER(*),                   INTENT(  OUT)   ::  ErrMsg              !< Error message if ErrStat /= ErrID_None
    
    ! local variables
    INTEGER(IntKi)                                  ::  TmpErrStat          !< Temporary error status
    CHARACTER(ErrMsgLen)                            ::  TmpErrMsg           !< temporary error message
    INTEGER(IntKi)                                  ::  UnitWind            ! Unit number for the InflowWind input file
    INTEGER(B2Ki)                                   ::  Dum_Int2            ! temporary variavle to store read in data from binary file
    CHARACTER( 1028 )                               ::  SumFile             ! length is LEN(ParamData%WindFileName) + the 4-character extension.
    CHARACTER( 1028 )                               ::  TwrFile             ! length is LEN(ParamData%WindFileName) + the 4-character extension.
    CHARACTER( 1028 )                               ::  binFile             ! length is LEN(ParamData%WindFileName) + the 4-character extension.

    REAL(ReKi)                                      ::  TI      (3)         ! turbulence intensities of the wind components as defined in the FF file, not necessarially the actual TI
    REAL(ReKi)                                      ::  Xvec    (30)
    REAL(ReKi)                                      ::  BinTI   (3)          ! turbulence intensities of the wind components as defined in the FF binary file, not necessarially the actual TI
    REAL(ReKi)                                      ::  UBar
    REAL(ReKi)                                      ::  ZCenter
    INTEGER(IntKi)                                  ::  I
    LOGICAL                                         ::  CWise
    LOGICAL                                         ::  Exists
    
    
    
    ErrStat =   0
    ErrMsg  =   ''
    

    TmpErrMsg   = ''
    TmpErrStat  = 0

    P%EvolutionFlag	= INPUTFILEDATA%EvolutionFlag	! store the flag in parameter
    !P%AvailabilityFlag	= INPUTFILEDATA%AvailabilityFlag	! store the flag in parameter
    !P%BladeBlockageFlag	= INPUTFILEDATA%BladeBlockageFlag	! store the flag in parameter
    
      ! Get a unit number to use

    CALL GetNewUnit(UnitWind, TmpErrStat, TmpErrMsg)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'Lidar wind evolution initialization')
    IF ( ErrStat >= AbortErrLev ) RETURN



    !...........................................................................................
    ! Create full-field binary file name from evo file root name.  Also get tower file
    ! name.
    !...........................................................................................

    CALL GetRoot(TRIM(InputFileData%EvolutionFilenameRoot), SumFile)
        
        binFile = TRIM(SumFile)//'.wnd'
        TwrFile = TRIM(SumFile)//'.twr'
        SumFile = TRIM(SumFile)//'.sum'
        
     
      !----------------------------------------------------------------------------------------------
      ! Open the binary file, read its "header" (first 2-byte integer) to determine what format
      ! binary file it is, and close it.
      !----------------------------------------------------------------------------------------------
    
    CALL OpenBInpFile (UnitWind, TRIM(binFile), TmpErrStat, TmpErrMsg)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'Lidar wind evolution initialization')
    IF ( ErrStat >= AbortErrLev ) RETURN

      ! Read the first binary integer from the file to get info on the type.
      ! Cannot use library read routines since this is a 2-byte integer.
    READ ( UnitWind, IOSTAT=TmpErrStat )  Dum_Int2
    CLOSE( UnitWind )

    IF (TmpErrStat /= 0) THEN
       CALL SetErrStat( ErrID_Fatal, ' Error reading first binary integer from file "'//TRIM(InputFileData%EvolutionFilenameRoot)//'."', &
            ErrStat, ErrMsg, 'Lidar wind evolution initialization')
       RETURN
    ENDIF
    



    !...........................................................................................
    ! Read the summary file to get necessary scaling information
    !...........................................................................................

    CALL Read_Summary_FF_evo (UnitWind, TRIM(SumFile), CWise, ZCenter, TI, TmpErrStat, TmpErrMsg )
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
        IF ( ErrStat >= AbortErrLev ) THEN
            CLOSE ( UnitWind )
            RETURN
        END IF
        
     UBar = p%MeanFFWS      ! temporary storage .... this is our only check to see if the summary and binary files "match"


         !...........................................................................................
         ! Open the binary file and read its header
         !...........................................................................................

            CALL OpenBInpFile (UnitWind, TRIM(binFile), TmpErrStat, TmpErrMsg )
               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
               IF ( ErrStat >= AbortErrLev ) THEN
                  CLOSE ( UnitWind )
                  RETURN
               END IF
            
            IF ( Dum_Int2 == -99 ) THEN                                                      ! Newer-style BLADED format
               CALL Read_Bladed_FF_Header1_evo (UnitWind, BinTI, TmpErrStat, TmpErrMsg)
                  CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
                  IF ( ErrStat >= AbortErrLev ) THEN
                     CLOSE ( UnitWind )
                     RETURN
                  END IF

                  ! If the TIs are also in the binary file (BinTI > 0),
                  ! use those numbers instead of ones from the summary file

               DO I =1,p%NFFComp
                  IF ( BinTI(I) > 0 ) TI(I) = BinTI(I)
               ENDDO

            ELSE
               CALL Read_Bladed_FF_Header0_evo (UnitWind, TmpErrStat, TmpErrMsg)     ! Older-style BLADED format
                  CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
                  IF ( ErrStat >= AbortErrLev ) THEN
                     CLOSE ( UnitWind )
                     RETURN
                  END IF

            ENDIF



         !...........................................................................................
         ! Let's see if the summary and binary FF wind files go together before continuing.
         !...........................................................................................

            IF ( ABS( UBar - p%MeanFFWS ) > 0.1 )  THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error: Incompatible mean hub-height wind speeds in FF wind files. '//&
                           '(Check that the .sum and .wnd files were generated together.)', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ENDIF


         !...........................................................................................
         ! Calculate the height of the bottom of the grid
         !...........................................................................................

            p%GridBase = ZCenter - p%FFZHWid         ! the location, in meters, of the bottom of the grid

          
            
            
         !...........................................................................................
         ! Close the wnd binary file, and now read the binary grids for evolution (converted to m/s) and close the file
         !...........................................................................................
            
            CLOSE ( UnitWind )
            
            CALL OpenBInpFile (UnitWind, TRIM(INPUTFILEDATA%EVOLUTIONFILENAMEROOT), TmpErrStat, TmpErrMsg )
               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
               IF ( ErrStat >= AbortErrLev ) THEN
                  CLOSE ( UnitWind )
                  RETURN
               END IF
            
            CALL Read_Bladed_Grids_evo( CWise, TI, Xvec,TmpErrStat, TmpErrMsg)
            CLOSE ( UnitWind )
            
               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
               IF ( ErrStat >= AbortErrLev ) RETURN
               
   
               !  calculate the simulation time and initial positions
   IF (p%Periodic) THEN
      p%InitXPosition = 0                ! start at the hub
      p%TotalTime     = p%NFFSteps*p%FFDTime
   ELSE
      p%InitXPosition = p%FFYHWid          ! start half the grid width ahead of the turbine
      p%TotalTime     = (p%NFFSteps-1)*p%FFDTime
   ENDIF    
      
    
    RETURN
    
    CONTAINS
    
  !====================================================================================================
   !>   Reads the binary headers from the turbulence files of the old Bladed variety.  Note that
   !!   because of the normalization, neither ParamData%NZGrids or ParamData%NYGrids are larger than 32 points.
   !!   21-Sep-2009 - B. Jonkman, NREL/NWTC.
   !!   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   !!   11-Nov-2020   
   SUBROUTINE Read_Bladed_FF_Header0_evo (UnitWind, ErrStat, ErrMsg)

      IMPLICIT                                              NONE

      CHARACTER(*),        PARAMETER                     :: RoutineName="Read_Bladed_FF_Header0_evo"

         ! Passed Variables:

      INTEGER(IntKi),                     INTENT(IN   )  :: UnitWind  !< unit number of already-opened wind file
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat   !< error status 
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg    !< error message 


         ! Local Variables:
      REAL(ReKi)                                         :: FFXDelt
      REAL(ReKi)                                         :: FFYDelt
      REAL(ReKi)                                         :: FFZDelt

      INTEGER(B2Ki)                                      :: Dum_Int2
      INTEGER(IntKi)                                     :: I


         ! Temporary Error Handling
      INTEGER(IntKi)                                     :: TmpErrStat     ! for checking the IOSTAT from a READ or Open statement
      CHARACTER(ErrMsgLen)                               :: TmpErrMsg      ! Temporary ErrMsg


      !-------------------------------------------------------------------------------------------------
      ! Initializations
      !-------------------------------------------------------------------------------------------------

      ErrStat  = ErrID_None
      ErrMsg   = ''


      !-------------------------------------------------------------------------------------------------
      ! Read the header (file has just been opened)
      !-------------------------------------------------------------------------------------------------

         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! -NFFC (file ID)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading number of wind components from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         p%NFFComp = -1*Dum_Int2


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! delta z (mm)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading dz from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         FFZDelt = 0.001*Dum_Int2
         p%InvFFZD = 1.0/FFZDelt


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! delta y (mm)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading dy from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         FFYDelt = 0.001*Dum_Int2
         p%InvFFYD = 1.0/FFYDelt


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! delta x (mm)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading dx from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         FFXDelt = 0.001*Dum_Int2


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! half the number of time steps

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading number of time steps from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         p%NFFSteps = 2*Dum_Int2


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! 10 times the mean full-field wind speed

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading mean full-field wind speed from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         p%MeanFFWS = 0.1*Dum_Int2
         p%InvMFFWS = 1.0/p%MeanFFWS
         p%FFDTime  = FFXDelt/p%MeanFFWS
         p%FFRate   = 1.0/p%FFDTime


      DO I = 1,5

         ! Read 2-byte integer. Can't use library routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                              ! unused variables: zLu, yLu, xLu, dummy, random seed

            IF (TmpErrStat /= 0) THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading 2-byte integers from binary FF file.', ErrStat, ErrMsg, RoutineName)
               RETURN
            ENDIF

      END DO


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! 1000*nz

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading nz from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         p%NZGrids  = Dum_Int2/1000
         p%FFZHWid  = 0.5*FFZDelt*( p%NZGrids - 1 )    ! half the vertical size of the grid


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! 1000*ny

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading ny from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         p%NYGrids  = Dum_Int2/1000
         p%FFYHWid  = 0.5*FFYDelt*( p%NYGrids - 1 )


      IF (p%NFFComp == 3) THEN

         DO I=1,6

               ! Read 2-byte integer. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                           ! unused variables: zLv, yLv, xLv, zLw, yLw, xLw

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 2-byte length scales from binary FF file.', ErrStat, ErrMsg, RoutineName)
                  RETURN
               ENDIF

         ENDDO !I

      ENDIF !NFFComp


      RETURN

   END SUBROUTINE Read_Bladed_FF_Header0_evo
   !====================================================================================================
   !>   Reads the binary headers from the turbulence files of the new Bladed variety.
   !!   16-May-2002 - Windward Engineering.
   !!   21-Sep-2009 - B. Jonkman, NREL.  updated to trap errors and add extra parameters for MANN model
   !!   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   SUBROUTINE Read_Bladed_FF_Header1_evo (UnitWind, TI, ErrStat, ErrMsg)

      IMPLICIT                                              NONE

      CHARACTER(*),        PARAMETER                     :: RoutineName="Read_Bladed_FF_Header1_evo"


         ! Passed Variables:

      INTEGER(IntKi),                     INTENT(IN   )  :: UnitWind  !< unit number of already-opened wind file
      REAL(ReKi),                         INTENT(  OUT)  :: TI(3)     !< turbulence intensity contained in file header 
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat   !< error status
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg    !< error message


         ! Local Variables:

      REAL(ReKi)                                         :: FFXDelt
      REAL(ReKi)                                         :: FFYDelt
      REAL(ReKi)                                         :: FFZDelt

      REAL(SiKi)                                         :: Dum_Real4
      INTEGER(B2Ki)                                      :: Dum_Int2
      INTEGER(B4Ki)                                      :: Dum_Int4

      INTEGER(IntKi)                                     :: I
      INTEGER(IntKi)                                     :: TurbType


         ! Temporary Error Handling
      INTEGER(IntKi)                                     :: TmpErrStat
      CHARACTER(ErrMsgLen)                               :: TmpErrMsg


      !-------------------------------------------------------------------------------------------------
      ! Initializations
      !-------------------------------------------------------------------------------------------------

      ErrStat  = ErrID_None
      ErrMsg   = ''

      TI(:) = -1                                                                                !Initialize to -1 (not all models contain TI)

      !-------------------------------------------------------------------------------------------------
      ! File reading
      !-------------------------------------------------------------------------------------------------

         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! -99 (file ID)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading integer from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF


         ! Read 2-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int2                                                 ! turbulence type

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading turbulence type from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         TurbType = Dum_Int2


      SELECT CASE (TurbType)
         CASE(1, 2)
            !----------------------------------------
            !1-component Von Karman (1) or Kaimal (2)
            !----------------------------------------
               p%NFFComp = 1

         CASE(3, 5)
            !----------------------------------------
            !3-component Von Karman (3) or IEC-2
            ! Kaimal (5)
            !----------------------------------------
               p%NFFComp = 3

         CASE(4)
            !----------------------------------------
            !improved Von Karman
            !----------------------------------------

                  ! Read 2-byte integer. Can't use library routines for this.
               READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                        ! number of components (should be 3)

                  IF (TmpErrStat /= 0) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading number of components from binary FF file.', ErrStat, ErrMsg, RoutineName)
                     RETURN
                  ENDIF
                  p%NFFComp = Dum_Int4

                  ! Read 4-byte real. Can't use library routines for this.
               READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                       ! Latitude (deg)

                  IF (TmpErrStat /= 0) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading latitude from binary FF file.', ErrStat, ErrMsg, RoutineName)
                     RETURN
                  ENDIF

                  ! Read 4-byte real. Can't use library routines for this.
               READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                       ! Roughness length (m)

                  IF (TmpErrStat /= 0) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading roughness length from binary FF file.', ErrStat, ErrMsg, RoutineName)
                     RETURN
                  ENDIF

                  ! Read 4-byte real. Can't use library routines for this.
               READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                       ! Reference height (m) = Z(1) + GridHeight / 2.0

                  IF (TmpErrStat /= 0) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading reference height from binary FF file.', ErrStat, ErrMsg, RoutineName)
                     RETURN
                  ENDIF


               DO I = 1,3
                     ! Read 4-byte real. Can't use library routines for this.
                  READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                    ! TI(u, v, w) (%)

                     IF (TmpErrStat /= 0) THEN
                        CALL SetErrStat( ErrID_Fatal, ' Error reading TI('//'TRIM(Num2LStr(I))'//') from binary FF file.', ErrStat, ErrMsg, RoutineName)
                        RETURN
                     ENDIF
                     TI(I) = Dum_Real4                                                          ! This overwrites the TI read in the summary file

               END DO !I


         CASE (7, 8)
            !----------------------------------------
            ! General Kaimal (7) or  Mann model (8)
            !----------------------------------------

                  ! Read 4-byte integer. Can't use library routines for this.
               READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                        ! number of bytes in header

                  IF (TmpErrStat /= 0) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading number of header records from binary FF file.', ErrStat, ErrMsg, RoutineName)
                     RETURN
                  ENDIF

                  ! Read 4-byte integer. Can't use library routines for this.
               READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                        ! number of components

                  IF (TmpErrStat /= 0) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading number of data from binary FF file.', ErrStat, ErrMsg, RoutineName)
                     RETURN
                  ENDIF
                  p%NFFComp = Dum_Int4


         CASE DEFAULT

            CALL SetErrStat( ErrID_Warn, ' InflowWind does not recognize the full-field turbulence file type ='// &
                        TRIM(Num2LStr(TurbType))//'.', ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN

      END SELECT !TurbType


         ! Read 4-byte real. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                                ! delta z (m)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading dz from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         FFZDelt = Dum_Real4
         p%InvFFZD = 1.0/FFZDelt


         ! Read 4-byte real. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                               ! delta y (m)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading dy from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         FFYDelt = Dum_Real4
         p%InvFFYD = 1.0/FFYDelt

         ! Read 4-byte real. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                               ! delta x (m)

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading dx from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         FFXDelt = Dum_Real4


         ! Read 4-byte integer. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                                ! half the number of time steps

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading number of time steps from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         p%NFFSteps = 2*Dum_Int4


         ! Read 4-byte real. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                               ! mean full-field wind speed

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading mean full-field wind speed from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         p%MeanFFWS = Dum_Real4
         p%InvMFFWS = 1.0/p%MeanFFWS
         p%FFDTime  = FFXDelt/p%MeanFFWS
         p%FFRate   = 1.0/p%FFDTime


      DO I = 1,3

            ! Read 4-byte real. Can't use library routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                            ! unused variables: zLu, yLu, xLu

            IF (TmpErrStat /= 0) THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte length scales from binary FF file.', ErrStat, ErrMsg, RoutineName)
               RETURN
            ENDIF

      END DO


      DO I = 1,2

         ! Read 4-byte integer. Can't use library routines for this.
         READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                             ! unused variables: dummy, random seed

            IF (TmpErrStat /= 0) THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte integers from binary FF file.', ErrStat, ErrMsg, RoutineName)
               RETURN
            ENDIF

      END DO


         ! Read 4-integer real. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                                ! nz

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading nz from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         p%NZGrids  = Dum_Int4
         p%FFZHWid  = 0.5*FFZDelt*( p%NZGrids - 1 )    ! half the vertical size of the grid


         ! Read 4-integer real. Can't use library routines for this.
      READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                                ! ny

         IF (TmpErrStat /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, ' Error reading ny from binary FF file.', ErrStat, ErrMsg, RoutineName)
            RETURN
         ENDIF
         p%NYGrids  = Dum_Int4
         p%FFYHWid  = 0.5*FFYDelt*( p%NYGrids - 1 )


      IF (p%NFFComp == 3) THEN

         DO I=1,6

               ! Read 4-real real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                         ! unused variables: zLv, yLv, xLv, zLw, yLw, xLw

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte length scales from binary FF file.', ErrStat, ErrMsg, RoutineName)
                  RETURN
               ENDIF

         ENDDO !I

      ENDIF !NFFComp



      IF ( TurbType == 7 ) THEN     ! General Kaimal model

               ! Read 4-real real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                         ! unused variable: coherence decay constant

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading coherence decay constant from binary FF file.', ErrStat, ErrMsg, RoutineName)
                  RETURN
               ENDIF

               ! Read 4-real real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                         ! unused variables: coherence scale parameter in m

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading coherence scale parameter from binary FF file.', ErrStat, ErrMsg, RoutineName)
                  RETURN
               ENDIF

      ELSE IF ( TurbType == 8 ) THEN     ! Mann model

         DO I=1,2

               ! Read 4-real real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                         ! unused variables: shear parameter (gamma), scale length

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte parameters from binary FF file.', ErrStat, ErrMsg, RoutineName)
                  RETURN
               ENDIF

         ENDDO !I

         DO I=1,4

               ! Read 4-real real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                         ! unused variables

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte parameters from binary FF file.', ErrStat, ErrMsg, RoutineName)
                  RETURN
               ENDIF

         ENDDO !I

         DO I=1,3

               ! Read 4-integer real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                          ! unused variables

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte parameters from binary FF file.', ErrStat, ErrMsg, RoutineName)
                  RETURN
               ENDIF

         ENDDO !I

         DO I=1,2

               ! Read 4-real real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                         ! unused variables

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte parameters from binary FF file.', ErrStat, ErrMsg, RoutineName)
                  RETURN
               ENDIF

         ENDDO !I

         DO I=1,3

               ! Read 4-integer real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Int4                                          ! unused variables

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte parameters from binary FF file.', ErrStat, ErrMsg, RoutineName)
                  RETURN
               ENDIF

         ENDDO !I

         DO I=1,2

               ! Read 4-real real. Can't use library routines for this.
            READ (UnitWind, IOSTAT=TmpErrStat)   Dum_Real4                                         ! unused variables

               IF (TmpErrStat /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading 4-byte parameters from binary FF file.', ErrStat, ErrMsg, RoutineName)
                  RETURN
               ENDIF

         ENDDO !I


      ENDIF !TurbType


      RETURN

   END SUBROUTINE Read_Bladed_FF_Header1_evo
   
   !====================================================================================================
   !> This subroutine continues reading UnitWind, starting after the headers have been read.
   !! It reads the Grids and converts the data to un-normalized wind speeds in m/s.
   !!   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   !   modified 14-Nov- 2020 F. Guo Flensburg University of Applied Sciences  some variables are not sotred in type name changed by adding _evo
   SUBROUTINE Read_Bladed_Grids_evo (CWise, TI, Xvec, ErrStat, ErrMsg )
      IMPLICIT                                              NONE

      CHARACTER(*),        PARAMETER                     :: RoutineName="Read_Bladed_Grids_evo"

         ! Passed variables


      LOGICAL,                            INTENT(IN   )  :: CWise          !< clockwise flag (determines if y is increasing or decreasing in file)
      REAL(ReKi),                         INTENT(IN   )  :: TI      (3)    !< turbulence intensities of the wind components as defined in the FF file, not necessarially the actual TI
      REAL(ReKi),                         INTENT(  OUT)  :: Xvec    (30)    !< a vector records the x position of unfrozen planes
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat        !< error status
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg         !< error message 

      REAL(ReKi),    PARAMETER                           :: FF_Offset(3) = (/ 1.0, 0.0, 0.0 /)  ! used for "un-normalizing" the data

      INTEGER(IntKi)                                     :: CFirst
      INTEGER(IntKi)                                     :: CLast
      INTEGER(IntKi)                                     :: CStep
      INTEGER(B2Ki)                                      :: Dum_Int2
      INTEGER(IntKi)                                     :: I
      INTEGER(IntKi)                                     :: IC
      INTEGER(IntKi)                                     :: IR
      INTEGER(IntKi)                                     :: IT
      INTEGER(IntKi)                                     :: NL    !number of unfrozen planes along x
      INTEGER(IntKi)                                     :: IL    !index of unfrozen planes along x
      INTEGER(IntKi)                                     :: TmpNumSteps

         ! Temporary variables for error handling

      INTEGER(IntKi)                                     :: TmpErrStat     ! for checking the result of IOSTAT on READ or Open statements
      CHARACTER(ErrMsgLen)                               :: TmpErrMsg

      ! read in the number of unfrozen planes and record the specific x positions
      READ (UnitWind,IOStat=TmpErrStat)  Dum_Int2
      NL = Dum_Int2
      
      
      !-------------------------------------------------------------------------------------------------
      ! Generate an informative message. Initialize the ErrStat.
      !-------------------------------------------------------------------------------------------------
         ! This could take a while, so we'll write a message to tell users what's going on:
         
      CALL WrScr( NewLine//'   Reading a '//TRIM(Num2LStr(NL))//'x'//TRIM( Num2LStr(p%NYGrids) )//'x'//TRIM( Num2LStr(p%NZGrids) )//  &
                  ' grid to include wind evolution. ' )
      ErrMsg   = ""
      ErrStat  =  ErrID_None
            
      !-------------------------------------------------------------------------------------------------
      ! read the x position of unfrozen planes and the corresponding number of planes
      !-------------------------------------------------------------------------------------------------
      
      
      CALL AllocAry(P%X_UNFROZEN, NL,'Position of unfrozen planes', TmpErrStat, TmpErrMsg )
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN
      
      
      DO IL=1,NL
          READ (UnitWind,IOStat=TmpErrStat)  Dum_Int2
          Xvec(IL) = Dum_Int2
          P%X_UNFROZEN(IL) = Dum_Int2
      end DO
      
      P%NLONGI = NL
      
      
      !-------------------------------------------------------------------------------------------------
      ! Allocate space for the FF array
      !-------------------------------------------------------------------------------------------------
    IF(p%NFFSteps/2.0==p%NFFSteps/2) Then
      TmpNumSteps = p%NFFSteps 
    ELSE
      TmpNumSteps = p%NFFSteps + 1       ! add another step, just in case there is an odd number of steps.
    ENDIF
    

		

    
   !bjj: should we reorganize this FFData array so we access the data faster?

      IF ( .NOT. ALLOCATED( p%FFData ) ) THEN
         CALL AllocAry( p%FFData, NL,p%NZGrids,p%NYGrids,p%NFFComp,TmpNumSteps, &
                  'Full-field wind data array.', TmpErrStat, TmpErrMsg )
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN

      ELSE
         IF (SIZE(p%FFData,1) /= p%NZGrids .OR. SIZE(p%FFData,2) /= p%NYGrids .OR. &
             SIZE(p%FFData,3) /= p%NFFComp .OR. SIZE(p%FFData,3) /= TmpNumSteps ) THEN

               ! Let's make the array the correct size (we should never get here, but you never know)

            DEALLOCATE( p%FFData )

            CALL AllocAry( p%FFData, NL,p%NZGrids,p%NYGrids,p%NFFComp,TmpNumSteps, &
                  'Full-field wind data array.', TmpErrStat, TmpErrMsg )
               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
               IF ( ErrStat >= AbortErrLev ) RETURN            

         ENDIF !Incorrect size
      ENDIF ! allocated

      !-------------------------------------------------------------------------------------------------
      ! Initialize the data and set column indexing to account for direction of turbine rotation (CWise)
      !-------------------------------------------------------------------------------------------------

      p%FFData(:,:,:,:,:) = 0.0                        ! we may have only one component

      IF ( CWise )  THEN
         CFirst    = p%NYGrids
         CLast     = 1
         CStep     = -1
      ELSE
         CFirst    = 1
         CLast     = p%NYGrids
         CStep     = 1
      ENDIF

    

          
      
      
      
      !-------------------------------------------------------------------------------------------------
      ! Loop through all the time steps, reading the data and converting to m/s
      !-------------------------------------------------------------------------------------------------
   !bjj: should we reorganize this FFData array so we access the data faster?

      p%NFFSteps = TmpNumSteps
      
     
  DO IL=1,NL
   TIME_LOOP:  DO IT=1,TmpNumSteps     ! time (add 1 to see if there is an odd number of grids)

         DO IR=1,p%NZGrids               ! the rows (vertical)

            DO IC=CFirst,CLast,CStep   ! the columns (lateral)

               DO I=1,p%NFFComp          ! wind components (U, V, W)

                     ! Get the next integer from the file.
                     ! This is a 2-byte integer, so we can't use the library read routines.
                  READ (UnitWind,IOStat=TmpErrStat)  Dum_Int2
                  IF (TmpErrStat /= 0) THEN
                     IF ( IT == TmpNumSteps ) THEN ! There really were an even number of steps
                        p%NFFSteps = TmpNumSteps - 1
                        ErrStat  = 0
                        EXIT TIME_LOOP
                     ELSE
                        CALL SetErrStat( ErrID_Fatal, ' Error reading binary data file. '// &
                                    'ic = '//TRIM(Num2LStr(ic))// &
                                    ', ir = '//TRIM(Num2LStr(ir))// &
                                    ', it = '//TRIM(Num2LStr(it))// &
                                    ', nffsteps = '//TRIM(Num2LStr(p%NFFSteps)), ErrStat, ErrMsg, RoutineName)
                        RETURN
                     ENDIF
                  ELSE
                     p%FFData(IL,IR,IC,I,IT) = p%MeanFFWS*(FF_Offset(I)+0.00001*TI(I)*Dum_Int2)
                  ENDIF

               END DO !I

            END DO !IC

         END DO !IR

   END DO TIME_LOOP !IT
   END DO ! IL

      IF ( p%Periodic ) THEN
         TmpErrMsg = '   Processed '//TRIM( Num2LStr( p%NFFSteps ) )//' time steps of '// &
                    TRIM( Num2LStr ( p%FFRate ) )//'-Hz full-field data (period of '// &
                    TRIM( Num2LStr( p%FFDTime*p%NFFSteps ) )//' seconds).'
         
      ELSE
         TmpErrMsg= '   Processed '//TRIM( Num2LStr( p%NFFSteps ) )//' time steps of '// &
                    TRIM( Num2LStr ( p%FFRate ) )//'-Hz full-field data ('// &
                    TRIM( Num2LStr( p%FFDTime*( p%NFFSteps - 1 ) ) )//' seconds).'
      ENDIF
      CALL WrScr( NewLine//TRIM(TmpErrMsg) )
      !CALL SetErrStat( ErrID_Info, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
      
      

   END SUBROUTINE Read_Bladed_Grids_evo
   
   
   !====================================================================================================
   !> This subroutine reads the text summary file to get normalizing parameters, the location of the
   !! grid, and the direction the grid was written to the binary file
   !!   modified from Read_Summary_FF by 16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
    !   modified 14-Nov- 2020 F. Guo Flensburg University of Applied Sciences  some variables are not sotred in type
   SUBROUTINE Read_Summary_FF_evo ( UnitWind, FileName, CWise, ZCenter, TI, ErrStat, ErrMsg )

      IMPLICIT                                              NONE

      CHARACTER(*),        PARAMETER                     :: RoutineName="Read_Summary_FF_evo"


         ! Passed variables
      INTEGER(IntKi),                     INTENT(IN   )  :: UnitWind       !< unit number for the file to open
      CHARACTER(*),                       INTENT(IN   )  :: FileName       !< name of the summary file
      LOGICAL,                            INTENT(  OUT)  :: CWise          !< rotation (for reading the order of the binary data)
      REAL(ReKi),                         INTENT(  OUT)  :: ZCenter        !< the height at the center of the grid
      REAL(ReKi),                         INTENT(  OUT)  :: TI      (3)    !< turbulence intensities of the wind components as defined in the FF file, not necessarially the actual TI
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat        !< returns 0 if no error encountered in the subroutine
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg         !< holds the error messages

        ! Local variables
      REAL(ReKi)                                         :: ZGOffset       ! The vertical offset of the turbine on rectangular grid (allows turbulence not centered on turbine hub)

      INTEGER, PARAMETER                                 :: NumStrings = 6 ! number of strings to be looking for in the file

      INTEGER(IntKi)                                     :: FirstIndx      ! The first character of a line where data is located
      INTEGER(IntKi)                                     :: I              ! A loop counter
      INTEGER(IntKi)                                     :: LastIndx       ! The last  character of a line where data is located
      INTEGER(IntKi)                                     :: LineCount      ! Number of lines that have been read in the file

      LOGICAL                                            :: StrNeeded(NumStrings)   ! if the string has been found

      CHARACTER(1024)                                    :: LINE           ! temporary storage for reading a line from the file

         ! Temporary variables for error handling
      INTEGER(IntKi)                                     :: TmpErrStat     ! temporary error status
      CHARACTER(ErrMsgLen)                               :: TmpErrMsg      ! temporary error message

         !----------------------------------------------------------------------------------------------
         ! Initialize some variables
         !----------------------------------------------------------------------------------------------

      ErrStat              = ErrID_None
      ErrMsg               = ''

      LineCount            = 0
      StrNeeded(:)         = .TRUE.
      ZGOffset             = 0.0
      p%RefHt      = 0.0
      p%Periodic   = .FALSE.

         !----------------------------------------------------------------------------------------------
         ! Open summary file.
         !----------------------------------------------------------------------------------------------

      CALL OpenFInpFile ( UnitWind, TRIM( FileName ), TmpErrStat, TmpErrMsg )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )                  
      IF ( ErrStat >= AbortErrLev ) RETURN


         !----------------------------------------------------------------------------------------------
         ! Read the summary file.
         !----------------------------------------------------------------------------------------------

      ! Here are the strings we're looking for, in this order:
      ! 1) 'CLOCKWISE'
      ! 2) 'HUB HEIGHT'
      ! 3)     (unused; decided we didn't need to read data also stored in the binary file)
      ! 4) 'UBAR'
      ! 5) 'HEIGHT OFFSET' (optional)
      ! 6) 'PERIODIC' (optional)
         
         
      DO WHILE ( ( ErrStat == ErrID_None ) .AND. StrNeeded(NumStrings) )

         LineCount = LineCount + 1

         READ ( UnitWind, '(A)', IOSTAT=TmpErrStat ) LINE
         IF ( TmpErrStat /= 0 ) THEN

               ! the "HEIGHT OFFSET" and "PERIODIC" parameters are not necessary.  We'll assume they are zero/false if we didn't find it.
            IF ( StrNeeded(1) .OR. StrNeeded(2) .OR. StrNeeded(4)  ) THEN
               CALL SetErrStat( ErrID_Fatal, ' Error reading line #'//TRIM(Num2LStr(LineCount))//' of the summary file, "'// &
                           TRIM(FileName)//'". Could not find all of the required parameters.', ErrStat, ErrMsg, RoutineName )                  
               RETURN
            ELSE
               EXIT
            ENDIF

         ENDIF

         CALL Conv2UC ( LINE )


         IF ( StrNeeded(1) ) THEN

            !-------------------------------------------------------------------------------------------
            ! #1: Get the rotation direction, using the string "CLOCKWISE"
            !-------------------------------------------------------------------------------------------

            IF ( INDEX( LINE, 'CLOCKWISE' ) > 0 ) THEN

               READ (LINE, *, IOSTAT = TmpErrStat)  CWise          ! Look for True/False values

               IF ( TmpErrStat /= 0 ) THEN                         ! Look for Yes/No values instead

                  LINE = ADJUSTL ( LINE )                      ! Remove leading spaces from input line

                  SELECT CASE (LINE(1:1) )
                     CASE ('Y')
                        CWise = .TRUE.
                     CASE ('N')
                        CWise = .FALSE.
                     CASE DEFAULT
                        CALL SetErrStat( ErrID_Fatal, ' Error reading rotation direction (CLOCKWISE) from FF summary file.', ErrStat, ErrMsg, RoutineName )                  
                        RETURN
                  END SELECT

               ENDIF ! TmpErrStat /= 0
               StrNeeded(1) = .FALSE.

            ENDIF   ! INDEX for "CLOCKWISE"

         ELSEIF ( StrNeeded(2) ) THEN

            !-------------------------------------------------------------------------------------------
            ! #2: Get the hub height, using the strings "HUB HEIGHT" or "ZHUB"
            !-------------------------------------------------------------------------------------------

            IF ( INDEX( LINE, 'HUB HEIGHT' ) > 0 .OR. INDEX( LINE, 'ZHUB' ) > 0 ) THEN

               READ (LINE, *, IOSTAT = TmpErrStat) p%RefHt

               IF ( TmpErrStat /= 0 ) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading hub height from FF summary file.', ErrStat, ErrMsg, RoutineName )                  
                  RETURN
               ENDIF
               StrNeeded(2) = .FALSE.

            ENDIF !INDEX for "HUB HEIGHT" or "ZHUB"

   !      ELSEIF ( StrNeeded(3) ) THEN
   !
   !         !-------------------------------------------------------------------------------------------
   !         ! #3: Get the grid width (& height, if available), using the strings "GRID WIDTH" or "RDIAM"
   !         !    If GRID HEIGHT is specified, use it, too. -- THIS IS UNNECESSARY AS IT'S STORED IN THE BINARY FILE
   !         !-------------------------------------------------------------------------------------------

         ELSEIF ( StrNeeded(4) ) THEN

            !-------------------------------------------------------------------------------------------
            ! #4: Get the mean wind speed "UBAR" and turbulence intensities from following lines for
            !     scaling Bladed-style FF binary files
            !-------------------------------------------------------------------------------------------

            IF ( INDEX( LINE, 'UBAR') > 0 ) THEN

               FirstIndx = INDEX( LINE, '=' ) + 1        ! Look for the equal siqn to find the number we're looking for

               READ ( LINE( FirstIndx:LEN(LINE) ), *, IOSTAT=TmpErrStat ) p%MeanFFWS

               IF ( TmpErrStat /= 0 ) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading UBar binary data normalizing parameter from FF summary file.', ErrStat, ErrMsg, RoutineName )                  
                  RETURN
               ENDIF

               DO I = 1,3

                  LineCount = LineCount + 1

                  READ ( UnitWind, '(A)', IOSTAT=TmpErrStat ) LINE
                  IF ( TmpErrStat /= 0 ) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading line #'//TRIM(Num2LStr(LineCount))//' of the summary file, "'//TRIM(FileName)//&
                                          '". Could not find all of the required parameters.', ErrStat, ErrMsg, RoutineName )                  
                     RETURN
                  ENDIF

                  FirstIndx = INDEX( LINE, '=' ) + 1     ! Read the number between the = and % signs
                  LastIndx  = INDEX( LINE, '%' ) - 1

                  IF ( LastIndx <= FirstIndx ) LastIndx = LEN( LINE )   ! If there's no % sign, read to the end of the line

                  READ ( LINE( FirstIndx:LastIndx ), *, IOSTAT=TmpErrStat ) TI(I)
                  IF ( TmpErrStat /= 0 ) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error reading TI('//TRIM(Num2LStr(I))// &
                                 ') binary data normalizing parameter from FF summary file.', ErrStat, ErrMsg, RoutineName )                  
                     RETURN
                  ENDIF

               ENDDO !I

               StrNeeded(4) = .FALSE.

             ENDIF

         ELSEIF ( StrNeeded(5) ) THEN

            !-------------------------------------------------------------------------------------------
            ! #5: Get the grid "HEIGHT OFFSET", if it exists (in TurbSim). Otherwise, assume it's zero
            !           ZGOffset = HH - GridBase - p%FFZHWid
            !-------------------------------------------------------------------------------------------
            IF ( INDEX( LINE, 'HEIGHT OFFSET' ) > 0  ) THEN

               FirstIndx = INDEX ( LINE, '=' ) + 1

               READ ( LINE( FirstIndx:LEN(LINE) ), *, IOSTAT=TmpErrStat ) ZGOffset

               IF ( TmpErrStat /= 0 ) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error reading height offset from FF summary file.', ErrStat, ErrMsg, RoutineName )                  
                  RETURN
               ENDIF

               StrNeeded(5) = .FALSE.

            ENDIF !INDEX for "HEIGHT OFFSET"

         ELSEIF ( StrNeeded(6) ) THEN

            !-------------------------------------------------------------------------------------------
            ! #5: Get the grid "PERIODIC", if it exists (in TurbSim). Otherwise, assume it's
            !        not a periodic file
            !-------------------------------------------------------------------------------------------
            IF ( INDEX( LINE, 'PERIODIC' ) > 0  ) THEN

               p%Periodic   = .TRUE.
               StrNeeded(6)         = .FALSE.

            ENDIF !INDEX for "PERIODIC"

         ENDIF ! StrNeeded


      ENDDO !WHILE


      !-------------------------------------------------------------------------------------------------
      ! Close the summary file
      !-------------------------------------------------------------------------------------------------

      CLOSE ( UnitWind )


      !!-------------------------------------------------------------------------------------------------
      !! calculate the height of the grid center
      !!-------------------------------------------------------------------------------------------------
      !
       zcenter  = p%refht - zgoffset


   END SUBROUTINE Read_Summary_FF_evo
    
    END SUBROUTINE LidarSim_InitializeWindEvolution 
    
    
!######################################################
    
    
    
     !this routine read the binary wind data of upstream positions, the code is mostly similarly to ifw_BladedFFWind module.
     !fot the initial implementing of wind evolution only bladed style file is supported when we generate the 4D turbulence
     !we also output the data as the bladed style
SUBROUTINE LidarSim_InitializeAvailability(p, InputFileData, ErrStat, ErrMsg)
    
    IMPLICIT                        NONE
    CHARACTER(*),                   PARAMETER       ::  RoutineName="LidarSim_InitializeAvailability"
    
    TYPE(LidarSim_ParameterType),   INTENT(INOUT)   ::  p           ! parameter data to write results in
    TYPE(LidarSim_InputFile),       INTENT(IN   )   ::  InputFileData       ! Inputdata from the input file
    INTEGER(IntKi),                 INTENT(  OUT)   ::  ErrStat             !< Error status of the operation
    CHARACTER(*),                   INTENT(  OUT)   ::  ErrMsg              !< Error message if ErrStat /= ErrID_None
    
    ! local variables
    INTEGER(IntKi)                                  ::  TmpErrStat          !< Temporary error status
    CHARACTER(ErrMsgLen)                            ::  TmpErrMsg           !< temporary error message
    INTEGER(IntKi)                                  ::  UnitWind            ! Unit number for the InflowWind input file
    INTEGER(B2Ki)                                   ::  Dum_Int2            ! temporary variavle to store read in data from binary file
    !CHARACTER( 1028 )                               ::  SumFile             ! length is LEN(ParamData%WindFileName) + the 4-character extension.
    !CHARACTER( 1028 )                               ::  TwrFile             ! length is LEN(ParamData%WindFileName) + the 4-character extension.
    CHARACTER( 1028 )                               ::  binFile             ! length is LEN(ParamData%WindFileName) + the 4-character extension.

    !REAL(ReKi)                                      ::  TI      (3)         ! turbulence intensities of the wind components as defined in the FF file, not necessarially the actual TI
   ! REAL(ReKi)                                      ::  Xvec    (30)
    !REAL(ReKi)                                      ::  BinTI   (3)          ! turbulence intensities of the wind components as defined in the FF binary file, not necessarially the actual TI
    INTEGER(IntKi)                                     ::  Ntime                ! number of time steps of LOS measurement
    INTEGER(IntKi)                                     ::  Nposition            ! number of lidar measurement positions
    REAL(ReKi)                                         ::  DataAvaiTimeEnd      ! end time of the lidar availability series
    REAL(ReKi)                                         ::  SampleTime           ! the unscalled sampling time of lidar
    INTEGER(B2Ki)                                      ::  SampleTimeScaled       ! the scalled sampling time of lidar
    INTEGER(IntKi)                                     ::  i_position             ! index number 
    INTEGER(IntKi)                                     ::  i_time                 ! index number
    !LOGICAL                                         ::  CWise
    !LOGICAL                                         ::  Exists
    
    
    
    ErrStat =   0
    ErrMsg  =   ''
    

    TmpErrMsg   = ''
    TmpErrStat  = 0

    P%AvailabilityFlag	= INPUTFILEDATA%AvailabilityFlag	! store the flag in parameter
    
    
    IF (INPUTFILEDATA%TRAJECTORYTYPE==0) THEN
        Nposition = 		INPUTFILEDATA%NUMBEROFPOINTS_CARTESIAN 
    ELSE IF (INPUTFILEDATA%TRAJECTORYTYPE==1) THEN
        Nposition = 		INPUTFILEDATA%NUMBEROFPOINTS_SPHERICAL*INPUTFILEDATA%GATESPERBEAM

    ELSE
    TmpErrStat = 1
    TmpErrMsg  = 'Measurement trajectory type is incorrect, 0 = cartesian coordinates; 1 = spherical coordinates'
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'Lidar measurement availability initialization')
    IF ( ErrStat >= AbortErrLev ) RETURN
    
    END IF

    
    
    
       !Get a unit number to use

    CALL GetNewUnit(UnitWind, TmpErrStat, TmpErrMsg)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'Lidar measurement availability initialization')
    IF ( ErrStat >= AbortErrLev ) RETURN
    
    
    
    !...........................................................................................
    ! Create full-field binary file name from evo file root name.  Also get tower file
    ! name.
    !...........................................................................................
    
    !CALL GetRoot(InputFileData%AvailabilityFileNameRoot, binFile)
        

      !----------------------------------------------------------------------------------------------
      ! Open the binary file
      !----------------------------------------------------------------------------------------------
    
    CALL OpenBInpFile (UnitWind, InputFileData%AvailabilityFileNameRoot, TmpErrStat, TmpErrMsg)
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'Lidar measurement availability initialization')
    IF ( ErrStat >= AbortErrLev ) RETURN
    
      ! Read the first binary integer from the file to get the info from the binary read in data (to cross check whether we are reading the correct binary data).
    
    READ ( UnitWind, IOSTAT=TmpErrStat )  SampleTimeScaled  ! the sample time scalled to integral
    READ ( UnitWind, IOSTAT=TmpErrStat )  Dum_Int2  ! end time of the lidar availability series
    DataAvaiTimeEnd         = Dum_Int2
    
    READ ( UnitWind, IOSTAT=TmpErrStat )  Dum_Int2    ! number of measurement positions
    
    IF (Dum_Int2 /= Nposition) THEN  
    TmpErrStat = 1
    TmpErrMsg  = 'The number of measurement positions in the availability binary data does not match with the lidar input, check the correct binary data is read!'
    END IF
    
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'Lidar measurement availability initialization')
    IF ( ErrStat >= AbortErrLev ) RETURN
    
    
    SampleTime = SampleTimeScaled*0.00001;
    
    IF (abs(SampleTime-INPUTFILEDATA%T_MEASUREMENT_INTERVAL)>1e-3) THEN  
    TmpErrStat = 1
    TmpErrMsg  = 'The time step in the availability binary data does not match with the lidar input (measurement interval), check the correct binary data is read!'
    END IF
    
    CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'Lidar measurement availability initialization')
    IF ( ErrStat >= AbortErrLev ) RETURN
    
    READ ( UnitWind, IOSTAT=TmpErrStat )  Dum_Int2      ! end time of the lidar availability series
    Ntime  = Dum_Int2
     ! allocate the availabiltiy data array in 

      IF ( .NOT. ALLOCATED(P%AVAIDATA)) THEN
      
       CALL AllocAry( P%AVAIDATA, NTIME,NPOSITION,'Lidar measurement availability data.', TmpErrStat, TmpErrMsg )
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN
      
      ELSE
      
      END IF ! allocated
      
      !-------------------------------------------------------------------------------------------------
      ! Initialize the data and set column indexing to account for direction of turbine rotation (CWise)
      !-------------------------------------------------------------------------------------------------
      
      P%AVAIDATA(:,:) = 1                        ! we may have only one component
      P%AvaiTimeEnd   = DataAvaiTimeEnd
      ! loop over to read in availability data
      DO i_position=1,Nposition
          TIME_LOOP: DO i_time=1,Ntime     ! time 
            ! Get the next integer from the file.
            ! This is a 2-byte integer, so we can't use the library read routines.
            READ (UnitWind,IOStat=TmpErrStat)  Dum_Int2
            IF (TmpErrStat /= 0) THEN
                IF ( i_time == Ntime) THEN ! we really reach the data end
                ErrStat  = 0
                EXIT TIME_LOOP
                ELSE
                CALL SetErrStat( ErrID_Fatal, ' Error reading binary data file. '// &
                            'i_position = '//TRIM(Num2LStr(i_position))// &
                            ', i_time = '//TRIM(Num2LStr(i_time)), ErrStat, ErrMsg, RoutineName)
                RETURN
                ENDIF
            ELSE
                !IF ( i_time == Ntime) THEN ! we really reach the data end, and there are still data read in so something is wrong
                !CALL SetErrStat( ErrID_Fatal, ' Error reading binary data file. The data length does not match with the specified value.'// &
                !            'i_position = '//TRIM(Num2LStr(i_position))// &
                !            ', i_time = '//TRIM(Num2LStr(i_time)), ErrStat, ErrMsg, RoutineName)
                !RETURN
                !ENDIF
                P%AVAIDATA(i_time,i_position) = Dum_Int2
            ENDIF

            END DO TIME_LOOP !i_time TIME_LOOP
       END DO ! i_position
      
      

    CLOSE( UnitWind )
    
    !IF (TmpErrStat /= 0) THEN
    !   CALL SetErrStat( ErrID_Fatal, ' Error reading first binary integer from file "'//TRIM(InputFileData%AvailabilityFileNameRoot)//'."', &
    !        ErrStat, ErrMsg, 'Lidar measurement availability initialization')
    !   RETURN
    !ENDIF
    
END  SUBROUTINE LidarSim_InitializeAvailability  
    
    
            
!######################################################
    
    SUBROUTINE LidarSim_CalculateVlos(p, UnitVector_I, Vlos, MeasuringPosition_I, LidarPosition_I,&
    Time,ErrStat, ErrMsg)
    
    IMPLICIT                                NONE
    CHARACTER(*),                           PARAMETER       ::  RoutineName="LidarSim_CalculateVlos"
    
    TYPE(LidarSim_ParameterType),           INTENT(INOUT)   ::  p                           !parameter data of the lidar module
    REAL(ReKi),                             INTENT(IN   )   ::  UnitVector_I(3)             !Line of Sight Unit Vector
    REAL(ReKi),                             INTENT(INOUT)   ::  MeasuringPosition_I(3)      !Position of the measuring point
    REAL(ReKi),                             INTENT(IN   )   ::  LidarPosition_I(3)      !Position of the measuring point
    REAL(ReKi),                             INTENT(  OUT)   ::  Vlos                        !Calculated speed in los direction
    REAL(DbKi),                             INTENT(IN   )   ::  Time                        !< Current simulation time in seconds
    
    
    !IfW Parameter
    !TYPE(InflowWind_ParameterType),         INTENT(IN   )   ::  IfW_p                       !< Parameters
    !TYPE(InflowWind_ContinuousStateType),   INTENT(IN   )   ::  IfW_ContStates              !< Continuous states at Time
    !TYPE(InflowWind_DiscreteStateType),     INTENT(IN   )   ::  IfW_DiscStates              !< Discrete states at Time
    !TYPE(InflowWind_ConstraintStateType),   INTENT(IN   )   ::  IfW_ConstrStates            !< Constraint states at Time
    !TYPE(InflowWind_OtherStateType),        INTENT(IN   )   ::  IfW_OtherStates             !< Other/optimization states at Time
    !TYPE(InflowWind_MiscVarType),           INTENT(INOUT)   ::  IfW_m                       !< Misc variables for optimization (not copied in glue code)
    
    !Error Variables
    INTEGER(IntKi),                         INTENT(  OUT)   ::  ErrStat                     !< Error status of the operation
    CHARACTER(*),                           INTENT(  OUT)   ::  ErrMsg                      !< Error message if ErrStat /= ErrID_None
    
    !Local Variables
    !TYPE(InflowWind_InputType)              ::  InputForCalculation                         ! input data field for the calculation in the InflowWind module
    !TYPE(InflowWind_OutputType)             ::  OutputForCalculation                        ! data field were the calculated speed is written from the InflowWind module
    INTEGER(IntKi)                          ::  Counter                                     ! Counter for the loop for the different weightings of the point
    REAL(ReKi),DIMENSION(:), ALLOCATABLE    ::  Vlos_tmp                                    !< Array with all temporary Vlos
    REAL(ReKi)                              ::  PositionXYZ(3,1)                                    !< Array with all temporary Vlos
    REAL(ReKi)                              ::  VelocityUVW(3,1)                                    !< Array with all temporary Vlos
    
    
    ! Temporary variables for error handling
    INTEGER(IntKi)                                          ::  ErrStat2        
    CHARACTER(ErrMsgLen)                                    ::  ErrMsg2
    
    !Initialize error values
    ErrStat        =  0
    ErrMsg         =  ""
    
    !CALL AllocAry(InputForCalculation%PositionXYZ, 3,1, 'InputForCalculation%PositionXYZ',ErrStat2, ErrMsg2)        !Allocating needed space for the input
    !CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
    !CALL AllocAry(OutputForCalculation%VelocityUVW, 3,1, 'OutputForCalculation%VelocityUVW',ErrStat2, ErrMsg2)      !Allocating needed space for the output
    !CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
    CALL AllocAry(Vlos_tmp ,SIZE(p%Weighting), 'Vlos_tmp%VelocityUVW',ErrStat2, ErrMsg2)                            !Allocating space for temporary windspeeds
    CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
     
    IF( P%WINDFILEFORMAT == 1 .OR. P%WINDFILEFORMAT == 2)Then !Uniform Wind 2 (und steady 1)
        MeasuringPosition_I(1) = MeasuringPosition_I(1)-LidarPosition_I(1)  !In the uniform wind case. the wind hits the turbine at the same time indepentend of the x shift
        DO Counter = 1, SIZE(p%Weighting)
            PositionXYZ(:,1) = MeasuringPosition_I + p%WeightingDistance(Counter) * UnitVector_I                                                    !position of the weighted measuring point
            !CALL  LidarSim_CalculateUVW(Time+DBLE(-PositionXYZ(1,1)/p%Uref),PositionXYZ, p, VelocityUVW, ErrStat2, ErrMsg2 )     !Calculation of the windspeed   still need to be updated
            !Vlos_tmp(Counter) = - DOT_PRODUCT(VelocityUVW(:,1),UnitVector_I)
            Vlos_tmp(Counter) = 10
        END DO
    ELSE IF(P%WINDFILEFORMAT ==  3 .OR. P%WINDFILEFORMAT == 4) THEN        !Bladed Turublent 4 ( und TurbSim 3)
        DO Counter = 1, SIZE(p%Weighting)
            
            PositionXYZ(:,1) = MeasuringPosition_I + p%WeightingDistance(Counter) * UnitVector_I   
            !position of the weighted measuring point
           
            CALL LidarSim_CalculateUVW(Time,PositionXYZ, p, VelocityUVW, ErrStat2, ErrMsg2 ) 
            Vlos_tmp(Counter) = - DOT_PRODUCT(VelocityUVW(:,1),UnitVector_I)
        END DO
        
    ELSE
        CALL SetErrStat( ErrID_Fatal, ' Error: The wind type '//TRIM(Num2LStr(P%WINDFILEFORMAT))//'. '// &
                      'is not currently support for lidar measurement simulation.', ErrStat, ErrMsg, RoutineName )
        
        
    END IF
    Vlos = DOT_PRODUCT(Vlos_tmp, p%Weighting)           !Calculation of the weighted windspeed
    
    !DEALLOCATE (InputForCalculation%PositionXYZ)        !Free Input Positions for the next measurement
    !DEALLOCATE (OutputForCalculation%VelocityUVW)       !Free Ouput Positions for the next measurement
    DEALLOCATE (Vlos_tmp)
    
    END SUBROUTINE LidarSim_CalculateVlos
    

    
    ! #######################################################
    ! Calculate the velocity based on the specific position
SUBROUTINE LidarSim_CalculateUVW( Time, PositionXYZ, p, VelocityUVW, ErrStat, ErrMsg )


      IMPLICIT                                                    NONE

      CHARACTER(*),              PARAMETER                     :: RoutineName="LidarSim_CalculateUVW"


         ! Inputs / Outputs

      REAL(DbKi),                               INTENT(IN   )  :: Time           !< Current simulation time in seconds
      REAL(ReKi),                               INTENT(IN   )  :: PositionXYZ(3)      !< Inputs at Time
      TYPE(LidarSim_ParameterType),           INTENT(INOUT)   ::  p                           !parameter data of the lidar module
      REAL(ReKi),                             INTENT(INOUT)   ::  VelocityUVW(3)      !Position of the measuring point
    
      
      INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat        !< Error status of the operation
      CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None


         ! Local variables
      !REAL(ReKi), ALLOCATABLE                                  :: PositionXYZprime(:,:)   !< PositionXYZ array in the prime (wind) coordinates
      REAL(ReKi)                                               :: DiskVel(3)     !< HACK for AD14: disk velocity output at Time

      INTEGER(IntKi)                                           :: I                   !< Generic counters


         ! Temporary variables for error handling
      INTEGER(IntKi)                                           :: TmpErrStat
      CHARACTER(ErrMsgLen)                                     :: TmpErrMsg            ! temporary error message



         ! Initialize ErrStat
      ErrStat  = ErrID_None
      ErrMsg   = ""


      !-----------------------------------------------------------------------
      !  Points coordinate transforms from to global to wind file coordinates
      !-----------------------------------------------------------------------
         !FG£º assume the wind always propogate in the positive x direction
      
      
      !!---------------------------------
      !!  
      !
      !   ! Compute the wind velocities by stepping through all the data points and calling the appropriate GetWindSpeed routine
      SELECT CASE (P%WINDFILEFORMAT)
         
      !  CASE (Steady_WindNumber, Uniform_WindNumber)
      !
      !         ! InputData only contains the Position array, so we can pass that directly.
      !      CALL  IfW_UniformWind_CalcOutput(  Time, PositionXYZprime, p%UniformWind, y%VelocityUVW, &
      !                                    DiskVel, m%UniformWind, TmpErrStat, TmpErrMsg)
      !
      !      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
      !      IF ( ErrStat >= AbortErrLev ) RETURN
      !
      !         ! Call IfW_UniformWind_CalcOutput again in order to get the values needed for the OutList -- note that we do not report errors from this
      !      IF ( p%NWindVel >= 1_IntKi .AND. FillWrOut ) THEN
      !            ! Move the arrays for the Velocity information
      !         CALL  IfW_UniformWind_CalcOutput(  Time, p%WindViXYZprime, p%UniformWind, m%WindViUVW, &
      !                                       DiskVel, m%UniformWind, TmpErrStat, TmpErrMsg)
      !      ENDIF
      !
         CASE (3,4) !Turbsim and Bladed style
      
               ! InputData only contains the Position array, so we can pass that directly.
            CALL  LidSim_TS_Bladed_FFWind_CalcOutput(Time, PositionXYZ, p, VelocityUVW, TmpErrStat, TmpErrMsg)
      
            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
            IF ( ErrStat >= AbortErrLev ) RETURN
      
      
      
      !   CASE (4) !also includes BladedFF_Shr_WindNumber
      !
      !         ! InputData only contains the Position array, so we can pass that directly.
      !      CALL  IfW_BladedFFWind_CalcOutput(  Time, PositionXYZprime, p%BladedFFWind, &
      !                                    y%VelocityUVW, DiskVel, m%BladedFFWind, TmpErrStat, TmpErrMsg)
      !
      !      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
      !      IF ( ErrStat >= AbortErrLev ) RETURN
      !
      !
      !         ! Call IfW_BladedFFWind_CalcOutput again in order to get the values needed for the OutList
      !      IF ( p%NWindVel >= 1_IntKi  .AND. FillWrOut ) THEN
      !            ! Move the arrays for the Velocity information
      !         CALL  IfW_BladedFFWind_CalcOutput(  Time, p%WindViXYZprime, p%BladedFFWind, &
      !                                       m%WindViUVW, DiskVel, m%BladedFFWind, TmpErrStat, TmpErrMsg)
      !
      !            ! Out of bounds errors will be ErrID_Severe, not ErrID_Fatal
      !         IF ( TmpErrStat >= ErrID_Fatal ) THEN
      !            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
      !            RETURN
      !         ENDIF
      !
      !      ENDIF
      !
      !
      !   CASE (User_WindNumber)
      !
      !         ! InputData only contains the Position array, so we can pass that directly.
      !      CALL  IfW_UserWind_CalcOutput(  Time, PositionXYZprime, p%UserWind, &
      !                                    y%VelocityUVW, DiskVel, m%UserWind, TmpErrStat, TmpErrMsg)
      !
      !      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
      !      IF ( ErrStat >= AbortErrLev ) RETURN
      !
      !
      !         ! Call IfW_UserWind_CalcOutput again in order to get the values needed for the OutList
      !      IF ( p%NWindVel >= 1_IntKi  .AND. FillWrOut ) THEN
      !            ! Move the arrays for the Velocity information
      !         CALL  IfW_UserWind_CalcOutput(  Time, p%WindViXYZprime, p%UserWind, &
      !                                       m%WindViUVW, DiskVel, m%UserWind, TmpErrStat, TmpErrMsg)
      !
      !            ! Out of bounds errors will be ErrID_Severe, not ErrID_Fatal
      !         IF ( TmpErrStat >= ErrID_Fatal ) THEN
      !            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
      !            RETURN
      !         ENDIF
      !
      !      ENDIF
      !
      !   CASE ( HAWC_WindNumber )
      !      
      !         ! InputData only contains the Position array, so we can pass that directly.
      !      CALL  IfW_HAWCWind_CalcOutput(  Time, PositionXYZprime, p%HAWCWind, &
      !                                    y%VelocityUVW, DiskVel, m%HAWCWind, TmpErrStat, TmpErrMsg)
      !
      !      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
      !      IF ( ErrStat >= AbortErrLev ) RETURN
      !
      !
      !         ! Call IfW_TSFFWind_CalcOutput again in order to get the values needed for the OutList
      !      IF ( p%NWindVel >= 1_IntKi  .AND. FillWrOut ) THEN
      !         CALL  IfW_HAWCWind_CalcOutput(  Time, p%WindViXYZprime, p%HAWCWind, &
      !                                       m%WindViUVW, DiskVel, m%HAWCWind, TmpErrStat, TmpErrMsg)
      !
      !            ! Out of bounds errors will be ErrID_Severe, not ErrID_Fatal
      !         IF ( TmpErrStat >= ErrID_Fatal ) THEN
      !            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
      !            RETURN
      !         ENDIF
      !      ENDIF
      !      
      !
      !   CASE ( FDext_WindNumber )
      !      
      !      CALL IfW_4Dext_CalcOutput(Time, PositionXYZprime, p%FDext, y%VelocityUVW,  m%FDext, TmpErrStat, TmpErrMsg) 
      !      DiskVel = 0.0_ReKi ! this is only for AD14, which we frankly don't care about in 4Dext wind
      !                  
      !      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
      !      IF ( ErrStat >= AbortErrLev ) RETURN
      !
      !
      !         ! Call IfW_4Dext_CalcOutput again in order to get the values needed for the OutList
      !      IF ( p%NWindVel >= 1_IntKi  .AND. FillWrOut ) THEN
      !            ! Move the arrays for the Velocity information
      !         CALL IfW_4Dext_CalcOutput(Time, p%WindViXYZprime, p%FDext, m%WindViUVW, m%FDext, TmpErrStat, TmpErrMsg)
      !         
      !            ! Out of bounds errors will be ErrID_Severe, not ErrID_Fatal
      !         IF ( TmpErrStat >= ErrID_Fatal ) THEN
      !            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
      !            RETURN
      !         ENDIF
      !
      !      ENDIF
      !      
      !      
      !      ! If it isn't one of the above cases, we have a problem and won't be able to continue
      !
         CASE DEFAULT
      
            CALL SetErrStat( ErrID_Fatal, ' Error: Undefined wind type '//TRIM(Num2LStr(P%WINDFILEFORMAT))//'. '// &
                      'Currently only turbsim and bladed wind types (3 and 4) are supported for lidar simulation.', ErrStat, ErrMsg, RoutineName )
      
            VelocityUVW(:) = 0.0
            RETURN
      
      END SELECT


            ! Add coherent turbulence to background wind

!!!         IF (p%CTTS_Flag) THEN
!!!
!!!            DO PointCounter = 1, SIZE(InputData%Position, 2)
!!!
!!!               TempWindSpeed = CTTS_GetWindSpeed(     Time, InputData%Position(:,PointCounter), ErrStat, ErrMsg )
!!!
!!!                  ! Error Handling -- move ErrMsg inside CTTS_GetWindSPeed and simplify
!!!               IF (ErrStat >= ErrID_Severe) THEN
!!!                  ErrMsg   = 'IfW_CalcOutput: Error in CTTS_GetWindSpeed for point number '//TRIM(Num2LStr(PointCounter))
!!!                  EXIT        ! Exit the loop
!!!               ENDIF
!!!
!!!               y%Velocity(:,PointCounter) = y%Velocity(:,PointCounter) + TempWindSpeed
!!!
!!!            ENDDO
!!!
!!!               ! If something went badly wrong, Return
!!!            IF (ErrStat >= ErrID_Severe ) RETURN
!!!
!!!         ENDIF
!!!
      !ENDIF





      !-----------------------------------------------------------------------
      !  Windspeed coordinate transforms from Wind file coordinates to global
      !-----------------------------------------------------------------------

         ! The VelocityUVW array data that has been returned from the sub-modules is in the wind file (X'Y'Z') coordinates at
         ! this point.  These must be rotated to the global XYZ coordinates.  So now we apply the coordinate transformation
         ! to the VelocityUVW(prime) coordinates (in wind X'Y'Z' coordinate frame) returned from the submodules to the XYZ
         ! coordinate frame, but only if PropagationDir is not zero.  This is only a rotation of the returned wind field, so
         ! UVW contains the direction components of the wind at XYZ after translation from the U'V'W' wind velocity components
         ! in the X'Y'Z' (wind file) coordinate frame.
         ! NOTE: rotations are about the hub at [ 0 0 H ].  See InflowWind_SetParameters for details.
      !IF ( p%RotateWindBox ) THEN
      !   DO I  = 1,SIZE(y%VelocityUVW,DIM=2)
      !      y%VelocityUVW(:,I)   =  MATMUL( p%RotFromWind, y%VelocityUVW(:,I) )
      !   ENDDO
      !ENDIF

         ! We also need to rotate the reference frame for the WindViUVW array
         ! NOTE: rotations are about the hub at [ 0 0 H ].  See InflowWind_SetParameters for details.
      !IF ( p%RotateWindBox .AND. FillWrOut ) THEN
      !   DO I  = 1,SIZE(m%WindViUVW,DIM=2)
      !      m%WindViUVW(:,I)   =  MATMUL( p%RotFromWind, m%WindViUVW(:,I) )
      !   ENDDO
      !ENDIF


         ! DiskVel values over to the output and apply the coordinate transformation
   !   y%DiskVel   =  MATMUL( p%RotFromWind, DiskVel )
   !
   !   
   !   ! Done with the prime coordinates for the XYZ position information that was passed in.
   !IF (ALLOCATED(PositionXYZprime)) DEALLOCATE(PositionXYZprime)



END SUBROUTINE LidarSim_CalculateUVW
                      

!##############################################################
    
    SUBROUTINE LidarSim_CalculateIMU(p,y,u)
    
    IMPLICIT                                   NONE
    CHARACTER(*),                              PARAMETER            ::  RoutineName="LidarSim_LidarSim_CalculateIMU"
    
    TYPE(LidarSim_ParameterType),              INTENT(INOUT)        ::  p
    TYPE(LidarSim_OutputType),                 INTENT(INOUT)        ::  y                       !Outputs computed at Time (IN for mesh reasons and data allocation)
    TYPE(LidarSim_InputType),                  INTENT(IN   )        ::  u   
    
    ! local variables
    REAL(ReKi)                                                      ::  CrossProduct(3)         !Variable for the crossproduct of the rotation and the lidar position ( in the nacelle coord. system)
    REAL(ReKi)                                                      ::  Rotation_L_I(3,3)
    REAL(ReKi)                                                      ::  Roll
    REAL(ReKi)                                                      ::  Pitch
    REAL(ReKi)                                                      ::  Yaw
    REAL(ReKi)                                                      ::  DisplacementNacelle(3)
    REAL(ReKi)                                                      ::  DisplacementLidar(3)
    REAL(ReKi)                                                      ::  LidarPosition_I(3)
    
    !FG: currently, OpenFAST does not write the velocity and acceleration for the hub node, we have to use the velocity and acceleration from the nacelle to calculate the IMU. This should be improved in the future
    ! However, we can use the rotational angles from the hub mode to determine the measurement position and eular angles for the lidar.
    
    IF (p%SpinnerMountedFlag==.True.) THEN
        
        Rotation_L_I = MATMUL(TRANSPOSE(u%NacelleMotion%Orientation(:,:,1)),p%LidarOrientation_N)
        IF(.NOT.(Rotation_L_I(3,1) == 1 .OR. Rotation_L_I(3,1) == -1)) THEN
            Pitch = -ASIN(Rotation_L_I(3,1))
            Roll = ATAN2(Rotation_L_I(3,2)/cos(Pitch),Rotation_L_I(3,3)/cos(Pitch))
            Yaw = ATAN2(Rotation_L_I(2,1)/cos(Pitch), Rotation_L_I(1,1)/cos(Pitch))  
        ELSE
            Yaw = 0
            IF(Rotation_L_I(3,1) == 1) THEN
                Pitch = PiBy2_D
                Roll = Yaw + ATAN2(Rotation_L_I(1,2),Rotation_L_I(1,3))
            ELSE
                Pitch = -PiBy2_D
                Roll = -Yaw + ATAN2(-Rotation_L_I(1,2),-Rotation_L_I(1,3))
            END IF
        END IF
    
		LidarPosition_I = MATMUL(TRANSPOSE(u%NacelleMotion%Orientation(:,:,1)),p%LidarPosition_N) !Rotates the position vector to the orientation of the inertial coord. system
    
    ELSE
        Rotation_L_I = MATMUL(TRANSPOSE(u%NacelleMotion%Orientation(:,:,1)),p%LidarOrientation_N)
        IF(.NOT.(Rotation_L_I(3,1) == 1 .OR. Rotation_L_I(3,1) == -1)) THEN
            Pitch = -ASIN(Rotation_L_I(3,1))
            Roll = ATAN2(Rotation_L_I(3,2)/cos(Pitch),Rotation_L_I(3,3)/cos(Pitch))
            Yaw = ATAN2(Rotation_L_I(2,1)/cos(Pitch), Rotation_L_I(1,1)/cos(Pitch))  
        ELSE
            Yaw = 0
            IF(Rotation_L_I(3,1) == 1) THEN
                Pitch = PiBy2_D
                Roll = Yaw + ATAN2(Rotation_L_I(1,2),Rotation_L_I(1,3))
            ELSE
                Pitch = -PiBy2_D
                Roll = -Yaw + ATAN2(-Rotation_L_I(1,2),-Rotation_L_I(1,3))
            END IF
        END IF
    
        LidarPosition_I = MATMUL(TRANSPOSE(u%NacelleMotion%Orientation(:,:,1)),p%LidarPosition_N) !Rotates the position vector to the orientation of the inertial coord. system
    
	END IF
	
	
        y%IMUOutputs(1)     = Roll                              !Roll Angle
        y%IMUOutputs(2)     = u%NacelleMotion%RotationVel(1,1)  !Roll Angle Velocity
        y%IMUOutputs(3)     = u%NacelleMotion%RotationAcc(1,1)  !Roll Angle Acceleration
        y%IMUOutputs(4)     = Pitch                             !Pitch Angle
        y%IMUOutputs(5)     = u%NacelleMotion%RotationVel(2,1)  !Pitch Angle Velocity
        y%IMUOutputs(6)     = u%NacelleMotion%RotationAcc(2,1)  !Pitch Angle Acceleration
        y%IMUOutputs(7)     = Yaw                               !Yaw Angle
        y%IMUOutputs(8)     = u%NacelleMotion%RotationVel(3,1)  !Yaw Angle Velocity
        y%IMUOutputs(9)     = u%NacelleMotion%RotationAcc(3,1)  !Yaw Angle Acceleration
    
        y%IMUOutputs(10)    = u%NacelleMotion%TranslationDisp(1,1)  !Displacement x 
        y%IMUOutputs(13)    = u%NacelleMotion%TranslationDisp(2,1)  !Displacement y
        y%IMUOutputs(16)    = u%NacelleMotion%TranslationDisp(3,1)  !Displacement z
    
        DisplacementNacelle(1) = u%NacelleMotion%TranslationDisp(1,1)
        DisplacementNacelle(2) = u%NacelleMotion%TranslationDisp(2,1)
        DisplacementNacelle(3) = u%NacelleMotion%TranslationDisp(3,1)
    
        DisplacementLidar   = LidarPosition_I + DisplacementNacelle
    
        y%IMUOutputs(10)    = DisplacementLidar(1)
        y%IMUOutputs(13)    = DisplacementLidar(2)
        y%IMUOutputs(16)    = DisplacementLidar(3)
    
        CrossProduct(1)     = u%NacelleMotion%RotationVel(2,1)*LidarPosition_I(3) - u%NacelleMotion%RotationVel(3,1)*LidarPosition_I(2)
        CrossProduct(2)     = u%NacelleMotion%RotationVel(3,1)*LidarPosition_I(1) - u%NacelleMotion%RotationVel(1,1)*LidarPosition_I(3)
        CrossProduct(3)     = u%NacelleMotion%RotationVel(1,1)*LidarPosition_I(2) - u%NacelleMotion%RotationVel(2,1)*LidarPosition_I(1)
    
        y%IMUOutputs(11)    = u%NacelleMotion%TranslationVel(1,1) + CrossProduct(1)    !Velocity x
        y%IMUOutputs(14)    = u%NacelleMotion%TranslationVel(2,1) + CrossProduct(2)    !Velocity y
        y%IMUOutputs(17)    = u%NacelleMotion%TranslationVel(3,1) + CrossProduct(3)    !Velocity z
    
        CrossProduct(1)     = u%NacelleMotion%RotationAcc(2,1)*LidarPosition_I(3) - u%NacelleMotion%RotationAcc(3,1)*LidarPosition_I(2)
        CrossProduct(2)     = u%NacelleMotion%RotationAcc(3,1)*LidarPosition_I(1) - u%NacelleMotion%RotationAcc(1,1)*LidarPosition_I(3)
        CrossProduct(3)     = u%NacelleMotion%RotationAcc(1,1)*LidarPosition_I(2) - u%NacelleMotion%RotationAcc(2,1)*LidarPosition_I(1)    
    
        y%IMUOutputs(12)    = u%NacelleMotion%TranslationAcc(1,1) + CrossProduct(1)     !Acceleration x
        y%IMUOutputs(15)    = u%NacelleMotion%TranslationAcc(2,1) + CrossProduct(2)     !Acceleration y
        y%IMUOutputs(18)    = u%NacelleMotion%TranslationAcc(3,1) + CrossProduct(3)     !Acceleration z
    

    END SUBROUTINE LidarSim_CalculateIMU
    
    
!####################################################
    
    SUBROUTINE LidarSim_InitializeOutputs(y,p, InitOutData, InputFileData, ErrStat, ErrMsg)
    
    IMPLICIT                               NONE
    CHARACTER(*),                          PARAMETER       ::  RoutineName='LidarSim_InitializeOutputs'
    
    TYPE(LidarSim_OutputType),             INTENT(  OUT)   ::  y                   ! Parameter for lidar outputs
    TYPE(LidarSim_ParameterType),          INTENT(INOUT)   ::  p                   ! Parameter data for the lidar module
    TYPE(LidarSim_InitOutputType),         INTENT(  OUT)   ::  InitOutData
    TYPE(LidarSim_InputFile),              INTENT(IN   )   ::  InputFileData
    INTEGER(IntKi),                        INTENT(  OUT)   ::  ErrStat              !< Temporary error status
    CHARACTER(ErrMsgLen),                  INTENT(  OUT)   ::  ErrMsg               !< temporary error message
    
    !Local error variables
    INTEGER(IntKi)          ::  TmpErrStat
    CHARACTER(ErrMsgLen)    ::  TmpErrMsg           !< temporary error message
    
    !Local variables
    INTEGER(IntKi)                          ::  SizeOutput
    INTEGER(IntKi)                          ::  LoopCounter
    INTEGER(IntKi)                          ::  NumberValidOutputs
    INTEGER(IntKi),DIMENSION(:),ALLOCATABLE ::  ValidOutputs
    CHARACTER(15) ,DIMENSION(:),ALLOCATABLE ::  ValidOutputNames
    CHARACTER(1024)                         ::  TmpIntegerToString
    CHARACTER(15)                           ::  OutListTmp
    
    INTEGER(IntKi), PARAMETER               ::  ROLLLI      =   1
    INTEGER(IntKi), PARAMETER               ::  ROLLDTLI    =   2
    INTEGER(IntKi), PARAMETER               ::  ROLLDTDTLI  =   3
    INTEGER(IntKi), PARAMETER               ::  PTCHLI      =   4
    INTEGER(IntKi), PARAMETER               ::  PTCHDTLI    =   5
    INTEGER(IntKi), PARAMETER               ::  PTCHDTDTLI  =   6
    INTEGER(IntKi), PARAMETER               ::  YAWLI       =   7
    INTEGER(IntKi), PARAMETER               ::  YAWDTLI     =   8
    INTEGER(IntKi), PARAMETER               ::  YAWDTDTLI   =   9
    INTEGER(IntKi), PARAMETER               ::  XLI         =   10
    INTEGER(IntKi), PARAMETER               ::  XDTLI       =   11
    INTEGER(IntKi), PARAMETER               ::  XDTDTLI     =   12
    INTEGER(IntKi), PARAMETER               ::  YLI         =   13
    INTEGER(IntKi), PARAMETER               ::  YDTLI       =   14
    INTEGER(IntKi), PARAMETER               ::  YDTDTLI     =   15
    INTEGER(IntKi), PARAMETER               ::  ZLI         =   16
    INTEGER(IntKi), PARAMETER               ::  ZDTLI       =   17
    INTEGER(IntKi), PARAMETER               ::  ZDTDTLI     =   18
    INTEGER(IntKi), PARAMETER               ::  AZIMUTHLI   =   19
    INTEGER(IntKi), PARAMETER               ::  ELEVATLI    =   20
    INTEGER(IntKi), PARAMETER               ::  MEASTIMELI  =   21
    INTEGER(IntKi), PARAMETER               ::  NEWDATALI   =   22
    INTEGER(IntKi), PARAMETER               ::  BEAMIDLI    =   23
    INTEGER(IntKi), PARAMETER               ::  BLDBLOCKLI    =   24
        
    CHARACTER(ChanLen), DIMENSION(:), ALLOCATABLE  :: ValidParamAry     ! This lists the names of the allowed parameters, which must be sorted alphabetically
    INTEGER(IntKi), DIMENSION(:), ALLOCATABLE  :: ParamIndexAry     ! This lists the index into AllOuts(:) of the allowed parameters ValidParamAry(:)
    CHARACTER(ChanLen), DIMENSION(:), ALLOCATABLE  :: ParamUnitsAry     ! This lists the units corresponding to the allowed parameters
    
    !Error Variables
    TmpErrStat  = 0
    TmpErrMsg   = ''
    
    SizeOutput = InputFileData%GatesPerBeam
        
    CALL AllocAry( ValidParamAry, 24+(2*SizeOutput), 'ValidParamAry', TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( ParamIndexAry, 24+(2*SizeOutput), 'ParamIndexAry', TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( ParamUnitsAry, 24+(2*SizeOutput), 'ParamUnitsAry', TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    
    ! Fill ValidParamAry and ParamIndexAry in alphabetical order
    ValidParamAry( 1 : 8 ) = (/ & 
        "AZIMUTHLI" ,"BEAMIDLI"  ,"ELEVATLI"  ,"MEASTIMELI","NEWDATALI" ,"PTCHDTDTLI", & 
        "PTCHDTLI"  ,"PTCHLI"    /)
    ParamIndexAry( 1 : 8 ) = (/ & 
        AZIMUTHLI   ,BEAMIDLI    ,ELEVATLI    ,MEASTIMELI  ,NEWDATALI   ,PTCHDTDTLI  , & 
        PTCHDTLI    ,PTCHLI      /)
    DO LoopCounter = 1,SizeOutput
        WRITE(UNIT=TmpIntegerToString,FMT='(I2.2)') LoopCounter
        ValidParamAry( 8 + LoopCounter ) = "RANGE"//TRIM(ADJUSTL(TmpIntegerToString))//"LI"
        ParamIndexAry( 8 + LoopCounter ) = 24 + LoopCounter
    ENDDO
    ValidParamAry( (9 + SizeOutput) : (11 + SizeOutput) ) = (/ & 
        "ROLLDTDTLI","ROLLDTLI"  ,"ROLLLI"    /)
    ParamIndexAry( (9 + SizeOutput) : (11 + SizeOutput) ) = (/ & 
        ROLLDTDTLI  ,ROLLDTLI    ,ROLLLI      /)
    DO LoopCounter = 1,SizeOutput
        WRITE(UNIT=TmpIntegerToString,FMT='(I2.2)') LoopCounter
        ValidParamAry( 11 + SizeOutput + LoopCounter ) = "VLOS"//TRIM(ADJUSTL(TmpIntegerToString))//"LI"
        ParamIndexAry( 11 + SizeOutput + LoopCounter ) = 24 + SizeOutput + LoopCounter
    ENDDO
    ValidParamAry( (12 + (2*SizeOutput)) : (24 + (2*SizeOutput)) ) = (/ & 
        "XDTDTLI"   ,"XDTLI"     ,"XLI"       ,"YAWDTDTLI" ,"YAWDTLI"   ,"YAWLI"     , & 
        "YDTDTLI"   ,"YDTLI"     ,"YLI"       ,"ZDTDTLI"   ,"ZDTLI"     ,"ZLI"       ,&
        "BLDBLOCKLI"       /)
    ParamIndexAry( (12 + (2*SizeOutput)) : (24 + (2*SizeOutput)) ) = (/ & 
        XDTDTLI     ,XDTLI       ,XLI         ,YAWDTDTLI   ,YAWDTLI     ,YAWLI       , & 
        YDTDTLI     ,YDTLI       ,YLI         ,ZDTDTLI     ,ZDTLI       ,ZLI         ,&
        BLDBLOCKLI        /)
    
    ! Fill ParamUnitsAry according to parameter order
    ParamUnitsAry( 1 : 24 ) = (/ & 
        "(rad)     ","(rad/s)   ","(rad/s^2) ","(rad)     ","(rad/s)   ","(rad/s^2) ", & 
        "(rad)     ","(rad/s)   ","(rad/s^2) ","(m)       ","(m/s)     ","(m/s^2)   ", & 
        "(m)       ","(m/s)     ","(m/s^2)   ","(m)       ","(m/s)     ","(m/s^2)   ", & 
        "(rad)     ","(rad)     ","(s)       ","()        ","()        ","()        "/)
    DO LoopCounter = 1,SizeOutput
        ParamUnitsAry( 24 + LoopCounter ) = "(m)       "
        ParamUnitsAry( 24 + SizeOutput + LoopCounter ) = "(m/s)     "
    ENDDO
    
    CALL AllocAry( y%IMUOutputs, 18, 'IMUOutputs', TmpErrStat, TmpErrMsg )                              !Allocate array for all IMU data ( rotation and position displacement, velocity, acceleration)
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( ValidOutputs, 24+(2*SizeOutput), 'ValidOutputs', TmpErrStat, TmpErrMsg )             !Allocate the maximal additional outputchannels
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( ValidOutputNames, 24+(2*SizeOutput), 'ValidOutputNames', TmpErrStat, TmpErrMsg )     !Allocate the maximal additional outputchannels
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    
    NumberValidOutputs=0
    DO LoopCounter = 1, InputFileData%NumOuts
        OutListTmp = InputFileData%OutList(LoopCounter)
        CALL Conv2UC(OutListTmp)
        ValidOutputs(NumberValidOutputs+1) = IndexCharAry(OutListTmp,ValidParamAry)
        IF(ValidOutputs(NumberValidOutputs+1) > 0) THEN
            ValidOutputs(NumberValidOutputs+1)      =   ParamIndexAry(ValidOutputs(NumberValidOutputs+1))
            ValidOutputNames(NumberValidOutputs+1)  =   InputFileData%OutList(LoopCounter)
            NumberValidOutputs = NumberValidOutputs + 1            
        ENDIF
    ENDDO
    
    IF(NumberValidOutputs>0) THEN
        CALL AllocAry( p%ValidOutputs, NumberValidOutputs, 'p%ValidOutputs', TmpErrStat, TmpErrMsg )                     !Allocate the fitting amount of outputchannels
        CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
        p%ValidOutputs=ValidOutputs(1:NumberValidOutputs)   
    ENDIF
    
    CALL AllocAry( y%WriteOutput, NumberValidOutputs, 'WriteOutput', TmpErrStat, TmpErrMsg )                     !Allocate the actual output array
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( InitOutData%WriteOutputHdr, NumberValidOutputs, 'WriteOutputHdr', TmpErrStat, TmpErrMsg )     !Name of the data output channels   (Size of the WriteOutputHdr array defines the number of outputs)
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    CALL AllocAry( InitOutData%WriteOutputUnt, NumberValidOutputs, 'WriteOutputUnt', TmpErrStat, TmpErrMsg )     !units of the output channels
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    
    DO LoopCounter = 1, NumberValidOutputs
         y%WriteOutput( LoopCounter ) = 0
         InitOutData%WriteOutputHdr( LoopCounter ) = ValidOutputNames(LoopCounter)
         InitOutData%WriteOutputUnt( LoopCounter ) = ParamUnitsAry(p%ValidOutputs(LoopCounter))
    ENDDO
    
    DEALLOCATE(ValidOutputNames)
    DEALLOCATE(ValidOutputs)
    
    ! Initialize SwapOutputs array
    CALL AllocAry( y%SwapOutputs, 2+SizeOutput+7, 'SwapOutputs', TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    
    y%SwapOutputs(1) = 0                        !NewData
    y%SwapOutputs(2) = 0                        !BeamID
    
    DO LoopCounter = 1,SizeOutput
         y%SwapOutputs(2 + LoopCounter) = 0     !Vlos
    ENDDO

    y%SwapOutputs(2 + SizeOutput + 1)   = 0     ! LdrRoll
    y%SwapOutputs(2 + SizeOutput + 2)   = 0     ! LdrPitch
    y%SwapOutputs(2 + SizeOutput + 3)   = 0     ! LdrYaw
    y%SwapOutputs(2 + SizeOutput + 4)   = 0     ! LdrXd
    y%SwapOutputs(2 + SizeOutput + 5)   = 0     ! LdrYd 
    y%SwapOutputs(2 + SizeOutput + 6)   = 0     ! LdrZd 
    y%SwapOutputs(2 + SizeOutput + 7)   = 0     ! LdrBladeBlock
    
    ! Initialize AllOutputs array
    CALL AllocAry( y%AllOutputs, 24+(2*SizeOutput), 'AllOutputs', TmpErrStat, TmpErrMsg )
    CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
    
    y%AllOutputs = 0
    
    END SUBROUTINE LidarSim_InitializeOutputs
    
    
!###########################################################
    
    SUBROUTINE LidarSim_SetOutputs(y,p,Vlos,BladeBlockageStatus,UnitVector,LoopGatesPerBeam,Time)
    
    IMPLICIT                                    NONE
    CHARACTER(*),                               PARAMETER       ::  RoutineName='LidarSim_SetOutputs'
    
    TYPE(LidarSim_OutputType),                  INTENT(INOUT)   ::  y
    TYPE(LidarSim_ParameterType),               INTENT(IN   )   ::  p
    REAL(ReKi),                                 INTENT(IN   )   ::  Vlos
    REAL(ReKi),                                 INTENT(IN   )   ::  BladeBlockageStatus
    REAL(ReKi),                                 INTENT(IN   )   ::  UnitVector(3)
    INTEGER(IntKi),                             INTENT(IN   )   ::  LoopGatesPerBeam                    !from 0 to p%GatesPerBeam-1
    REAL(DbKi),                                 INTENT(IN   )   ::  Time
    
    !Local variables
    INTEGER(IntKi)                                              ::  LoopCounter
    REAL(ReKi)                                                  ::  Dot_LidarPosition_I(3)
    
    Dot_LidarPosition_I(1) = y%IMUOutputs(11)
    Dot_LidarPosition_I(2) = y%IMUOutputs(14)
    Dot_LidarPosition_I(3) = y%IMUOutputs(17)
    
    y%AllOutputs( 1 : 18 )  = y%IMUOutputs ( 1 : 18 )
    y%AllOutputs( 19 )      = p%MeasuringPoints_Spherical_L(2,p%LastMeasuringPoint+LoopGatesPerBeam)        !Azimuth
    y%AllOutputs( 20 )      = p%MeasuringPoints_Spherical_L(3,p%LastMeasuringPoint+LoopGatesPerBeam)        !Elevation
    y%AllOutputs( 21 )      = Time                                                                          !MeasTime
    y%AllOutputs( 22 )      = 1                                                                             !NewData
    y%AllOutputs( 23 )      = REAL(p%NextBeamID)                                                            !BeamID
    y%AllOutputs( 24 )      = BladeBlockageStatus                                                           !BladeBlockStatus

    y%AllOutputs( 25 + LoopGatesPerBeam ) = p%MeasuringPoints_Spherical_L(1,p%LastMeasuringPoint+LoopGatesPerBeam)   !Rangegates
    y%AllOutputs( 25 + p%GatesPerBeam + LoopGatesPerBeam ) = Vlos + (  DOT_PRODUCT(Dot_LidarPosition_I,UnitVector))  !Output the measured V_los. Consiting of the windspeed and the movement of the measuring system itself note here the LOS towards lidar is positive
    
    DO LoopCounter = 1,SIZE(p%ValidOutputs)
        y%WriteOutput( LoopCounter ) = y%AllOutputs( p%ValidOutputs(LoopCounter) )
    END DO
    
    y%SwapOutputs( 1                        )   = y%AllOutputs( 22 )                                        !NewData
    y%SwapOutputs( 2                        )   = y%AllOutputs( 23 )                                        !BeamID
    y%SwapOutputs( 3 + LoopGatesPerBeam     )   = y%AllOutputs( 25 + p%GatesPerBeam + LoopGatesPerBeam )    !Vlos

    y%SwapOutputs( 2 + p%GatesPerBeam   + 1 )   = y%IMUOutputs(  1 )    ! LdrRoll
    y%SwapOutputs( 2 + p%GatesPerBeam   + 2 )   = y%IMUOutputs(  4 )    ! LdrPitch
    y%SwapOutputs( 2 + p%GatesPerBeam   + 3 )   = y%IMUOutputs(  7 )    ! LdrYaw
    y%SwapOutputs( 2 + p%GatesPerBeam   + 4 )   = y%IMUOutputs( 11 )    ! LdrXd
    y%SwapOutputs( 2 + p%GatesPerBeam   + 5 )   = y%IMUOutputs( 14 )    ! LdrYd
    y%SwapOutputs( 2 + p%GatesPerBeam   + 6 )   = y%IMUOutputs( 17 )    ! LdrZd
    y%SwapOutputs( 2 + p%GatesPerBeam   + 7 )   = y%AllOutputs( 24 )    ! LdrBladeBlock
    END SUBROUTINE
    

!###################################################
    
    SUBROUTINE LidarSim_CheckBladeBlockage_14(p,y,u,LidarPosition_I,MeasuringPosition_I,BladeBlockageStatus)
    
    
    IMPLICIT                                   NONE
    CHARACTER(*),                              PARAMETER            ::  RoutineName="LidarSim_CheckBladeBlockage_14"
    
    TYPE(LidarSim_ParameterType),              INTENT(INOUT)        ::  p
    TYPE(LidarSim_OutputType),                 INTENT(INOUT)        ::  y                       !Outputs computed at Time (IN for mesh reasons and data allocation)
    TYPE(LidarSim_InputType),                  INTENT(IN   )        ::  u   

    REAL(ReKi),                                INTENT(IN   )        ::  LidarPosition_I(3)             ! the lidar position, origin of the lidar beam
    REAL(ReKi),                                INTENT(IN   )        ::  MeasuringPosition_I(3)          !Transformed Measuring Position
    REAL(ReKi),                                INTENT(  OUT)        ::  BladeBlockageStatus            ! 1 block 0 not block
   
    
    ! local variables
    
    
    REAL(ReKi)                                                      ::  AbsPitch                     ! the absolute pitch, the angle between chord line and the tangential diraction of rotation
    REAL(ReKi)                                                      ::  CL
    REAL(ReKi)                                                      ::  AeroCenterCoordinate_inner_I(3)     !coordinate of one aerodynamic center position
    REAL(ReKi)                                                      ::  AeroCenterCoordinate_outer_I(3)     !coordinate of one aerodynamic center position
    REAL(ReKi)                                                      ::  LeadEdgeCoordinate_inner_I(3)       !coordinate of one leading edge position for one blade element
    REAL(ReKi)                                                      ::  TrailEdgeCoordinate_inner_I(3)      !coordinate of one trailing edge position for one blade element
    REAL(ReKi)                                                      ::  LeadEdgeCoordinate_outer_I(3)       !coordinate of one leading edge position for one blade element
    REAL(ReKi)                                                      ::  TrailEdgeCoordinate_outer_I(3)      !coordinate of one trailing edge position for one blade element
    REAL(ReKi)                                                      ::  IElement                      ! index number of the blade elements
    REAL(ReKi)                                                      ::  NELM                          ! number of blade elements
    REAL(ReKi)                                                      ::  CLprojectionLength            ! the length of chord projected to yz plane
    REAL(ReKi)                                                      ::  IBlade                        ! index number of the blade 
    REAL(ReKi)                                                      ::  NBlade                        ! number of the blade 
    REAL(ReKi)                                                      ::  Beam_at_LeadEdge(3)           ! the  coordinate of the beam at the x plane where leading edge located 
    REAL(ReKi)                                                      ::  Beam_at_TrailEdge(3)          ! the coordinate of the beam at the x plane where trailing edge located 
    !REAL(ReKi)                                                      ::  Orientation(3,3)              ! the coordinate of the beam at the x plane where trailing edge located 
    REAL(ReKi)                                                      ::  vpt                           ! the dot product between lidar beam vector and the normal vector of rotor plane
    REAL(ReKi)                                                      ::  linevector(3)                 ! the vector of lidar beam
    REAL(ReKi)                                                      ::  LidarBeamInRotor_I(3)         ! the coordinate of the beam point that locates in the rotor plane
    REAL(ReKi)                                                      ::  TEMP                          ! a temporary variable   
    REAL(ReKi)                                                      ::  r_lidarBeamInRotor_I          ! the radius of lidar beam in rotor point relative to the rotor center
    REAL(ReKi)                                                      ::  r_node_inner                  ! the radius of blade element node relative to the rotor center (closer to center)
    REAL(ReKi)                                                      ::  r_node_outer                  ! the radius of blade element node relative to the rotor center furter to center)
    REAL(ReKi)                                                      ::  edge1(3)                      ! the first vector in the triangle 
    REAL(ReKi)                                                      ::  edge2(3)                      ! the second vector in the triangl, share same point with first vector
    REAL(ReKi)                                                      ::  tvec(3)                       ! the vector from the sharing point of the triangle to the ray origin
    REAL(ReKi)                                                      ::  pvec(3)                       ! intermediate vector
    REAL(ReKi)                                                      ::  qvec(3)                       ! intermediate vector
    REAL(ReKi)                                                      ::  det                           ! determinant of the matrix M, M = dot(edge1,pvec)
    REAL(ReKi)                                                      ::  tb                             ! 
    REAL(ReKi)                                                      ::  vb                             ! 
    REAL(ReKi)                                                      ::  ub                             ! 
    REAL(ReKi)                                                      ::  TriangleCoordinate(3,3)
    REAL(ReKi)                                                      ::  HubOrientationX(3)             ! 
    REAL(ReKi)                                                      ::  HubPosition(3)             ! 
    
    NBlade       =      P%NBLADE	! number of blade 
    NELM         = 		P%NELM   ! number of blade element
    CL           =      1
    IBlade       =      1 
    !   the absolute pitch which is the angle between chordline and tangential direction of rotation 
    !Orientation    = AD14_u%TurbineComponents%Blade(1)%Orientation
    BladeBlockageStatus    = 0     ! by default there is no blockage
    
    HubOrientationX = (/U%HUBMOTION%ORIENTATION(1,1,1),U%HUBMOTION%ORIENTATION(1,2,1),U%HUBMOTION%ORIENTATION(1,3,1)/)
    HubPosition     = (/U%HUBMOTION%POSITION(1,1),U%HUBMOTION%POSITION(2,1),U%HUBMOTION%POSITION(3,1)/)
    
    !  firt find the xyz coordinate in inertial system of the lidar beam in the rotor plane, we assume the rotor plane pass the hubcenter node, and the first orientation vector is normal to the rotor plane
    linevector   =       MeasuringPosition_I-LidarPosition_I
    vpt          =      DOT_PRODUCT(linevector,HubOrientationX)    ! the cos of the angle between lidar beam line vector and the rotor plane norm vector
    
    IF (vpt /= 0) THEN  ! the lidar beam is not in parallel with the rotor swept plane
    
       TEMP     = DOT_PRODUCT(HubPosition-LidarPosition_I,HubOrientationX)/vpt   ! temporary vector 
       LidarBeamInRotor_I(1)     = LidarPosition_I(1)+linevector(1)*TEMP   ! the intersection of lidar beam and the rotor plane
       LidarBeamInRotor_I(2)     = LidarPosition_I(2)+linevector(2)*TEMP
       LidarBeamInRotor_I(3)     = LidarPosition_I(3)+linevector(3)*TEMP
       
       
       Do IBlade = 1,NBlade
           Do IElement = 1,NELM-1
        !+		U%BLADEMOTION(1)%POSITION	(1:3,1:19)	REAL(4) 

        
        IF  (BladeBlockageStatus == 0) THEN      ! if Blockage Status is 1 which means it is not blocked by previous blade element, then we should enter the loop
               
                    AeroCenterCoordinate_inner_I   = U%BLADEMOTION(IBlade)%POSITION(:,IElement)      ! coordinate of the aerodynamic center of a blade element node
                    r_node_inner                   = NORM2(AeroCenterCoordinate_inner_I -HubPosition)  ! assume the aero center is in the swept plane
                    AeroCenterCoordinate_outer_I   = U%BLADEMOTION(IBlade)%POSITION(:,IElement+1)      ! coordinate of the aerodynamic center of a blade element node
                    r_node_outer                   = NORM2(AeroCenterCoordinate_outer_I-HubPosition)  ! assume the aero center is in the swept plane
                    r_lidarBeamInRotor_I           = NORM2(LidarBeamInRotor_I-HubPosition)
                    
                    IF ((r_lidarBeamInRotor_I>r_node_inner).AND.(r_lidarBeamInRotor_I<=r_node_outer)) THEN        
                        CL                       = P%BLADEELMCL(IElement) 	 !AD14_P%BLADE%C(IElement)   ! chord length of the node
                    
                        ! two edges for the inner node
                        LeadEdgeCoordinate_inner_I (1)  = AeroCenterCoordinate_inner_I(1)-U%BLADEMOTION(IBlade)%Orientation(2,1,IElement)*CL *0.25
                        LeadEdgeCoordinate_inner_I (2)  = AeroCenterCoordinate_inner_I(2)-U%BLADEMOTION(IBlade)%Orientation(2,2,IElement)*CL *0.25
                        LeadEdgeCoordinate_inner_I (3)  = AeroCenterCoordinate_inner_I(3)-U%BLADEMOTION(IBlade)%Orientation(2,3,IElement)*CL *0.25
                        TrailEdgeCoordinate_inner_I (1) = AeroCenterCoordinate_inner_I(1)+U%BLADEMOTION(IBlade)%Orientation(2,1,IElement)*CL *0.75
                        TrailEdgeCoordinate_inner_I (2) = AeroCenterCoordinate_inner_I(2)+U%BLADEMOTION(IBlade)%Orientation(2,2,IElement)*CL *0.75    
                        TrailEdgeCoordinate_inner_I (3) = AeroCenterCoordinate_inner_I(3)+U%BLADEMOTION(IBlade)%Orientation(2,3,IElement)*CL *0.75    
                    
                        ! two edges for the outer node
                        LeadEdgeCoordinate_outer_I(1)  = AeroCenterCoordinate_outer_I(1)-U%BLADEMOTION(IBlade)%Orientation(2,1,IElement+1)*CL *0.25
                        LeadEdgeCoordinate_outer_I(2)  = AeroCenterCoordinate_outer_I(2)-U%BLADEMOTION(IBlade)%Orientation(2,2,IElement+1)*CL *0.25
                        LeadEdgeCoordinate_outer_I(3)  = AeroCenterCoordinate_outer_I(3)-U%BLADEMOTION(IBlade)%Orientation(2,3,IElement+1)*CL *0.25
                        TrailEdgeCoordinate_outer_I(1) = AeroCenterCoordinate_outer_I(1)+U%BLADEMOTION(IBlade)%Orientation(2,1,IElement+1)*CL *0.75
                        TrailEdgeCoordinate_outer_I(2) = AeroCenterCoordinate_outer_I(2)+U%BLADEMOTION(IBlade)%Orientation(2,2,IElement+1)*CL *0.75    
                        TrailEdgeCoordinate_outer_I(3) = AeroCenterCoordinate_outer_I(3)+U%BLADEMOTION(IBlade)%Orientation(2,3,IElement+1)*CL *0.75  
                    
                        ! formulate first triangle
                        TriangleCoordinate(1,:)        = LeadEdgeCoordinate_inner_I
                        TriangleCoordinate(2,:)        = TrailEdgeCoordinate_inner_I
                        TriangleCoordinate(3,:)        = LeadEdgeCoordinate_outer_I
                    
                         !check weather the first triangle is cross with the beam or not
                        edge1                          = TriangleCoordinate(2,:)-TriangleCoordinate(1,:)         
                        edge2                          = TriangleCoordinate(3,:)-TriangleCoordinate(1,:)         
                        tvec                           = LidarPosition_I-TriangleCoordinate(1,:)   
                        pvec                           = CROSS_PRODUCT(linevector, edge2)
                        det                            = DOT_PRODUCT(edge1,pvec)
                    
                        IF (ABS(det)>0.00001) THEN     ! there is intersection in the triangle plane, then we check whether the point is inside the triangle
                                
                                ub    = DOT_PRODUCT(tvec,pvec)/det;       ! 1st barycentric coordinate
                                qvec  = CROSS_PRODUCT(tvec, edge1);             ! prepare to test V parameter
                                vb    = DOT_PRODUCT(linevector,qvec)/det;       ! 2nd barycentric coordinate
                                tb    = DOT_PRODUCT(edge2,qvec)/det;            ! 'position on the line' coordinate
                                 !test if line/plane intersection is within the triangle
                                IF   ((ub>=-0.0) .AND. (vb>=-0.0) .AND. (ub+vb<=1.0)) THEN 
                                    BladeBlockageStatus = 1 
                                ELSE 
                                END IF  ! the point is inside the triangle the beam is blocked
                                
                    
                        ELSE   ! there is no intersection for this triangle, we will check next triangle
                    
                                !formulate second triangle
                                TriangleCoordinate(1,:)        = TrailEdgeCoordinate_inner_I
                                TriangleCoordinate(2,:)        = LeadEdgeCoordinate_outer_I
                                TriangleCoordinate(3,:)        = TrailEdgeCoordinate_outer_I
                                !check weather the second triangle is cross with the beam or not
                                edge1                          = TriangleCoordinate(2,:)-TriangleCoordinate(1,:)         
                                edge2                          = TriangleCoordinate(3,:)-TriangleCoordinate(1,:)         
                                tvec                           = LidarPosition_I-TriangleCoordinate(1,:)   
                                pvec                           = CROSS_PRODUCT(linevector, edge2)
                                det                            = DOT_PRODUCT(edge1,pvec)
                                
                                IF (ABS(det)>0.00001) THEN     ! there is intersection in the second triangle plane, then we check whether the point is inside the triangle
                               
                                qvec = CROSS_PRODUCT(tvec, edge1);             ! prepare to test V parameter
                                vb    = DOT_PRODUCT(linevector,qvec)/det;       ! 2nd barycentric coordinate
                                tb    = DOT_PRODUCT(edge2,qvec)/det;            ! 'position on the line' coordinate
                                 !test if line/plane intersection is within the triangle
                                   IF   ((ub>=-0.0) .AND. (vb>=-0.0) .AND. (ub+vb<=1.0)) THEN 
                                       BladeBlockageStatus = 1 
                                   ELSE  
                                   END IF  ! the point is inside the triangle the beam is blocked
                                
                                ELSE
                                END IF
                    
                        END IF
                    ELSE
                    END IF ! if the radius is still between inner and outer
                    !
                  ELSE   ! the loop should be existed 
                  END IF  !IF  BlockFlag == 0 
                  
           END DO !IElement = 1,NELM-1
        END DO !IBlade = 1,NBlade
        
    ELSE    ! (vpt == 0)
 
    END IF

    END SUBROUTINE LidarSim_CheckBladeBlockage_14

!###################################################
    
    SUBROUTINE LidarSim_CheckBladeBlockage_15(p,y,u,LidarPosition_I,MeasuringPosition_I,BladeBlockageStatus)
    
    
    IMPLICIT                                   NONE
    CHARACTER(*),                              PARAMETER            ::  RoutineName="LidarSim_CheckBladeBlockage_15"
    
    TYPE(LidarSim_ParameterType),              INTENT(INOUT)        ::  p
    TYPE(LidarSim_OutputType),                 INTENT(INOUT)        ::  y                       !Outputs computed at Time (IN for mesh reasons and data allocation)
    TYPE(LidarSim_InputType),                  INTENT(IN   )        ::  u   
    REAL(ReKi),                                INTENT(IN   )        ::  LidarPosition_I(3)             ! the lidar position, origin of the lidar beam
    REAL(ReKi),                                INTENT(IN   )        ::  MeasuringPosition_I(3)          !Transformed Measuring Position
    REAL(ReKi),                                INTENT( OUT)        ::  BladeBlockageStatus            ! 1 block 0 not block
   
    
    ! local variables
    
    
    REAL(ReKi)                                                      ::  AbsPitch                     ! the absolute pitch, the angle between chord line and the tangential diraction of rotation
    REAL(ReKi)                                                      ::  CL
    REAL(ReKi)                                                      ::  AeroCenterCoordinate_inner_I(3)     !coordinate of one aerodynamic center position
    REAL(ReKi)                                                      ::  AeroCenterCoordinate_outer_I(3)     !coordinate of one aerodynamic center position
    REAL(ReKi)                                                      ::  LeadEdgeCoordinate_inner_I(3)       !coordinate of one leading edge position for one blade element
    REAL(ReKi)                                                      ::  TrailEdgeCoordinate_inner_I(3)      !coordinate of one trailing edge position for one blade element
    REAL(ReKi)                                                      ::  LeadEdgeCoordinate_outer_I(3)       !coordinate of one leading edge position for one blade element
    REAL(ReKi)                                                      ::  TrailEdgeCoordinate_outer_I(3)      !coordinate of one trailing edge position for one blade element
    REAL(ReKi)                                                      ::  IElement                      ! index number of the blade elements
    REAL(ReKi)                                                      ::  NELM                          ! number of blade elements
    REAL(ReKi)                                                      ::  CLprojectionLength            ! the length of chord projected to yz plane
    REAL(ReKi)                                                      ::  IBlade                        ! index number of the blade 
    REAL(ReKi)                                                      ::  NBlade                        ! number of the blade 
    REAL(ReKi)                                                      ::  Beam_at_LeadEdge(3)           ! the  coordinate of the beam at the x plane where leading edge located 
    REAL(ReKi)                                                      ::  Beam_at_TrailEdge(3)          ! the coordinate of the beam at the x plane where trailing edge located 
    !REAL(ReKi)                                                      ::  Orientation(3,3)              ! the coordinate of the beam at the x plane where trailing edge located 
    REAL(ReKi)                                                      ::  vpt                           ! the dot product between lidar beam vector and the normal vector of rotor plane
    REAL(ReKi)                                                      ::  linevector(3)                 ! the vector of lidar beam
    REAL(ReKi)                                                      ::  LidarBeamInRotor_I(3)         ! the coordinate of the beam point that locates in the rotor plane
    REAL(ReKi)                                                      ::  TEMP                          ! a temporary variable   
    REAL(ReKi)                                                      ::  r_lidarBeamInRotor_I          ! the radius of lidar beam in rotor point relative to the rotor center
    REAL(ReKi)                                                      ::  r_node_inner                  ! the radius of blade element node relative to the rotor center (closer to center)
    REAL(ReKi)                                                      ::  r_node_outer                  ! the radius of blade element node relative to the rotor center furter to center)
    REAL(ReKi)                                                      ::  edge1(3)                      ! the first vector in the triangle 
    REAL(ReKi)                                                      ::  edge2(3)                      ! the second vector in the triangl, share same point with first vector
    REAL(ReKi)                                                      ::  tvec(3)                       ! the vector from the sharing point of the triangle to the ray origin
    REAL(ReKi)                                                      ::  pvec(3)                       ! intermediate vector
    REAL(ReKi)                                                      ::  qvec(3)                       ! intermediate vector
    REAL(ReKi)                                                      ::  det                           ! determinant of the matrix M, M = dot(edge1,pvec)
    REAL(ReKi)                                                      ::  tb                             ! 
    REAL(ReKi)                                                      ::  vb                             ! 
    REAL(ReKi)                                                      ::  ub                             ! 
    REAL(ReKi)                                                      ::  TriangleCoordinate(3,3)
    REAL(ReKi)                                                      ::  HubOrientationX(3)             ! 
    REAL(ReKi)                                                      ::  HubPosition(3)             ! 
    
    NBlade       =      P%NBLADE	! number of blade 
    NELM         = 		P%NELM   ! number of blade element
    CL           =      1
    IBlade       =      1 
    !   the absolute pitch which is the angle between chordline and tangential direction of rotation 
    !Orientation    = AD14_u%TurbineComponents%Blade(1)%Orientation
    BladeBlockageStatus    = 0     ! by default there is no blockage
    
    HubOrientationX = (/U%HUBMOTION%ORIENTATION(1,1,1),U%HUBMOTION%ORIENTATION(1,2,1),U%HUBMOTION%ORIENTATION(1,3,1)/)
    HubPosition     = (/U%HUBMOTION%POSITION(1,1),U%HUBMOTION%POSITION(2,1),U%HUBMOTION%POSITION(3,1)/)
    
    !  firt find the xyz coordinate in inertial system of the lidar beam in the rotor plane, we assume the rotor plane pass the hubcenter node, and the first orientation vector is normal to the rotor plane
    linevector   =       MeasuringPosition_I-LidarPosition_I
    vpt          =      DOT_PRODUCT(linevector,HubOrientationX)    ! the cos of the angle between lidar beam line vector and the rotor plane norm vector
    
    IF (vpt /= 0) THEN  ! the lidar beam is not in parallel with the rotor swept plane
    
       TEMP     = DOT_PRODUCT(HubPosition-LidarPosition_I,HubOrientationX)/vpt   ! temporary vector 
       LidarBeamInRotor_I(1)     = LidarPosition_I(1)+linevector(1)*TEMP   ! the intersection of lidar beam and the rotor plane
       LidarBeamInRotor_I(2)     = LidarPosition_I(2)+linevector(2)*TEMP
       LidarBeamInRotor_I(3)     = LidarPosition_I(3)+linevector(3)*TEMP
       
       
       Do IBlade = 1,NBlade
           Do IElement = 1,NELM-1
        !+		U%BLADEMOTION(1)%POSITION	(1:3,1:19)	REAL(4) 

        
        IF  (BladeBlockageStatus == 0) THEN      ! if Blockage Status is already 1 by previous blade element, then we should not enter the loop
               
                    AeroCenterCoordinate_inner_I   = U%BLADEMOTION(IBlade)%POSITION(:,IElement)+U%BLADEMOTION(IBlade)%TRANSLATIONDISP(:,IElement)! coordinate of the aerodynamic center of a blade element node
                    r_node_inner                   = NORM2(AeroCenterCoordinate_inner_I -HubPosition)  ! assume the aero center is in the swept plane
                    AeroCenterCoordinate_outer_I   = U%BLADEMOTION(IBlade)%POSITION(:,IElement+1)+U%BLADEMOTION(IBlade)%TRANSLATIONDISP(:,IElement+1)      ! coordinate of the aerodynamic center of a blade element node
                    r_node_outer                   = NORM2(AeroCenterCoordinate_outer_I-HubPosition)  ! assume the aero center is in the swept plane
                    r_lidarBeamInRotor_I           = NORM2(LidarBeamInRotor_I-HubPosition)
                    
                    IF ((r_lidarBeamInRotor_I>r_node_inner).AND.(r_lidarBeamInRotor_I<=r_node_outer)) THEN        
                        CL                       = P%BLADEELMCL(IElement) 	 !AD14_P%BLADE%C(IElement)   ! chord length of the node
                    
                        ! two edges for the inner node
                        LeadEdgeCoordinate_inner_I (1)  = AeroCenterCoordinate_inner_I(1)-U%BLADEMOTION(IBlade)%Orientation(2,1,IElement)*CL *0.25
                        LeadEdgeCoordinate_inner_I (2)  = AeroCenterCoordinate_inner_I(2)-U%BLADEMOTION(IBlade)%Orientation(2,2,IElement)*CL *0.25
                        LeadEdgeCoordinate_inner_I (3)  = AeroCenterCoordinate_inner_I(3)-U%BLADEMOTION(IBlade)%Orientation(2,3,IElement)*CL *0.25
                        TrailEdgeCoordinate_inner_I (1) = AeroCenterCoordinate_inner_I(1)+U%BLADEMOTION(IBlade)%Orientation(2,1,IElement)*CL *0.75
                        TrailEdgeCoordinate_inner_I (2) = AeroCenterCoordinate_inner_I(2)+U%BLADEMOTION(IBlade)%Orientation(2,2,IElement)*CL *0.75    
                        TrailEdgeCoordinate_inner_I (3) = AeroCenterCoordinate_inner_I(3)+U%BLADEMOTION(IBlade)%Orientation(2,3,IElement)*CL *0.75    
                    
                        ! two edges for the outer node
                        LeadEdgeCoordinate_outer_I(1)  = AeroCenterCoordinate_outer_I(1)-U%BLADEMOTION(IBlade)%Orientation(2,1,IElement+1)*CL *0.25
                        LeadEdgeCoordinate_outer_I(2)  = AeroCenterCoordinate_outer_I(2)-U%BLADEMOTION(IBlade)%Orientation(2,2,IElement+1)*CL *0.25
                        LeadEdgeCoordinate_outer_I(3)  = AeroCenterCoordinate_outer_I(3)-U%BLADEMOTION(IBlade)%Orientation(2,3,IElement+1)*CL *0.25
                        TrailEdgeCoordinate_outer_I(1) = AeroCenterCoordinate_outer_I(1)+U%BLADEMOTION(IBlade)%Orientation(2,1,IElement+1)*CL *0.75
                        TrailEdgeCoordinate_outer_I(2) = AeroCenterCoordinate_outer_I(2)+U%BLADEMOTION(IBlade)%Orientation(2,2,IElement+1)*CL *0.75    
                        TrailEdgeCoordinate_outer_I(3) = AeroCenterCoordinate_outer_I(3)+U%BLADEMOTION(IBlade)%Orientation(2,3,IElement+1)*CL *0.75  
                    
                        ! formulate first triangle
                        TriangleCoordinate(1,:)        = LeadEdgeCoordinate_inner_I
                        TriangleCoordinate(2,:)        = TrailEdgeCoordinate_inner_I
                        TriangleCoordinate(3,:)        = LeadEdgeCoordinate_outer_I
                    
                         !check weather the first triangle is cross with the beam or not
                        edge1                          = TriangleCoordinate(2,:)-TriangleCoordinate(1,:)         
                        edge2                          = TriangleCoordinate(3,:)-TriangleCoordinate(1,:)         
                        tvec                           = LidarPosition_I-TriangleCoordinate(1,:)   
                        pvec                           = CROSS_PRODUCT(linevector, edge2)
                        det                            = DOT_PRODUCT(edge1,pvec)
                    
                        IF (ABS(det)>0.00001) THEN     ! there is intersection in the triangle plane, then we check whether the point is inside the triangle
                                
                                ub    = DOT_PRODUCT(tvec,pvec)/det;       ! 1st barycentric coordinate
                                qvec  = CROSS_PRODUCT(tvec, edge1);             ! prepare to test V parameter
                                vb    = DOT_PRODUCT(linevector,qvec)/det;       ! 2nd barycentric coordinate
                                tb    = DOT_PRODUCT(edge2,qvec)/det;            ! 'position on the line' coordinate
                                 !test if line/plane intersection is within the triangle
                                IF   ((ub>=-0.0) .AND. (vb>=-0.0) .AND. (ub+vb<=1.0)) THEN 
                                    BladeBlockageStatus = 1 
                                ELSE 
                                END IF  ! the point is inside the triangle the beam is blocked
                                
                    
                        ELSE   ! there is no intersection for this triangle, we will check next triangle
                    
                                !formulate second triangle
                                TriangleCoordinate(1,:)        = TrailEdgeCoordinate_inner_I
                                TriangleCoordinate(2,:)        = LeadEdgeCoordinate_outer_I
                                TriangleCoordinate(3,:)        = TrailEdgeCoordinate_outer_I
                                !check weather the second triangle is cross with the beam or not
                                edge1                          = TriangleCoordinate(2,:)-TriangleCoordinate(1,:)         
                                edge2                          = TriangleCoordinate(3,:)-TriangleCoordinate(1,:)         
                                tvec                           = LidarPosition_I-TriangleCoordinate(1,:)   
                                pvec                           = CROSS_PRODUCT(linevector, edge2)
                                det                            = DOT_PRODUCT(edge1,pvec)
                                
                                IF (ABS(det)>0.00001) THEN     ! there is intersection in the second triangle plane, then we check whether the point is inside the triangle
                               
                                qvec = CROSS_PRODUCT(tvec, edge1);             ! prepare to test V parameter
                                vb    = DOT_PRODUCT(linevector,qvec)/det;       ! 2nd barycentric coordinate
                                tb    = DOT_PRODUCT(edge2,qvec)/det;            ! 'position on the line' coordinate
                                 !test if line/plane intersection is within the triangle
                                   IF   ((ub>=-0.0) .AND. (vb>=-0.0) .AND. (ub+vb<=1.0)) THEN 
                                       BladeBlockageStatus = 1 
                                   ELSE  
                                   END IF  ! the point is inside the triangle the beam is blocked
                                
                                ELSE
                                END IF
                    
                        END IF
                    ELSE
                    END IF ! if the radius is still between inner and outer
                    !
                  ELSE   ! the loop should be existed 
                  END IF  !IF  BlockFlag == 0 
                  
           END DO !IElement = 1,NELM-1
        END DO !IBlade = 1,NBlade
        
    ELSE    ! (vpt == 0)
 
    END IF

    END SUBROUTINE LidarSim_CheckBladeBlockage_15
    
 ! ###################################################   
    
SUBROUTINE LidarSim_InitializeFrozenflow(p, Ifw_p,InputFileData, ErrStat, ErrMsg)
    
    IMPLICIT                        NONE
    CHARACTER(*),                   PARAMETER       ::  RoutineName="LidarSim_InitializeFrozenflow"
    
    TYPE(LidarSim_ParameterType),   INTENT(INOUT)   ::  p           ! parameter data to write results in
    TYPE(InflowWind_ParameterType), INTENT(INOUT)   ::  Ifw_p               !< Inflow Parameters deliver to lidar
    TYPE(LidarSim_InputFile),       INTENT(IN   )   ::  InputFileData       ! Inputdata from the input file
    INTEGER(IntKi),                 INTENT(  OUT)   ::  ErrStat             !< Error status of the operation
    CHARACTER(*),                   INTENT(  OUT)   ::  ErrMsg              !< Error message if ErrStat /= ErrID_None
    
    ! local variables
    INTEGER(IntKi)                                  ::  TmpErrStat          !< Temporary error status
    CHARACTER(ErrMsgLen)                            ::  TmpErrMsg           !< temporary error message
    
    ErrStat =   0
    ErrMsg  =   ''
    
    TmpErrMsg   = ''
    TmpErrStat  = 0
    
    
    P%WINDFILEFORMAT	= 		IFW_P%WINDTYPE
    
    
    IF(IFW_P%WINDTYPE ==  3) THEN	!turbsim format
        
		    P%NLONGI	     =      1	
		    P%PERIODIC	     = 		IFW_P%TSFFWIND%FF%PERIODIC
		    P%TOWERDATAEXIST =  	'FALSE'
		    P%DT	         =      IFW_P%TSFFWIND%FF%FFDTIME
		    !P%FFTOWER	  FG not needed for lidar sim, might be updated later
		    P%FFDTIME	     =      IFW_P%TSFFWIND%FF%FFDTIME
		    P%FFRATE	     =      IFW_P%TSFFWIND%FF%FFRATE
		    P%FFYHWID	     =      IFW_P%TSFFWIND%FF%FFYHWID
		    P%FFZHWID        =      IFW_P%TSFFWIND%FF%FFZHWID
		    P%REFHT	         =		IFW_P%TSFFWIND%FF%REFHT	
		    P%GRIDBASE	     =      IFW_P%TSFFWIND%FF%GRIDBASE
		    P%INITXPOSITION  =      IFW_P%TSFFWIND%FF%INITXPOSITION
		    P%INVFFYD	     =      IFW_P%TSFFWIND%FF%INVFFYD
		    P%INVFFZD	     =      IFW_P%TSFFWIND%FF%INVFFZD
		    P%INVMFFWS	     =      IFW_P%TSFFWIND%FF%INVMFFWS
		    P%MEANFFWS	     =      IFW_P%TSFFWIND%FF%MEANFFWS
		    P%TOTALTIME	     =      IFW_P%TSFFWIND%FF%TOTALTIME
		    P%NFFCOMP	     =      IFW_P%TSFFWIND%FF%NFFCOMP
		    P%NFFSTEPS	     =      IFW_P%TSFFWIND%FF%NFFSTEPS 
		    P%NYGRIDS	     =      IFW_P%TSFFWIND%FF%NYGRIDS
		    P%NZGRIDS	     =      IFW_P%TSFFWIND%FF%NZGRIDS
		    P%NTGRIDS	     =      IFW_P%TSFFWIND%FF%NTGRIDS
		

         

          CALL AllocAry(P%X_UNFROZEN, 1,'Position of frozen planes', TmpErrStat, TmpErrMsg )
             CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
             IF ( ErrStat >= AbortErrLev ) RETURN
          P%X_UNFROZEN = 0
      
      
      
          !-------------------------------------------------------------------------------------------------
          ! Allocate space for the FF array
          !-------------------------------------------------------------------------------------------------

          IF ( .NOT. ALLOCATED( p%FFData ) ) THEN
             CALL AllocAry( p%FFData, 1,p%NZGrids,p%NYGrids,p%NFFComp,P%NFFSTEPS, &
                      'Full-field wind data array.', TmpErrStat, TmpErrMsg )
             CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
             IF ( ErrStat >= AbortErrLev ) RETURN

          ELSE
             IF (SIZE(p%FFData,1) /= p%NZGrids .OR. SIZE(p%FFData,2) /= p%NYGrids .OR. &
                 SIZE(p%FFData,3) /= p%NFFComp .OR. SIZE(p%FFData,3) /= P%NFFSTEPS ) THEN

                   ! Let's make the array the correct size (we should never get here, but you never know)

                DEALLOCATE( p%FFData )

                CALL AllocAry( p%FFData, 1,p%NZGrids,p%NYGrids,p%NFFComp,P%NFFSTEPS, &
                      'Full-field wind data array.', TmpErrStat, TmpErrMsg )
                   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
                   IF ( ErrStat >= AbortErrLev ) RETURN            

             ENDIF !Incorrect size
          ENDIF ! allocated
          p%FFData(1,:,:,:,:)  =IFW_P%TSFFWIND%FF%FFDATA
      
      
      ELSEIF (IFW_P%WINDTYPE ==  4) THEN	!bladed format
            P%NLONGI	     =      1	
		    P%PERIODIC	     = 		IFW_P%BLADEDFFWIND%FF%PERIODIC
		    P%TOWERDATAEXIST =  	'FALSE'
		    P%DT	         =      IFW_P%BLADEDFFWIND%FF%FFDTIME
		    !P%FFTOWER	  FG not needed for lidar sim, might be updated later
		    P%FFDTIME	     =      IFW_P%BLADEDFFWIND%FF%FFDTIME
		    P%FFRATE	     =      IFW_P%BLADEDFFWIND%FF%FFRATE
		    P%FFYHWID	     =      IFW_P%BLADEDFFWIND%FF%FFYHWID
		    P%FFZHWID        =      IFW_P%BLADEDFFWIND%FF%FFZHWID
		    P%REFHT	         =		IFW_P%BLADEDFFWIND%FF%REFHT	
		    P%GRIDBASE	     =      IFW_P%BLADEDFFWIND%FF%GRIDBASE
		    P%INITXPOSITION  =      IFW_P%BLADEDFFWIND%FF%INITXPOSITION
		    P%INVFFYD	     =      IFW_P%BLADEDFFWIND%FF%INVFFYD
		    P%INVFFZD	     =      IFW_P%BLADEDFFWIND%FF%INVFFZD
		    P%INVMFFWS	     =      IFW_P%BLADEDFFWIND%FF%INVMFFWS
		    P%MEANFFWS	     =      IFW_P%BLADEDFFWIND%FF%MEANFFWS
		    P%TOTALTIME	     =      IFW_P%BLADEDFFWIND%FF%TOTALTIME
		    P%NFFCOMP	     =      IFW_P%BLADEDFFWIND%FF%NFFCOMP
		    P%NFFSTEPS	     =      IFW_P%BLADEDFFWIND%FF%NFFSTEPS 
		    P%NYGRIDS	     =      IFW_P%BLADEDFFWIND%FF%NYGRIDS
		    P%NZGRIDS	     =      IFW_P%BLADEDFFWIND%FF%NZGRIDS
		    P%NTGRIDS	     =      IFW_P%BLADEDFFWIND%FF%NTGRIDS
		

         

          CALL AllocAry(P%X_UNFROZEN, 1,'Position of frozen planes', TmpErrStat, TmpErrMsg )
             CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
             IF ( ErrStat >= AbortErrLev ) RETURN
          P%X_UNFROZEN = 0
      
      
      
          !-------------------------------------------------------------------------------------------------
          ! Allocate space for the FF array
          !-------------------------------------------------------------------------------------------------

          IF ( .NOT. ALLOCATED( p%FFData ) ) THEN
             CALL AllocAry( p%FFData, 1,p%NZGrids,p%NYGrids,p%NFFComp,P%NFFSTEPS, &
                      'Full-field wind data array.', TmpErrStat, TmpErrMsg )
             CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
             IF ( ErrStat >= AbortErrLev ) RETURN

          ELSE
             IF (SIZE(p%FFData,1) /= p%NZGrids .OR. SIZE(p%FFData,2) /= p%NYGrids .OR. &
                 SIZE(p%FFData,3) /= p%NFFComp .OR. SIZE(p%FFData,3) /= P%NFFSTEPS ) THEN

                   ! Let's make the array the correct size (we should never get here, but you never know)

                DEALLOCATE( p%FFData )

                CALL AllocAry( p%FFData, 1,p%NZGrids,p%NYGrids,p%NFFComp,P%NFFSTEPS, &
                      'Full-field wind data array.', TmpErrStat, TmpErrMsg )
                   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
                   IF ( ErrStat >= AbortErrLev ) RETURN            

             ENDIF !Incorrect size
          ENDIF ! allocated
          p%FFData(1,:,:,:,:)  =IFW_P%BLADEDFFWIND%FF%FFDATA
      
      
    END IF
    
END SUBROUTINE  LidarSim_InitializeFrozenflow  
    

!----------------------------------------the subroutine for both forzen and evolving turbulence case, get the UVW
SUBROUTINE LidSim_TS_Bladed_FFWind_CalcOutput(Time, PositionXYZ, ParamData, Velocity, ErrStat, ErrMsg)

   IMPLICIT                                                       NONE

   CHARACTER(*),        PARAMETER                              :: RoutineName="LidSim_TS_Bladed_FFWind_CalcOutput"

      ! Passed Variables
   REAL(DbKi),                                  INTENT(IN   )  :: Time              !< time from the start of the simulation
   REAL(ReKi),                                  INTENT(IN   )  :: PositionXYZ(3,1)  !< Array of XYZ coordinates, 3xN
   TYPE(LidarSim_ParameterType),                INTENT(IN   )  :: ParamData         !< Parameters
   REAL(ReKi),                                  INTENT(INOUT)  :: Velocity(3,1)     !< Velocity output at Time    (Set to INOUT so that array does not get deallocated)
   !REAL(ReKi),                                  INTENT(  OUT)  :: DiskVel(3)        !< HACK for AD14: disk velocity output at Time
   !TYPE(IfW_BladedFFWind_MiscVarType),          INTENT(INOUT)  :: MiscVars          !< misc/optimization data (storage for the main data)

      ! Error handling
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< error status
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< The error message

      ! local variables
   INTEGER(IntKi)                                              :: NumPoints         ! Number of points specified by the PositionXYZ array

      ! local counters
   INTEGER(IntKi)                                              :: PointNum          ! a loop counter for the current point

      ! temporary variables
   INTEGER(IntKi)                                              :: TmpErrStat        ! temporary error status
   CHARACTER(ErrMsgLen)                                        :: TmpErrMsg         ! temporary error message


      !-------------------------------------------------------------------------------------------------
      ! Check that the module has been initialized.
      !-------------------------------------------------------------------------------------------------

   ErrStat     = ErrID_None
   ErrMsg      = ''

      !-------------------------------------------------------------------------------------------------
      ! Initialize some things
      !-------------------------------------------------------------------------------------------------


      ! The array is transposed so that the number of points is the second index, x/y/z is the first.
      ! This is just in case we only have a single point, the SIZE command returns the correct number of points.
   NumPoints   =  SIZE(PositionXYZ,2)


      ! Step through all the positions and get the velocities
   DO PointNum = 1, NumPoints

         ! Calculate the velocity for the position
      Velocity(:,PointNum) = FF_Interp_evo(Time,PositionXYZ(:,PointNum),ParamData,TmpErrStat,TmpErrMsg)

         ! Error handling
      IF (TmpErrStat /= ErrID_None) THEN  !  adding this so we don't have to convert numbers to strings every time
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, "LidSim_TS_Bladed_FFWind_CalcOutput [position=("//   &
                                                      TRIM(Num2LStr(PositionXYZ(1,PointNum)))//", "// &
                                                      TRIM(Num2LStr(PositionXYZ(2,PointNum)))//", "// &
                                                      TRIM(Num2LStr(PositionXYZ(3,PointNum)))//")]" )
         IF (ErrStat >= AbortErrLev) RETURN
      END IF

   ENDDO


   RETURN

END SUBROUTINE LidSim_TS_Bladed_FFWind_CalcOutput


!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
   !>    This function is used to interpolate into the full-field wind array or tower array if it has
   !!    been defined and is necessary for the given inputs.  It receives X, Y, Z and
   !!    TIME from the calling routine.  It then computes a time shift due to a nonzero X based upon
   !!    the average windspeed.  The modified time is used to decide which pair of time slices to interpolate
   !!    within and between.  After finding the two time slices, it decides which four grid points bound the
   !!    (Y,Z) pair.  It does a bilinear interpolation for each time slice. Linear interpolation is then used
   !!    to interpolate between time slices.  This routine assumes that X is downwind, Y is to the left when
   !!    looking downwind and Z is up.  It also assumes that no extrapolation will be needed.
   !!
   !!    If tower points are used, it assumes the velocity at the ground is 0.  It interpolates between
   !!    heights and between time slices, but ignores the Y input.
   !!
   !!    11/07/1994 - Created by M. Buhl from the original TURBINT.
   !!    09/25/1997 - Modified by M. Buhl to use f90 constructs and new variable names.  Renamed to FF_Interp.
   !!    09/23/2009 - Modified by B. Jonkman to use arguments instead of modules to determine time and position.
   !!                 Height is now relative to the ground
   !!   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   !!   23-Nov-2020   Feng Guo Flensburg University of Applied Sciences, adjust to consider evolving wind,
   !!                 Remove the 3D linear interpolation for lidar LOS simulation, because this will alter the spectral properity, see: https://wes.copernicus.org/preprints/wes-2021-91/

   FUNCTION FF_Interp_evo(Time, Position, ParamData, ErrStat, ErrMsg)

      IMPLICIT                                              NONE

      CHARACTER(*),           PARAMETER                     :: RoutineName="FF_Interp_evo"

      REAL(DbKi),                            INTENT(IN   )  :: Time           !< time (s)
      REAL(ReKi),                            INTENT(IN   )  :: Position(3)    !< takes the place of XGrnd, YGrnd, ZGrnd
      TYPE(LidarSim_ParameterType),  INTENT(IN   )          :: ParamData      !< Parameters
      REAL(ReKi)                                            :: FF_Interp_evo(3)   !< The U, V, W velocities

      INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat        !< error status
      CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg         !< error message 

         ! Local Variables:

      REAL(ReKi)                                            :: TimeShifted
      REAL(ReKi),PARAMETER                                  :: Tol = 1.0E-3   ! a tolerance for determining if two reals are the same (for extrapolation)
      REAL(ReKi)                                            :: T
      REAL(ReKi)                                            :: TGRID
      REAL(ReKi)                                            :: Y
      REAL(ReKi)                                            :: YGRID
      REAL(ReKi)                                            :: Z
      REAL(ReKi)                                            :: ZGRID
      REAL(ReKi)                                            :: N(8)           ! array for holding scaling factors for the interpolation algorithm
      REAL(ReKi)                                            :: u(8)           ! array for holding the corner values for the interpolation algorithm across a cubic volume
      REAL(ReKi)                                            :: M(4)           ! array for holding scaling factors for the interpolation algorithm
      REAL(ReKi)                                            :: v(4)           ! array for holding the corner values for the interpolation algorithm across an area

      REAL(DbKi)                                            :: TimeToCalculate           !< time (s)
      
      
      INTEGER(IntKi)                                        :: IDIM
      INTEGER(IntKi)                                        :: ITHI         ! index of higher integral in t dimention
      INTEGER(IntKi)                                        :: ITLO         ! index of lower integral in t dimention
      INTEGER(IntKi)                                        :: IYHI         ! index of higher integral in y dimention
      INTEGER(IntKi)                                        :: IYLO         ! index of lower integral in t dimention
      INTEGER(IntKi)                                        :: IZHI         ! index of higher integral in y dimention
      INTEGER(IntKi)                                        :: IZLO         ! index of lower integral in t dimention
      INTEGER(IntKi)                                        :: IXCL           ! index of closest integral in x this determines which unfrozen y-z grid we read
      INTEGER(IntKi)                                        :: IYCL           ! index of closest integral in y dimention
      INTEGER(IntKi)                                        :: IZCL           ! index of closest integral in z dimention
      REAL(ReKi)                                            :: IXIND(PARAMDATA%NLONGI)           ! array of distance
      INTEGER(IntKi)                                        :: ITCL          ! index of closest integral in t dimention
      
      LOGICAL                                               :: OnGrid

      !-------------------------------------------------------------------------------------------------
      ! Initialize variables
      !-------------------------------------------------------------------------------------------------
      
      FF_Interp_evo(:)        = 0.0_ReKi                         ! the output velocities (in case ParamData%NFFComp /= 3)

      ErrStat              = ErrID_None
      ErrMsg               = ""
      
      
      TimeToCalculate      = Time
      !-------------------------------------------------------------------------------------------------
      ! Find the bounding time slices.
      !-------------------------------------------------------------------------------------------------

      
      ! Feng: since the inported upstream wind field is not time shifted (they only have the decorrelation due to evolution), so here we first check which upstream
      ! field has nearest distance in x direction (time) to the point we are asking. 
      ! for example: we have two 3D field:  t = 2.5s, and t = 5s   if lidar is measuring upsream at 40m, and mean wind speed is 20m/s.
      ! the pre-time is 2s, then we chose the first wind field at t = 2.5s. After that, in this t = 2.5 wind field we find the point at x = 40m and y z at specifed value to get the wind speed.
      
      
      ! To simulate lidar, we define several 3D(y,z,t) turbulence at different x positions,
      ! Specially, if we only have 1 3D turbulence, then it will follow the Taylor Frozen theory
      DO IDIM =1,PARAMDATA%NLONGI       ! all the x positions of the 3D wind fields
          ! we will find the index to determine which 3D(y,z,t) turbulence we will use for the current x position
         IXIND   =  ABS(PARAMDATA%X_UNFROZEN - (-Position(1)))    ! not need ParamData%InitXPosition
      END DO
      IXCL = INT(MINLOC(IXIND,1))
      
      
      ! Check if the field is periodic and we have enough time length
      
      IF (TimeToCalculate > PARAMDATA%TOTALTIME) THEN
          
          IF (PARAMDATA%PERIODIC==.True.) THEN
              !write (*,*) 'Turbulence time boundary is violated, the field is read periodically to simulate lidar measurement.'
                  Do while (TimeToCalculate > PARAMDATA%TOTALTIME) 
                    TimeToCalculate = TimeToCalculate- PARAMDATA%TOTALTIME
                  end do
              
            ELSE 
         
            ErrMsg   = ' Turbulence wind array boundaries violated. Turbulence Grid too small in time dimension '// &
                       '(time (t='//TRIM(Num2LStr(Time))//' m) is longer than the turbulence field, the lidar measurement can not be calculated). Extend the field or use periodic field.'
            ErrStat  = ErrID_Fatal
      
            ENDIF
     ENDIF 
      
      ! Perform the time shift.  At time=0, a point half the grid width downstream (ParamData%FFYHWid) will index into the zero time slice.
      ! If we did not do this, any point downstream of the tower at the beginning of the run would index outside of the array.
      ! This all assumes the grid width is at least as large as the rotor.  If it isn't, then the interpolation will not work.

      
      TimeShifted = TimeToCalculate + ( ParamData%InitXPosition - Position(1) )*ParamData%InvMFFWS    ! in distance, X: InputInfo%Position(1) - ParamData%InitXPosition - TIME*ParamData%MeanFFWS


      IF ( ParamData%Periodic ) THEN ! translate TimeShifted to ( 0 <= TimeShifted < ParamData%TotalTime )

         TimeShifted = MODULO( TimeShifted, ParamData%TotalTime )
             ! If TimeShifted is a very small negative number, modulo returns the incorrect value due to internal rounding errors.
             ! See bug report #471
         IF (TimeShifted == ParamData%TotalTime) TimeShifted = 0.0_ReKi

         TGRID = TimeShifted*ParamData%FFRate
         ITLO  = INT( TGRID )             ! convert REAL to INTEGER (add 1 later because our grids start at 1, not 0)
         
         ITCL  = NINT(TGRID) +1             ! the closest t index
         
         
         T     = 2.0_ReKi * ( TGRID - REAL(ITLO, ReKi) ) - 1.0_ReKi     ! a value between -1 and 1 that indicates a relative position between ITLO and ITHI

         ITLO = ITLO + 1
         IF ( ITLO == ParamData%NFFSteps ) THEN
            ITHI = 1
         ELSE
            ITHI = ITLO + 1
         ENDIF


      ELSE

         TGRID = TimeShifted*ParamData%FFRate
         ITLO  = INT( TGRID )             ! convert REAL to INTEGER (add 1 later because our grids start at 1, not 0)
         T     = 2.0_ReKi * ( TGRID - REAL(ITLO, ReKi) ) - 1.0_ReKi     ! a value between -1 and 1 that indicates a relative position between ITLO and ITHI
         
         ITCL  = NINT(TGRID)              ! the closest t index
         
         ITLO = ITLO + 1                  ! add one since our grids start at 1, not 0
         ITHI = ITLO + 1

         IF ( ITLO >= ParamData%NFFSteps .OR. ITLO < 1 ) THEN
            IF ( ITLO == ParamData%NFFSteps  ) THEN
               ITHI = ITLO
               IF ( T <= TOL ) THEN ! we're on the last point
                  T = -1.0_ReKi
               ELSE  ! We'll extrapolate one dt past the last value in the file
                  ITLO = ITHI - 1
               ENDIF
            ELSE
               ErrMsg   = ' Error: FF wind array was exhausted at '//TRIM( Num2LStr( REAL( TIME,   ReKi ) ) )// &
                          ' seconds (trying to access data at '//TRIM( Num2LStr( REAL( TimeShifted, ReKi ) ) )//' seconds).'
               ErrStat  = ErrID_Fatal
               RETURN
            ENDIF
         ENDIF

      ENDIF


      !-------------------------------------------------------------------------------------------------
      ! Find the bounding rows for the Z position. [The lower-left corner is (1,1) when looking upwind.]
      !-------------------------------------------------------------------------------------------------

      ZGRID = ( Position(3) - ParamData%GridBase )*ParamData%InvFFZD

      IF (ZGRID > -1*TOL) THEN
         OnGrid = .TRUE.

            ! Index for start and end slices
         IZLO = INT( ZGRID ) + 1             ! convert REAL to INTEGER, then add one since our grids start at 1, not 0
         IZHI = IZLO + 1
         
         
         IZCL = MAX(NINT(ZGRID)+1,1)     ! the closet Z index
         
         
            ! Set Z as a value between -1 and 1 for the relative location between IZLO and IZHI.
            ! Subtract 1_IntKi from Z since the indices are starting at 1, not 0
         Z = 2.0_ReKi * (ZGRID - REAL(IZLO - 1_IntKi, ReKi)) - 1.0_ReKi

         IF ( IZLO < 1 ) THEN
            IF ( IZLO == 0 .AND. Z >= 1.0-TOL ) THEN
               Z    = -1.0_ReKi
               IZLO = 1
            ELSE
               ErrMsg   = ' FF wind array boundaries violated. Grid too small in Z direction (Z='//&
                          TRIM(Num2LStr(Position(3)))//' m is below the grid).'
               ErrStat  = ErrID_Fatal
               RETURN
            ENDIF
         ELSEIF ( IZLO >= ParamData%NZGrids ) THEN
            IF ( IZLO == ParamData%NZGrids .AND. Z <= TOL ) THEN
               Z    = -1.0_ReKi
               IZHI = IZLO                   ! We're right on the last point, which is still okay
            ELSE
               ErrMsg   = ' FF wind array boundaries violated. Grid too small in Z direction (Z='//&
                          TRIM(Num2LStr(Position(3)))//' m is above the grid).'
               ErrStat  = ErrID_Fatal
               RETURN
            ENDIF
         ENDIF

      ELSE

         OnGrid = .FALSE.  ! this is on the tower    this should cause an error
         
         

         IF ( ParamData%NTGrids < 1 ) THEN
            ErrMsg   = ' Turbulence wind array boundaries violated. Evolution Turbulence Grid too small in Z direction '// &
                       '(height (Z='//TRIM(Num2LStr(Position(3)))//' m) is below the lidar measurement can not be calculated).'
            ErrStat  = ErrID_Fatal
            RETURN
         ENDIF

         IZLO = INT( -1.0*ZGRID ) + 1            ! convert REAL to INTEGER, then add one since our grids start at 1, not 0


         IF ( IZLO >= ParamData%NTGrids ) THEN  !our dz is the difference between the bottom tower point and the ground
            IZLO  = ParamData%NTGrids

               ! Check that this isn't zero.  Value between -1 and 1 corresponding to the relative position.
            Z = 1.0_ReKi - 2.0_ReKi * (Position(3) / (ParamData%GridBase - REAL(IZLO - 1_IntKi, ReKi)/ParamData%InvFFZD))

         ELSE

               ! Set Z as a value between -1 and 1 for the relative location between IZLO and IZHI.  Used in the interpolation.
            Z = 2.0_ReKi * (ABS(ZGRID) - REAL(IZLO - 1_IntKi, ReKi)) - 1.0_ReKi

         ENDIF
         IZHI = IZLO + 1

      ENDIF


      IF ( OnGrid ) THEN      ! The tower points don't use this

         !-------------------------------------------------------------------------------------------------
         ! Find the bounding columns for the Y position. [The lower-left corner is (1,1) when looking upwind.]
         !-------------------------------------------------------------------------------------------------

            YGRID = ( Position(2) + ParamData%FFYHWid )*ParamData%InvFFYD    ! really, it's (Position(2) - -1.0*ParamData%FFYHWid)

            IYLO = INT( YGRID ) + 1             ! convert REAL to INTEGER, then add one since our grids start at 1, not 0
            IYHI = IYLO + 1
            
            IYCL = MAX(NINT(YGRID)+1,1)     ! the closet Y index
            
            
               ! Set Y as a value between -1 and 1 for the relative location between IYLO and IYHI.  Used in the interpolation.
               ! Subtract 1_IntKi from IYLO since grids start at index 1, not 0
            Y = 2.0_ReKi * (YGRID - REAL(IYLO - 1_IntKi, ReKi)) - 1.0_ReKi

            IF ( IYLO >= ParamData%NYGrids .OR. IYLO < 1 ) THEN
               IF ( IYLO == 0 .AND. Y >= 1.0-TOL ) THEN
                  Y    = -1.0_ReKi
                  IYLO = 1
               ELSE IF ( IYLO == ParamData%NYGrids .AND. Y <= TOL ) THEN
                  Y    = -1.0_ReKi
                  IYHI = IYLO                   ! We're right on the last point, which is still okay
               ELSE
                  ErrMsg   = ' FF wind array boundaries violated: Grid too small in Y direction. Y='// &
                             TRIM(Num2LStr(Position(2)))//'; Y boundaries = ['//TRIM(Num2LStr(-1.0*ParamData%FFYHWid))// &
                             ', '//TRIM(Num2LStr(ParamData%FFYHWid))//']'
                  ErrStat = ErrID_Fatal         ! we don't return anything
                  RETURN
               ENDIF
            ENDIF

         !-------------------------------------------------------------------------------------------------
         ! Interpolate on the grid    the code below uses 3d linear interpolation
         !-------------------------------------------------------------------------------------------------

         DO IDIM=1,ParamData%NFFComp       ! all the components

            IF (ParamData%NearestInterpFlag==.True.) THEN

            !the code below uses nearest point, this is used to avoid filtering effect due to interpolation
             ITCL = MIN(ITCL,PARAMDATA%NFFSTEPS)
            FF_Interp_evo(IDIM)  = ParamData%FFData(IXCL, IZCL, IYCL, IDIM, ITCL )
            
            ELSE
              !New Algorithm here,   the code below uses 3d linear interpolation
                N(1)  = ( 1.0_ReKi + Z )*( 1.0_ReKi - Y )*( 1.0_ReKi - T )
                N(2)  = ( 1.0_ReKi + Z )*( 1.0_ReKi + Y )*( 1.0_ReKi - T )
                N(3)  = ( 1.0_ReKi - Z )*( 1.0_ReKi + Y )*( 1.0_ReKi - T )
                N(4)  = ( 1.0_ReKi - Z )*( 1.0_ReKi - Y )*( 1.0_ReKi - T )
                N(5)  = ( 1.0_ReKi + Z )*( 1.0_ReKi - Y )*( 1.0_ReKi + T )
                N(6)  = ( 1.0_ReKi + Z )*( 1.0_ReKi + Y )*( 1.0_ReKi + T )
                N(7)  = ( 1.0_ReKi - Z )*( 1.0_ReKi + Y )*( 1.0_ReKi + T )
                N(8)  = ( 1.0_ReKi - Z )*( 1.0_ReKi - Y )*( 1.0_ReKi + T )
                N     = N / REAL( SIZE(N), ReKi )  ! normalize
            
            
                u(1)  = ParamData%FFData(IXCL, IZHI, IYLO, IDIM, ITLO )
                u(2)  = ParamData%FFData(IXCL, IZHI, IYHI, IDIM, ITLO )
                u(3)  = ParamData%FFData(IXCL, IZLO, IYHI, IDIM, ITLO )
                u(4)  = ParamData%FFData(IXCL, IZLO, IYLO, IDIM, ITLO )
                u(5)  = ParamData%FFData(IXCL, IZHI, IYLO, IDIM, ITHI )
                u(6)  = ParamData%FFData(IXCL, IZHI, IYHI, IDIM, ITHI )
                u(7)  = ParamData%FFData(IXCL, IZLO, IYHI, IDIM, ITHI )
                u(8)  = ParamData%FFData(IXCL, IZLO, IYLO, IDIM, ITHI )
            
                FF_Interp_evo(IDIM)  =  SUM ( N * u  )
             END IF

            
         END DO !IDIM

      ELSE

      !-------------------------------------------------------------------------------------------------
      ! Interpolate on the tower array   should not be reached
      !-------------------------------------------------------------------------------------------------
         
         !DO IDIM=1,ParamData%NFFComp    ! all the components
         !
         !   !----------------------------------------------------------------------------------------------
         !   ! Interpolate between the two times using an area interpolation.
         !   !----------------------------------------------------------------------------------------------
         !
         !      ! Setup the scaling factors.  Set the unused portion of the array to zero
         !   M(1)  =  ( 1.0_ReKi + Z )*( 1.0_ReKi - T )
         !   M(2)  =  ( 1.0_ReKi + Z )*( 1.0_ReKi + T )
         !   M(3)  =  ( 1.0_ReKi - Z )*( 1.0_ReKi - T )
         !   M(4)  =  ( 1.0_ReKi - Z )*( 1.0_ReKi + T )
         !   M     =  M / 4.0_ReKi               ! normalize
         !
         !   IF (IZHI > ParamData%NTGrids) THEN
         !      v(1)  =  0.0_ReKi  ! on the ground
         !      v(2)  =  0.0_ReKi  ! on the ground
         !   ELSE
         !      v(1)  =  ParamData%FFTower( IDIM, IZHI, ITLO )
         !      v(2)  =  ParamData%FFTower( IDIM, IZHI, ITHI )
         !   END IF
         !   
         !   v(3)  =  ParamData%FFTower( IDIM, IZLO, ITLO )
         !   v(4)  =  ParamData%FFTower( IDIM, IZLO, ITHI )
         !   
         !   FF_Interp(IDIM)  =  SUM ( M * v ) 
         
         
         !END DO !IDIM

      ENDIF ! OnGrid
      RETURN

   END FUNCTION FF_Interp_evo

!############################################    
    END MODULE LidarSim_Subs