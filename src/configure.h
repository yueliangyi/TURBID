/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVBASE_CONFIGURE_H_
#define YLY_PS_LVBASE_CONFIGURE_H_



//
// Configuation for the basic Armadillo matrix library.
// For better performance, do NOT comment out the last three lines.
// The definition of ARMA_NO_DEBUG is for debugging.
// Comment out for a release version.
//
//#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
#define ARMA_USE_ARPACK



//#define PS_DEBUG_TRACE
#define PS_DEBUG_CHECK



//
// Directory and name for system file.
// Note that CONFIGURATION_FILE_NAME has a default suffix ".ini".
// Others should include the suffix with name if needed.
#define IO_DEFAULT_DIRECTORY              "./"
#define IO_DEFAULT_FILE_NAME              "default"
#define IO_FILE_CLASS_BOUNDARY            "boundary.dat"
#define IO_FILE_CLASS_VARIABLE            "variable.dat"
#define IO_FILE_CLASS_DOMAIN              "domain.dat"
#define IO_FILE_CLASS_STEPPER_TIME        "steppertime.dat"
#define IO_FILE_CLASS_STEPPER_SPACE       "stepperspace.dat"
#define IO_FILE_CLASS_PHASE               "phase.dat"
#define IO_FILE_LOG_PREPROCESSING         "logpreprocessing"
#define IO_FILE_LOG_POSTPROCESSING        "logpostprocessing"


#define IO_FILE_SEPERATOR_LINE string(120,'-')
#define IO_FILE_SEPERATOR_STAR string(120,'*')
#define IO_FILE_HEAD_COLUMN_LENGTH        30
#define IO_FILE_LOG_WRITE_INTERVAL        10


//------------------------------------------------------------------------------
#define SECTION_SYSTEM                      "System"
#define KEY_SS_PRINT_LOG_WHEN_RUNNING       "print_log_when_running"
#define KEY_SS_WRITE_LOG_WHEN_RUNNING       "write_log_when_running"
#define KEY_SS_TOTAL_PHASE_NUMBER           "total_phase_number"
#define KEY_SD_MPI_GRID_POINT_NUMBER        "mpi_grid_point_number"
#define KEY_SD_THREAD_NUMBER_OPENMP         "thread_number_openmp"
//------------------------------------------------------------------------------
#define SECTION_PRE_PROCESSOR               "PreProcessor"
//------------------------------------------------------------------------------
#define SECTION_POST_PROCESSOR              "PostProcessor"
#define KEY_SPOST_WORK_INDEX_RANGE          "work_index_range"
//------------------------------------------------------------------------------
#define SECTION_DOMAIN                      "Domain"
#define KEY_SD_APPLIED_DOMAIN_TYPE          "applied_domain_type"
#define KEY_SD_DIMENSIONLESS_LENGTH         "dimensionless_length"
#define KEY_SD_GRID_POINT_NUMBER            "grid_point_number"
#define KEY_SD_TOP_BOTTOM_SLOPE             "top_bottom_slope"
#define KEY_SD_TOP_BOTTOM_SLOPE_DIRECTION   "slope_direction"
#define KEY_SD_READ_PROFILE_FROM_FILE       "read_profile_from_file"
#define KEY_SD_RIPPLE_BED_AMPLITUDE         "ripple_bed_amplitude"
//------------------------------------------------------------------------------
#define SECTION_STEPPER_SPACE               "Space Stepper"
#define KEY_SSS_APPLIED_STEPPER_TYPE        "applied_stepper_type"
#define KEY_SSS_MAPPING_FUNCTION            "mapping_function"
#define KEY_SSS_MAPPING_COEFFICIENT         "mapping_coefficient"
//------------------------------------------------------------------------------
#define SECTION_STEPPER_TIME                "Time Stepper"
#define KEY_SST_APPLIED_STEPPER_TYPE        "applied_stepper_type"
#define KEY_SST_MAXIMUM_CFL_NUMBER          "maximum_cfl_number"
#define KEY_SST_OUTPUT_BASED_ON_STEP        "output_based_on_step"

#define KEY_SST_MAXIMUM_TASK_STEP           "maximum_task_step"
#define KEY_SST_RESULT_WRITE_INTERVAL       "result_write_interval"
#define KEY_SST_MARCHING_TIME_STEP          "marching_time_step"
#define KEY_SST_TASK_START_STEP             "task_start_step"

#define KEY_SST_MAX_MARCHING_TIME           "max_marching_time"
#define KEY_SST_OUTPUT_INTERVAL             "output_interval"
#define KEY_SST_INITIAL_TIME_STEP           "initial_time_step"
#define KEY_SST_TASK_START_TIME_INDEX       "task_start_time_index"
#define KEY_SST_CFL_ADJUST_SCALE            "cfl_adjust_scale"

#define KEY_SST_RESULT_WRITE_INITIAL        "result_write_initial"
//------------------------------------------------------------------------------
#define SECTION_PHASE_CARRIER               "Phase Carrier"
#define KEY_SPC_APPLIED_FORCE_TYPE          "applied_force_type"
#define KEY_SPC_FORCE_DIRECTION             "force_direction"
#define KEY_SPC_FORCE_AMPLITUDE             "force_amplitude"
#define KEY_SPC_REYNOLDS_NUMBER             "reynolds_number"
#define KEY_SPC_FROUDE_NUMBER               "froude_number"
#define KEY_SPC_AVERAGED_CONCENTRATION      "averaged_concentration"
#define KEY_SPC_VERTICAL_TYPE_U1            "vertical_type_u1"
#define KEY_SPC_VERTICAL_TYPE_U2            "vertical_type_u2"
#define KEY_SPC_VERTICAL_TYPE_U3            "vertical_type_u3"
#define KEY_SPC_VERTICAL_TYPE_P             "vertical_type_p"
#define KEY_SPC_SPECIFIC_GRAVITY            "specific_gravity"
#define KEY_SPC_ROSSBY_NUMBER               "rossby_number"
#define KEY_SPC_AFFECTED_BY_SCALAR_PHASE    "affected_by_scalar_phase"
#define KEY_SPC_DAMPING_COEFFICIENT         "damping_coefficient"
#define KEY_SPC_ITERATION_COEFFICIENT       "iteration_coefficient"
//------------------------------------------------------------------------------
#define SECTION_PHASE_SCALAR                "Phase Scalar"
#define KEY_SPS_SCHMIDT_NUMBER              "schmidt_number"
#define KEY_SPS_SETTLING_VELOCITY           "settling_velocity"
#define KEY_SPS_VERTICAL_BOUNDARY_TYPE      "vertical_boundary_type"
#define KEY_SPS_CRITICAL_SHEAR_STRESS       "critical_shear_stress"
#define KEY_SPS_EROSION_RATE                "erosion_rate"
//------------------------------------------------------------------------------






#endif // YLY_PS_LVBASE_CONFIGURE_H_
