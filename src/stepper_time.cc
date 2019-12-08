/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include "stepper_time.h"
#include <iomanip>

namespace yueliangyi {
namespace pseudospectral {
namespace levelone {
  
  
using namespace levelbase;
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool TimeStepper::Configuration::Load(IniFile& fileini) {
  PS_DEBUG_TRACE_ENTER("TimeStepper::Configuration::Load(IniFile)")
  
  maxcfl = fileini.GetValueDouble(SECTION_STEPPER_TIME,KEY_SST_MAXIMUM_CFL_NUMBER);
  output_based_on_step = fileini.GetValueInt(SECTION_STEPPER_TIME,KEY_SST_OUTPUT_BASED_ON_STEP);
  
  maxstep = fileini.GetValueInt(SECTION_STEPPER_TIME,KEY_SST_MAXIMUM_TASK_STEP);
  writeinterval = fileini.GetValueInt(SECTION_STEPPER_TIME,KEY_SST_RESULT_WRITE_INTERVAL);
  timestep = fileini.GetValueDouble(SECTION_STEPPER_TIME,KEY_SST_MARCHING_TIME_STEP);
  startstep = fileini.GetValueInt(SECTION_STEPPER_TIME,KEY_SST_TASK_START_STEP);
  
  max_marching_time = fileini.GetValueDouble(SECTION_STEPPER_TIME,KEY_SST_MAX_MARCHING_TIME);
  output_interval = fileini.GetValueDouble(SECTION_STEPPER_TIME,KEY_SST_OUTPUT_INTERVAL);
  initial_time_step = fileini.GetValueDouble(SECTION_STEPPER_TIME,KEY_SST_INITIAL_TIME_STEP);
  start_time_index = fileini.GetValueInt(SECTION_STEPPER_TIME,KEY_SST_TASK_START_TIME_INDEX);
  cfl_adjust_scale  = fileini.GetVectorDouble(SECTION_STEPPER_TIME,KEY_SST_CFL_ADJUST_SCALE);
  
  writeinitial = fileini.GetValueBool(SECTION_STEPPER_TIME,KEY_SST_RESULT_WRITE_INITIAL);

  PS_DEBUG_TRACE_LEAVE("TimeStepper::Configuration::Load(IniFile)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
TimeStepper::TimeStepper(const mpi_t& mpi, const dcp_t& dcp, const Base::cfg_t* cfg) :
  Base(mpi,dcp,cfg) {
  PS_DEBUG_TRACE_ENTER(name_+"TimeStepper(mpi_t,dcp_t,cfg_t)")
  
    
#ifdef PS_DEBUG_CHECK
    if (cfg_ == NULL) {
      name_ += "TimeStepper(mpi_t,dcp_t,cfg_t) > Can not Initialize TimeStepper Without Configuration!";
      throw std::runtime_error(name_);
    }
#endif
    
    // Set the name of output file.
    filename_ = string("phase_00_")+IO_FILE_CLASS_STEPPER_TIME;
    
    
  PS_DEBUG_TRACE_LEAVE(name_+"TimeStepper(mpi_t,dcp_t,cfg_t)")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool TimeStepper::Update(const double maxulen) {
  PS_DEBUG_TRACE_ENTER(name_+"Update(void)")
  
  // If the computation flag is flase, do nothing.
  if (!isforward_) { return true; }
  
  
  cfg_t* pcfg = (cfg_t*)cfg_;
  
  
  // Update time and step interval.
  step_ += 1;
//  time_  = step_*dt_;
  time_  = time_+dt_;
  
  // Update the computation and output flag.
  if (step_ < pcfg->maxstep) {
    isforward_ = true;
  } else {
    isforward_ = false;
  }
  if (step_%pcfg->writeinterval == 0) {
    isoutput_ = true;
  } else {
    isoutput_ = false;
  }
  
  // Compute CFL number in a rigorous way.
  // Note that the computed CFL number belongs to previous stage, usually.
//  cfl_ = (maxu/mindx+maxv/mindy+maxw/mindz)*dt_;
  cfl_ = maxulen*dt_;
  
  // If the CFL number does not meet the condition requirement, turn off the
  // computation flag while turn on the output flag for record.
  if (cfl_ > pcfg->maxcfl) {
    isforward_ = false;
    isoutput_ = true;
  }
  
  PS_DEBUG_TRACE_LEAVE(name_+"Update(void)")
  return true;
}
  
  
  
//
//
//
bool TimeStepper::Update2(const double maxulen) {
  PS_DEBUG_TRACE_ENTER(name_+"Update(void)")
  
  // If the computation flag is flase, do nothing.
  if (!isforward_) { return true; }
  
  
  cfg_t* pcfg = (cfg_t*)cfg_;
  
  
  // Update time and step interval.
  step_ += 1;
  time_ = time_+dt_;
  
  // Update the computation and output flag.
  if (time_ < pcfg->max_marching_time) {
    isforward_ = true;
  } else {
    isforward_ = false;
  }
//  if (abs(time_/pcfg->output_interval-output_time_index_-1) < 0.00001/100) {
//    isoutput_ = true;
//    output_time_index_ = output_time_index_+1;
//  } else {
//    isoutput_ = false;
//  }
  
  // Compute CFL number in a rigorous way.
  // Note that the computed CFL number belongs to previous stage, usually.
  cfl_ = maxulen*dt_;
  
  
  // adjust dt
  // dt will be changed every time step (not every stage)!
  // since maxulen will be updated every time step
  if (mpi_.isMain() && cfl_>0) {
//    dt_ = 0.95*pcfg->maxcfl/maxulen;
//    const double adjscaletop = 0.90;
//    const double adjscalebot = 0.85;
    const double adjscaletop = pcfg->cfl_adjust_scale(0);
    const double adjscalebot = pcfg->cfl_adjust_scale(1);
    if (cfl_>adjscaletop*pcfg->maxcfl) {
      dt_ = dt_*adjscaletop;
      while (maxulen*dt_>adjscaletop*pcfg->maxcfl) {
        dt_ = dt_*adjscaletop;
      }
    } else if (cfl_<adjscalebot*pcfg->maxcfl) {
//      dt_ = dt_*(2-adjscaletop);
      dt_ = dt_*1.05;
//      while (maxulen*dt_<adjscalebot*pcfg->maxcfl) {
//        dt_ = dt_*1.05;
//      }
    }
    
    double dttmp = 0;
    if ((time_+dt_)/pcfg->output_interval-output_time_index_-1 >=0) {
      dttmp = (output_time_index_+1)*pcfg->output_interval-time_;
      if (dttmp < pcfg->initial_time_step/100) {
        // keep the adjusted dt!
        dt_ = dt_/20;
      } else {
        dt_ = dttmp;
      }
      output_time_index_ = output_time_index_+1;
      isoutput_ = true;
    } else {
      isoutput_ = false;
    }
    
//    std::cout << dt_ << std::endl;
  }
  mpi_.BcastDouble(time_);
  mpi_.BcastDouble(dt_);
  mpi_.BcastInt(output_time_index_);
  mpi_.BcastBool(isoutput_);

  
  
  // If the CFL number does not meet the condition requirement, turn off the
  // computation flag while turn on the output flag for record.
  const double relaxscalecfl = 1.1;
  if (cfl_>pcfg->maxcfl*relaxscalecfl || dt_<pcfg->initial_time_step/1000) {
    isforward_ = false;
    isoutput_ = true;
  }
  
  PS_DEBUG_TRACE_LEAVE(name_+"Update(void)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
bool TimeStepper::Configure(void) {
  PS_DEBUG_TRACE_ENTER(name_+"Configure(void)")
  
  cfg_t* pcfg = (cfg_t*)cfg_;
  
  
  if (pcfg->output_based_on_step) {
  // Initialize step counter with default value 0.
  step_ = std::max(0,pcfg->startstep);
  // Load time and step interval.
  dt_   = pcfg->timestep;
  time_ = step_*dt_;
  // Update flags at the initialization stage.
  if (step_ <= pcfg->maxstep) { isforward_ = true; }
  } else {
    step_ = 0;
    dt_   = pcfg->initial_time_step;
    if (pcfg->start_time_index > 0) {
      time_ = pcfg->start_time_index*pcfg->output_interval;
    } else {
      time_ = 0;
    }
    output_time_index_ = pcfg->start_time_index;
    if (time_ <= pcfg->max_marching_time) { isforward_ = true; }
  }
  isoutput_ = pcfg->writeinitial;
  
  isconfigured_ = true;
  
  PS_DEBUG_TRACE_LEAVE(name_+"Configure(void)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void TimeStepper::Print(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER(name_+"Print(ostringstream)")
  
  cfg_t* pcfg = (cfg_t*)cfg_;
  
  stream
  << "Maximum Task Step "
  << pcfg->maxstep
  << std::endl
  << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << ""
  << "Output by Step "
  << pcfg->writeinterval
  << std::endl
  << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << ""
  << "Fixed Time Step "
  << pcfg->timestep
  << std::endl
  << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << ""
  << "Start at "
  << pcfg->startstep
  << std::endl
  << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << ""
  << "Maximum CFL Number "
  << pcfg->maxcfl
  << std::endl
  << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << ""
  << "Write Initial "
  << pcfg->writeinitial;
  
  
  PS_DEBUG_TRACE_LEAVE(name_+"Print(ostringstream)")
}
  
  
  
  
  
  
  
} // end namespace levelone
} // end namespace pseudospectral
} // end namespace yueliangyi
