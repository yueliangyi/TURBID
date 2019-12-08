/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include <iomanip>
#include <omp.h>
#include "system.h"

namespace yueliangyi {
namespace pseudospectral {
namespace leveltwo {
  
  
  using namespace levelbase;  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
System::System(int& argc, char**& argv) :
  size_(3),mpigrid_(2),
  mpi_(argc,argv),dcp_(mpi_,mpigrid_,size_),
  filelog_(mpi_,dcp_),
  fileini_(mpi_,dcp_){
  PS_DEBUG_TRACE_ENTER("System::System()")
    
    
    IO_FILE_CONFIGURATION = "configuration.ini";
    IO_FILE_LOG_COMPUTATION = "logcomputation";
    if (argc>1) {
      IO_FILE_CONFIGURATION = string(argv[1]);
    }
    if (argc>2) {
      IO_FILE_LOG_COMPUTATION = string(argv[2]);
    }

    
    
    fileini_.SetFileName(IO_FILE_CONFIGURATION);
    fileini_.SetSensitive(0);
    fileini_.Load();
    

    printlog_ = fileini_.GetValueBool(SECTION_SYSTEM,KEY_SS_PRINT_LOG_WHEN_RUNNING);
    writelog_ = fileini_.GetValueBool(SECTION_SYSTEM,KEY_SS_WRITE_LOG_WHEN_RUNNING);
    phasenum_ = fileini_.GetValueInt(SECTION_SYSTEM,KEY_SS_TOTAL_PHASE_NUMBER);
        
    
    ivec mpigrid = fileini_.GetVectorInt(SECTION_SYSTEM,KEY_SD_MPI_GRID_POINT_NUMBER);
    mpigrid_(0) = mpigrid(0);
    mpigrid_(1) = mpigrid(1);
    if (mpi_.GetSize() != mpigrid_(0)*mpigrid_(1)) {
      mpigrid_(1) = mpi_.GetSize()/mpigrid_(0);
    }
    
    int threads = fileini_.GetValueInt(SECTION_SYSTEM,KEY_SD_THREAD_NUMBER_OPENMP);
    omp_set_num_threads(threads);
    
    size_ = fileini_.GetValueString(SECTION_DOMAIN,KEY_SD_GRID_POINT_NUMBER);
    
    dcp_.Configure();
    
    
  PS_DEBUG_TRACE_LEAVE("System::System()")
}
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void System::Print(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER("System::Print()")
  
  
  time_t time_head_;
  struct tm* tm_head_ = NULL;
  
  time(&time_head_);
  tm_head_ = localtime(&time_head_);
  
  
  stream
  << "Task Start Time "
  << tm_head_->tm_mon+1     << '/'
  << tm_head_->tm_mday      << '/'
  << tm_head_->tm_year+1900 << '-'
  << tm_head_->tm_hour      << ':'
  << tm_head_->tm_min       << ':'
  << tm_head_->tm_sec
  << std::endl
  << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << ""
  << "Total Phase Number "
  << phasenum_;
  
  
  
  PS_DEBUG_TRACE_LEAVE("System::Print()")
}
  
  
  
  
  
  
  
} // end namespace leveltwo
} // end namespace pseudospectral
} // end namespace yueliangyi
