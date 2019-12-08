/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/

#include <stdexcept>
#include <iomanip>
#include <omp.h>
#include "base.h"
#include "p3dfft.h"

namespace yueliangyi {
namespace pseudospectral {
  
  
  
// Define the global output message container here.
// Not a very good design for debugging output.
// But it makes the collection of information easy.
ostringstream _message_;
  
string IO_FILE_CONFIGURATION;
string IO_FILE_LOG_COMPUTATION;
  
  
  
namespace levelbase {
  
  
  
//------------------------------------------------------------------------------
// Function : Constructor
//            Initialize the parallelization scheme for system. The system is
//            based on MPI, meaning that there is at least one MPI thread.
//            Always trying to implement OpenMP. Set the number of threads
//            equals to one if there is no require of OpenMP at all. It is
//            recommended to check the MPI and OpenMP settings before running.
//------------------------------------------------------------------------------
// Parameter: IN  - argc      > argument count
//            IN  - argv      > argument value
//------------------------------------------------------------------------------
MPIGroup::MPIGroup(int& argc, char**& argv) {
  PS_DEBUG_TRACE_ENTER(name_+"MPIGroup(int,char**)")
  
//  int hostlength = 0;
//  char* host = new char(128);
  int errflag = MPI_SUCCESS;
  
  errflag += MPI_Init(&argc,&argv);
  errflag += MPI_Comm_size(MPI_COMM_WORLD,&size_);
  errflag += MPI_Comm_rank(MPI_COMM_WORLD,&rank_);
//  errflag += MPI_Get_processor_name(host,&hostlength);
//  host_ = string(host);
  
  if (errflag != MPI_SUCCESS) {
    name_ += "MPIGroup(int,char**) > Failed to Initialize MPI Environment!";
    throw std::runtime_error(name_);
  }
  
  // Be careful for the usage of OpenMP.
  // For wrong setting of OMP_NUM_THREADS, the behavior could be bad.
  threads_ = omp_get_max_threads();
  
//  delete [] host;
  PS_DEBUG_TRACE_LEAVE(name_+"MPIGroup(int,char**)")
}
  
  
  
//------------------------------------------------------------------------------
// Function : Broadcast string
//------------------------------------------------------------------------------
// Parameter: None
//------------------------------------------------------------------------------
void MPIGroup::BcastString(string& target, const int root) const {
  PS_DEBUG_TRACE_ENTER(name_+"BcastString(string,int)")
  
  int errflag = MPI_SUCCESS;
  int strlength = target.length();
  
  errflag += MPI_Bcast(&strlength,1,MPI_UNSIGNED,root,MPI_COMM_WORLD);
  if (strlength == 0) { target = ""; return; }
  if (GetRank() != root) { target.resize(strlength); }
  
  errflag += MPI_Bcast(&target.at(0),
                       strlength,
                       MPI_CHAR,
                       root,
                       MPI_COMM_WORLD);
  
  if (errflag != MPI_SUCCESS) {
    throw std::runtime_error("MPIGroup::BcastString(string,int) > Failed to Broadcast!");
  }
  
  PS_DEBUG_TRACE_LEAVE(name_+"BcastString(string,int)")
}
  
  
  
//------------------------------------------------------------------------------
// Function : Broadcast int
//------------------------------------------------------------------------------
// Parameter: None
//------------------------------------------------------------------------------
void MPIGroup::BcastInt(int& target, const int root) const {
  PS_DEBUG_TRACE_ENTER(name_+"BcastInt(int,int)")
  
  int errflag = MPI_SUCCESS;
  errflag += MPI_Bcast(&target,1,MPI_INT,root,MPI_COMM_WORLD);
  
  if (errflag != MPI_SUCCESS) {
    throw std::runtime_error("MPIGroup::BcastInt(int,int) > Failed to Broadcast!");
  }
  
  PS_DEBUG_TRACE_LEAVE(name_+"BcastInt(int,int)")
}
  
  
  
//------------------------------------------------------------------------------
// Function : Broadcast bool
//------------------------------------------------------------------------------
// Parameter: None
//------------------------------------------------------------------------------
void MPIGroup::BcastBool(bool& target, const int root) const {
  PS_DEBUG_TRACE_ENTER(name_+"BcastBool(bool,int)")
  
  int errflag = MPI_SUCCESS;
  errflag += MPI_Bcast(&target,1,MPI_C_BOOL,root,MPI_COMM_WORLD);
  
  if (errflag != MPI_SUCCESS) {
    throw std::runtime_error("MPIGroup::BcastBool(bool,int) > Failed to Broadcast!");
  }
  
  PS_DEBUG_TRACE_LEAVE(name_+"BcastBool(bool,int)")
}
  
  
  
//------------------------------------------------------------------------------
// Function : Broadcast double
//------------------------------------------------------------------------------
// Parameter: None
//------------------------------------------------------------------------------
void MPIGroup::BcastDouble(double& target, const int root) const {
  PS_DEBUG_TRACE_ENTER(name_+"BcastDouble(double,int)")
  
  int errflag = MPI_SUCCESS;
  errflag += MPI_Bcast(&target,1,MPI_DOUBLE,root,MPI_COMM_WORLD);
  
  if (errflag != MPI_SUCCESS) {
    throw std::runtime_error("MPIGroup::BcastDouble(double,int) > Failed to Broadcast!");
  }
  
  PS_DEBUG_TRACE_LEAVE(name_+"BcastDouble(double,int)")
}
  
  
  
//------------------------------------------------------------------------------
// Function : Print current status
//------------------------------------------------------------------------------
// Parameter: IN  - stream    > target stream as container
//------------------------------------------------------------------------------
void MPIGroup::Print(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER(name_+"Print(ostringstream)")
  
  stream
  << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << "MPIGroup > "
  << "Current Rank " << GetRank() << std::endl
  << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << ""
  << "Name of Holder \"" << GetHost() << "\"" << std::endl
  << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << ""
  << "Total MPI Number " << GetSize() << std::endl
  << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << ""
  << "Thread Per MPI " << GetThreads() << std::endl;
  
  PS_DEBUG_TRACE_LEAVE(name_+"Print(ostringstream)")
}
  
  
  
  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
DomainDecomposition::~DomainDecomposition(void) {
  Cp3dfft_clean();
}
  
  
  
  
//------------------------------------------------------------------------------
// Function : Configure the data distribution with specific MPI grid and data
//            size. Note that, the library P3DFFT should be compiled with flag
//            "--enable-stride1" to make data arranged in z-direction first.
//            What's more, the data for wavenumber domain is complex.
//------------------------------------------------------------------------------
// Parameter: None
//------------------------------------------------------------------------------
bool DomainDecomposition::Configure(void) {
  PS_DEBUG_TRACE_ENTER(name_+"Configure(void)")
  
#ifdef PS_DEBUG_CHECK
  if (mpigrid_.n_elem < 2) {
    name_ += "Configure(void) > MPI Grid Size Should Have at Least Two Dimensions!";
    throw std::runtime_error(name_);
  }
  if (size_.n_elem < 3) {
    name_ += "Configure(void) > Data Size Should Have no Less Than Three Elements!";
    throw std::runtime_error(name_);
  }
  if (size_(0)==0 || size_(1)==0 || size_(2)==0) {
    name_ += "Configure(void) > Empty Data Dimension is not Supported!";
    throw std::runtime_error(name_);
  }
  if (size_(0)%2!=0 && size_(1)%2!=0 && size_(2)%2==0) {
    name_ += "Configure(void) > Data Size Should be Even on Plane and Odd in Z!";
    throw std::runtime_error(name_);
  }
#endif
  
  // Initialize the P3DFFT library before using it.
  // Note that there is no pruned dimensions at all.
  ivec mpigrid = mpigrid_;
  ivec sizememory(3);
  Cp3dfft_setup(mpigrid.memptr(),
                size_(0),
                size_(1),
                size_(2),
                MPI_Comm_c2f(MPI_COMM_WORLD),
                size_(0)/3*2,
                size_(1)/3*2,
                size_(2),
                1,
                sizememory.memptr());
  
  
  const int sizedimension = 50;
  ivec dimensionlocal(sizedimension,arma::fill::zeros);
  
  // 1. Fill in the grid of MPI.
  dimensionlocal(10*0+0) = mpigrid_(0);
  dimensionlocal(10*0+1) = mpigrid_(1);
  dimensionlocal(10*0+2) = mpigrid_(0)*mpigrid_(1);
  
  // 2. Fill in the global data sizes.
  dimensionlocal(10*1+0) = size_(0);
  dimensionlocal(10*1+1) = size_(1);
  dimensionlocal(10*1+2) = size_(2);
  dimensionlocal(10*1+3) = size_(0)*size_(1);
  dimensionlocal(10*1+4) = size_(0)*size_(1)*size_(2);
  dimensionlocal(10*1+5) = size_(2);
  dimensionlocal(10*1+6) = size_(1)/3*2;
  dimensionlocal(10*1+7) = size_(0)/3*2/2+1;
  dimensionlocal(10*1+8) = size_(2)*size_(1)/3*2;
  dimensionlocal(10*1+9) = size_(2)*size_(1)/3*2*(size_(0)/3*2/2+1);
  
  // 3. Get dimensions for data in x-pencil.
  Cp3dfft_get_dims(dimensionlocal.memptr()+10*2+0,
                   dimensionlocal.memptr()+10*2+3,
                   dimensionlocal.memptr()+10*2+6,
                   1);
  dimensionlocal(10*2+9) = dimensionlocal(10*2+6)*
                           dimensionlocal(10*2+7)*
                           dimensionlocal(10*2+8);
  
  // 4. Get dimensions for data in z-pencil.
  Cp3dfft_get_dims(dimensionlocal.memptr()+10*3+0,
                   dimensionlocal.memptr()+10*3+3,
                   dimensionlocal.memptr()+10*3+6,
                   2);
  dimensionlocal(10*3+9) = dimensionlocal(10*3+6)*
                           dimensionlocal(10*3+7)*
                           dimensionlocal(10*3+8);
  
  
  
  
  
  
  //
  //
  double* memoryrealx = new double [dimensionlocal(10*2+9)];
  double* memoryrealz = new double [dimensionlocal(10*2+9)*2];
  
  Cp3dfft_rtran_x2z(memoryrealx,
                    memoryrealz,
                    memoryrealx,
                    memoryrealz,
                    dimensionlocal.memptr()+10*4+0,
                    dimensionlocal.memptr()+10*4+3,
                    dimensionlocal.memptr()+10*4+6);
  
  dimensionlocal(10*4+9) = dimensionlocal(10*4+6)*
                           dimensionlocal(10*4+7)*
                           dimensionlocal(10*4+8);
  
  delete [] memoryrealx;
  delete [] memoryrealz;
  
  
  
  
  
  
  
  
  
  const int sizempi = mpi_.GetSize();
  dimension_.set_size(sizedimension,sizempi);
  // 5. Synchronize the dimension table for every MPI process.
  MPI_Allgather(dimensionlocal.memptr(),
                sizedimension,
                MPI_INT,
                dimension_.memptr(),
                sizedimension,
                MPI_INT,
                MPI_COMM_WORLD);
  
  
  int offsetx = 0;
  int offsetz = 0;
  mpibuffercountsreal_.set_size(sizempi); mpibuffercountsreal_.fill(0);
  mpibufferdisplsreal_.set_size(sizempi); mpibufferdisplsreal_.fill(0);
  mpibuffercountscplx_.set_size(sizempi); mpibuffercountscplx_.fill(0);
  mpibufferdisplscplx_.set_size(sizempi); mpibufferdisplscplx_.fill(0);
  
  // Compute the table for data transformation.
  for (int id=0; id<sizempi; ++id) {
    mpibuffercountsreal_(id) = dimension_(10*2+9,id);
    mpibuffercountscplx_(id) = dimension_(10*3+9,id);
    mpibufferdisplsreal_(id) = offsetx;
    mpibufferdisplscplx_(id) = offsetz;
    offsetx += mpibuffercountsreal_(id);
    offsetz += mpibuffercountscplx_(id);
  }
  
  PS_DEBUG_TRACE_LEAVE(name_+"Configure(void)")
  return true;
}
  
  
  
//------------------------------------------------------------------------------
// Function : Print current status
//------------------------------------------------------------------------------
// Parameter: IN  - stream    > target stream as container
//------------------------------------------------------------------------------
void DomainDecomposition::Print(ostringstream& stream) const {
  PS_DEBUG_TRACE_ENTER(name_+"Print(ostringstream)")
  
  const int rank = mpi_.GetRank();
  
  stream
  << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << "DomainDecomposition > "
  << "Size of Domain ["
  << GetDimension(10*1+0,rank) << " "
  << GetDimension(10*1+1,rank) << " "
  << GetDimension(10*1+2,rank) << "]"
  << std::endl
  << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << ""
  << "MPI Grid ["
  << GetDimension(10*0+0,rank) << " "
  << GetDimension(10*0+1,rank) << "]"
  << std::endl
  << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << ""
  << "Local X-Pencil Data ["
  << GetDimension(10*2+0,rank) << "-" << GetDimension(10*2+3,rank) << " "
  << GetDimension(10*2+1,rank) << "-" << GetDimension(10*2+4,rank) << " "
  << GetDimension(10*2+2,rank) << "-" << GetDimension(10*2+5,rank) << "]"
  << std::endl
  << std::setw(IO_FILE_HEAD_COLUMN_LENGTH) << ""
  << "Local Z-Pencil Data ["
  << GetDimension(10*3+0,rank) << "-" << GetDimension(10*3+3,rank) << " "
  << GetDimension(10*3+1,rank) << "-" << GetDimension(10*3+4,rank) << " "
  << GetDimension(10*3+2,rank) << "-" << GetDimension(10*3+5,rank) << "]"
  << std::endl;
  
  PS_DEBUG_TRACE_LEAVE(name_+"Print(ostringstream)")
}



} // end namespace levelbase
} // end namespace pseudospectral
} // end namespace yueliangyi
