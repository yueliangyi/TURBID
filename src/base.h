/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#ifndef YLY_PS_LVBASE_BASE_H_
#define YLY_PS_LVBASE_BASE_H_

#include <fftw3-mpi.h>
#include "common.h"

namespace yueliangyi {
namespace pseudospectral {
namespace levelbase {
  
  
  
// The class MPIGroup processes all the prerequisites of MPI for the system. For
// different thread group, specific MPIGroup information should be created. All
// objects created for one thread group will share the same MPIGroup. Remember
// to create the MPIGroup in the beginning.
class MPIGroup {
  
public:
  
  MPIGroup(int& argc, char**& argv);
  ~MPIGroup(void) { MPI_Finalize(); }
  
  int GetSize   (void) const { return size_   ; }
  int GetRank   (void) const { return rank_   ; }
  int GetThreads(void) const { return threads_; }
  string GetHost(void) const { return host_   ; }
  
  bool isMain  (void          ) const { return rank_==0    ? true : false; }
  bool isRank  (const int rank) const { return rank_==rank ? true : false; }
  bool isThread(void          ) const { return threads_; }
  
  void BcastString(string& target, const int root = 0) const;
  void BcastInt   (int   & target, const int root = 0) const;
  void BcastBool  (bool  & target, const int root = 0) const;
  void BcastDouble(double& target, const int root = 0) const;
  
  void Print(ostringstream& stream) const;
  
protected:
  
private:
  
  // Disable copy constructor with feature in C++11.
  MPIGroup(const MPIGroup& other) = delete;
  
  int size_    = 0;       // size of group
  int rank_    = 0;       // rank in group
  int threads_ = 0;       // number of threads (OpenMP)
  string host_;           // name of host
  
  string name_ = "MPIGroup::";
  
}; typedef MPIGroup mpi_t;
  
  
  
// The class DomainDecomposition splits the whole domain into x and z-pencils,
// depending on the number of MPI and using the P3DFFT splitting strategy. Data
// in physical domain is in x-pencil while in wavenumber domain is in z-pencil.
class DomainDecomposition {
  
public:
  
  DomainDecomposition(const mpi_t& mpi, const ivec& mpigrid, const ivec& size) :
    mpi_(mpi),mpigrid_(mpigrid),size_(size) { }
  
  ~DomainDecomposition(void);
  
  const ivec& GetSize   (void           ) const { return size_       ; }
  int         GetSize   (const size_t id) const { return size_(id)   ; }
  const ivec& GetMpiGrid(void           ) const { return mpigrid_    ; }
  int         GetMpiGrid(const size_t id) const { return mpigrid_(id); }
  
  const imat& GetDimension(void                         ) const { return dimension_          ; }
  int         GetDimension(const int loc, const int rank) const { return dimension_(loc,rank); }
  
  const ivec& GetMpiBufferCountsReal(void) const { return mpibuffercountsreal_; }
  const ivec& GetMpiBufferDisplsReal(void) const { return mpibufferdisplsreal_; }
  const ivec& GetMpiBufferCountsCplx(void) const { return mpibuffercountscplx_; }
  const ivec& GetMpiBufferDisplsCplx(void) const { return mpibufferdisplscplx_; }
  
  bool isLoadX(void) const { return GetDimension(10*2+9,mpi_.GetRank()); }
  bool isLoadZ(void) const { return GetDimension(10*3+9,mpi_.GetRank()); }

  // Check whether the thread has the mean component.
  bool isMeanZ(void) const {
    return GetDimension(10*3+1,mpi_.GetRank())==1 &&
           GetDimension(10*3+2,mpi_.GetRank())==1;
  }
  
  //
  bool isRealBottom(void) const {
    return GetDimension(10*2+2,mpi_.GetRank()) == 1;
  }
  bool isRealTop(void) const {
    return GetDimension(10*2+5,mpi_.GetRank()) == size_(2);
  }
  
  bool Configure(void);
  void Print(ostringstream& stream) const;
  
protected:
  
private:
  
  // Disable copy constructor with feature in C++11.
  DomainDecomposition(const DomainDecomposition& other) = delete;
  
  const mpi_t& mpi_    ;
  const ivec & mpigrid_;
  const ivec & size_   ;
  
  // A table contains dimensions from all MPI processes.
  imat dimension_;
  
  // Table created for data transformation between MPI processes.
  ivec mpibuffercountsreal_;
  ivec mpibufferdisplsreal_;
  ivec mpibuffercountscplx_;
  ivec mpibufferdisplscplx_;
  
  string name_ = "DomainDecomposition::";
  
}; typedef DomainDecomposition dcp_t;
  
  
  
// The class File is designed for file IO.
class File {
  
public:
  
  //
  // Definition of I/O operation level.
  // Used for different I/O requirement.
  //
  enum IOControlLevel {
    IOControlLevelOne   = 1,
    IOControlLevelTwo   = 2,
    IOControlLevelThree = 3,
  }; typedef IOControlLevel ioctrl_t;
  
  File(const mpi_t& mpi,
       const dcp_t& dcp,
       const string filename = IO_DEFAULT_FILE_NAME,
       const string directory = IO_DEFAULT_DIRECTORY) :
    mpi_(mpi),dcp_(dcp),directory_(directory),filename_(filename) { }
  
  const mpi_t& GetMpi(void) const { return mpi_; }
  const dcp_t& GetDcp(void) const { return dcp_; }
  
  const string& GetDirectory(void) const { return directory_; }
  const string& GetFileName (void) const { return filename_ ; }
  
  void SetDirectory(const string& target) { directory_ = target; }
  void SetFileName (const string& target) { filename_  = target; }
  
  virtual bool Save(const string subfolder = "",
                    const ioctrl_t level = IOControlLevelOne) = 0;
  virtual bool Load(const string subfolder = "",
                    const ioctrl_t level = IOControlLevelOne) = 0;
  
protected:
  
  const mpi_t & mpi_ ;
  const dcp_t& dcp_;
  
  string directory_;
  string filename_;
  fstream filestream_;
  
  string name_ = "File::";
  
private:
  
};
  
  
  
// The class Base works as a template for upper classes, such as the domain and
// steppers. Derived from this class, the upper ones will share the same
// operations and data members. However, try to reload the configuration if
// the derived class needs special treatment.
class Base : public File {
  
public:
  
  struct Configuration { }; typedef Configuration cfg_t;
  
  Base(const mpi_t& mpi, const dcp_t& dcp, const cfg_t* cfg = NULL) :
    File(mpi,dcp),cfg_(cfg) { }
  virtual ~Base(void) { }
  
  const cfg_t * GetCfg (void) const { return cfg_ ; }
  
  bool isConfigured(void) { return isconfigured_; }
  
  virtual bool Configure(void) = 0;
  virtual void Print(ostringstream& stream) const = 0;
  
protected:
  
  // Disable copy constructor with feature in C++11.
  Base(const Base& other) = delete;
  
  const cfg_t * cfg_  = NULL;
  bool  isconfigured_ = false;
  
  string name_ = "Base::";
  
private:
  
};


} // end namespace levelbase
} // end namespace pseudospectral
} // end namespace yueliangyi

#endif // YLY_PS_LVBASE_BASE_H_
