/*******************************************************************************
 *------------------------------------------------------------------------------
 *------------------------------------------------------------------------------
 ******************************************************************************/
#include <stdexcept>
#include "solver_channel_flat.h"

using namespace std;
using namespace yueliangyi::pseudospectral;


int main(int argc, char **argv) {
  
  
  try {
    leveltwo::FlatChannelSolver solver(argc,argv);
    
    
    solver.Configure();
    if (!solver.isConfigured()) {
    }
    
    
    solver.Execute();
    
    
  } catch (runtime_error& err) {
    cout << err.what() << endl;
    throw err;
  }
  
  
  return 0;
  
}
