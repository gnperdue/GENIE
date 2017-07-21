#include <TPrimary.h>

#include "Fragmentation/PythiaSingleton.h"

#include "Conventions/Constants.h"
#include "Conventions/GBuild.h"
#include "Messenger/Messenger.h"

using namespace genie;
using namespace genie::constants;

ClassImp(PythiaSingleton)

PythiaSingleton * PythiaSingleton::fgInstance = 0;

//____________________________________________________________________________
PythiaSingleton::PythiaSingleton()
{
    // Constructor
    if (fgInstance) {
      // TODO raise assertion?
      LOG("PythiaSingleton", pERROR) <<
          "Instance of PythiaSingleton already exists";
      return;
    }

    fPythia = new Pythia8::Pythia();
}
//____________________________________________________________________________
PythiaSingleton::~PythiaSingleton()
{
    // Destructor
    if (fgInstance) {
        fgInstance->Delete();
        delete fgInstance;
        fgInstance = 0;
    }
    delete fPythia;
}
//____________________________________________________________________________
PythiaSingleton* PythiaSingleton::Instance() 
{
    return fgInstance ? fgInstance : (fgInstance = new PythiaSingleton());
}

