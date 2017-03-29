// @(#)root/pythia6:$Id$
// Author: Piotr Golonka   17/09/97
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  class GMCParticle                                                         //
//                                                                            //
// This class serves as a data storage for description of one particle.       //
// It is especially convenient to store information taken from LUJETS common, //
// which is done by interface class TPythia.                                  //
////////////////////////////////////////////////////////////////////////////////

#include "Fragmentation/GMCParticle.h"
#include <TPrimary.h>
#include "PDG/PDGLibrary.h"

ClassImp(GMCParticle)

//______________________________________________________________________________
void GMCParticle::ls(Option_t *) const
{
   printf("(%2i,%4i) <-%3i, =>[%3i,%3i]",fKS,fKF,fParent,
          fFirstChild,fLastChild);
   printf(":  p=(%7.3f,%7.3f,%9.3f) ;",fPx,fPy,fPz);

   printf(" E=%8.3f ; m=%7.3f ; V=(%g,%g,%g); t=%g, tau=%g\n",
          fEnergy,fMass,fVx,fVy,fVz,fTime,fLifetime);
}

//______________________________________________________________________________
const char *GMCParticle::GetName() const
{
   const char * name;
   genie::PDGLibrary * pdg_database = genie::PDGLibrary::Instance();
   name = pdg_database->Find(fKF)->GetName();

   return name;
}
