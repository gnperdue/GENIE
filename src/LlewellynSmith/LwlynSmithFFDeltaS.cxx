//____________________________________________________________________________
/*
 Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Hugh Gallagher  <hugh.gallagher \at tufts.edu>
         Tufts University

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Conventions/Constants.h"
#include "LlewellynSmith/LwlynSmithFFDeltaS.h"
#include "Messenger/Messenger.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
LwlynSmithFFDeltaS::LwlynSmithFFDeltaS() :
LwlynSmithFF("genie::LwlynSmithFFDeltaS")
{

}
//____________________________________________________________________________
LwlynSmithFFDeltaS::LwlynSmithFFDeltaS(string config) :
LwlynSmithFF("genie::LwlynSmithFFDeltaS", config)
{

}
//____________________________________________________________________________
LwlynSmithFFDeltaS::~LwlynSmithFFDeltaS()
{

}
//____________________________________________________________________________
double LwlynSmithFFDeltaS::F1V(const Interaction * interaction) const
{
  LOG("LwlynSmith", pDEBUG) << "Calling the Strange F1V";
  return LwlynSmithFF::StrangeF1V(interaction);
}
//____________________________________________________________________________
double LwlynSmithFFDeltaS::xiF2V(const Interaction * interaction) const
{
  LOG("LwlynSmith", pDEBUG) << "Calling the Strange xiF2V";
  return LwlynSmithFF::StrangexiF2V(interaction);
}
//____________________________________________________________________________
double LwlynSmithFFDeltaS::FA(const Interaction * interaction) const
{
  LOG("LwlynSmith", pDEBUG) << "Calling the Strange FA";
  return LwlynSmithFF::StrangeFA(interaction);
}
//____________________________________________________________________________
double LwlynSmithFFDeltaS::Fp(const Interaction * interaction) const
{
  return LwlynSmithFF::Fp(interaction);
}
//____________________________________________________________________________
