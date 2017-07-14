//____________________________________________________________________________
/*
 Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 14, 2009 - CA
   Was first added in v2.5.1 
 @ Apr 24, 2010 - CA
   Add code to decay the off-the-mass-shell W- using PYTHIA6. 
   First complete version of the GLRES event thread.
*/
//____________________________________________________________________________

#include <cstring>

#include <TClonesArray.h>
#include <TMath.h>

#include "Conventions/Constants.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "VHE/GLRESGenerator.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
GLRESGenerator::GLRESGenerator() :
EventRecordVisitorI("genie::GLRESGenerator")
{

}
//___________________________________________________________________________
GLRESGenerator::GLRESGenerator(string config) :
EventRecordVisitorI("genie::GLRESGenerator", config)
{

}
//___________________________________________________________________________
GLRESGenerator::~GLRESGenerator()
{

}
//___________________________________________________________________________
void GLRESGenerator::ProcessEventRecord(GHepRecord * event) const
{
  GHepParticle * nu = event -> Probe();
  GHepParticle * el = event -> HitElectron();
  assert(nu);
  assert(el);

  // Pick electron velocity
  //
  //... neglect for now

  //
  // Add remnant nucleus 
  //
  GHepParticle * target = event -> TargetNucleus();
  if(target) {
    int pdgc = target->Pdg();
    TLorentzVector p4 ( * target->P4() );
    TLorentzVector x4 ( * target->X4() );
    event->AddParticle(pdgc, kIStStableFinalState, 1,-1,-1,-1, p4, x4);
  }

  //
  // Add resonance
  //

  TLorentzVector p4_nu (* nu->P4());
  TLorentzVector p4_el (* el->P4());

  TLorentzVector x4(* el->X4());

  TLorentzVector p4_W = p4_nu + p4_el;

  event->AddParticle(kPdgWM, kIStDecayedState, 0,-1,-1,-1, p4_W, x4);

  //
  // Decay resonance and add decay products
  //

  double mass = p4_W.M();
  fPythia8->Pythia8()->readString("Print:quiet           = on");
  fPythia8->Pythia8()->readString("PDF:lepton            = off");
  fPythia8->Pythia8()->readString("WeakBosonExchange:all = on");

  fPythia8->Pythia8()->settings.mode("Beams:idA", -12); // nu_ebar
  fPythia8->Pythia8()->settings.mode("Beams:idB", 11); // e-
  fPythia8->Pythia8()->settings.mode("Beams:frameType", 1);
  fPythia8->Pythia8()->settings.parm("Beams:eCM", mass);

  fPythia8->Pythia8()->next();
  fPythia8->Pythia8()->event.list();

  // copy PYTHIA container to a new TClonesArray so as to transfer ownership
  // of the container and of its elements to the calling method
  Pythia8::Event &fEvent = fPythia8->Pythia8()->event;
  int numpart = fEvent.size();
  assert(numpart>0);

  // Vector defining rotation from LAB to LAB' (z:= \vec{resonance momentum})
  TVector3 unitvq = p4_W.Vect().Unit();

  // Boost velocity LAB' -> Resonance rest frame
  TVector3 beta(0,0,p4_W.P()/p4_W.Energy());

  for (int i = 1; i < numpart; ++i) {
    int pdgc = fEvent[i].id();
    int ist  = fEvent[i].status();
    if(ist > 0) {
        TLorentzVector p4o(fEvent[i].px(), fEvent[i].py(), fEvent[i].pz(), fEvent[i].e());
        p4o.Boost(beta); 
        TVector3 p3 = p4o.Vect();
        p3.RotateUz(unitvq); 
        TLorentzVector p4(p3,p4o.Energy());
        event->AddParticle(pdgc, kIStStableFinalState, 4,-1,-1,-1, p4, x4);
     }
  }

}
//___________________________________________________________________________
void GLRESGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESGenerator::LoadConfig(void)
{
  fPythia8 = PythiaSingleton::Instance();

  // sync GENIE/PYTHIA8 seed number
  RandomGen::Instance();
}
//____________________________________________________________________________

