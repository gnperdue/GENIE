// @(#)root/pythia6:$Id$
// Author: Piotr Golonka   17/09/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_GMCParticle
#define ROOT_GMCParticle

#ifndef ROOT_TObject
#include "TObject.h"
#endif
#ifndef ROOT_TAttLine
#include "TAttLine.h"
#endif
#ifndef ROOT_TPrimary
#include "TPrimary.h"
#endif


class GMCParticle : public TObject, public TAttLine {

private:

   Int_t    fKS;            // status of particle
   Int_t    fKF;            // KF flavour code
   Int_t    fParent;        // parrent's id
   Int_t    fFirstChild;    // id of first child
   Int_t    fLastChild;     // id of last  child

   Float_t  fPx;            // X momenta [GeV/c]
   Float_t  fPy;            // Y momenta [GeV/c]
   Float_t  fPz;            // Z momenta [GeV/c]
   Float_t  fEnergy;        // Energy    [GeV]
   Float_t  fMass;          // Mass      [Gev/c^2]

   Float_t  fVx;            // X vertex  [mm]
   Float_t  fVy;            // Y vertex  [mm]
   Float_t  fVz;            // Z vertex  [mm]
   Float_t  fTime;          // time of procuction [mm/c]
   Float_t  fLifetime;      // proper lifetime [mm/c]

   // new Pythia8 parameters
   // color and anticolor tags, Les Houches Accord style
   // so far only used in test/gtestHadronization
   Int_t    fColor;         // color tag for Pythia8
   Int_t    fAntiColor;     // anticolor tag for Pythia8

public:
   GMCParticle() : fKS(0), fKF(0), fParent(0), fFirstChild(0),
     fLastChild(0), fPx(0), fPy(0), fPz(0), fEnergy(0), fMass(0),
     fVx(0), fVy(0), fVz(0), fTime(0), fLifetime(0), fColor(0), fAntiColor(0) {}

            GMCParticle(Int_t kS, Int_t kF, Int_t parent,
                        Int_t firstchild, Int_t lastchild,
                        Float_t px, Float_t py, Float_t pz,
                        Float_t energy, Float_t mass,
                        Float_t vx, Float_t vy, Float_t vz,
                        Float_t time, Float_t lifetime,
                        Int_t color=0, Int_t anticolor=0) :

               fKS(kS),
               fKF(kF),
               fParent(parent),
               fFirstChild(firstchild),
               fLastChild(lastchild),
               fPx(px),
               fPy(py),
               fPz(pz),
               fEnergy(energy),
               fMass(mass),
               fVx(vx),
               fVy(vy),
               fVz(vz),
               fTime(time),
               fLifetime(lifetime),
               fColor(color),
               fAntiColor(anticolor) { }


   virtual             ~GMCParticle() { }

   Int_t       GetKS() const {return fKS;}
   Int_t       GetKF() const {return fKF;}
   Int_t       GetParent() const {return fParent;}
   Int_t       GetFirstChild() const {return fFirstChild;}
   Int_t       GetLastChild() const {return fLastChild;}

   Float_t     GetPx() const {return fPx;}
   Float_t     GetPy() const {return fPy;}
   Float_t     GetPz() const {return fPz;}
   Float_t     GetEnergy() const {return fEnergy;}
   Float_t     GetMass() const {return fMass;}

   Float_t     GetVx() const {return fVx;}
   Float_t     GetVy() const {return fVy;}
   Float_t     GetVz() const {return fVz;}
   Float_t     GetTime() const {return fTime;}
   Float_t     GetLifetime() const {return fLifetime;}

   Int_t       GetColor() const {return fColor;}
   Int_t       GetAntiColor() const {return fAntiColor;}
   virtual const char     *GetName() const;

   virtual void        SetKS(Int_t kS) {fKS=kS;}
   virtual void        SetKF(Int_t kF) {fKF=kF;}
   virtual void        SetParent(Int_t parent) {fParent=parent;}
   virtual void        SetFirstChild(Int_t first) {fFirstChild=first;}
   virtual void        SetLastChild(Int_t last) {fLastChild=last;}

   virtual void        SetPx(Float_t px) {fPx=px;}
   virtual void        SetPy(Float_t py) {fPy=py;}
   virtual void        SetPz(Float_t pz) {fPz=pz;}
   virtual void        SetEnergy(Float_t energy) {fEnergy=energy;}
   virtual void        SetMass(Float_t mass) {fMass=mass;}

   virtual void        SetVx(Float_t vx) {fVx=vx;}
   virtual void        SetVy(Float_t vy) {fVy=vy;}
   virtual void        SetVz(Float_t vz) {fVz=vz;}
   virtual void        SetTime(Float_t time) {fTime=time;}
   virtual void        SetLifetime(Float_t lifetime) {fLifetime=lifetime;}

   virtual void        SetColor(Int_t color) {fColor=color;}
   virtual void        SetAntiColor(Int_t anticolor) {fAntiColor=anticolor;}

   virtual void        ls(Option_t* option) const;

   ClassDef(GMCParticle,1)  // LUJETS particles data record.
};

#endif
