//____________________________________________________________________________
/*!

\class    genie::PDFLIB

\brief    LHAPDF/PDFLIB library interface.
          Concrete implementation of the PDFModelI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  June 06, 2004

\cpright  Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PDFLIB_H_
#define _PDFLIB_H_

#include "PDF/PDFModelI.h"

namespace genie {

class PDFLIB : public PDFModelI {

public:

  PDFLIB();
  PDFLIB(string config);
  virtual ~PDFLIB();

  // Implement PDFModelI interface

  double UpValence   (double x, double Q2) const;
  double DownValence (double x, double Q2) const;
  double UpSea       (double x, double Q2) const;
  double DownSea     (double x, double Q2) const;
  double Strange     (double x, double Q2) const;
  double Charm       (double x, double Q2) const;
  double Bottom      (double x, double Q2) const;
  double Top         (double x, double Q2) const;
  double Gluon       (double x, double Q2) const;
  PDF_t  AllPDFs     (double x, double Q2) const;

  // Override the default "Confugure" implementation 
  // of the Algorithm interface

  void Configure (const Registry & config);
  void Configure (string config);

private:

  void   Initialize          (void) const;
  void   SetPDFSetFromConfig (void) const;
};

}         // genie namespace

#endif    // _PDFLIB_H_
