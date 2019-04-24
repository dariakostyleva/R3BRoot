// -------------------------------------------------------------------------
// -----                      R3BTraPoint source file                  -----
// -------------------------------------------------------------------------

#include "R3BTraPoint.h"

#include <iostream>

using std::cout;
using std::endl;
using std::flush;


// -----   Default constructor   -------------------------------------------
R3BTraPoint::R3BTraPoint() : FairMCPoint() {
  fX_out      = fY_out  = fZ_out  = fEloss = 0.;
  fPx_out     = fPy_out = fPz_out = 0.; fPdi=0;
}
// -------------------------------------------------------------------------



// -----   Standard constructor   ------------------------------------------
R3BTraPoint::R3BTraPoint(Int_t trackID, Int_t detID, Int_t detCopyID,
			 TVector3 posIn, 
			 TVector3 posOut, TVector3 momIn, TVector3 momOut,
			 Double_t tof, Double_t length, Double_t eLoss, int PDI) 
  : FairMCPoint(trackID, detID, posIn, momIn, tof, length, eLoss) {
  fDetCopyID = detCopyID;  // added by Marc
  fX_out  = posOut.X();
  fY_out  = posOut.Y();
  fZ_out  = posOut.Z();
  fEloss  = eLoss;
  fPdi = PDI;
  fPx_out = momOut.Px();
  fPy_out = momOut.Py();
  fPz_out = momOut.Pz();
}
// -------------------------------------------------------------------------



// -----   Destructor   ----------------------------------------------------
R3BTraPoint::~R3BTraPoint() { }
// -------------------------------------------------------------------------




// -----   Public method Print   -------------------------------------------
void R3BTraPoint::Print(const Option_t* opt) const {
  cout << "-I- R3BTraPoint: STS Point for track " << fTrackID 
       << " in detector " << fDetectorID << endl;	
  cout << "    Position (" << fX << ", " << fY << ", " << fZ
       << ") cm" << endl;
  cout << "    Momentum (" << fPx << ", " << fPy << ", " << fPz
       << ") GeV" << endl;
  cout << "    Time " << fTime << " ns,  Length " << fLength 
       << " cm,  Energy loss " << fELoss*1.0e06 << " keV" << endl;
}
// -------------------------------------------------------------------------



// -----   Point x coordinate from linear extrapolation   ------------------
Double_t R3BTraPoint::GetX(Double_t z) const {
  //  cout << fZ << " " << z << " " << fZ_out << endl;
  if ( (fZ_out-z)*(fZ-z) >= 0. ) return (fX_out+fX)/2.;
  Double_t dz = fZ_out - fZ;
  return ( fX + (z-fZ) / dz * (fX_out-fX) );
}
// -------------------------------------------------------------------------



// -----   Point y coordinate from linear extrapolation   ------------------
Double_t R3BTraPoint::GetY(Double_t z) const {
  if ( (fZ_out-z)*(fZ-z) >= 0. ) return (fY_out+fY)/2.;
  Double_t dz = fZ_out - fZ;
  //  if ( TMath::Abs(dz) < 1.e-3 ) return (fY_out+fY)/2.;
  return ( fY + (z-fZ) / dz * (fY_out-fY) );
}
// -------------------------------------------------------------------------



// -----   Public method IsUsable   ----------------------------------------
Bool_t R3BTraPoint::IsUsable() const {
  Double_t dz = fZ_out - fZ;
  if ( TMath::Abs(dz) < 1.e-4 ) return kFALSE;
  return kTRUE;
}
// -------------------------------------------------------------------------


ClassImp(R3BTraPoint)
