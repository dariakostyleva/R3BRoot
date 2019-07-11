// -------------------------------------------------------------------------
// -----                R3BAsciiGenerator source file                 -----
// -------------------------------------------------------------------------
#include "R3BAsciiGenerator.h"

#include "FairPrimaryGenerator.h"
#include "FairIon.h"
#include "FairRunSim.h"
#include "FairLogger.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TRandom.h"

#include <iostream>

using namespace std;

using std::cout;
using std::endl;
using std::map;
using std::ifstream;

// -----   Default constructor   ------------------------------------------
R3BAsciiGenerator::R3BAsciiGenerator()
  : fInputFile(NULL), fFileName(""),
    fPDG(NULL), fIonMap(), 
    fX(0.), fY(0.), fZ(0.), fPointVtxIsSet(kFALSE),
    fDX(0.), fDY(0.), fDZ(0.), fBoxVtxIsSet(kFALSE)
{
}
// ------------------------------------------------------------------------



// -----   Standard constructor   -----------------------------------------
R3BAsciiGenerator::R3BAsciiGenerator(const char* fileName)
  : fInputFile(NULL), fFileName(fileName),
    fPDG(TDatabasePDG::Instance()), fIonMap(), 
    fX(0.), fY(0.), fZ(0.), fPointVtxIsSet(kFALSE),
    fDX(0.), fDY(0.), fDZ(0.), fBoxVtxIsSet(kFALSE)
{
  cout << "-I- R3BAsciiGenerator: Opening input file " << fileName << endl;
  // Open first the file to register all new ions.
  fInputFile = new ifstream(fFileName);
  if ( ! fInputFile->is_open() ) 
    LOG(fatal) << "R3BAsciiGenerator: Cannot open input file.";
  cout << "-I- R3BAsciiGenerator: Looking for ions..." << endl;

  Int_t nIons = RegisterIons();
  cout << "-I- R3BAsciiGenerator: " << nIons << " ions registered." 
       << endl;
  CloseInput();

  // Re-Open the file for standard streaming ...
  fInputFile = new ifstream(fFileName);
}
// ------------------------------------------------------------------------


R3BAsciiGenerator::R3BAsciiGenerator(const R3BAsciiGenerator& right)
  : fInputFile(right.fInputFile), fFileName(right.fFileName),
    fPDG(right.fPDG), fIonMap(right.fIonMap), 
    fX(right.fX), fY(right.fY), fZ(right.fZ),
    fPointVtxIsSet(right.fPointVtxIsSet),
    fDX(right.fDX), fDY(right.fDY), fDZ(right.fDZ),
    fBoxVtxIsSet(right.fBoxVtxIsSet)
{
}


// -----   Destructor   ---------------------------------------------------
R3BAsciiGenerator::~R3BAsciiGenerator()
{
  CloseInput();
}
// ------------------------------------------------------------------------



// -----   Public method ReadEvent   --------------------------------------
Bool_t R3BAsciiGenerator::ReadEvent(FairPrimaryGenerator* primGen) {

  // Check for input file
  if ( ! fInputFile->is_open() ) {
    cout << "-E- R3BAsciiGenerator: Input file not open!" << endl;
    return kFALSE;
  }

  // Define event variable to be read from file
  Int_t    eventId = 0;
  Int_t    nTracks = 1;
  Double_t pBeam   = 0.;
  Double_t b       = 0.;


  // Define track variables to be read from file
  Int_t    iPid   = -1;
  Int_t    iA      = 30;
  Int_t    iZ      = 17;
  Double_t px      = 0.;
  Double_t py      = 0.;
  Double_t pz      = 0.;
  Double_t vx      = 0.;
  Double_t vy      = 0.;
  Double_t vz      = 0.;
  Double_t iMass      = 0.;


  //******* Defining new particles, registered in evt_gen3.dat ****
  Int_t pdgType_29S = 1000160290;
  //***************************************************************

  //******* Defining angles for isotropic distribution ************
  Double_t costheta = 2.*gRandom->Uniform(0,1)-1.0;
  Double_t sintheta = std::sqrt((1.0 - costheta)*(1.0 + costheta));
  Double_t phi  = 2*TMath::Pi()*gRandom->Uniform(0,1);
  //cout << "cos theta = " << costheta << endl;
  //cout << "theta = " << TMath::ASin(sintheta) << endl;
  //cout << "sintheta = " << sintheta << endl;
  //cout << "phi = " << phi << endl;
  //***************************************************************

  //******* Defining masses, mother energy, tlife **********
  
  Double_t AMS29=29.*0.93149406-0.003156; //mass of 29sulphur
  Double_t QDECAY2 = 0.002;
  Double_t AMASS1 = 0.938272297; //mass of proton
  iMass = AMS29 + AMASS1 + QDECAY2; //mass of 30 chlorine
  Double_t Ekin_30Cl = 0.618*30; //kin energy of chlorine
  //Double_t Ekin_30Cl = 0.002;
  Double_t tlife = 6.58e-19; //tlife of 30Cl in sec
  Double_t iMass1 = iMass + 3.291086E-25/tlife*std::tan(TMath::Pi());
  Double_t Pmom_30Cl = sqrt(Ekin_30Cl*(Ekin_30Cl+2.*iMass1));
  Double_t qvalue_out = 0.; // to check the reproduction of the q-value
  Double_t P_hi = 0.0, Pp = 0.0;         // to check the reproduction of the q-value

  //***************************************************************
  //Double_t P_30Cl[3];
  //P_30Cl[1] = Pmom_30Cl*sintheta*std::cos(phi);
  //P_30Cl[2] = Pmom_30Cl*sintheta*std::sin(phi);
  //P_30Cl[3] = Pmom_30Cl*costheta;

  //******* Two-body decay ****************************************
  Double_t Etot_29S, Etot_29S_fort;
  Double_t Etot_p, Etot_p_fort;
  Double_t Etot_30Cl;
  Etot_30Cl = sqrt(Pmom_30Cl*Pmom_30Cl+iMass1*iMass1);
  Double_t daughtermass[2];
  Double_t daughtermomentum;
  daughtermass[0] = AMS29;
  daughtermass[1] = AMASS1;
  daughtermomentum = Pcm(iMass1,daughtermass[0],daughtermass[1]); 
  cout << "Etot_30Cl = " << Etot_30Cl << endl;
  printf("daughtermomentum = %f\n", daughtermomentum);
  //energies in c.m.s of decay products
  Etot_29S= sqrt(daughtermass[0]*daughtermass[0] + daughtermomentum*daughtermomentum);
  Etot_p= sqrt(daughtermass[1]*daughtermass[1] + daughtermomentum*daughtermomentum); 

  //another method to calc energies, like it was in G3, same result
  /*
  Etot_29S_fort= (iMass1*iMass1 + AMS29*AMS29 - AMASS1*AMASS1)/(2*iMass1);
  Etot_p_fort = sqrt(Pcm_daughter[0]*Pcm_daughter[0] + Pcm_daughter[1]*Pcm_daughter[1]+Pcm_daughter[2]*Pcm_daughter[2]+AMASS1*AMASS1);
  cout << "Etot_p_fort = " << Etot_p_fort <<endl;
  cout << "Etot_p = " << Etot_p <<endl;
  cout << "Etot_29S_fort = " << Etot_29S_fort <<endl;
  */
  cout << "Etot_29S in cm before Lorentz = " << Etot_29S <<endl;
  
  //at this point we have decay products in the rest frame of mother
  //****************************************************************


  //******** Transfromations from REST to LAB **********************
  // from mother we take beta (velocity) and gamma
  // Lorentz boost for daughters in the moving frame of mother (z only)
  Double_t gamma, beta;
  gamma = Etot_30Cl/iMass1;
  beta = - Pmom_30Cl/Etot_30Cl;

 // cout << "beta = " << beta << endl;
 // cout << "gamma = "<< gamma << endl;
  
  //transformation of pcm from spherical coordinates to cartesian (three components)
  Double_t Pcm_daughter[3];
  Pcm_daughter[0] = daughtermomentum*sintheta*std::cos(phi);
  Pcm_daughter[1] = daughtermomentum*sintheta*std::sin(phi);
  Pcm_daughter[2] = daughtermomentum*costheta;

  Double_t Plab_29S[3];
  Double_t Plab_p[3];
  Double_t Elab_29S;

  for (Int_t i=0; i<2; i++) {
    Plab_29S[i] = Pcm_daughter[i];
    Plab_p[i] = (-1)*Pcm_daughter[i];
  }
  cout << "Pcm_daughter[2] " << Pcm_daughter[2] << endl;
  
  //now Lorentz boost on z
  Plab_29S[2] = Pcm_daughter[2] + gamma*beta*((gamma*beta*Pcm_daughter[2])/(gamma+1) - Etot_29S);
  Plab_p[2] = (-1)*Pcm_daughter[2] + gamma*beta*((gamma*beta*(-1)*Pcm_daughter[2])/(gamma+1) - Etot_p);

  Elab_29S = gamma*(Etot_29S - beta*Pcm_daughter[2]);
  cout << "Elab_29S in lab after Lorentz = " << Elab_29S <<endl;


  //********************* check qvalue**********************************************
  Pp = sqrt(Plab_p[2]*Plab_p[2] + Plab_p[1]*Plab_p[1] + Plab_p[0]*Plab_p[0]);
  P_hi = sqrt(Plab_29S[2]*Plab_29S[2] + Plab_29S[1]*Plab_29S[1] + Plab_29S[0]*Plab_29S[0]);
  //this gives exacltly the decay energy
  qvalue_out = sqrt(Pp*Pp+AMASS1*AMASS1) + sqrt(P_hi*P_hi + AMS29*AMS29) - sqrt(Pmom_30Cl*Pmom_30Cl + iMass1*iMass1) - (AMASS1 + AMS29 - iMass1);


    //some stupid calculations 
   //    pcm =  Phi + g*b*((g*b*Phi)/(g+1) + e_hi_lab); //here I calculate pcm from plab
  //qvalue_out = 4*iMass1*iMass1*daughtermomentum*daughtermomentum/((iMass1 + AMASS1 - AMS29)*(iMass1 - AMASS1 + AMS29)*(iMass1 + AMASS1 + AMS29));
  printf("Qvalue from pcm from plab and masses = %f\n",qvalue_out);
 

 // printf("T kin of S before detectors %f\n", sqrt(P_hi*P_hi + AMS29*AMS29) - AMS29);
 // printf("Qvalue from print = %f\n",qvalue_out);

  //qvalue_out = 4*iMass1*Pcm_daughter[2]*Pcm_daughter[2]/((iMass1+AMASS1-AMS29)*(iMass1-AMASS1-AMS29)*(iMass1+AMASS1+AMS29));
 // printf("Qvalue from pcm and masses = %f\n",qvalue_out);

  //********************************************************************************

  //just for "cheking" recoil effect
  // theta determines where the particle in rest frame goes - along mother vector or against
  //if(costheta > 0.9 && costheta <= 1) std::cout<<"Pz of sulphur " << Plab_29S[2] << ", pz of proton " << Plab_p[2] <<" when theta is close to 0 degrees with Pcm_daughter[2] = " << Pcm_daughter[2] << " and theta is "<<TMath::ACos(costheta) << std::endl;
  //if(costheta >= -1. && costheta < 0.) std::cout<<"Pz of sulphur " << Plab_29S[2] << ", pz of proton " << Plab_p[2] <<" when theta is close to 180 degrees with Pcm_daughter[2] = " << Pcm_daughter[2] << " and theta is "<<TMath::ACos(costheta) << std::endl;

  //PRINTOUT OF MOMENTA
  /*
  for (Int_t i=0; i<3; i++) {
    cout << "Pcm_daughter[i] = " << Pcm_daughter[i] << " for i = " << i << endl;
    cout << "Plab_29S[i] = " << Plab_29S[i] << " for i = " << i << endl;
    cout << "Plab_p[i] = " << Plab_p[i] << " for i = " << i << endl;
  }
  */

  //*****************************************************************

  // Read event header line from input file
  //*fInputFile >> eventId >> nTracks >> pBeam >> b;

  if ( fInputFile->eof() ) {
      cout << "-I- R3BAsciiGenerator: End of input file reached " << endl;

      //*****go back to the beggining! i.e. reads the same lines over and over again
      fInputFile->clear();
      fInputFile->seekg(0, ios::beg);
      *fInputFile >> eventId >> nTracks >> pBeam >> b;
      //*****

      //normal way to reach file end and close it
     // CloseInput();
     // return kFALSE;
  }

  cout << "-I- R3BAsciiGenerator: Reading Event: " << eventId << ",  pBeam = "
      << pBeam << "GeV, b = " << b << " fm, multiplicity " << nTracks
      << endl;


  // Loop over tracks in the current event
  for (Int_t itrack=0; itrack<nTracks; itrack++) {

   //*fInputFile >> iPid  >> iZ >> iA >> px >> py >> pz >> vx >> vy >> vz >> iMass;
     
    /*      cout << "-I- R3BAsciiGenerator: iPid: " << iPid <<
	  ",   A = " << iA << " Z = " << iZ <<
	  " px = "  << px <<
	  " py = "  << py <<
	  " pz = "  << pz <<
	  " vx = "  << vx <<
	  " vy = "  << vy <<
	  " vz = " << vz << endl;
    */  

      Int_t pdgType=0;

      // Ion case ( iPid = -1 )
    if ( iPid < 0 ) {
	    char ionName[20];
	    if(1 == iZ && 2 == iA) {
	      sprintf(ionName, "Deuteron");
	    } else if(1 == iZ && 3 == iA) {
	      sprintf(ionName, "Triton");
	    } else if(2 == iZ && 3 == iA) {
	      sprintf(ionName, "HE3");
	    } else if (2 == iZ && 4 == iA) {
	      sprintf(ionName, "Alpha");
	    } else {
	      sprintf(ionName, "Ion_%d_%d", iA, iZ);
	    }
	    TParticlePDG* part = fPDG->GetParticle(ionName);
	    if ( ! part ) {
	        cout << "-W- R3BAsciiGenerator::ReadEvent: Cannot find "
		    << ionName << " in database!" << endl;
	      continue;
	    }
	    pdgType = part->PdgCode();
    }
    //else pdgType = iPid;  // "normal" particle
    else pdgType = iA;  // "normal" particle

      // Give track to PrimaryGenerator
      //cout << "PDG : " << pdgType << endl;

    if (fPointVtxIsSet){ 
      vx = fX;
      vy = fY;
      vz = fZ;
      if (fBoxVtxIsSet){
        vx = gRandom->Gaus(fX,fDX);
        vy = gRandom->Gaus(fY,fDY);
        vz = gRandom->Gaus(fZ,fDZ); 
      }         	
    }
    //cout << "pdgType = " << pdgType << endl;
   // primGen->AddTrack(pdgType, px, py, pz, vx, vy, vz);
    //primGen->AddTrack(1000170300, px, py, pz, vx, vy, vz);

    //Put here the right values for the momentum of the two decay particles/ions: proton and 29S
    primGen->AddTrack(2212, Plab_p[0], Plab_p[1], Plab_p[2], vx, vy, vz);
    primGen->AddTrack(pdgType_29S, Plab_29S[0], Plab_29S[1], Plab_29S[2], vx, vy, vz);

  }//! tracks

  return kTRUE;
}
// ------------------------------------------------------------------------
// -----   Private method CloseInput   ------------------------------------
void R3BAsciiGenerator::CloseInput() {
  if ( fInputFile ) {
    if ( fInputFile->is_open() ) {
       cout << "-I- R3BAsciiGenerator: Closing input file " 
	    << fFileName << endl;
       fInputFile->close();
    }
    delete fInputFile;
    fInputFile = NULL;
  }
}
// ------------------------------------------------------------------------

// -----   Private method RegisterIons   ----------------------------------
Int_t R3BAsciiGenerator::RegisterIons() {

  Int_t nIons = 0;
  Int_t eventId, nTracks;
  Double_t pBeam,b;

  // Define track variables to be read from file
  Int_t    iPid   = -1;
  Int_t    iA      = 0;
  Int_t    iZ      = 0;
  Double_t px      = 0.;
  Double_t py      = 0.;
  Double_t pz      = 0.;
  Double_t vx      = 0.;
  Double_t vy      = 0.;
  Double_t vz      = 0.;
  Double_t iMass      = 0.;

  fIonMap.clear();

  while ( ! fInputFile->eof()) {
    *fInputFile >> eventId >> nTracks >> pBeam >> b;

    if ( fInputFile->eof() ) continue;

    for (Int_t iTrack=0; iTrack<nTracks; iTrack++) {
      *fInputFile >> iPid >> iZ >> iA >> px >> py >> pz >> vx >> vy >> vz >> iMass;

      // Ion Case
      if ( iPid < 0 ) {
	char buffer[20];
 	sprintf(buffer, "Ion_%d_%d", iA, iZ);
	TString ionName(buffer);
	if ( fIonMap.find(ionName) == fIonMap.end() ) {
	  //FairIon* ion = new FairIon(ionName, iZ, iA, iZ);
	  FairIon* ion = new FairIon(ionName, iZ, iA, iZ, 0., iMass);
	  fIonMap[ionName] = ion;
	  nIons++;
	}  // new ion
      } // ion
    }// !tracks

  }//!events

  FairRunSim* run = FairRunSim::Instance();
  map<TString, FairIon*>::iterator mapIt;
  for (mapIt=fIonMap.begin(); mapIt!=fIonMap.end(); mapIt++) {
    FairIon* ion = (*mapIt).second;
    run->AddNewIon(ion);
  }

  return nIons;
}

Double_t R3BAsciiGenerator::Pcm(Double_t e, Double_t p1, Double_t p2) {
   // calcurate momentum of daughter particles in two-body decay
   Double_t ppp = (e+p1+p2)*(e+p1-p2)*(e-p1+p2)*(e-p1-p2)/(4.0*e*e);
   if (ppp>0) return std::sqrt(ppp);
   else       return -1.;
   // mom of decay products is equal in cms and inverse in direction
}
// ------------------------------------------------------------------------

ClassImp(R3BAsciiGenerator)

