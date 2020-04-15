#ifndef Constants_h
#define Constants_h

#include <cstring>
#include <vector>
#include <TObject.h>
#include <TMath.h>

// Set these compiler directives to adjust the logic flow of this file
// #define USEDEBUG // Uncomment to print more debug information
#define FINALDESIGN

#ifdef USEDEBUG
#define showDebug(x) std::cout << x
#else
#define showDebug(x)
#endif

#define NX 9000
#define DX 0.030
#define DY 0.030
#define NY 5500
#define DZ 0.435

// When the final design is used, add the air gap and Al spacer in mid detector:
#ifdef FINALDESIGN
#undef DZ
#define DZ 2.00
#endif

// -------------------------------------

Bool_t   kHelium = true;
Bool_t   kSpotScanning = false;
Int_t    kRotation = 90; // degrees
Bool_t   kPhantom = false;
Bool_t   kDoTracking = true; 
Bool_t   kFilterNuclearInteractions = false; 
Bool_t   kDoDiffusion = true;
Int_t    kEventsPerRun = 100;
Int_t    kSkipTracks = 0; // during readout
const Int_t    kEnergy = 917; // 917 MeV_Helium ~= 230 MeV_proton // 600 HeC phantom
Bool_t   kUseDegrader = true; 
Bool_t   kUseAlpide = true;
Bool_t   kDoTrackerPropagation = true;

const Int_t sizeOfEventID = 25;
const Int_t nChildrenInNode = 2; // max concurrent track segments to follow
Float_t     kMaxTrackScore = 0.3; // This number is a placeholder, optimized through MC scans in Classes/Clusters/findTracks.C
Float_t     kMaxTrackAngle = 0.05; // allow for consecutive 50 mrad changes
Bool_t      kConcatenateHits = true;

// natural unit is mm
const Float_t cm = 0.1;
const Float_t um = 1000;
const Float_t kRad = 3.14159265/180.;

// Some general run parameters
const Int_t nx = NX;
const Int_t ny = NY;
const Int_t nTrackers = 4;
const Float_t kAbsorberThickness = 3.5;

// nLayers are loaded in MaterialConstants.C according to the detector geometry
const Float_t dx = DX; // mm
const Float_t dy = DY; // mm
const Float_t dz = DZ + kAbsorberThickness;
const Float_t dz2 = 52.4;

#ifdef FINALDESIGN
const Bool_t kFinalDesign = true;
#else
const Bool_t kFinalDesign = false;
#endif

// Used for treatment of available experimental data files
const Int_t nEnergies = 6;
Int_t energies[nEnergies] = {122, 140, 150, 170, 180, 188};


Float_t kSpotX = 0;

enum eFrameType {kCalorimeter, kTracker};
enum eDataType {kMC, kData};
enum eMaterial {kTungsten, kAluminum, kPMMA, kWater, kCarbon};
enum eOutputUnit {kPhysical, kWEPL, kUnitEnergy};
const Int_t kMaterial = kAluminum;

Int_t kOutputUnit = kWEPL;
Bool_t kUseCSDA = false; // Use CSDA for range calculations and MC truth input

// Use experimental ALPIDE data clustering model -- empirical model with updated parameters
// Otherwise use "old" gaussian model from FOCAL data
const Bool_t kUseExperimentalClustering = true;
const Bool_t kUseExperimentalClusterPainting = false;
const Bool_t kUseRefinedClustering = true;

// Tracking parameters
Float_t initialSearchRadius = 50 * dx; // 20 if dimensionless
Float_t searchRadius = 40 * dx; // 20 if dimensionless
Float_t kMCSFactor = 3;
const Double_t SpreadNumber = 0.1; // number of iterations in gaussian spread
const Double_t SpreadSigma  = 0.7; // sigma (in pixels) to spread

Bool_t kSaveCWT = false;

// Minimum track length (in mm) to account for the track in the bragg peak fitting analysis.
const Int_t kMinimumTracklength = 5;

// How much above the average edep must the bragg peak (last two layers) be?
const Float_t kBPFactorAboveAverage = 1.3;

Int_t    GlobalLayerID = 0;
Float_t  run_energy = 0;
Float_t  run_degraderThickness = 0;
Int_t    kDataType = kMC;
Long64_t lastJentry_ = 0;

#endif
