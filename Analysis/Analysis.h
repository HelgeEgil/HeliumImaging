#ifndef Analysis_h
#define Analysis_h

class TCanvas;
class TH1F;
class TF1;

#include "Classes/Track/Tracks.h"
#include "Classes/Hit/Hits.h"

using namespace DTC;

void		getTracksReconstructionEfficiency(Int_t dataType, Float_t energy = 250, Float_t degraderThickness = 0);
void   		drawTracksRangeHistogram(Int_t Runs, Int_t dataType = kMC, Bool_t recreate = 0, Float_t energy = 250, Float_t degraderThickness = 0, Int_t eventsPerRun = kEventsPerRun, Int_t outputFileIdx = -1, Bool_t drawFitResults = true, Bool_t doTracking = kDoTracking, Bool_t excludeNuclearInteractions = kFilterNuclearInteractions, Int_t skipTracks = 0);
void   		drawTracksDepthDose(Int_t Runs, Int_t tracksperrun = kEventsPerRun, Float_t degraderThickness = 160, Bool_t recreate = true, Bool_t doTracking = kDoTracking);
void		drawTracks3D(Int_t Runs, Int_t tracksperrun = kEventsPerRun, Float_t degraderThickness = 160, Bool_t recreate = true, Bool_t doTracking = kDoTracking);
void        analyseSecondaryFilters(Int_t Runs = 10, Int_t tracksPerRun = 100, Float_t degradetThickness = 100, Bool_t doTracking = true);
void        makeOutputFileForImageReconstructionRad(Int_t Runs, Int_t tracksperrun, Int_t rotation, Int_t useSpotX, TString phantomName);
void        makeOutputFileForImageReconstructionCT(Int_t Runs, Int_t tracksperrun, Int_t rotation, TString phantomName);


#endif
