#ifndef Analysis_cxx
#define Analysis_cxx

#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <TH2.h>
#include <TH3.h>
#include <TPolyLine.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TRandom3.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TAxis3D.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TEllipse.h>
#include <TLegend.h>
#include <TStopwatch.h>
#include <TPaveStats.h>
#include <TView.h>
#include <TLeaf.h>
#include <TArrow.h>
#include <TF1.h>
#include <Math/ProbFunc.h>

#include "Analysis/Analysis.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "GlobalConstants/RangeAndEnergyCalculations.h"
#include "GlobalConstants/Misalign.h"
#include "Classes/Track/Tracks.h"
#include "Classes/Hit/Hits.h"
#include "Classes/DataInterface/DataInterface.h"
#include "HelperFunctions/Tools.h"
#include "HelperFunctions/getTracks.h"

using namespace std;
using namespace DTC;

void getTracksReconstructionEfficiency(Int_t dataType, Float_t energy, Float_t degraderThickness) {
   Int_t nRuns = 0;

   run_energy = energy;
   run_degraderThickness = degraderThickness;

   Int_t nRunArray[13] = {4,8,16,25,32,64,85,100,128,181,256,350,512};

   for (Int_t i=0; i<13; i++) { // 1 -> 30
      nRuns = nRunArray[i];

      kEventsPerRun = nRuns;
      Float_t factor = 2;

      Int_t totalNumberOfRuns = 10000 / kEventsPerRun;
      if (totalNumberOfRuns < 1) totalNumberOfRuns = 1;
      if (totalNumberOfRuns > 10000) totalNumberOfRuns = 10000;

      Tracks * tracks = loadOrCreateTracks(1, totalNumberOfRuns, dataType, energy);

      Track *thisTrack;
      Int_t EID, thisEID;
      Int_t nTotal = 0;
      Int_t nFirstAndLastAllTracks = 0;

      for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
         thisTrack = tracks->At(j);
         if (!thisTrack) continue;
         nTotal++;

         if (thisTrack->isFirstAndLastEventIDEqual() && tracks->getNMissingClustersWithEventID(thisTrack->getEventID(0), thisTrack->Last()->getLayer(), thisTrack->At(0)->getLayer()) == 0) {
            nFirstAndLastAllTracks++;
         }
      }

      Float_t ratioFirstAndLastAllTracks = (float) nFirstAndLastAllTracks / nTotal;
      Float_t readoutAbsorber = (roundf(kAbsorberThickness) == kAbsorberThickness) ? kAbsorberThickness : kAbsorberThickness*10;
   
      Float_t  cutAngle = 60; // 2.5 sigma
      Float_t  cutEdep = 12;
   
//      tracks->removeHighAngularChangeTracks(cutMaxAngle); // mrad
      tracks->doTrackFit();
      tracks->removeShortTracks(4);
      tracks->removeNuclearInteractions();
      tracks->removeTracksWithMinWEPL(60);
      tracks->removeThreeSigmaShortTracks();
      tracks->removeHighAngleTracks(cutAngle); // mrad
      
      Int_t nTotalAfterFilter = 0;
      Int_t nFirstAndLastAllTracksAfterFilter = 0;

      for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
         thisTrack = tracks->At(j);
         if (!thisTrack) continue;
         nTotalAfterFilter++;

         if (thisTrack->isFirstAndLastEventIDEqual() && tracks->getNMissingClustersWithEventID(thisTrack->getEventID(0), thisTrack->Last()->getLayer(), thisTrack->At(0)->getLayer()) == 0) {
            nFirstAndLastAllTracksAfterFilter++;
         }
      }
      
      Float_t ratioFirstAndLastAllTracksAfterFilter = (float) nFirstAndLastAllTracksAfterFilter / nTotalAfterFilter;

      ofstream file2(Form("OutputFiles/lastLayerCorrect_different_nRuns_Helium_5_proton.csv"), ofstream::out | ofstream::app);
      file2 << readoutAbsorber << " " << nRuns << " " << " " << ratioFirstAndLastAllTracks << " " << ratioFirstAndLastAllTracksAfterFilter << endl;
      file2.close();
      
      delete tracks;
   }
}

void drawTracksDepthDose(Int_t Runs, Int_t tracksperrun /* kEventsPerRun */, Float_t degraderThickness /* 160 */, Bool_t recreate /* 1 */, Bool_t doTracking /* 1 */) {
   run_degraderThickness = degraderThickness;
   run_energy = float(kEnergy);
   kEventsPerRun = tracksperrun;
   Int_t dataType = kMC;
  
   if (kUseDegrader) {
      run_energy = getEnergyFromDegraderThickness(degraderThickness);
      printf("Using degrader, expecting nominal residual energy %.2f MeV\n", run_energy);
   }
   
   Bool_t         removeHighAngleTracks = true;
   Bool_t         removeNuclearInteractions = true;
   Bool_t         removeShortTracks = false;
   Float_t        fitRange, fitScale, fitError, fitSigma;
   Int_t          nCutDueToTrackEndingAbruptly = 0;
   Int_t          nPlotX = 3, nPlotY = 3;
   Int_t          skipPlot = 0;
   Int_t          fitIdx = 0, plotSize = nPlotX*nPlotY;
   TGraphErrors * outputGraph;
   char         * sDataType = getDataTypeChar(dataType);
   char         * sMaterial = getMaterialChar();
   TCanvas      * cGraph = new TCanvas("cGraph", "Fitted data points", nPlotX*500, nPlotY*500);

   cGraph->Divide(nPlotX,nPlotY, 0.000001, 0.000001, 0);
   cGraph->cd();
   gPad->SetBorderMode(0); gStyle->SetFrameBorderMode(0);
   gPad->SetTickx(1); gPad->SetTicky(1);
   gPad->SetTopMargin(0.05); gPad->SetRightMargin(0.05);
   gPad->SetBottomMargin(0.05);
   gPad->SetLeftMargin(0.15);

   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, run_energy,15,15);
   tracks->removeEmptyTracks();
//   tracks->fillOutIncompleteTracks(0.4);
   
   Float_t  cutHalo = 12;
   Float_t  cutMaxAngle = 60;
   Float_t  cutAngle = 40;
   Float_t  cutEdep = 10;
  
   Track *thisTrack = nullptr;

   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);
      if (!thisTrack) continue;
      if (tracks->isTrackIncomplete(thisTrack)) {
         thisTrack->Last()->setSecondary(true);
      }
   }

   tracks->doTrackFit(true); 

   for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
//      if (j < skipPlot) continue;

      Track *thisTrack = tracks->At(j);
      if (!thisTrack) continue;

      if (thisTrack->Last()->isSecondary()) continue;
//      if (thisTrack->getRangemm() >50) continue;
      //if (thisTrack->Last()->getDepositedEnergy() > thisTrack->At(thisTrack->GetEntriesFast()-2)->getDepositedEnergy()) continue;
      if (thisTrack->Last()->getDepositedEnergy() > 5 || thisTrack->Last()->isSecondary()) continue;

      // Do track fit, extract all parameters for this track
      outputGraph = (TGraphErrors*) thisTrack->doTrackFit(true, kUseCSDA); // (bool isScaleVariable, bool useTrackLength (~ CSDA))
      if (!outputGraph) continue;

      fitRange = thisTrack->getFitParameterRange();
      fitScale = thisTrack->getFitParameterScale();
      fitError = quadratureAdd(thisTrack->getFitParameterError(), dz*0.28867); // latter term from error on layer position
      fitSigma = thisTrack->getFitParameterSigma();
      Float_t rf = thisTrack->getRiseFactor(5);

//      if (thisTrack->Last()->getDepositedEnergy() > 3.5) continue;

//      if (fitRange > 50) continue;

      if (fitIdx < plotSize) {
         drawIndividualGraphs(cGraph, outputGraph, fitRange, fitScale, fitError, fitIdx++);
         if (thisTrack->Last()->isSecondary()) {
            outputGraph->SetTitle(Form("Secondary track RF %.1f sigma %.1f comb %.1f", rf, fitSigma, rf/fitSigma));
         }
         else {
            outputGraph->SetTitle(Form("Primary track RF %.1f sigma %.1f comb %.1f", rf, fitSigma, rf/fitSigma));
         }
         printf("Drawing plot number %d.\n", fitIdx);
      }
   
      else break;
   }
}

void drawTracksRangeHistogram(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy, Float_t degraderThickness, Int_t eventsPerRun, Int_t outputFileIdx, Bool_t drawFitResults, Bool_t doTracking, Bool_t excludeNuclearInteractions, Int_t skipTracks) {
   run_degraderThickness = degraderThickness;
   run_energy = energy;
   kEventsPerRun = eventsPerRun;
   kDoTracking = doTracking;
   kFilterNuclearInteractions = excludeNuclearInteractions;
  
   if (kUseDegrader) {
      run_energy = getEnergyFromDegraderThickness(degraderThickness);
   }

   printf("Using water degrader of thickness %.0f mm, the initial energy of %.0f MeV is reduced to %.1f MeV.\n", degraderThickness, energy, run_energy);

   kDataType = dataType;
   Bool_t         removeHighAngleTracks = true;
   Bool_t         removeNuclearInteractions = true;
   Bool_t         drawVerticalLayerLines = false;

   Float_t        finalEnergy = 0;
   Float_t        fitRange, fitScale, fitError;
   Int_t          nCutDueToTrackEndingAbruptly = 0;
   TGraphErrors * outputGraph;
   char         * sDataType = getDataTypeChar(dataType);
   char         * sMaterial = getMaterialChar();
   char         * hTitle = Form("Fitted energy of a %.2f MeV beam in %s (%s)", run_energy, sMaterial, sDataType);

   Int_t nEnergyBins = getUnitFromEnergy(run_energy);
   gStyle->SetOptTitle(0);

   if (kUseDegrader) {
      hTitle = Form("Fitted energy of a %.0f MeV nominal beam on %s DTC w/%.1f mm water degrader", energy, sMaterial, degraderThickness);
   }

   Float_t lowHistogramLimit = getUnitFromEnergy(0);
   Float_t highHistogramLimit = getUnitFromEnergy(run_energy)*1.2 + 10;
   if (isnan(highHistogramLimit)) highHistogramLimit = getUnitFromEnergy(run_energy) + 30;
   TH1F * hFitResults = new TH1F("fitResult", hTitle, fmax(nEnergyBins,100), lowHistogramLimit, highHistogramLimit);
   TH1F * hLastLayer = new TH1F("hLastLayer", "Last Layer;Layer;Entries", 50, 0, 50);
 
   printf("Using material: %s\n", sMaterial);
   printf("Histogram limits: %.2f to %.2f.\n", lowHistogramLimit, highHistogramLimit);
   printf("At energy %.0f, expecting range %.2f mm and WEPL %.2f mm.\n", run_energy, getTLFromEnergy(run_energy), getWEPLFromEnergy(run_energy));
   printf("This corresponds to a WEPL factor of %.2f.\n", getWEPLFactorFromEnergy(run_energy));
   printf("(Using direct calculation, nominal WET is %.2f mm)\n", getWETFromDegrader(run_degraderThickness));

   hFitResults->SetLineColor(kBlack); hFitResults->SetFillColor(kGreen-5);

   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
   
   Float_t  cutHalo = 12;
   Float_t  cutMaxAngle = 60;
   Float_t  cutAngle = 50;
   Float_t  cutEdep = 10;

   tracks->removeThreeSigmaShortTracks();


   for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
      Track *thisTrack = tracks->At(j);
      if (!thisTrack) continue;
      if (j < skipTracks) continue;
    
      if (thisTrack->doesTrackEndAbruptly()) {
         nCutDueToTrackEndingAbruptly++;
      }

      // Do track fit, extract all parameters for this track
      outputGraph = (TGraphErrors*) thisTrack->doTrackFit(false, false); // kUseCSDA); // (bool isScaleVariable, bool useTrackLength (~ CSDA))
      if (!outputGraph) continue;
   
      delete outputGraph;

      fitRange = thisTrack->getFitParameterRange();
//      hFitResults->Fill(getUnitFromTL(fitRange));
      hFitResults->Fill(fitRange);
      hLastLayer->Fill(thisTrack->Last()->getLayer());
   }
  
   // Draw expected gaussian distribution of results from initial energy
   Float_t expectedStraggling = 0, expectedMean = 0, dlayer_down = 0, dlayer = 0;
   Float_t separationFactor = 0.9, nullStraggling = 0;
   Float_t sigma_energy = getSigmaEnergy(run_energy);
   
   expectedMean = getUnitFromEnergy(run_energy);
   expectedStraggling = getUnitStragglingFromEnergy(run_energy, sigma_energy);
   nullStraggling = getUnitStragglingFromEnergy(run_energy, 0);
  
   expectedMean = getWETFromDegrader(run_degraderThickness);

   cout << "OutputUnit is " << kOutputUnit << " and the expected mean value is " << expectedMean 
       << ". The straggling including / excluding energy variation is " << expectedStraggling << " / " << nullStraggling << ".\n";
       
   Float_t means[10] = {};
   Float_t sigmas[10] = {};

   TF1 *gauss = doSimpleGaussianFit(hFitResults, means, sigmas, outputFileIdx);
   Float_t empiricalMean = means[9];
   Float_t empiricalSigma = sigmas[9];
   
   Float_t energySigma = getEnergyFromUnit(empiricalMean +  empiricalSigma/ 2) - getEnergyFromUnit(empiricalMean - empiricalSigma / 2);

   if (drawFitResults) {
      TCanvas * cFitResults = new TCanvas("cFitResults", hTitle, 1000, 1000);
      
      if       (kOutputUnit == kPhysical) hFitResults->SetXTitle("Physical range [mm]");
      else if  (kOutputUnit == kWEPL)     hFitResults->SetXTitle("Range in Water Equivalent Path Length [mm]");
      else if  (kOutputUnit == kUnitEnergy)   hFitResults->SetXTitle("Energy [MeV]");

      hFitResults->SetYTitle("Number of protons");
      hFitResults->GetXaxis()->SetTitleFont(22);
      hFitResults->GetXaxis()->SetLabelFont(22);
      hFitResults->GetYaxis()->SetTitleFont(22);
      hFitResults->GetYaxis()->SetLabelFont(22);
      hFitResults->GetYaxis()->SetTitleOffset(1.5);

      hFitResults->SetTitle("");   
      hFitResults->Draw();
      gPad->Update();

      if (drawVerticalLayerLines) {
         TLine *l = nullptr;
         Float_t line_z = 0;
         for (Int_t i=0; i<65; i++) {
            line_z = getUnitFromTL(getLayerPositionmm(i));
            if (line_z > gPad->GetUxmax()) break;
            l = new TLine(line_z, 0, line_z, hFitResults->GetMaximum()*1.05);
            l->SetLineColor(kBlack); l->SetLineWidth(2); l->Draw();
         }
      }

      Float_t calibratedMean = empiricalMean * 1.0036 + 2.93; 

      gPad->Update();
      gStyle->SetOptStat(11);
      TPaveStats *ps = (TPaveStats*) cFitResults->GetPrimitive("stats");
      hFitResults->SetBit(TH1::kNoStats);
      ps->SetX1NDC(0.4176); ps->SetX2NDC(0.9257);
      ps->SetY1NDC(0.7415); ps->SetY2NDC(0.9712);
      ps->SetTextFont(22);
      ps->AddText(Form("Nominal WEPL = %.2f #pm %.2f", expectedMean, expectedStraggling));
      ps->AddText(Form("Calculated WEPL = %.2f #pm %.2f", empiricalMean, empiricalSigma));
      ps->AddText(Form("WEPL deviation = %.2f #pm %.2f", empiricalMean - expectedMean, sqrt(pow(empiricalSigma, 2) - pow(expectedStraggling, 2))));
      ps->AddText(Form("Cal. WEPL deviation = %.2f #pm %.2f", calibratedMean - expectedMean, sqrt(pow(empiricalSigma, 2) - pow(expectedStraggling, 2))));
      cFitResults->Modified();

      cFitResults->SaveAs(Form("OutputFiles/RangeHistogram/%.0f_%.0f.png", kAbsorberThickness, degraderThickness));

      TCanvas *c2 = new TCanvas();
      hLastLayer->Draw();
   }

   delete tracks;
}

void analyseHelium(Int_t Runs, Int_t tracksperrun /* kEventsPerRun */, Float_t degraderThickness /* 160 */, Bool_t recreate /* 1 */, Bool_t doTracking /* 1 */) {
   Int_t dataType = kMC;
   Float_t energy = float(kEnergy);
   run_energy = energy;
   run_degraderThickness = degraderThickness;
   kEventsPerRun = tracksperrun;
   kDoTracking = doTracking;
   
   Clusters * savedClusters = new Clusters();

   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy,0,0, savedClusters);
   tracks->removeEmptyTracks();
   printf("Loaded tracks, checking for incomplete tracks...\n");
   
   Track * thisTrack;
   Int_t n=0;
   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);
      if (!thisTrack)   continue;
      if (thisTrack->isIncomplete()) {
         thisTrack->Last()->setSecondary(true);
         n++;
      }
   }
   printf("Done! Set %d (MC tagged) too short tracks to 2nd.\n", n);
   
   Cluster *c = nullptr;
   // propagate secondary status
   Int_t nSecFirst = 0;

   Int_t nSecondariesEnd = 0, nPrimariesEnd = 0, nTotal = tracks->GetEntriesFast();
   Float_t nPrimaries = kEventsPerRun * Runs;
   for (Int_t i=0; i<nTotal; i++) {
      if (!tracks->At(i)) continue;
      nSecondariesEnd += tracks->At(i)->Last()->isSecondary();
      nPrimariesEnd += (!tracks->At(i)->Last()->isSecondary());
   }

   Float_t  cutHalo = 15;
   Float_t  cutMaxAngle = 70;
   Float_t  cutAngle = 50;
   Float_t  cutEdep = 3.5;
   Float_t  cutChi2 = 500;
   Float_t  cutEdepPlateau = 1.5;
   if (kHelium) {
      cutHalo = 15;
      cutMaxAngle = 70;
      cutAngle = 45; // 2.5 sigma 
      cutEdep = 8;
      cutChi2 = 500;
      cutEdepPlateau = 3.5;
   }

   nSecondariesEnd = 87625; // Since CS filter removes many
   nPrimariesEnd = 49958; // Since CS filter removes many

   printf("Found in total %d tracks. %d primaries, %d secondaries\n", nTotal, nPrimariesEnd, nSecondariesEnd);

//   tracks->doTrackFit(true);


// Use
   tracks->removeEmptyTracks();
   tracks->removeShortTracks(4);
   tracks->removeNuclearInteractions();
   tracks->removeThreeSigmaShortTracks();
//   tracks->removeHighAngleTracks(cutAngle); // mrad
  
// Don't use
//   track->removeHighChiSquare(210);
//   tracks->removeTracksEndingInHalo(30);
//   tracks->removeHighAngularChangeTracks(cutMaxAngle); // mrad

   Int_t nfSecondariesEnd = 0, nfPrimariesEnd = 0, nfTotal = 0;
   Int_t nfHe = 0, nfHe3 = 0;
   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      c = tracks->At(i)->Last();
      nfSecondariesEnd += c->isSecondary();
      nfPrimariesEnd += (!c->isSecondary());
      nfHe += (c->isSecondary() && c->getPDG() == 1000020040);
      nfHe3 += (c->isSecondary() && c->getPDG() == 1000020030);
   }
   nfTotal = tracks->GetEntries();

   printf("After filtering: %d tracks. %d primaries, %d secondaries\n", nfTotal, nfPrimariesEnd, nfSecondariesEnd);
   printf("%d He, %d He3\n", nfHe, nfHe3);

   printf("True positives (filtered secondaries): %.2f %%\n", float(nSecondariesEnd-nfSecondariesEnd)/nPrimaries*100);
   printf("False positives (filtered primaries): %.2f %%\n", float(nPrimariesEnd-nfPrimariesEnd)/nPrimaries*100);
   printf("True negatives (unfiltered primaries): %.2f %%\n", (nfPrimariesEnd)/float(nPrimaries) * 100);
   printf("False negatives (unfiltered secondaries): %.2f %%\n", (nfSecondariesEnd)/float(nPrimaries) * 100);

   Cluster *a, *b;
   Float_t angle, range, edep, bragg;

   TCanvas *cPosition = new TCanvas("cPosition", "cPosition", 1200, 600);
   cPosition->Divide(2,1,1e-5,1e-5);
   TH1F *hPositionP = new TH1F("hPositionP", "Primary particle;Radial profile first layer;Entries",100,0,100);
   TH1F *hPositionS = new TH1F("hPositionS", "Secondary particle;Radial profile first layer;Entries",100,0,100);
   
   TCanvas *cRiseFactor = new TCanvas("cRiseFactor", "cRiseFactor", 1200, 600);
   cRiseFactor->Divide(2,1,1e-5,1e-5);
   TH2F *hRiseFactorP = new TH2F("hRiseFactorP", "Primary;Edep;line/BK sigma;",100,0,50,100,0,5);
   TH2F *hRiseFactorS = new TH2F("hRiseFactorS", "Secondary;Edep;line/BK sigma;",100,0,50,100,0,5);

   TCanvas *cEdep = new TCanvas("cEdep","cEdep",1200,600);
   cEdep->Divide(2,1,1e-5,1e-5);
   TH1F * edepP = new TH1F("edepP", "Primary particle;E_{dep} last layer [kev/#mum];Entries", 100, 0, 100);
   TH1F * edepS = new TH1F("edepS", "Secondary particle;E_{dep} last layer [kev/#mum];Entries", 100, 0, 100);

   TCanvas *cEdep2D = new TCanvas("cEdep2D","cEdep2D",1200,600);
   cEdep2D->Divide(2,1,1e-5,1e-5);
   TH2F *edep2DP = new TH2F("edep2DP", "Primary particle;Edep plateau [keV/#mum];Edep last layer [keV/#mum]", 60, 0, 30, 60, 0, 30);
   TH2F *edep2DS = new TH2F("edep2DS", "Secondary particle;Edep plateau [keV/#mum];Edep last layer [keV/#mum];", 60, 0, 30,60, 0, 30);

   TCanvas *cAllEdep = new TCanvas("cAllEdep","cAllEdep",1200,600);
   cAllEdep->Divide(2,1,1e-5,1e-5);
   TH1F * cAllEdepP = new TH1F("allEdepP", "Primary particle;Average E_{dep} (plateau) [keV/#mum];Entries",100,0,20);
   TH1F * cAllEdepS = new TH1F("allEdepS", "Secondary particle;Average E_{dep} (plateau) [kev/#mum];Entries",100,0,20);
   
   TCanvas *cAllVarEdep = new TCanvas("cAllVarEdep","cAllVarEdep",1200,600);
   cAllVarEdep->Divide(2,1,1e-5,1e-5);
   TH1F * cAllVarEdepP = new TH1F("allVarEdepP", "Primary particle;E_{dep} all clusters;Entries",100,0,75);
   TH1F * cAllVarEdepS = new TH1F("allVarEdepS", "Secondary particle;E_{dep} all clusters;Entries",100,0,75);
   
   TCanvas *cAngle = new TCanvas("cAngle","cAngle",1200,600);
   cAngle->Divide(2,1,1e-5,1e-5);
   TH1F * angleP = new TH1F("angleP", "Primary particle;Incoming angle [mrad];Entries;", 100, 0, 200);
   TH1F * angleS = new TH1F("angleS", "Secondary particle;Incoming angle [mrad];Entries;", 100, 0, 200);

   TCanvas *cAngleAll = new TCanvas("cAngleAll", "cAngleAll", 700, 600);
   TH1F * angleAll = new TH1F("angleAll", "Primary + secondary;Incoming angle [mrad];Entries;", 100, 0, 200);
   
   TCanvas *cRange = new TCanvas("cRange","cRange",1200,600);
   cRange->Divide(2,1,1e-5,1e-5);
   TH1F * rangeP = new TH1F("rangeP", "Primary particle;Residual range [WEPL mm];Entries;", 100, 0, 300);
   TH1F * rangeS = new TH1F("rangeS", "Secondary particle;Residual range [WEPL mm];Entries;", 100, 0, 300);

   TCanvas *cRange2 = new TCanvas("cRange2", "cRange", 800, 600);
   TH1F * rangePS = new TH1F("rangePS", "All;Residual range [WEPL mm];Entries;", 100, 0, 300);

   TCanvas *cBragg = new TCanvas("cBragg","cBragg",1200,600);
   cBragg->Divide(2,1,1e-5,1e-5);
   TH1F * braggP = new TH1F("braggP", "Primary particle;Depth-dose fit log_{10} #chi^{2};Entries;", 100, 1, 2000);
   TH1F * braggS = new TH1F("braggS", "Secondary particle;Depth-dose fit log_{10} #chi^{2};Entries;", 100, 0, 2000);
   
   TCanvas *cSecondary = new TCanvas();
   TH1F *hSec = new TH1F("hSec", "Secondary prouction depth;Detector layer;Entries", 46, -0.5, 45.5);
   TH1F *hPrim = new TH1F("hPrim", "Primary beam;layer;freq", 50, 0, 50);

   TCanvas *cMaxDeltaTheta = new TCanvas("cMaxDeltaTheta","cMaxDeltaTheta", 1200,600);
   cMaxDeltaTheta->Divide(2,1,1e-5,1e-5);
   TH1F *hMaxDeltaThetaP = new TH1F("hMaxDeltaThetaP", "Primary particle;Max layer-wise angular change [mrad];Entries",100,0,300);
   TH1F *hMaxDeltaThetaS = new TH1F("hMaxDeltaThetaS", "Secondary particle;Max layer-wise angular change [mrad];Entries",100,0,300);

   TCanvas *cCS = new TCanvas("cCS", "Cluster sizes", 1200,600);
   cCS->Divide(2,1,1e-5,1e-5);
   TH1F *hCSP = new TH1F("hCSP", "Primary particle;Cluster size;Entries", 50, 0, 50);
   TH1F *hCSS = new TH1F("hCSS", "Secondary particle;Cluster size;Entries", 50, 0, 50);

   TCanvas *cHits = new TCanvas("cHits", "cHits", 1200,600);
   cHits->Divide(2,1,1e-5,1e-5);
   TH1F *cHitsP = new TH1F("cHitsP", "Primary particle;Beam profile at entrance layer #sqrt{x^{2}+y^{2}} [mm];Entries", 100, 0, 100);
   TH1F *cHitsS = new TH1F("cHitsS", "Secondary particle;Beam profile at entrance layer #sqrt{x^{2}+y^{2}} [mm];Entries", 100, 0, 100);

   for (int i=0; i<tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);
      if (!thisTrack) continue;
      if (!thisTrack->Last()) continue;
//      printf("Track length is GEF %d / GE %d\n", thisTrack->GetEntriesFast(), thisTrack->GetEntries());

      bool was = false;
      float biggestChange = 0;
      for (int j=0; j<thisTrack->GetEntriesFast(); j++) {
         if (!thisTrack->At(j)) continue;
         if (!was) {
            if (thisTrack->At(j)->isSecondary()) {
               hSec->Fill(thisTrack->getLayer(j));
               was = true;
            }
         }
      }
      int last = thisTrack->GetEntriesFast() - 1;
      if (last>=5) {
         Float_t meanEdep = 0, sdEdep = 0;
         for (Int_t j=0; j<thisTrack->GetEntriesFast()-5; j++) {
            if (!thisTrack->At(j)) continue;
            meanEdep += thisTrack->At(j)->getDepositedEnergy();
         }

         if (thisTrack->GetEntriesFast() < 5 ) meanEdep = thisTrack->At(0)->getDepositedEnergy();
         else meanEdep /= (thisTrack->GetEntriesFast()-5);

         if (!thisTrack->Last()->isSecondary()) {
               cAllEdepP->Fill(meanEdep);
               for (Int_t k=0; k<thisTrack->GetEntriesFast(); k++) {
                  if (!thisTrack->At(k)) continue;
                  cAllVarEdepP->Fill(thisTrack->At(k)->getDepositedEnergy());
               }
               edep2DP->Fill(thisTrack->Last()->getDepositedEnergy(), meanEdep);
//               hRiseFactorP->Fill(thisTrack->Last()->getDepositedEnergy(), thisTrack->getRiseFactor(5)/thisTrack->getFitParameterSigma());
         }
         else {
               cAllEdepS->Fill(meanEdep);
               for (Int_t k=0; k<thisTrack->GetEntriesFast(); k++) {
                  if (!thisTrack->At(k)) continue;
                  cAllVarEdepS->Fill(thisTrack->At(k)->getDepositedEnergy());
               }
               edep2DS->Fill(thisTrack->Last()->getDepositedEnergy(), meanEdep);
//               hRiseFactorS->Fill(thisTrack->Last()->getDepositedEnergy(), thisTrack->getRiseFactor(5)/thisTrack->getFitParameterSigma());
         }
      }
      
      Cluster *a = thisTrack->At(0);
      Cluster *b = thisTrack->At(1);
      if (!a || !b) continue;
      angle = getDotProductAngle(a, a, b) * 1000;
      if (isnan(thisTrack->getFitParameterRange())) continue;

      range = getUnitFromTL(thisTrack->getRangemm());
      bragg = thisTrack->getFitParameterChiSquare();
      edep = thisTrack->getDepositedEnergy(last);

      rangePS->Fill(range);
      Float_t maxRadius = 0;
      for (Int_t j=0; j<thisTrack->GetEntriesFast(); j++) {
         if (!thisTrack->At(j)) continue;
         maxRadius = max(maxRadius, sqrt(pow(thisTrack->At(j)->getXmm(),2) + pow(thisTrack->At(j)->getYmm(),2)));
      }

      angleAll->Fill(angle);
      if (!thisTrack->Last()->isSecondary()) {
         edepP->Fill(edep);
         angleP->Fill(angle);
         rangeP->Fill(range);
         braggP->Fill(bragg);
         hPositionP->Fill(maxRadius);
         hCSP->Fill(thisTrack->Last()->getSize());
         hMaxDeltaThetaP->Fill(thisTrack->getMaximumSlopeAngleChange() * 3.1415 / 180 * 1000);
      }
      else {
         edepS->Fill(edep);
         angleS->Fill(angle);
         rangeS->Fill(range);
         braggS->Fill(bragg);
         hCSS->Fill(thisTrack->Last()->getSize());
         hPositionS->Fill(maxRadius);
         hMaxDeltaThetaS->Fill(thisTrack->getMaximumSlopeAngleChange() * 3.1415 / 180 * 1000);
      }
   }

   rangePS->SetFillColor(kRed);
   edepP->SetFillColor(kGray);
   edepS->SetFillColor(kGreen-3);
   angleP->SetFillColor(kGray);
   angleS->SetFillColor(kGreen-3);
   rangeP->SetFillColor(kGray);
   rangeS->SetFillColor(kGreen-3);
   braggP->SetFillColor(kGray);
   braggS->SetFillColor(kGreen-3);
   cAllEdepP->SetFillColor(kGray);
   cAllEdepS->SetFillColor(kGreen-3);
   cAllVarEdepP->SetFillColor(kGray);
   cAllVarEdepS->SetFillColor(kGreen-3);
   hSec->SetFillColor(kGreen-3);
   hMaxDeltaThetaP->SetFillColor(kGray);
   hMaxDeltaThetaS->SetFillColor(kGreen-3);
   cHitsP->SetFillColor(kGray);
   cHitsS->SetFillColor(kGreen-3);
   hCSP->SetFillColor(kGray);
   hCSS->SetFillColor(kGreen-3);
   hPositionP->SetFillColor(kGray);
   hPositionS->SetFillColor(kGreen-3);
//   hRiseFactorP->SetFillColor(kGray);
//   hRiseFactorS->SetFillColor(kGreen-3);

   // Please tell me how to do this properly with gStyle
   rangePS->SetLineColor(kBlack);
   edepP->SetLineColor(kBlack);
   edepS->SetLineColor(kBlack);
   angleP->SetLineColor(kBlack);
   angleS->SetLineColor(kBlack);
   rangeP->SetLineColor(kBlack);
   rangeS->SetLineColor(kBlack);
   braggP->SetLineColor(kBlack);
   braggS->SetLineColor(kBlack);
   cAllEdepP->SetLineColor(kBlack);
   cAllEdepS->SetLineColor(kBlack);
   cAllVarEdepP->SetLineColor(kBlack);
   cAllVarEdepS->SetLineColor(kBlack);
   hSec->SetLineColor(kBlack);
   hMaxDeltaThetaP->SetLineColor(kBlack);
   hMaxDeltaThetaS->SetLineColor(kBlack);
   cHitsP->SetLineColor(kBlack);
   cHitsS->SetLineColor(kBlack);
   hCSP->SetLineColor(kBlack);
   hCSS->SetLineStyle(kBlack);
   hPositionP->SetLineColor(kBlack);
   hPositionS->SetLineColor(kBlack);
//  hRiseFactorP->SetLineColor(kBlack);
//  hRiseFactorS->SetLineColor(kBlack);

   TLine *l1, *l2;

   cEdep->cd(1);
   edepP->Draw();
   l1 = new TLine(cutEdep, 0, cutEdep, edepP->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();
   cEdep->cd(2);
   edepS->Draw();
   l1 = new TLine(cutEdep, 0, cutEdep, edepS->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();

   cAngle->cd(1);
   angleP->Draw();
   l1 = new TLine(cutAngle, 0, cutAngle, angleP->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();
   cAngle->cd(2);
   angleS->Draw();
   l1 = new TLine(cutAngle, 0, cutAngle, angleS->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();

   cAngleAll->cd();
   angleAll->Draw();

   Float_t mu = rangePS->GetXaxis()->GetBinCenter(rangePS->GetMaximumBin());
   Float_t sigma = 10; // as expected from ~homogenous area // WEPL 

   cRange->cd(1);
   rangeP->Draw();
   l1 = new TLine(mu-3*sigma, 0, mu-3*sigma, rangeP->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();
   l2 = new TLine(mu+3*sigma, 0, mu+3*sigma, rangeP->GetMaximum()*1.05);
   l2->SetLineStyle(7);
   l2->SetLineColor(kBlack);
   l2->Draw();

   cRange->cd(2);
   rangeS->Draw();
   l1 = new TLine(mu-3*sigma, 0, mu-3*sigma, rangeS->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();
   l2 = new TLine(mu+1.5*sigma, 0, mu+1.5*sigma, rangeS->GetMaximum()*1.05);
   l2->SetLineStyle(7);
   l2->SetLineColor(kBlack);
   l2->Draw();

   cRange2->cd();
   rangePS->Draw();
   l1 = new TLine(mu-3*sigma, 0, mu-3*sigma, rangePS->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();
   l2 = new TLine(mu+1.5*sigma, 0, mu+1.5*sigma, rangePS->GetMaximum()*1.05);
   l2->SetLineStyle(7);
   l2->SetLineColor(kBlack);
   l2->Draw();

   cBragg->cd(1);
   braggP->Draw();
   l1 = new TLine(cutChi2, 0, cutChi2, braggP->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();
   cBragg->cd(2);
   braggS->Draw();
   l1 = new TLine(cutChi2, 0, cutChi2, braggS->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();

   cAllEdep->cd(1);
   cAllEdepP->Draw();
   l1 = new TLine(cutEdepPlateau, 0, cutEdepPlateau, cAllEdepP->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();
   cAllEdep->cd(2);
   cAllEdepS->Draw();
   l1 = new TLine(cutEdepPlateau, 0, cutEdepPlateau, cAllEdepS->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();
   
   cAllVarEdep->cd(1);
   cAllVarEdepP->Draw();
   cAllVarEdep->cd(2);
   cAllVarEdepS->Draw();

   cEdep2D->cd(1);
   edep2DP->Draw("colz");
   cEdep2D->cd(2);
   edep2DS->Draw("colz");

   cPosition->cd(1);
   hPositionP->Draw();
   cPosition->cd(2);
   hPositionS->Draw();
   
   cRiseFactor->cd(1);
   hRiseFactorP->Draw("colz");
   cRiseFactor->cd(2);
   hRiseFactorS->Draw("colz");

   cHits->cd(1);
   cHitsP->Draw();
   l1 = new TLine(cutHalo, 0, cutHalo, cHitsP->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();
   cHits->cd(2);
   cHitsS->Draw();
   l1 = new TLine(cutHalo, 0, cutHalo, cHitsS->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();

   cCS->cd(1);
   hCSP->Draw();
   cCS->cd(2);
   hCSS->Draw();

   // True positive: Secondaries removed by filter
   // True negative: Primaries not removed filter
   // False positive: Primaries removed by filter
   // False negative: Secondaries not removed by filter
   // Normalize by the total number of particles

   int bmin, bmax;
   int ntracks = tracks->GetEntries();

   bmin = cHitsP->GetXaxis()->FindBin(cutHalo);
   int totalHalo = cHitsS->Integral() + cHitsP->Integral();
   int removedByHalo = cHitsS->Integral(bmin,100) + cHitsP->Integral(bmin,100);
   int removedByHalo2nd = cHitsS->Integral(bmin,100);
   
   float haloTP = cHitsS->Integral(bmin,100) / float(totalHalo);
   float haloTN = (cHitsP->Integral() - cHitsP->Integral(bmin,100)) / float(totalHalo);
   float haloFP = cHitsP->Integral(bmin,100) / float(totalHalo);
   float haloFN = (cHitsS->Integral() - cHitsS->Integral(bmin,100)) / float(totalHalo);

   bmin = angleS->GetXaxis()->FindBin(40);
   int removedByAngle = angleS->Integral(bmin,100) + angleP->Integral(bmin,100);
   int removedByAngle2nd = angleS->Integral(bmin,100);
   
   float angleTP = angleS->Integral(bmin,100) / float(ntracks);
   float angleTN = (angleP->Integral() - angleP->Integral(bmin,100)) / float(ntracks);
   float angleFP = angleP->Integral(bmin,100) / float(ntracks);
   float angleFN = (angleS->Integral() - angleS->Integral(bmin,100)) / float(ntracks);

   bmin = hMaxDeltaThetaS->GetXaxis()->FindBin(50);
   int removedByMaxAngle = hMaxDeltaThetaS->Integral(bmin,100) + hMaxDeltaThetaP->Integral(bmin,100);
   int removedByMaxAngle2nd = hMaxDeltaThetaS->Integral(bmin,100);
   
   float maxAngleTP = hMaxDeltaThetaS->Integral(bmin,100) / float(ntracks);
   float maxAngleTN = (hMaxDeltaThetaP->Integral() - hMaxDeltaThetaP->Integral(bmin,100)) / float(ntracks);
   float maxAngleFP = hMaxDeltaThetaP->Integral(bmin,100) / float(ntracks);
   float maxAngleFN = (hMaxDeltaThetaS->Integral() - hMaxDeltaThetaS->Integral(bmin,100)) / float(ntracks);

   bmin = rangeS->GetXaxis()->FindBin(mu-3*sigma);
   bmax = rangeS->GetXaxis()->FindBin(mu+3*sigma);
   int removedByR = rangeS->Integral(0,bmin-1) + rangeS->Integral(bmax+1,100) + rangeP->Integral(0,bmin-1) + rangeP->Integral(bmax+1,100);
   int removedByR2nd = rangeS->Integral(0,bmin-1) + rangeS->Integral(bmax+1,100);
   
   float rangeTP = (rangeS->Integral() - rangeS->Integral(bmin-1,bmax+1)) / float(ntracks);
   float rangeTN = (rangeP->Integral(bmin,bmax)) / float(ntracks);
   float rangeFP = (rangeP->Integral() - rangeP->Integral(bmin-1,bmax+1)) / float(ntracks);
   float rangeFN = (rangeS->Integral(bmin,bmax)) / float(ntracks);

   bmin = edepP->GetXaxis()->FindBin(cutEdep);
   int removedByED = edepS->Integral(0,bmin) + edepP->Integral(0,bmin);
   int removedByED2nd = edepS->Integral(0,bmin);
   
   float edepTP = edepS->Integral(0,bmin) / float(ntracks);
   float edepTN = (edepP->Integral() - edepP->Integral(0,bmin)) / float(ntracks);
   float edepFP = edepP->Integral(0,bmin) / float(ntracks);
   float edepFN = (edepS->Integral() - edepS->Integral(0,bmin)) / float(ntracks);

   bmin = braggP->GetXaxis()->FindBin(cutChi2);
   int removedByChi2 = braggS->Integral(bmin,100) + braggP->Integral(bmin,100);
   int removedByChi22nd = braggS->Integral(bmin,100);
   
   float chi2TP = braggS->Integral(bmin,100) / float(ntracks);
   float chi2TN = (braggP->Integral() - braggP->Integral(bmin,100)) / float(ntracks);
   float chi2FP = braggP->Integral(bmin,100) / float(ntracks);
   float chi2FN = (braggS->Integral() - braggS->Integral(bmin,100)) / float(ntracks);

   bmin = cAllEdepP->GetXaxis()->FindBin(20);
   int removedByAllEdep = cAllEdepS->Integral(0,bmin) + cAllEdepP->Integral(0, bmin);
   int removedByAllEdep2nd = cAllEdepS->Integral(0,bmin);
   
   float allEdepTP = cAllEdepS->Integral(0,bmin) / float(ntracks);
   float allEdepTN = (cAllEdepP->Integral() - cAllEdepP->Integral(0,bmin)) / float(ntracks);
   float allEdepFP = cAllEdepP->Integral(0,bmin) / float(ntracks);
   float allEdepFN = (cAllEdepS->Integral() - cAllEdepS->Integral(0,bmin)) / float(ntracks);
      
   printf("Total number of tracks = %d\n", ntracks);
   printf("Tracks removed by Halo: %d (%d secondaries) (%.1f%% purity) (%.1f%% removed)\n", removedByHalo, removedByHalo2nd, 100.*removedByHalo2nd/removedByHalo, 100.*removedByHalo/totalHalo);
   printf("Tracks removed by Angle: %d (%d secondaries) (%.1f%% purity) (%.1f%% removed)\n", removedByAngle, removedByAngle2nd, 100.*removedByAngle2nd/removedByAngle, 100.*removedByAngle/ntracks);
   printf("Tracks removed by MaxAngle: %d (%d secondaries) (%.1f%% purity) (%.1f%% removed)\n", removedByMaxAngle, removedByMaxAngle2nd, 100.*removedByMaxAngle2nd/removedByMaxAngle, 100.*removedByMaxAngle/ntracks);
   printf("Tracks removed by R: %d (%d secondaries) (%.1f%% purity) (%.1f%% removed)\n", removedByR, removedByR2nd, 100.*removedByR2nd/removedByR, 100.*removedByR/ntracks);
   printf("Tracks removed by last ED: %d (%d secondaries) (%.1f%% purity) (%.1f%% removed)\n", removedByED, removedByED2nd, 100.*removedByED2nd/removedByED, 100.*removedByED/ntracks);
   printf("Tracks removed by last2 ED: %d (%d secondaries) (%.1f%% purity) (%.1f%% removed)\n", removedByAllEdep, removedByAllEdep2nd, 100.*removedByAllEdep2nd/removedByAllEdep, 100.*removedByAllEdep/ntracks);
   printf("Tracks removed by Chi2: %d (%d secondaries) (%.1f%% purity) (%.1f%% removed)\n", removedByChi2, removedByChi22nd, 100.*removedByChi22nd/removedByChi2, 100.*removedByChi2/ntracks);

   printf("True/False Positive/Negative: \n\n");
   printf("Halo: TP = %.1f, FP = %.1f, TN = %.1f, FN = %.1f\n", haloTP*100, haloFP*100, haloTN*100, haloFN*100);
   printf("Angle: TP = %.1f, FP = %.1f, TN = %.1f, FN = %.1f\n", angleTP*100, angleFP*100, angleTN*100, angleFN*100);
   printf("MaxAngle: TP = %.1f, FP = %.1f, TN = %.1f, FN = %.1f\n", maxAngleTP*100, maxAngleFP*100, maxAngleTN*100, maxAngleFN*100);
   printf("Range: TP = %.1f, FP = %.1f, TN = %.1f, FN = %.1f\n", rangeTP*100, rangeFP*100, rangeTN*100, rangeFN*100);
   printf("Edep: TP = %.1f, FP = %.1f, TN = %.1f, FN = %.1f\n", edepTP*100, edepFP*100, edepTN*100, edepFN*100);
   printf("Chi2: TP = %.1f, FP = %.1f, TN = %.1f, FN = %.1f\n", chi2TP*100, chi2FP*100, chi2TN*100, chi2FN*100);
   printf("AllEdep: TP = %.1f, FP = %.1f, TN = %.1f, FN = %.1f\n", allEdepTP*100, allEdepFP*100, allEdepTN*100, allEdepFN*100);

   cSecondary->cd();
//   hPrim->Draw();
   hSec->Draw();

   cMaxDeltaTheta->cd(1);
   hMaxDeltaThetaP->Draw();
   l1 = new TLine(cutMaxAngle, 0, cutMaxAngle, hMaxDeltaThetaP->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();
   cMaxDeltaTheta->cd(2);
   hMaxDeltaThetaS->Draw();
   l1 = new TLine(cutMaxAngle, 0, cutMaxAngle, hMaxDeltaThetaS->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();

   Int_t particleListIndex[100] = {};
   Int_t particleListSort[100] = {};
   Int_t particleList[100] = {};
   Int_t idxPDG = 0, pdg;
   Bool_t inList = false;

   // last particle
   for (Int_t i=0; i<=tracks->GetEntriesFast(); i++) {
      Track *thisTrack = tracks->At(i);
      if (!thisTrack) continue;

      c = thisTrack->Last();
      if (!c) continue;
      pdg = c->getPDG();
      if (!c->isSecondary()) continue;
      
      inList = false;
      for (int j=0; j<idxPDG; j++) {
         if (particleListIndex[j] == pdg) {
            inList = true;
            particleList[j]++;
            break;
         }
      }
      if (!inList) {
         particleListIndex[idxPDG] = pdg;
         particleList[idxPDG++] = 1;
      }
   }

   // Make indexed sorting list
   Int_t currentMax = 0;
   Int_t maxIdx = 0;
   Int_t sortIdx = 0;
   while (sortIdx < idxPDG) {
      currentMax = 0;
      for (Int_t i=0; i<idxPDG; i++) {
         if (particleList[i] >= currentMax) {
            inList = false;
            for (Int_t j=0; j<sortIdx; j++) inList += (particleListSort[j] == i);
            if (inList) continue;

            currentMax = particleList[i];
            maxIdx = i;
         }
      }
      particleListSort[sortIdx++] = maxIdx;
   }

   string name;
   Int_t p;

   printf("Secondary particle types: \n");
   for (Int_t i=0; i<idxPDG; i++) {
      p = particleListIndex[particleListSort[i]];
      if (p == 11) name = "Electron";
      else if (p == -11) name = "Positron";
      else if (p == 2212) name = "Proton";
      else if (p == 22) name = "Gamma";
      else if (p == 1000030040) name = "Litium4"; 
      else if (p == 1000030060) name = "Litium";
      else if (p == 1000020040) name = "Helium";
      else if (p == 1000020030) name = "Helium3";
      else if (p == 1000010030) name = "Tritium";
      else if (p == 1000010020) name = "Deuterium";
      else if (p == 1000130270) name = "Aluminium";
      else if (p == 1000140280) name = "Silicon";
      else if (p == 1000140260) name = "Silicon26";
      else if (p == 1000140299) name = "Silicon29";
      else if (p == 1000120240) name = "Magnesium";
      else name = Form("%d", p);

      cout << name << ": " << particleList[particleListSort[i]] << " tracks\n";
   }

   TCanvas *cPDG = new TCanvas("cPDG", "PDG(z)");
   cPDG->cd();
   float xfrom = -0.5;
   float xto = 40.5;
   int xlen = 41;
   
   gStyle->SetOptStat(0);

   TH1F *hPDGe = new TH1F("hPDGe", "Electrons", xlen, xfrom, xto);
   TH1F *hPDGpos = new TH1F("hPDGpos", "Positrons", xlen, xfrom, xto);
   TH1F *hPDGpro = new TH1F("hPDGpro", ";Depth in detector [layer number];Fraction of secondary species", xlen, xfrom, xto);
   TH1F *hPDGneu = new TH1F("hPDGneu", ";Depth in detector [layer number];Fraction of secondary species", xlen, xfrom, xto);
   TH1F *hPDGhe = new TH1F("hPDGhe", "Helium", xlen, xfrom, xto);
   TH1F *hPDGhe3 = new TH1F("hPDGhe3", ";Depth in detector [layer number];Fraction of secondary species", xlen, xfrom, xto);
   TH1F *hPDGgam = new TH1F("hPDGgam", "Gamma", xlen, xfrom, xto);
   TH1F *hPDGli = new TH1F("hPDGli", "Litium", xlen, xfrom, xto);
   TH1F *hPDGdeu = new TH1F("hPDGdeu", "Deuterium", xlen, xfrom, xto);
   TH1F *hPDGtri = new TH1F("hPDGtri", "Tritium", xlen, xfrom, xto);
   TH1F *hPDGsi = new TH1F("hPDGsi", "Silicon", xlen, xfrom, xto);
   TH1F *hPDGmg = new TH1F("hPDGmg", "Magnesium", xlen, xfrom, xto);
   TH1F *hPDGal = new TH1F("hPDGal", "Aluminum", xlen, xfrom, xto);

   Int_t l;
   for (Int_t i=0; i<=tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);
      if (!thisTrack) continue;
      for (Int_t j=0; j<=thisTrack->GetEntriesFast(); j++) {
         if (!thisTrack->At(j)) continue;
         l = thisTrack->getLayer(j);
         p = thisTrack->At(j)->getPDG();
         if (!thisTrack->At(j)->isSecondary()) continue;

         if (p == 11) hPDGe->Fill(l);
         else if (p == -11) hPDGpos->Fill(l);
         else if (p == 2212) hPDGpro->Fill(l);
         else if (p == 2112) hPDGneu->Fill(l);
         else if (p == 22) hPDGgam->Fill(l);
         else if (p == 1000030040) hPDGli->Fill(l);
         else if (p == 1000030060) hPDGli->Fill(l);
         else if (p == 1000020040) hPDGhe->Fill(l);
         else if (p == 1000020030) hPDGhe3->Fill(l);
         else if (p == 1000010030) hPDGtri->Fill(l);
         else if (p == 1000010020) hPDGdeu->Fill(l);
         else if (p == 1000130270) hPDGal->Fill(l);
         else if (p == 1000140280) hPDGsi->Fill(l);
         else if (p == 1000140260) hPDGsi->Fill(l);
         else if (p == 1000140299) hPDGsi->Fill(l);
         else if (p == 1000120240) hPDGmg->Fill(l);
      }
   }

   hPDGhe->SetLineColor(kOrange); hPDGpro->SetLineColor(1); hPDGdeu->SetLineColor(2); hPDGhe3->SetLineColor(3);
   hPDGtri->SetLineColor(4); hPDGe->SetLineColor(kMagenta); hPDGgam->SetLineColor(6); hPDGpos->SetLineColor(7);
   hPDGal->SetLineColor(14); hPDGsi->SetLineColor(16); hPDGmg->SetLineColor(18); hPDGli->SetLineColor(12);
   hPDGneu->SetLineColor(28); 

   hPDGhe->SetLineWidth(3); hPDGpro->SetLineWidth(3); hPDGdeu->SetLineWidth(3); hPDGhe3->SetLineWidth(3);
   hPDGtri->SetLineWidth(3); hPDGe->SetLineWidth(3); hPDGli->SetLineWidth(3); hPDGgam->SetLineWidth(3);
   hPDGal->SetLineWidth(3); hPDGsi->SetLineWidth(3); hPDGmg->SetLineWidth(3); hPDGpos->SetLineWidth(3);
   hPDGneu->SetLineWidth(3);

   // NORMALIZE|
   hPDGhe->Scale(100 / nPrimaries);
   hPDGpro->Scale(100 / nPrimaries);
   hPDGdeu->Scale(100 / nPrimaries);
   hPDGhe3->Scale(100 / nPrimaries);
   hPDGtri->Scale(100 / nPrimaries);
   hPDGe->Scale(100 / nPrimaries);
   hPDGgam->Scale(100 / nPrimaries);
   hPDGpos->Scale(100 / nPrimaries);
   hPDGal->Scale(100 / nPrimaries);
   hPDGsi->Scale(100 / nPrimaries);
   hPDGmg->Scale(100 / nPrimaries);
   hPDGli->Scale(100 / nPrimaries);
   hPDGneu->Scale(100 / nPrimaries);

   hPDGpro->GetYaxis()->SetTitleOffset(1.5);
   hPDGhe3->GetYaxis()->SetTitleOffset(1.5);

   TLegend *leg = new TLegend(0.8,0.4,0.95,0.95);
   leg->AddEntry(hPDGpro, "Protons", "L");
   leg->AddEntry(hPDGdeu, "Deuterium", "L");
   leg->AddEntry(hPDGe, "Electrons", "L");
   leg->AddEntry(hPDGhe3, "Helium3", "L");
   leg->AddEntry(hPDGtri, "Tritium", "L");
   leg->AddEntry(hPDGhe, "Helium4 (2^{nd})", "L"); 

   hPDGhe3->Draw("hist"); 
   hPDGpro->Draw("same hist");
   
   hPDGhe->Draw("same hist"); 

   if (hPDGe->Integral()) {
      hPDGe->Draw("same hist"); 
   }
   
   hPDGpro->GetYaxis()->SetRangeUser(0, 1.1*max(hPDGhe->GetMaximum(), max(hPDGpro->GetMaximum(), hPDGhe3->GetMaximum())));
   hPDGhe3->GetYaxis()->SetRangeUser(0, 1.1*max(hPDGhe->GetMaximum(), max(hPDGpro->GetMaximum(), hPDGhe3->GetMaximum())));
   
   if (hPDGdeu->Integral()) { 
      hPDGdeu->Draw("same hist"); 
   }
   
   
   if (hPDGtri->Integral()) { 
      hPDGtri->Draw("same hist"); 
   }

   if (hPDGneu->Integral()) {
   //   hPDGneu->Draw("same hist");
   }

   if (hPDGpos->Integral() && false) {
      hPDGpos->Draw("same hist");
      leg->AddEntry(hPDGpos, "Positrons", "L");
   }

   if (hPDGgam->Integral() && false) { 
      hPDGgam->Draw("same hist");
      leg->AddEntry(hPDGgam, "Gamma", "L");
   }

   /*
   if (hPDGli->Integral()) { 
      hPDGli->Draw("same hist"); 
      leg->AddEntry(hPDGli, "Litium", "L");
   }
*/
/*
   if (hPDGsi->Integral()) { 
      hPDGsi->Draw("same hist"); 
      leg->AddEntry(hPDGsi, "Silicon", "L");
   }

   if (hPDGmg->Integral()) { 
      hPDGmg->Draw("same hist"); 
      leg->AddEntry(hPDGmg, "Magnesium", "L");
   }

   if (hPDGal->Integral()) { 
      hPDGal->Draw("same hist"); 
      leg->AddEntry(hPDGal, "Aluminum", "L");
   }
*/
   leg->SetTextFont(22);
   leg->SetTextSize(0.045);
   leg->Draw();
 
   // ADD % TO AXIS
   TText *t = new TText();
   hPDGpro->GetYaxis()->SetLabelOffset(5);
   hPDGhe3->GetYaxis()->SetLabelOffset(5);
   t->SetTextAlign(32);
   t->SetTextSize(0.05);
   t->SetTextFont(22);
   Float_t at;
   for (Int_t i=0; i<=4;i++) {
      // max > i*6
      // i < max/6
      
      at = i/4. * 20; // 10*round(max(hPDGhe->GetMaximum(), max(hPDGpro->GetMaximum(), hPDGhe3->GetMaximum())/10);
      t->DrawText(-0.42, at, Form("%.2f%%", at));
   }

}

void makeOutputFileForImageReconstructionRad(Int_t Runs, Int_t tracksperrun, Int_t rotation, Int_t useSpotX, TString phantomName) {
//   run_energy = 760;
   kSplitSpotColumnsPerRotation = true;
   kPhantom = true;
   kSpotScanning = true;
   kSaveCWT = false;
   kPhantomName = phantomName;

   run_degraderThickness = 180; // This number is useless I hope
   run_energy = 230;
   if (kHelium) run_energy = 917;

   kEventsPerRun = tracksperrun;
   kDoTracking = true;
   kSpotX = useSpotX; // Use this to open correct file
   kRotation = rotation;

   Bool_t  kDraw = false;

   Float_t lastEdep, ultEdep, penUltEdep, antePenUltEdep;
   Int_t ultLayer, penUltLayer, antePenUltLayer;
   Bool_t isSingleEventID, isLastHitSecondary, isAngleTooHigh;
   Float_t spotXlim, spotYlim, spotStep;

   spotStep = 7;

   if (kPhantomName == "linePair" || kPhantomName == "CTP404") {
      spotYlim = 28;
   }

   else if (kPhantomName == "wedge") {
      spotYlim = 4;
   }

   else { // Headphantom
     spotYlim = 88;
   }  

   TH2F *hSimpleImage = nullptr;
   TH2F *hSimpleImageNorm = nullptr;

   if (kDraw) {
      hSimpleImage = new TH2F("hSimpleImage", "Simple Image;X [mm];Y [mm]", 100, -125, -75, 150, -75,75);
      hSimpleImageNorm = new TH2F("hSimpleImageNorm", "Simple Image;X [mm];Y [mm]", 100, -125, -75, 150, -75,75);
   }

   Float_t outSpotX, outSpotY, outWEPL, outX2x, outX2y, outP2x, outP2y, residualRange, wepl, wepl_calibrated, outWEPL_uncalibrated;
   Track * thisTrack = nullptr;
   Float_t trackLateralDistance, trackOutgoingAngle;

   TFile *fOut = nullptr;
   if (!kHelium) {
      fOut = new TFile(Form("OutputFiles/%s/%s_rotation%03ddeg_spotx%04.f.root", kPhantomName.Data(), kPhantomName.Data(), kRotation, kSpotX), "recreate");
   }
   else {
      fOut = new TFile(Form("OutputFiles/%s_he/%s_rotation%03ddeg_spotx%04.f.root", kPhantomName.Data(), kPhantomName.Data(), kRotation, kSpotX), "recreate");
   }

   TTree *tOut = new TTree("WEPLData", "proton beam");

   tOut->Branch("spotX", &outSpotX, "spotX/F");
   tOut->Branch("spotY", &outSpotY, "spotY/F");
   tOut->Branch("WEPL", &outWEPL, "WEPL/F");
   tOut->Branch("X2x", &outX2x, "X2x/F");
   tOut->Branch("X2y", &outX2y, "X2y/F");
   tOut->Branch("P2x", &outP2x, "P2x/F");
   tOut->Branch("P2y", &outP2y, "P2y/F");
   tOut->Branch("ultEdep", &ultEdep, "ultEdep/F");
   tOut->Branch("penUltEdep", &penUltEdep, "penUltEdep/F");
   tOut->Branch("antePenUltEdep", &antePenUltEdep, "antePenUltEdep/F");
   tOut->Branch("ultLayer", &ultLayer, "ultLayer/I");
   tOut->Branch("penUltLayer", &penUltLayer, "penUltLayer/I");
   tOut->Branch("antePenUltLayer", &antePenUltLayer, "antePenUltLayer/I");
//   tOut->Branch("isSingleEventID", &isSingleEventID, "isSingleEventID/O");
//   tOut->Branch("isLastHitSecondary", &isLastHitSecondary, "isLastHitSecondary/O");
//   tOut->Branch("isAngleTooHigh", &isAngleTooHigh, "isAngleTooHigh/O");
//   tOut->Branch("trackLateralDifference", &trackLateralDistance, "trackLateralDifference/F");
//   tOut->Branch("trackOutgoingAngle", &trackOutgoingAngle, "trackOutgoingAngle/F");

   Int_t    nruns = 0;
   Float_t  cutMaxAngle = 60;
   Float_t  cutAngle = 50;
   Float_t  cutEdep = 12;

   Float_t spotX = float(useSpotX);
   kSpotX = spotX;
   Float_t angleXmrad, angleYmrad;

   // Find the angles based on the geometry
   angleXmrad = getAngleAtSpot(spotX) * 1000;
   Float_t maxWEPL = 333.7;
   if (kHelium) maxWEPL = 332.3;

   Cluster *ultCluster, *penUltCluster, *antePenUltCluster;

   for (float spotY = -spotYlim; spotY <= spotYlim; spotY += spotStep) {
      kSpotY = spotY;
      angleYmrad = getAngleAtSpot(spotY) * 1000;
      printf("Spot (%.0f,%.0f)\n", spotX, spotY);
      Tracks * tracks = loadOrCreateTracks(true, Runs, false, run_energy, spotX, spotY);
      printf("lastJentry = %lld\n", lastJentry_);
      tracks->removeHighAngleTracksRelativeToSpot(45, angleXmrad, angleYmrad);
      tracks->removeTracksWithMinWEPL(30);
//      tracks->removeThreeSigmaShortTracks();
  
      for (Int_t i=0; i<=tracks->GetEntriesFast(); i++) {
         thisTrack = tracks->At(i);
         if (!thisTrack) continue;

         outSpotX = spotX;
         outSpotY = spotY;
         wepl = thisTrack->getFitParameterRange();
         wepl_calibrated =  1.0024 * wepl - 0.41; // pol1 calibration for protons
//         outWEPL = maxWEPL - wepl_calibrated; // 333.7 mm: 230 MeV proton @ 78 eV H2O (extrapolated from data...)
         outWEPL = maxWEPL - wepl; // 333.7 mm: 230 MeV proton @ 78 eV H2O (extrapolated from data...)
         if (isnan(outWEPL)) continue;
         if (!thisTrack->At(0) || !thisTrack->At(1)) continue;

         outX2x = thisTrack->At(0)->getXmm();
         outX2y = thisTrack->At(0)->getYmm();
         outP2x = (thisTrack->At(1)->getXmm() - outX2x) / dz2;
         outP2y = (thisTrack->At(1)->getYmm() - outX2y) / dz2;

         ultCluster = thisTrack->Last();
         penUltCluster = thisTrack->At(thisTrack->GetEntriesFast()-2);
         antePenUltCluster = thisTrack->At(thisTrack->GetEntriesFast()-3);

         if (!ultCluster || !penUltCluster || !antePenUltCluster) continue;

         ultEdep = ultCluster->getDepositedEnergy();
         penUltEdep = penUltCluster->getDepositedEnergy();
         antePenUltEdep = antePenUltCluster->getDepositedEnergy();

         ultLayer = ultCluster->getLayer();
         penUltLayer = penUltCluster->getLayer();
         antePenUltLayer = antePenUltCluster->getLayer();

         /*
         lastEdep = thisTrack->Last()->getDepositedEnergy();
         isSingleEventID = thisTrack->isFirstAndLastEventIDEqual();
         isLastHitSecondary = thisTrack->Last()->isSecondary();
         
         Int_t last = thisTrack->GetEntriesFast()-1;
         Cluster *c0 = thisTrack->At(0);
         Cluster *c1 = thisTrack->At(last-1);
         Cluster *c2 = thisTrack->At(last);
         
         if (c1&&c2) {
            trackOutgoingAngle = getDotProductAngle(c1, c1, c2);
         }
         else {
            trackOutgoingAngle = 0;
         }

         if (c0&&c2) {
            trackLateralDistance = sqrt(pow(c2->getXmm() - c0->getXmm(),2) + pow(c2->getYmm() - c0->getYmm(),2));
         }
         else {
            trackLateralDistance = 0;
         }
*/
         tOut->Fill();

         if (kDraw) {
            hSimpleImage->Fill(outX2x, outX2y, outWEPL);
            hSimpleImageNorm->Fill(outX2x, outX2y);
         }
      }
       nruns++;
       delete tracks;
   }

   fOut->Write();
   fOut->Close();
   printf("Image Reconstruction Input data written to Output/%s/%s_rotation%03d_spotx%04.f.root.\n", kPhantomName.Data(), kPhantomName.Data(), kRotation, kSpotX);
   
   if (kDraw) {
      hSimpleImage->Divide(hSimpleImageNorm);
      hSimpleImage->Draw("COLZ");
   }
}


void makeOutputFileForImageReconstructionCT(Int_t Runs, Int_t tracksperrun, Int_t rotation, TString phantomName) {
//   run_energy = 760;
   kSplitSpotColumnsPerRotation = false;
   kPhantomName = phantomName;
   kPhantom = true;
   kSpotScanning = true;
   kSaveCWT = false;

   run_degraderThickness = 180; // This number is useless I hope
   run_energy = 230;
   kEventsPerRun = tracksperrun;
   kDoTracking = true;
   kRotation = rotation;

   Bool_t  kDraw = false;
   
   Float_t spotXlim, spotYlim, spotStep;

   spotStep = 7;

   if (kPhantomName == "linePair" || kPhantomName == "CTP404") {
      spotXlim = 84;
      spotYlim = 28;
   }

   else if (kPhantomName == "wedge") {
      spotXlim = 84;
      spotYlim = 7;
   }

   else { // Headphantom
     spotXlim = 98;
     spotYlim = 88;
   }  

   TH2F *hSimpleImage = nullptr;
   TH2F *hSimpleImageNorm = nullptr;

   if (kDraw) {
      hSimpleImage = new TH2F("hSimpleImage", "Simple Image;X [mm];Y [mm]", 100, -125, -75, 150, -75,75);
      hSimpleImageNorm = new TH2F("hSimpleImageNorm", "Simple Image;X [mm];Y [mm]", 100, -125, -75, 150, -75,75);
   }

   Float_t outSpotX, outSpotY, outWEPL, outX2x, outX2y, outP2x, outP2y, residualRange, wepl, wepl_calibrated, outWEPL_uncalibrated;
   Float_t lastEdep;
   Bool_t isSingleEventID, isLastHitSecondary;
   Track * thisTrack = nullptr;

   TFile *fOut = new TFile(Form("OutputFiles/%s/%s_rotation%03ddeg.root", kPhantomName.Data(), kPhantomName.Data(), kRotation), "recreate");
   TTree *tOut = new TTree("WEPLData", "proton beam");

   tOut->Branch("spotX", &outSpotX, "spotX/F");
   tOut->Branch("spotY", &outSpotY, "spotY/F");
   tOut->Branch("WEPL", &outWEPL, "WEPL/F");
   tOut->Branch("X2x", &outX2x, "X2x/F");
   tOut->Branch("X2y", &outX2y, "X2y/F");
   tOut->Branch("P2x", &outP2x, "P2x/F");
   tOut->Branch("P2y", &outP2y, "P2y/F");
//    tOut->Branch("lastEdep", &lastEdep, "lastEdep/F");
//   tOut->Branch("isSingleEventID", &isSingleEventID, "isSingleEventID/O");
//   tOut->Branch("isLastHitSecondary", &isLastHitSecondary, "isLastHitSecondary/O");

   Int_t nruns = 0;
   Float_t  cutMaxAngle = 60;
   Float_t  cutAngle = 50;
   Float_t  cutEdep = 12;

   Float_t angleXmrad, angleYmrad;

   for (float spotX = -spotXlim; spotX <= spotXlim; spotX += spotStep) {
      kSpotX = spotX;
      angleXmrad = getAngleAtSpot(spotX) * 1000;
      for (float spotY = -spotYlim; spotY <= spotYlim; spotY += spotStep) {
         kSpotY = spotY;
   //   Float_t spotY = 0;
         angleYmrad = getAngleAtSpot(spotY) * 1000;

         printf("Spot (%.0f,%.0f)\n", spotX, spotY);
         Tracks * tracks = loadOrCreateTracks(true, Runs, false, run_energy, spotX, spotY);
         printf("lastJentry = %lld\n", lastJentry_);
//         tracks->removeHighAngleTracksRelativeToSpot(65, angleXmrad, angleYmrad);
//         tracks->removeTracksWithMinWEPL(100);
//         tracks->removeNuclearInteractions();
//         tracks->removeThreeSigmaShortTracks();

         for (Int_t i=0; i<=tracks->GetEntriesFast(); i++) {
            thisTrack = tracks->At(i);
            if (!thisTrack) continue;

            outSpotX = spotX;
            outSpotY = spotY;
            wepl = thisTrack->getFitParameterRange();
            // wepl_calibrated + 0.41 = 1.0024 * wepl;
            wepl_calibrated =  1.0024 * wepl - 0.41; // pol1 calibration for protons
            outWEPL = 333.7 - wepl_calibrated; // 333.7 mm: 230 MeV proton @ 78 eV H2O (extrapolated from data...)
            if (isnan(outWEPL)) continue;

   //            outWEPL = 238.8 - residualRange; // 190 MeV/u Helium
            if (!thisTrack->At(0) || !thisTrack->At(1)) continue;

            outX2x = thisTrack->At(0)->getXmm();
            outX2y = thisTrack->At(0)->getYmm();
            outP2x = (thisTrack->At(1)->getXmm() - outX2x) / dz2;
            outP2y = (thisTrack->At(1)->getYmm() - outX2y) / dz2;

            lastEdep = thisTrack->Last()->getDepositedEnergy();
            isSingleEventID = thisTrack->isFirstAndLastEventIDEqual();
            isLastHitSecondary = thisTrack->Last()->isSecondary();

            tOut->Fill();

            if (kDraw) {
               hSimpleImage->Fill(outX2x, outX2y, outWEPL);
               hSimpleImageNorm->Fill(outX2x, outX2y);
            }
         }
          nruns++;
          delete tracks;
      }
   }

   fOut->cd();
   fOut->Write();
   fOut->Close();
   printf("Image Reconstruction Input data written to Output/%s/%s_rotation%03d.root.\n", kPhantomName.Data(), kPhantomName.Data(), kRotation);
   
   if (kDraw) {
      hSimpleImage->Divide(hSimpleImageNorm);
      hSimpleImage->Draw("COLZ");
   }
}

void drawTracks3D(Int_t Runs, Int_t tracksperrun /* kEventsPerRun */, Float_t degraderThickness /* 160 */, Bool_t recreate /* 1 */, Bool_t doTracking /* 1 */) {
   Int_t dataType = 0;
   Int_t switchLayer = 100;
   Float_t energy = float(kEnergy);

   run_energy = energy;
   run_degraderThickness = degraderThickness;
   kEventsPerRun = tracksperrun;
   kDoTracking = doTracking;
   
   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy, 0, 0);
//   tracks->removeEmptyTracks();

   tracks->removeShortTracks(4);
//   tracks->sortTracksByLength();
//   tracks->fillOutIncompleteTracks(0.4);
   // Complete incomplete tracks

   printf("Found %d tracks before filtering.\n", tracks->GetEntries());

   Int_t numberOfPrimaries = 0;
   Int_t numberOfSecondaries = 0;
   for (int i=0; i<tracks->GetEntriesFast(); i++) {
      if (!tracks->At(i)) continue;
      if (!tracks->At(i)->Last()->isSecondary()) numberOfPrimaries++;
      else numberOfSecondaries++;
   }
   cout << "Number of primaries = " << numberOfPrimaries << ", number of secondaries = " << numberOfSecondaries << endl;
   
   Float_t  cutHalo = 15;
   Float_t  cutMaxAngle = 60;
   Float_t  cutAngle = 45;
   Float_t  cutEdep = 12;
   
   tracks->removeHighAngleTracks(cutAngle); // mrad
//   tracks->removeHighAngularChangeTracks(cutMaxAngle);
//   tracks->removeNuclearInteractions();
   tracks->removeThreeSigmaShortTracks();

   Bool_t   kDraw = false;

   TH1I  *hWEPLCorrect = new TH1I("hWEPLCorrect", ";Range in detector [mm WEPL];Frequency", 400, 0, 250);
   TH1I  *hWEPLSecondary = new TH1I("hWEPLSecondary", "", 400, 0, 250);
   TH1I  *hWEPLConfused = new TH1I("hWEPLConfused", "", 400, 0, 250);
   TH1I  *hProjConfused = new TH1I("hProjConfused", ";Projected error on phantom 10 cm from tracker;Frequency", 100, 0, 50);

   Float_t means[10];
   Float_t sigmas[10];

   TCanvas *c1 = nullptr;
   if (kDraw) {
      c1 = new TCanvas("c1", "c1", 1000, 800);
      c1->SetTitle(Form("Tracks from %.2f MeV protons on %s", energy, getMaterialChar()));
   }

   TView *view = nullptr; 
   if (kDraw) view = TView::CreateView(1);
   float fromx = 0.1 * nx;
   float tox = 0.9 * nx;
   float fromy = 0.1 * ny;
   float toy = 0.9 * ny;

   /*
   fromy = 0, toy = ny;
   fromx = 0, tox = nx;
   */

   Int_t zoom = 750; // 750

   fromx = nx/2 - zoom*2;
   fromy = ny/2 - zoom*2;
   tox = nx/2 + zoom*2;
   toy = ny/2 + zoom*2;

   Int_t iret;
   Float_t theta = 337;
   Float_t phi = 76;

   if (kDraw) {
      view->SetRange(fromx, 0, fromy, tox, 350, toy);
      view->SetView(theta, phi, 0, iret);
   }

//   TClonesArray *restPoints = tracks->getClustersWithoutTrack();
   Clusters * conflictClusters = nullptr;

   Int_t nClusters = 0;
   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      if (!tracks->At(i)) continue;
      nClusters += tracks->GetEntriesFast(i);
   }
   
   Int_t restPrimary = 0;
   Int_t restSecondary = 0;
   
//   for (Int_t i=0; i<restPoints->GetEntriesFast(); i++) {
   for (Int_t i=0; i<tracks->GetEntriesFastCWT(); i++) {
      if (!tracks->AtCWT(i)) 
         continue;

      Cluster *thisCluster = (Cluster*) tracks->AtCWT(i);

      if (thisCluster->isSecondary()) restSecondary++;
      else restPrimary++;
   }
   printf("Removed clusters: S = %d, P = %d", restSecondary, restPrimary);

   TPolyMarker3D *pMarker = new TPolyMarker3D(restPrimary, 7);
   TPolyMarker3D *EIDMarker = new TPolyMarker3D(restSecondary, 7);
   TPolyMarker3D *conflictMarker = new TPolyMarker3D(nClusters, 7);
   pMarker->SetMarkerColor(kBlue); // Missing cluster
   EIDMarker->SetMarkerColor(kGreen);
   conflictMarker->SetMarkerColor(kRed); // Conflicting cluster
  

   Int_t iPrim = 0, iSec = 0;

   for (Int_t i=0; i<tracks->GetEntriesFastCWT(); i++) {
      if (!tracks->AtCWT(i)) 
         continue;

      Cluster *thisCluster = (Cluster*) tracks->AtCWT(i);
      Float_t x = thisCluster->getX();
      Float_t z = thisCluster->getY();
      Float_t y = thisCluster->getLayermm();
   
      if (thisCluster->isSecondary() && thisCluster->getPDG() > 1000) {
         EIDMarker->SetPoint(iSec++,x,y,z);
      }
      else if (!thisCluster->isSecondary()) {
         pMarker->SetPoint(iPrim++, x, y, z);
      }
   }

   Int_t ntracks = tracks->GetEntriesFast();
   Int_t EIDidx = 0;
   Int_t conflictIdx = 0;

   Int_t medianEventID = -1;

   Int_t nTrueTracks = 0;
   Int_t nOKTracks = 0;
   Int_t nOKTracksAllClusters = 0;
   Int_t nOKTracksAllClustersOK2nd = 0;
   Int_t nOKMinusTracks = 0;
   Int_t nOKLastLayers = 0;
   Int_t nOneWrong = 0;
   Int_t nMissingEID;
   Int_t nOkPrimary = 0;
   Int_t nPrimaryIncomplete = 0;
   Int_t nPrimaryConfused = 0;
   Int_t nPrimaryConfusedIncomplete = 0;
   Int_t nOkSecondary = 0;
   Int_t nSecondaryIncomplete = 0;
   Int_t nSecondaryConfused = 0;
   Int_t nSecondaryConfusedIncomplete = 0;
   numberOfPrimaries = 0;
   numberOfSecondaries = 0;
   Float_t correctPosx;
   Float_t thisPosx;
   Float_t correctPosy;
   Float_t thisPosy;
   Float_t distToPhantom = 100;
   Cluster *a = nullptr;
   Cluster *b = nullptr;
   Cluster *aTrue = nullptr;
   Cluster *bTrue = nullptr;

   Track * thisTrack = nullptr;
//   tracks->createEIDSortList();
   
   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);
      if (!thisTrack) continue;
      if (isnan(thisTrack->getFitParameterRange())) continue;

      nMissingEID = tracks->getNMissingClustersWithEventID(thisTrack->getEventID(0), thisTrack->Last()->getLayer(), thisTrack->At(0)->getLayer()); 

      if (!thisTrack->Last()->isSecondary() && !thisTrack->isSecondary(0)) { // primary
         numberOfPrimaries++;

         if (thisTrack->isFirstAndLastEventIDEqual()) { // No confused tracks
            if (nMissingEID == 0) { //  && !tracks->isTrackIncomplete(thisTrack)) { // No missing clusters
               nOkPrimary++;
               hWEPLCorrect->Fill(getUnitFromTL(thisTrack->getFitParameterRange()));
            }
            
            else { // Missing clusters
               nPrimaryIncomplete++;
//               cout << "Incomplete track: " << *thisTrack << endl << endl;
               hWEPLConfused->Fill(getUnitFromTL(thisTrack->getFitParameterRange()));
            }

         }

         else { // Confused tracks
            a = thisTrack->At(0);
            b = thisTrack->At(1);
            int eid = thisTrack->Last()->getEventID();
            // Track *corTrack = tracks->getTrackWithEID(eid);
            Track * corTrack = nullptr;
            if (corTrack) {
               aTrue = corTrack->At(0);
               bTrue = corTrack->At(1);

               if (a && b && aTrue && bTrue) {
                  float deltax = b->getXmm() - a->getXmm();
                  float deltay = b->getYmm() - a->getYmm();
                  float deltaz = b->getLayermm() - a->getLayermm();
                  thisPosx = a->getXmm() - (deltax/deltaz) * distToPhantom;
                  thisPosy = a->getYmm() - (deltay/deltaz) * distToPhantom;

                  deltax = bTrue->getXmm() - aTrue->getXmm();
                  deltay = bTrue->getYmm() - aTrue->getYmm();
                  deltaz = bTrue->getLayermm() - aTrue->getLayermm();
                  correctPosx = aTrue->getXmm() - (deltax/deltaz) * distToPhantom;
                  correctPosy = aTrue->getYmm() - (deltay/deltaz) * distToPhantom;

                  float delta = sqrt(pow(thisPosx - correctPosx, 2) + pow(thisPosy - correctPosy, 2));
                  hProjConfused->Fill(delta);
               }
            }
            if (nMissingEID == 0) { // && !tracks->isTrackIncomplete(thisTrack)) { // No missing clusters
               nPrimaryConfused++;
            }

            else { // Confused and missing clusters
               nPrimaryConfusedIncomplete++;
            }
         }
      }

      else { // secondary
         numberOfSecondaries++;
         hWEPLSecondary->Fill(getUnitFromTL(thisTrack->getFitParameterRange()));

         if (thisTrack->isFirstAndLastEventIDEqual()) { // No confused tracks
            if (nMissingEID == 0) { // No missing clusters
               nOkSecondary++;
            }

            else { // Missing clusters
               nSecondaryIncomplete++;
            }
         }

         else { // Confused tracks
            if (nMissingEID == 0) { //  && !tracks->isTrackIncomplete(thisTrack)) { // No missing clusters
               nSecondaryConfused++;
            }
            else { // Confused AND Missing clusters
               nSecondaryConfusedIncomplete++;
            }
         }
      }

      if (thisTrack->isOneEventID()) nTrueTracks++;
      
      
      if (thisTrack->isFirstAndLastEventIDEqual()) nOKTracks++;
      else {
         if (!thisTrack->Last()) continue;
         Int_t lastEID = thisTrack->Last()->getEventID();   
         Bool_t lastSecondary = thisTrack->Last()->isSecondary();
         if (lastSecondary) continue;

         Int_t trackID = tracks->getTrackIdxFromFirstLayerEID(lastEID);

         if (trackID < 0) continue;
         if (!tracks->At(trackID)->At(0)) continue;
   
         Float_t delta = diffmmXY(tracks->At(trackID)->At(0), thisTrack->At(0));
         Float_t phi0 = thisTrack->getSlopeAngleBetweenLayers(1);
         Float_t phi1 = tracks->At(trackID)->getSlopeAngleBetweenLayers(1);

         Float_t deltaphi = fabs(phi0 - phi1);

         if (delta < 0.5 && deltaphi < 1) {
            nOKMinusTracks++;
            nOKLastLayers++;
         }
         else if (thisTrack->getWEPL() < 0.2 * getWEPLFromEnergy(run_energy)) {
            // Bad track ends early. OK...
            nOKLastLayers++;
         }
      }
   }

   nOKMinusTracks += nOKTracks;
   nOKLastLayers += nOKTracks;

   Int_t numberOfTracks = tracks->GetEntries();

   cout << endl << "Total number of tracks = " << numberOfPrimaries + numberOfSecondaries << endl;
   cout << "Number of primaries = " << numberOfPrimaries << ", number of secondaries = " << numberOfSecondaries << endl;
   cout << "Number of primaries tracked correctly = " << nOkPrimary << " (" << 100 * (float) nOkPrimary / numberOfPrimaries << "%)\n";
   cout << "Number of primaries confused = " << nPrimaryConfused << " (" << 100 * (float) nPrimaryConfused / numberOfPrimaries << "%)\n";
   cout << "Number of primaries incompletely tracked = " << nPrimaryIncomplete << " (" << 100 * (float) nPrimaryIncomplete / numberOfPrimaries << "%)\n";
   cout << "Number of primaries confused + incompletely tracked = " << nPrimaryConfusedIncomplete << " (" << 100 * (float) nPrimaryConfusedIncomplete / numberOfPrimaries << "%)\n";
   
   cout << endl << endl;
   cout << "Total number of primaries = " << numberOfPrimaries << endl;
   cout << "Primaries tracked correctly = " << nOkPrimary << " (" << 100 * (float) nOkPrimary / numberOfPrimaries;
   cout << "%),  Secondaries tracked correctly = " << nOkSecondary << " (" << 100 * (float) nOkSecondary / numberOfSecondaries << "%)\n";

   /*
   cout << endl << endl;
   cout << "Total number of tracks = " << numberOfPrimaries + numberOfSecondaries << endl;
   cout << "Tracks tracked correctly = " << nOkPrimary + nOkSecondary << " (" << 100 * (float) ( nOkPrimary + nOkSecondary) / (numberOfPrimaries + numberOfSecondaries) << "%), of which " << nOkPrimary << " (" << 100 * (float) nOkPrimary / (nOkPrimary + nOkSecondary);
   cout << "%) are primaries, and " << nOkSecondary << " (" << 100 * (float) nOkSecondary / (nOkSecondary + nOkPrimary) << "%) are secondaries.\n";
   cout << "Tracks confused = " << nPrimaryConfused + nSecondaryConfused << " (" << 100 * (float) (nPrimaryConfused + nSecondaryConfused) / (numberOfPrimaries + numberOfSecondaries) << "%), of which " << nPrimaryConfused << " (" << 100 * (float) nPrimaryConfused / (nPrimaryConfused + nSecondaryConfused);
   cout << "%) are primaries and " << nSecondaryConfused << " (" << 100 * (float) nSecondaryConfused / (nPrimaryConfused + nSecondaryConfused);
   cout << "%) are secondaries.\n";
   cout << "Tracks incompletely tracked = " << nPrimaryIncomplete + nSecondaryIncomplete << " (" << 100 * (float) (nPrimaryIncomplete + nSecondaryIncomplete) / (numberOfPrimaries + numberOfSecondaries) << "%), of which " << nPrimaryIncomplete << " (" << 100 * (float) nPrimaryIncomplete / (nPrimaryIncomplete + nSecondaryIncomplete);
   cout << "%) are primaries and " << nSecondaryIncomplete << " (" << 100 * (float) nSecondaryIncomplete/ (nPrimaryIncomplete+ nSecondaryIncomplete);
   cout << "%) are secondaries.\n";
   cout << "Tracks confused and incompletely tracked = " << nPrimaryConfusedIncomplete + nSecondaryConfusedIncomplete << " (" << 100 * (float) (nPrimaryConfusedIncomplete + nSecondaryConfusedIncomplete) / (numberOfPrimaries + numberOfSecondaries) << "%), of which " << nPrimaryConfusedIncomplete << " (" << 100 * (float) nPrimaryConfusedIncomplete / (nPrimaryConfusedIncomplete + nSecondaryConfusedIncomplete);
   cout << "%) are primaries and " << nSecondaryConfusedIncomplete << " (" << 100 * (float) nSecondaryConfusedIncomplete/ (nPrimaryConfusedIncomplete+ nSecondaryConfusedIncomplete);
   cout << "%) are secondaries.\n";

   Float_t factorEIDOK = 100 * ((float) nOKTracks / numberOfTracks);
   Float_t factorEIDOKAllClusters = 100 * ((float) nOKTracksAllClusters / numberOfTracks);
   Float_t factorEIDOKAllClustersOK2nd = 100 * ((float) nOkPrimary / numberOfTracks);
   Float_t factorEIDOKMinus = 100 * ((float) nOKMinusTracks / numberOfTracks);
   Float_t factorLastLayers = 100 * ((float) nOKLastLayers / numberOfTracks);

   cout << endl << endl;
   cout << nOKTracks << " of total " << numberOfTracks << " tracks has the same first/last ID (" << factorEIDOK << "%)\n";
   cout << nOKTracksAllClustersOK2nd << " of total " << numberOfTracks << " track has first/last event ID + no missing clusters (" << factorEIDOKAllClustersOK2nd << "%)\n";
   cout << nOKMinusTracks << " of total " << numberOfTracks << " tracks has a close match (0.5 mm, 1 degree) on first / last cluster (" << factorEIDOKMinus << "%)\n";
   cout << nOKLastLayers << " of total " << numberOfTracks << " tracks has a close match (0.5 mm, 1 degree) or is a very short track (" << factorLastLayers << "%)\n";
   */

   Int_t badSecondary = 0;
   Int_t badPrimary = 0;
   Int_t okPrimary = 0;
   Int_t badShort = 0;
   Int_t incompletePrimary = 0;

   if (kDraw) {
      for (Int_t i=0; i<ntracks; i++) {
         Track *thisTrack = tracks->At(i);
         if (!thisTrack) continue;


         if (thisTrack->Last()->getDepositedEnergy() > 5 || thisTrack->Last()->isSecondary()) continue;

         Int_t n = thisTrack->GetEntriesFast();

         TPolyLine3D *l = new TPolyLine3D(n);
         l->SetLineWidth(2);
         TPolyMarker3D *trackPoints = new TPolyMarker3D(nClusters, 7);
         nMissingEID = tracks->getNMissingClustersWithEventID(thisTrack->getEventID(0), thisTrack->Last()->getLayer(), thisTrack->At(0)->getLayer());
         if (!thisTrack->isFirstAndLastEventIDEqual()) {
            l->SetLineColor(kRed);
         }

         else if (nMissingEID>0) { // || tracks->isTrackIncomplete(thisTrack)) {
            l->SetLineColor(kGray);
         }
        
         if (thisTrack->isSecondary(0) || thisTrack->Last()->isSecondary()) {
            l->SetLineColor(kGreen);
         }

         Int_t lineElementNumber = 0;
         Int_t pointNumber = 0;
         for (Int_t j=0; j<n; j++) {
            if (!thisTrack->At(j)) continue;

            Float_t x = thisTrack->getX(j);
            Float_t z = thisTrack->getY(j);
            Float_t y = thisTrack->getLayermm(j);
            
            if (thisTrack->getLayer(j) < switchLayer) {
               l->SetPoint(lineElementNumber++,x,y,z);
            }
            else {
               trackPoints->SetPoint(pointNumber++, x, y, z);
            }
         }

         conflictClusters = (Clusters*) thisTrack->getConflictClusters();
         for (Int_t j=0; j<conflictClusters->GetEntriesFast(); j++) {
            if (!conflictClusters->At(j)) continue;
            
            Float_t x = conflictClusters->getX(j);
            Float_t z = conflictClusters->getY(j);
            Float_t y = conflictClusters->getLayer(j);
            
            conflictMarker->SetPoint(conflictIdx++, x,y,z);
         }
         if (kDraw) {
            l->SetLineWidth(3);
            if (l->GetLineColor() == kRed) l->Draw();
            if (l->GetLineColor() == kGreen) {
               l->SetLineWidth(2);
//               l->Draw();
            }
            if (l->GetLineColor() == kGray) l->Draw();
            if (l->GetLineColor() == kBlack) l->Draw();
         }

         if (l->GetLineColor() == kGreen) badSecondary++;
         if (l->GetLineColor() == kRed)   badPrimary++;
         else if (l->GetLineColor() == kGray) incompletePrimary++;
         else okPrimary++;
         
         if (kDraw) {
//            trackPoints->Draw();
            EIDMarker->Draw();
//            conflictMarker->Draw();
            pMarker->Draw();
         }
      }
      
      view->ShowAxis(); // comment for pure display
      c1->Update();

      TAxis3D *axis = TAxis3D::GetPadAxis();
      axis->SetTitle("3D view of tracks and clusters");
      axis->SetLabelColor(kBlack);
      axis->SetAxisColor(kBlack);
      axis->SetXTitle("Pixels in X");
      axis->SetYTitle("Layer number");
      axis->SetZTitle("Pixels in Y");
      axis->SetLabelSize(0.035);
      axis->SetTitleOffset(1.3);

   }

   TCanvas *c2 = new TCanvas();
   TH1I *changeEID = new TH1I("changeEID", "Change EID;Layer;EID changes", 100, 0, 100);
   for (int i=0; i<tracks->GetEntriesFast(); i++) {
      Track * thisTrack = tracks->At(i);
      if (!thisTrack) continue;
      int firstEID = thisTrack->getEventID(0);
      for (int j=0; j<thisTrack->GetEntriesFast(); j++) {
         if (thisTrack->getEventID(j) != firstEID) {
            changeEID->Fill(thisTrack->getLayer(j));
            break;
         }
      }
   }
   changeEID->Draw();

   Int_t badTotal = badPrimary + badSecondary;
   printf("Of %d bad tracks, %d (%.2f %%) are primaries and %d (%.2f %%) are secondaries.\n", badTotal, badPrimary, 100 * float(badPrimary) / badTotal, badSecondary, 100 * float(badSecondary) / badTotal);
   printf("okPrimary = %d\n", okPrimary);
   printf("Incomplete primary = %d\n", incompletePrimary);

   /*
   vector<Int_t> * conflictTracks = tracks->getTracksWithConflictClusters();
   vector<Int_t> * oneConflictPair = nullptr;
   vector<Int_t> * allConflictPairs = new vector<Int_t>;
   
   for (UInt_t i=0; i<conflictTracks->size(); i++) {
      oneConflictPair = tracks->getConflictingTracksFromTrack(conflictTracks->at(i));

      Int_t idx0 = oneConflictPair->at(0);
      Int_t idx1 = -1;
      if (oneConflictPair->size() > 1) {
         idx1 = oneConflictPair->at(1);
      }

      allConflictPairs->push_back(idx0);
      allConflictPairs->push_back(idx1);
   }

   // PRINTING
   
   cout << "Found the following tracks with conflicting clusters: ";
   for (UInt_t i=0; i<conflictTracks->size(); i++) {
      cout << conflictTracks->at(i) << " (eventID " << tracks->At(conflictTracks->at(i))->getEventID(0) << "), ";
   }
   cout << "\n";

   for (UInt_t i=0; i<allConflictPairs->size() / 2; i++) {
      if (allConflictPairs->at(2*i+1) < 0) continue;
      
      Track * trackA = tracks->At(allConflictPairs->at(2*i));
      Track * trackB = tracks->At(allConflictPairs->at(2*i+1));

      cout << "Track pair number " << i+1 << " found is: \n\tTRACK A: ";
      for (Int_t j=0; j<trackA->GetEntriesFast(); j++) { 
         if ( ! trackA->At(j) ) continue;
         cout << *trackA->At(j) << ", ";
      }

      cout << "\n\tTRACK B: ";
      for (Int_t j=0; j<trackB->GetEntriesFast(); j++) { 
         if ( ! trackB->At(j) ) continue;
         cout << *trackB->At(j) << ", ";
      }
      cout << endl;
   }
   */

   if (kDraw) c1->SaveAs(Form("OutputFiles/figures/testOutput_switchLayer%d.png", switchLayer));
/*
   TCanvas *c2 = new TCanvas();
   c2->Divide(2,1,1e-5,1e-5);
   c2->cd(1);
   hWEPLCorrect->SetLineColor(kBlack);
   hWEPLCorrect->Draw();
   hWEPLConfused->SetLineColor(kRed);
   hWEPLConfused->Draw("same");
   hWEPLSecondary->SetLineColor(kGreen);
   hWEPLSecondary->Draw("same");

   TLegend *l = new TLegend(0.3,0.3,0.4,0.4);
   l->AddEntry(hWEPLCorrect, "Good tracks", "L");
   l->AddEntry(hWEPLSecondary, "Nuclear interacting tracks");
   l->AddEntry(hWEPLConfused, "Tracks confused during recon");
   l->Draw();

   c2->cd(2);
   hProjConfused->SetLineColor(kRed);
   hProjConfused->Draw();
*/
   delete tracks;
//   delete conflictTracks;
//   delete allConflictPairs;
}

#endif
