#include <TSystem.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>
#include <TNtuple.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TText.h>
#include <iostream>
#include <fstream>
#include <TLatex.h>
#include <string.h>
#include <TString.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSpline.h>
#include <THStack.h>
#include <TRandom3.h>
#include <TPad.h>
#include <TMath.h>

using namespace std;

Bool_t kUseCarbon = true;
const Int_t arraySize = 1500;
const Int_t xFrom = 40;

void plotEnergyVsRange() {
   TCanvas *c1 = new TCanvas("c1", "Range accuracy", 1200, 800);
   
   TPaveLabel *Xtitle = new TPaveLabel(0.01, 0.1, 0.03, 0.9, "Range deviation [mm WEPL]");
   Xtitle->SetBorderSize(0);
   Xtitle->SetFillColor(kWhite);
   Xtitle->SetTextAngle(90);
   Xtitle->SetTextFont(22);
   Xtitle->SetTextSize(0.05);
   Xtitle->Draw();

   TPad *graphPad = new TPad("Graphs", "Graphs", 0.05, 0.05, 0.95, 0.95);
   graphPad->Draw();
   graphPad->cd();
   graphPad->Divide(1,5,0,0);

   /*
   Double_t W = 0.3;
   Int_t Ny = 3;
   Double_t Ym = (1-(Ny*W))/2;
   Double_t dw = (W*0.1)/4;

   TPad *p1 = new TPad("p1", "p1", 0.2, Ym, 0.8, Ym+W+dw, 0, 0, 0);
   p1->SetBottomMargin(0);
   p1->Draw();

   TPad *p2 = new TPad("p1", "p1", 0.2, Ym+W+dw, 0.8, Ym+2*W-dw, 0, 0, 0);
   p2->SetTopMargin(0);
   p2->SetBottomMargin(0);
   p2->Draw();
   
   TPad *p3 = new TPad("p1", "p1", 0.2, Ym+2*W-dw, 0.8, Ym+3*W, 0, 0, 0);
   p3->SetTopMargin(0);
   p3->Draw();
*/

   Float_t  arrayE2[arraySize] = {0}; // energy MC
   Float_t  arrayE3[arraySize] = {0}; // energy MC
   Float_t  arrayE4[arraySize] = {0}; // energy MC
   Float_t  arrayE5[arraySize] = {0}; // energy MC
   Float_t  arrayE6[arraySize] = {0}; // energy MC
   Float_t  arrayMC2[arraySize] = {0}; // range MC
   Float_t  arrayMC3[arraySize] = {0}; // range MC
   Float_t  arrayMC4[arraySize] = {0}; // range MC
   Float_t  arrayMC5[arraySize] = {0}; // range MC
   Float_t  arrayMC6[arraySize] = {0}; // range MC

   Double_t energies[arraySize] = {0};
   Double_t thicknesses[arraySize] = {0};

   gStyle->SetOptStat(0);

   Float_t nomrange_, estrange_, sigmaRange_, lastRange_, nomsigma_, waterphantomthickness_, dummy0;
   Int_t energy_, thickness_;
   Float_t estimatedStraggling;

   ifstream in1;
   in1.open("../../Data/Ranges/EnergyAfterDegrader.csv");
   Int_t thick, n=0;
   Double_t energy;
   while (1) {
      in1 >> thick >> energy >> dummy0;
      if (!in1.good()) break;   
      thicknesses[n] = thick;
      energies[n++] = energy;
   }
   in1.close();

   TSpline3 *energySpline = new TSpline3("energySpline", thicknesses, energies, n);

   ifstream in;
   if (!kUseCarbon) {
      in.open("../../OutputFiles/result_makebraggpeakfit.csv");
   }
   else {
      in.open("../../OutputFiles/result_makebraggpeakfitCarbon.csv");
   }

   Int_t nlines6 = 0, nlines2 = 0, nlines3 = 0, nlines4 = 0, nlines5 = 0;
   
   while (1) {
      in >> thickness_ >> energy_ >> nomrange_ >> estrange_ >> nomsigma_ >> sigmaRange_;

      if (!in.good()) {
         break;
      }

      energy = energySpline->Eval(energy_);

      if (thickness_ == 2) arrayE2[nlines2] = energy;
      if (thickness_ == 3) arrayE3[nlines3] = energy;
      if (thickness_ == 4) arrayE4[nlines4] = energy;
      if (thickness_ == 5) arrayE5[nlines5] = energy;
      if (thickness_ == 6) arrayE6[nlines6] = energy;

      if (thickness_ == 2) arrayMC2[nlines2++] = -nomrange_ + estrange_;
      if (thickness_ == 3) arrayMC3[nlines3++] = estrange_ - nomrange_;
      if (thickness_ == 4) arrayMC4[nlines4++] = -nomrange_ + estrange_;
      if (thickness_ == 5) arrayMC5[nlines5++] = estrange_ - nomrange_;
      if (thickness_ == 6) arrayMC6[nlines6++] = -nomrange_ + estrange_;
   }
   
   in.close();

   TGraph *hMC2 = new TGraph(nlines2, arrayE2, arrayMC2);
   TGraph *hMC3 = new TGraph(nlines3, arrayE3, arrayMC3);
   TGraph *hMC4 = new TGraph(nlines4, arrayE4, arrayMC4);
   TGraph *hMC5 = new TGraph(nlines5, arrayE5, arrayMC5);
   TGraph *hMC6 = new TGraph(nlines6, arrayE6, arrayMC6);

   hMC2->SetTitle(";Energy [MeV];");
   hMC3->SetTitle(";Energy [MeV];");
   hMC4->SetTitle(";Energy [MeV];");
   hMC5->SetTitle(";Energy [MeV];");
   hMC6->SetTitle(";Energy [MeV];");

   Float_t ysize = 0.11;

   hMC6->GetXaxis()->SetTitleFont(22);
   hMC6->GetXaxis()->SetTitleSize(0.11);
   hMC6->GetXaxis()->SetLabelFont(22);
   hMC6->GetXaxis()->SetLabelSize(0.11);
   hMC4->GetYaxis()->SetLabelSize(ysize);
   hMC4->GetYaxis()->SetLabelFont(22);
   hMC2->GetYaxis()->SetLabelSize(ysize);
   hMC2->GetYaxis()->SetLabelFont(22);
   hMC6->GetYaxis()->SetLabelSize(ysize);
   hMC6->GetYaxis()->SetLabelFont(22);
   hMC3->GetYaxis()->SetLabelSize(ysize);
   hMC3->GetYaxis()->SetLabelFont(22);
   hMC5->GetYaxis()->SetLabelSize(ysize);
   hMC5->GetYaxis()->SetLabelFont(22);

   hMC2->SetLineColor(kRed+4);
   hMC3->SetLineColor(kRed+3);
   hMC4->SetLineColor(kRed+2);
   hMC5->SetLineColor(kRed+1);
   hMC6->SetLineColor(kRed);
   hMC2->SetLineWidth(3);
   hMC3->SetLineWidth(3);
   hMC4->SetLineWidth(3);
   hMC5->SetLineWidth(3);
   hMC6->SetLineWidth(3);

   Float_t yfrom = -3.5;
   Float_t yto = 0.5;

   hMC2->GetXaxis()->SetRangeUser(60, 240);
   hMC3->GetXaxis()->SetRangeUser(60, 240);
   hMC4->GetXaxis()->SetRangeUser(60, 240);
   hMC5->GetXaxis()->SetRangeUser(60, 240);
   hMC6->GetXaxis()->SetRangeUser(60, 240);
   hMC2->GetYaxis()->SetRangeUser(yfrom, yto);
   hMC3->GetYaxis()->SetRangeUser(yfrom, yto);
   hMC4->GetYaxis()->SetRangeUser(yfrom, yto);
   hMC5->GetYaxis()->SetRangeUser(yfrom, yto);
   hMC6->GetYaxis()->SetRangeUser(yfrom, yto);
   hMC2->GetYaxis()->SetNdivisions(404);
   hMC3->GetYaxis()->SetNdivisions(404);
   hMC4->GetYaxis()->SetNdivisions(404);
   hMC5->GetYaxis()->SetNdivisions(404);
   hMC6->GetYaxis()->SetNdivisions(404);

   graphPad->cd(1);
   graphPad->SetGridy();
   hMC2->Draw("LA");
   TText *t2 = new TText();
   t2->SetTextSize(0.14);
   t2->SetTextFont(22);
   t2->DrawText(100, -0.5, "2 mm Al absorber"); 

   graphPad->cd(2);
   hMC3->Draw("LA");
   t2->DrawText(100, -0.5, "3 mm Al absorber");

   graphPad->cd(3);
   hMC4->Draw("LA");
   t2->DrawText(100, -0.5, "4 mm Al absorber");
   
   graphPad->cd(4);
   hMC5->Draw("LA");
   t2->DrawText(100, -0.5, "5 mm Al absorber");
   
   graphPad->cd(5);
   hMC6->Draw("LA");
   t2->DrawText(100, -0.5, "6 mm Al absorber");


}