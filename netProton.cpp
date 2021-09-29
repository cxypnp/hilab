#include <cmath>
#include <fstream>
#include <iostream>

#include "TH2F.h"
#include "TH1I.h"
#include "TH3F.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include "f14.h"
#include "TMultiGraph.h"
using namespace std;

void centralityBound(double *centraclityNchg, TH1D *refmult3)
{
   TH1 *cumul = refmult3->GetCumulative(false);
   cumul->Scale(1.0 / refmult3->Integral());
   for (int k = 0; k < refmult3->GetXaxis()->GetNbins(); k++)
   {
      for (int l = 0; l < 9; l++)
      {
         if (l == 0)
         {
            if (cumul->GetBinContent(k) > 0.05)
            {
               centraclityNchg[0] = cumul->GetBinLowEdge(k);
            }
         }
         else
         {
            if (cumul->GetBinContent(k) > l * 0.1)
            {
               centraclityNchg[l] = cumul->GetBinLowEdge(k);
            }
         }
      }
   }
}

double getPhiAbs(double x, double y)
{
   return atan2(abs(y), abs(x));
}

double getPhi(double x, double y)
{
   return atan(y / x);
}

void fillHistograms(TH2F *h2NcgNpp, TH2F **h2NcgNppPhi, bool reconstruction)
{
   TH1F *h1EventPlane = new TH1F("h1EventPlane", "h1EventPlane", 1000, -0.5 * TMath::Pi(), 0.5 * TMath::Pi());
   TTree *event = NULL;

   f14 reader(event);

   Long64_t nentries = reader.fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry = 0; jentry < nentries; jentry++)
   {
      if (jentry % 50000 == 0)
      {
         cout << "Progress:" << jentry * 1.0 / nentries * 100 << endl;
      }
      Long64_t ientry = reader.LoadTree(jentry);
      if (ientry < 0)
         break;
      int _nNetProton = 0;
      int _nCharged = 0;
      int _nTotalCharged = 0;
      double sumQx = 0, sumQy = 0;
      nb = reader.fChain->GetEntry(jentry);
      nbytes += nb;

      // if (Cut(ientry) < 0) continue;
      if (reader.fColHdr_timestep != 200 || reader.fColHdr_Ntrack == 0)
      {
         continue;
      }

      // Read a Event, Determine the Centrality and Event plane of this event
      for (Int_t j = 0; j < reader.fColHdr_Ntrack; j++)
      {
         // Charged Particle
         if ((reader.fTracksOut_chg[j] != 0))
         {
            float pT = TMath::Sqrt(reader.fTracksOut_px[j] * reader.fTracksOut_px[j] + reader.fTracksOut_py[j] * reader.fTracksOut_py[j]);
            // Is it (anti-)proton?
            if ((reader.fTracksOut_ityp[j] == 1) || (reader.fTracksOut_ityp[j] == -1))
            {
               // (anti-)proton

               float y = 0.5 * log((reader.fTracksOut_p0[j] + reader.fTracksOut_pz[j]) / (reader.fTracksOut_p0[j] - reader.fTracksOut_pz[j]));

               if (TMath::Abs(y) < 0.5 && pT < 2.0 && pT > 0.4)
               {
                  if (reader.fTracksOut_ityp[j] == 1)
                     _nNetProton++;
                  else if (reader.fTracksOut_ityp[j] == -1)
                     _nNetProton--;
               }
            }
            else
            {
               // Other Particle
               float p = TMath::Sqrt(reader.fTracksOut_px[j] * reader.fTracksOut_px[j] + reader.fTracksOut_py[j] * reader.fTracksOut_py[j] + reader.fTracksOut_pz[j] * reader.fTracksOut_pz[j]);
               float eta = 0.5 * log((p + reader.fTracksOut_pz[j]) / (p - reader.fTracksOut_pz[j]));
               if (TMath::Abs(eta) < 1 && pT > 0.1)
               {
                  _nCharged++;
               }
            }
            sumQx += pT * cos(2 * getPhi(reader.fTracksOut_rx[j], reader.fTracksOut_ry[j]));
            sumQy += pT * sin(2 * getPhi(reader.fTracksOut_rx[j], reader.fTracksOut_ry[j]));
            _nTotalCharged++;
         }
      }
      double Psi = 1.0 / 2.0 * atan(sumQy / sumQx);
      h1EventPlane->Fill(Psi);
      h2NcgNpp->Fill(_nNetProton, _nCharged);
      // Have known the Ncg and Npp, fill the h2 with the phi of every particle
      int _nNetProtonPhi[6] = {0, 0, 0, 0, 0, 0};
      for (Int_t j = 0; j < reader.fColHdr_Ntrack; j++)
      {
         // Charged Particle
         if ((reader.fTracksOut_chg[j] != 0))
         {
            // Is it (anti-)proton?
            if ((reader.fTracksOut_ityp[j] == 1) || (reader.fTracksOut_ityp[j] == -1))
            {
               // (anti-)proton
               float pT = TMath::Sqrt(reader.fTracksOut_px[j] * reader.fTracksOut_px[j] + reader.fTracksOut_py[j] * reader.fTracksOut_py[j]);
               float y = 0.5 * log((reader.fTracksOut_p0[j] + reader.fTracksOut_pz[j]) / (reader.fTracksOut_p0[j] - reader.fTracksOut_pz[j]));
               if (TMath::Abs(y) < 0.5 && pT < 2.0 && pT > 0.4)
               {
                  int index = 0;
                  if (reconstruction)
                  {
                     index = getPhiAbs(reader.fTracksOut_rx[j] * cos(Psi) + reader.fTracksOut_ry[j] * sin(Psi), reader.fTracksOut_ry[j] * cos(Psi) - reader.fTracksOut_rx[j] * sin(Psi)) / (TMath::Pi() * 1.0 / 12);
                  }
                  else
                  {
                     index = getPhiAbs(reader.fTracksOut_rx[j], reader.fTracksOut_ry[j]) / (TMath::Pi() * 1.0 / 12);
                  }
                  if (reader.fTracksOut_ityp[j] == 1)
                     _nNetProtonPhi[index]++;
                  else if (reader.fTracksOut_ityp[j] == -1)
                     _nNetProtonPhi[index]--;
               }
            }
            else
            {
               //          // // Other Particle
            }
         }
      }
      for (int j = 0; j < 6; j++)
      {
         h2NcgNppPhi[j]->Fill(_nNetProtonPhi[j], _nCharged);
      }
   }
   TCanvas *cTemp = new TCanvas();
   h1EventPlane->Draw();
   cTemp->SaveAs("h1EventPlane.png");
}

TH3I *drawEvent()
{
   TH3I *particles = new TH3I("particles", "particles", 100, -50, 50, 100, -50, 50, 100, -50, 50);
   TTree *event = NULL;

   f14 reader(event);

   Long64_t nentries = reader.fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry = 0; jentry < 1000; jentry++)
   {
      if (jentry % 10000 == 0)
      {
         cout << jentry * 1.0 / nentries * 100 << endl;
      }
      Long64_t ientry = reader.LoadTree(jentry);
      if (ientry < 0)
         break;
      nb = reader.fChain->GetEntry(jentry);
      nbytes += nb;
      for (Int_t j = 0; j < reader.fColHdr_Ntrack; j++)
      {
         particles->Fill(reader.fTracksOut_rx[j], reader.fTracksOut_ry[j], reader.fTracksOut_rz[j]);
      }
   }
   return particles;
}

void getCumulants(TH2F *h2NcgNpp, double *C1, double *C2, double *C3, double *C4, double *centraclityNchg)
{
   int binMax, binMin;
   int sum;
   double mean[9] = {0}, sigma[9] = {0}, fS[9] = {0}, kappa[9] = {0};
   // For centain phi, calculate Ci
   // CBWM included
   for (int l = 0; l < 9; l++)
   {
      sum = 0;
      if (l == 0)
      {
         binMin = h2NcgNpp->ProjectionY()->GetBin(centraclityNchg[l]);
         binMax = h2NcgNpp->ProjectionY()->GetMinimumBin();
      }
      else
      {
         binMin = h2NcgNpp->ProjectionY()->GetBin(centraclityNchg[l]);
         binMax = h2NcgNpp->ProjectionY()->GetBin(centraclityNchg[l - 1]);
      }
      // Loop on every centrality bin and do the modification
      for (int j = binMin; j < binMax; j++)
      {
         TH1D *_netProtonDist = h2NcgNpp->ProjectionX("", j, j);
         int nEvents = _netProtonDist->Integral();

         if (nEvents == 0 or nEvents == 1)
            continue;

         float _mean = 0, _sigma = 0, _fS = 0, _kappa = 0;
         _mean = _netProtonDist->GetMean();
         _sigma = _netProtonDist->GetStdDev();
         _fS = _netProtonDist->GetSkewness();
         if (isnan(_fS))
         {
            _fS = 0;
            cout << "omitEvent: " << nEvents << endl;
         }
         _kappa = _netProtonDist->GetKurtosis();
         if (isnan(_kappa))
         {
            _kappa = 0;
            cout << "omitEvent: " << nEvents << endl;
         }

         mean[l] += nEvents * _mean;
         sigma[l] += nEvents * _sigma;
         fS[l] += nEvents * _fS;
         kappa[l] += nEvents * _kappa;
         sum += nEvents;
      }

      mean[l] = mean[l] / sum;
      sigma[l] = sigma[l] / sum;
      fS[l] = fS[l] / sum;
      kappa[l] = kappa[l] / sum;

      C1[l] = mean[l];
      C2[l] = sigma[l] * sigma[l];
      C3[l] = fS[l] * C2[l] * sigma[l];
      C4[l] = kappa[l] * sigma[l] * sigma[l] * sigma[l] * sigma[l];
   }
}

void netProton(TString outputplotsfolder = "outputFolder/")
{
   system("rm -rf " + outputplotsfolder);
   system("mkdir -p " + outputplotsfolder);
   for (int i = 0; i < 9; i++)
   {
      system(Form("mkdir -p " + outputplotsfolder + "/%d", i));
      system(Form("mkdir -p " + outputplotsfolder + "/%d", i));
      system(Form("mkdir -p " + outputplotsfolder + "/%d", i));
      system(Form("mkdir -p " + outputplotsfolder + "/%d", i));
   }

   TCanvas *c1 = new TCanvas();
   c1->cd();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0101);

   TFile f("outputFile.root");
   TH2F *h2NcgNpp = NULL;
   TH2F *h2NcgNppPhi[6];
   f.GetObject("h2NcgNpp", h2NcgNpp);
   for (int i = 0; i < 6; i++)
   {
      f.GetObject(Form("h2NcgNppPhi%d", i), h2NcgNppPhi[i]);
   }

   if (!h2NcgNpp)
   {
      h2NcgNpp = new TH2F("h2NcgNpp", "h2NcgNpp", 100, -5.5, 94.5, 450, -0.5, 449.5);
      h2NcgNpp->GetYaxis()->SetTitle("Charged Multiplicity");
      h2NcgNpp->GetXaxis()->SetTitle("Net-Proton");
      for (int i = 0; i < 6; i++)
      {
         h2NcgNppPhi[i] = new TH2F(Form("h2NcgNppPhi%d", i), Form("h2NcgNppPhi%d", i), 100, -5.5, 94.5, 450, -0.5, 449.5);
         h2NcgNppPhi[i]->GetYaxis()->SetTitle("Charged Multiplicity");
         h2NcgNppPhi[i]->GetXaxis()->SetTitle("Net-Proton");
         h2NcgNppPhi[i]->GetZaxis()->SetTitle("#Phi");
      }
      fillHistograms(h2NcgNpp, h2NcgNppPhi, true);
   }

   double centraclityNchg[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
   centralityBound(centraclityNchg, h2NcgNpp->ProjectionY());
   TFile *hfile = TFile::Open("outputFile.root", "RECREATE");
   for (int i = 0; i < 6; i++)
   {
      h2NcgNppPhi[i]->Write();
      h2NcgNppPhi[i]->Draw();
      c1->SaveAs(Form(outputplotsfolder + "h2NcgNppPhi%d.png", i));
   }

   // ofstream logOut;
   // logOut.open(outputplotsfolder + "run.log", ios_base::app);
   // for (int l = 0; l < 9; l++)
   // {
   //    if (l == 0)
   //    {
   //       logOut << 5 << "%-10%:" << centraclityNchg[l] << endl;
   //    }
   //    else
   //    {
   //       logOut << l * 10 << "%" << centraclityNchg[l] << endl;
   //    }
   // }

   // Loop on every angle

   h2NcgNpp->Write();

   double x1Bin = 0.5 * TMath::Pi() / 2;
   double x1err = 0.5 * TMath::Pi() / 2;
   double x2Bin[2] = {0.125 * TMath::Pi(), 0.375 * TMath::Pi()};
   double x2err[2] = {0.125 * TMath::Pi(), 0.125 * TMath::Pi()};
   double x3Bin[3] = {1.0 * TMath::Pi() / 12, 3.0 * TMath::Pi() / 12, 5.0 * TMath::Pi() / 12};
   double x3err[3] = {1.0 * TMath::Pi() / 12, 1.0 * TMath::Pi() / 12, 1.0 * TMath::Pi() / 12};
   double x6Bin[6] = {1.0 * TMath::Pi() / 24, 3.0 * TMath::Pi() / 24, 5.0 * TMath::Pi() / 24, 7.0 * TMath::Pi() / 24, 9.0 * TMath::Pi() / 24, 11.0 * TMath::Pi() / 24};
   double x6err[6] = {1.0 * TMath::Pi() / 24, 1.0 * TMath::Pi() / 24, 1.0 * TMath::Pi() / 24, 1.0 * TMath::Pi() / 24, 1.0 * TMath::Pi() / 24, 1.0 * TMath::Pi() / 24};

   double C2oC1VsPhi6Bin[9][6];
   double C3oC2VsPhi6Bin[9][6];
   double C4oC2VsPhi6Bin[9][6];

   double C2oC1VsPhi3Bin[9][3];
   double C3oC2VsPhi3Bin[9][3];
   double C4oC2VsPhi3Bin[9][3];

   double C2oC1VsPhi2Bin[9][2];
   double C3oC2VsPhi2Bin[9][2];
   double C4oC2VsPhi2Bin[9][2];

   double C2oC1VsPhi1Bin[9];
   double C3oC2VsPhi1Bin[9];
   double C4oC2VsPhi1Bin[9];
   // 6 Bins
   for (int i = 0; i < 6; i++)
   {
      double C1[9] = {0}, C2[9] = {0}, C3[9] = {0}, C4[9] = {0};
      double C2oC1[9] = {0}, C3oC2[9] = {0}, C4oC2[9] = {0};
      getCumulants(h2NcgNppPhi[i], C1, C2, C3, C4, centraclityNchg);
      for (int l = 0; l < 9; l++)
      {
         C2oC1[l] = C2[l] / C1[l];
         C3oC2[l] = C3[l] / C2[l];
         C4oC2[l] = C4[l] / C2[l];
      }
      for (int j = 0; j < 9; j++)
      {
         C2oC1VsPhi6Bin[j][i] = C2oC1[j];
         C3oC2VsPhi6Bin[j][i] = C3oC2[j];
         C4oC2VsPhi6Bin[j][i] = C4oC2[j];
      }
   }

   for (int i = 0; i < 6;)
   {
      double C1[9] = {0}, C2[9] = {0}, C3[9] = {0}, C4[9] = {0};
      double C2oC1[9] = {0}, C3oC2[9] = {0}, C4oC2[9] = {0};
      TH2F temp = (*h2NcgNppPhi[i] + *h2NcgNppPhi[i + 1]);
      getCumulants(&temp, C1, C2, C3, C4, centraclityNchg);
      for (int l = 0; l < 9; l++)
      {
         C2oC1[l] = C2[l] / C1[l];
         C3oC2[l] = C3[l] / C2[l];
         C4oC2[l] = C4[l] / C2[l];
      }
      for (int j = 0; j < 9; j++)
      {
         C2oC1VsPhi3Bin[j][i / 2] = C2oC1[j];
         C3oC2VsPhi3Bin[j][i / 2] = C3oC2[j];
         C4oC2VsPhi3Bin[j][i / 2] = C4oC2[j];
      }
      i += 2;
   }

   for (int i = 0; i < 6;)
   {
      double C1[9] = {0}, C2[9] = {0}, C3[9] = {0}, C4[9] = {0};
      double C2oC1[9] = {0}, C3oC2[9] = {0}, C4oC2[9] = {0};
      TH2F temp = (*h2NcgNppPhi[i] + *h2NcgNppPhi[i + 1]);
      temp = temp + *h2NcgNppPhi[i + 2];
      getCumulants(&temp, C1, C2, C3, C4, centraclityNchg);
      for (int l = 0; l < 9; l++)
      {
         C2oC1[l] = C2[l] / C1[l];
         C3oC2[l] = C3[l] / C2[l];
         C4oC2[l] = C4[l] / C2[l];
      }
      for (int j = 0; j < 9; j++)
      {
         C2oC1VsPhi2Bin[j][i / 3] = C2oC1[j];
         C3oC2VsPhi2Bin[j][i / 3] = C3oC2[j];
         C4oC2VsPhi2Bin[j][i / 3] = C4oC2[j];
      }
      i += 3;
   }

   TH2F temp = *h2NcgNppPhi[0];
   for (int i = 1; i < 6; i++)
   {
      temp = temp + *h2NcgNppPhi[i];
   }
   temp.Write();
   double C1[9] = {0}, C2[9] = {0}, C3[9] = {0}, C4[9] = {0};
   double C2oC1[9] = {0}, C3oC2[9] = {0}, C4oC2[9] = {0};
   getCumulants(&temp, C1, C2, C3, C4, centraclityNchg);
   for (int l = 0; l < 9; l++)
   {
      C2oC1[l] = C2[l] / C1[l];
      C3oC2[l] = C3[l] / C2[l];
      C4oC2[l] = C4[l] / C2[l];
   }
   for (int j = 0; j < 9; j++)
   {
      C2oC1VsPhi1Bin[j] = C2oC1[j];
      C3oC2VsPhi1Bin[j] = C3oC2[j];
      C4oC2VsPhi1Bin[j] = C4oC2[j];
      cout << C4oC2[j] << endl;
   }

   for (int i = 0; i < 9; i++)
   {

      c1->Clear();
      TGraphErrors *h1C2oC1VsPhi1Bin = new TGraphErrors(1, &x1Bin, &C2oC1VsPhi1Bin[i], &x1err, NULL);
      h1C2oC1VsPhi1Bin->SetLineColor(2);
      h1C2oC1VsPhi1Bin->SetMarkerSize(1);
      h1C2oC1VsPhi1Bin->SetMarkerStyle(20);
      TGraphErrors *h1C2oC1VsPhi2Bin = new TGraphErrors(2, x2Bin, C2oC1VsPhi2Bin[i], x2err, NULL);
      h1C2oC1VsPhi2Bin->SetLineColor(3);
      h1C2oC1VsPhi2Bin->SetMarkerSize(1);
      h1C2oC1VsPhi2Bin->SetMarkerStyle(21);
      TGraphErrors *h1C2oC1VsPhi3Bin = new TGraphErrors(3, x3Bin, C2oC1VsPhi3Bin[i], x3err, NULL);
      h1C2oC1VsPhi3Bin->SetLineColor(4);
      h1C2oC1VsPhi3Bin->SetMarkerSize(1);
      h1C2oC1VsPhi3Bin->SetMarkerStyle(22);
      TGraphErrors *h1C2oC1VsPhi6Bin = new TGraphErrors(6, x6Bin, C2oC1VsPhi6Bin[i], x6err, NULL);
      h1C2oC1VsPhi6Bin->SetLineColor(7);
      h1C2oC1VsPhi6Bin->SetMarkerSize(1);
      h1C2oC1VsPhi6Bin->SetMarkerStyle(23);
      TMultiGraph *mg = new TMultiGraph();
      mg->Add(h1C2oC1VsPhi1Bin);
      mg->Add(h1C2oC1VsPhi2Bin);
      mg->Add(h1C2oC1VsPhi3Bin);
      mg->Add(h1C2oC1VsPhi6Bin);
      mg->Draw("zpA");
      c1->SaveAs(Form(outputplotsfolder + "%d/h1C2oC1VsPhiCentral.png", i));
      c1->SetTitle(Form("h1C2oC1VsPhiCentral%d", i));

      TGraphErrors *h1C3oC2VsPhi1Bin = new TGraphErrors(1, &x1Bin, &C3oC2VsPhi1Bin[i], &x1err, NULL);
      h1C3oC2VsPhi1Bin->SetLineColor(2);
      h1C3oC2VsPhi1Bin->SetMarkerSize(1);
      h1C2oC1VsPhi1Bin->SetMarkerStyle(20);
      TGraphErrors *h1C3oC2VsPhi2Bin = new TGraphErrors(2, x2Bin, C3oC2VsPhi2Bin[i], x2err, NULL);
      h1C3oC2VsPhi2Bin->SetLineColor(3);
      h1C3oC2VsPhi2Bin->SetMarkerSize(1);
      h1C3oC2VsPhi2Bin->SetMarkerStyle(21);
      TGraphErrors *h1C3oC2VsPhi3Bin = new TGraphErrors(3, x3Bin, C3oC2VsPhi3Bin[i], x3err, NULL);
      h1C3oC2VsPhi3Bin->SetLineColor(4);
      h1C3oC2VsPhi3Bin->SetMarkerSize(1);
      h1C3oC2VsPhi3Bin->SetMarkerStyle(22);
      TGraphErrors *h1C3oC2VsPhi6Bin = new TGraphErrors(6, x6Bin, C3oC2VsPhi6Bin[i], x6err, NULL);
      h1C3oC2VsPhi6Bin->SetLineColor(7);
      h1C3oC2VsPhi6Bin->SetMarkerSize(1);
      h1C3oC2VsPhi6Bin->SetMarkerStyle(23);
      mg = new TMultiGraph();
      mg->Add(h1C3oC2VsPhi1Bin);
      mg->Add(h1C3oC2VsPhi2Bin);
      mg->Add(h1C3oC2VsPhi3Bin);
      mg->Add(h1C3oC2VsPhi6Bin);
      mg->Draw("zpA");
      c1->Size(0, 0);
      c1->SaveAs(Form(outputplotsfolder + "%d/h1C3oC2VsPhiCentral.png", i));
      c1->SetTitle(Form("h1vVsPhiCentral%d", i));

      c1->Clear();
      TGraphErrors *h1C4oC2VsPhi1Bin = new TGraphErrors(1, &x1Bin, &C4oC2VsPhi1Bin[i], &x1err, NULL);
      h1C4oC2VsPhi1Bin->SetLineColor(2);
      h1C4oC2VsPhi1Bin->SetMarkerSize(1);
      h1C4oC2VsPhi1Bin->SetMarkerStyle(20);
      TGraphErrors *h1C4oC2VsPhi2Bin = new TGraphErrors(2, x2Bin, C4oC2VsPhi2Bin[i], x2err, NULL);
      h1C4oC2VsPhi2Bin->SetLineColor(3);
      h1C4oC2VsPhi2Bin->SetMarkerSize(1);
      h1C4oC2VsPhi2Bin->SetMarkerStyle(21);
      TGraphErrors *h1C4oC2VsPhi3Bin = new TGraphErrors(3, x3Bin, C4oC2VsPhi3Bin[i], x3err, NULL);
      h1C4oC2VsPhi3Bin->SetLineColor(4);
      h1C4oC2VsPhi3Bin->SetMarkerSize(1);
      h1C4oC2VsPhi3Bin->SetMarkerStyle(22);
      TGraphErrors *h1C4oC2VsPhi6Bin = new TGraphErrors(6, x6Bin, C4oC2VsPhi6Bin[i], x6err, NULL);
      h1C4oC2VsPhi6Bin->SetLineColor(7);
      h1C4oC2VsPhi6Bin->SetMarkerSize(1);
      h1C3oC2VsPhi6Bin->SetMarkerStyle(23);
      mg = new TMultiGraph();
      mg->Add(h1C4oC2VsPhi1Bin);
      mg->Add(h1C4oC2VsPhi2Bin);
      mg->Add(h1C4oC2VsPhi3Bin);
      mg->Add(h1C4oC2VsPhi6Bin);
      mg->Draw("zpA");
      c1->Size(0, 0);
      c1->SaveAs(Form(outputplotsfolder + "%d/h1C4oC2VsPhiCentral.png", i));
      c1->SetTitle(Form("h1vVsPhiCentral%d", i));
   }

   c1->SetLogy(1);
   h2NcgNpp->ProjectionX()->DrawNormalized("E");
   h2NcgNpp->SetTitle("Net-proton production vs. N_{chg}");
   c1->SaveAs(outputplotsfolder + "netPro.pdf");
   c1->SaveAs(outputplotsfolder + "netPro.png");

   h2NcgNpp->ProjectionY()->DrawNormalized("E");
   c1->SaveAs(outputplotsfolder + "refmult3.pdf");
   c1->SaveAs(outputplotsfolder + "refmult3.png");

   TH1 *cumul = h2NcgNpp->ProjectionY()->GetCumulative(false);
   cumul->SetTitle("Cumulative Charged Particle Multiplicity");
   cumul->Scale(1.0 / h2NcgNpp->ProjectionY()->Integral());
   cumul->Draw("E");
   c1->SaveAs(outputplotsfolder + "refmult3Cumu.pdf");
   c1->SaveAs(outputplotsfolder + "refmult3Cumu.png");

   c1->SetLogz(1);
   c1->SetLogy(0);
   h2NcgNpp->Draw("COLZ");
   c1->SaveAs(outputplotsfolder + "h2NcgNpp.pdf");
   c1->SaveAs(outputplotsfolder + "h2NcgNpp.png");

   int xbins = h2NcgNpp->GetXaxis()->GetNbins();
   int ybins = h2NcgNpp->GetYaxis()->GetNbins();

   // drawEvent()->Write();
   hfile->Write();
   hfile->Close();
}
