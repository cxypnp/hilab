#include <cmath>
#include <fstream>
#include <iostream>

#include "TH2F.h"
#include "TH3F.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include "f14.h"
#include "TMultiGraph.h"
#include "TLorentzVector.h"
#include "TRandom.h"

using namespace std;

void centralityBound(double *centralityNchg, TH1D *refmult3)
{
   // This function is for determine the centrality bins
   TH1 *cumul = refmult3->GetCumulative();
   cumul->Scale(1.0 / refmult3->Integral());
   for (int k = 0; k < refmult3->GetXaxis()->GetNbins(); k++)
   {
      for (int l = 0; l < 9; l++)
      {
         if (l == 0)
         {
            // 0-5% bin
            if (cumul->GetBinContent(k) <= 0.05)
            {
               centralityNchg[0] = cumul->GetBinLowEdge(k);
            }
         }
         else
         {
            // l*10% bins
            if (cumul->GetBinContent(k) <= l * 0.1)
            {
               centralityNchg[l] = cumul->GetBinLowEdge(k);
            }
         }
      }
   }
}

bool netProtonCut(TLorentzVector *particle)
{
   // y < 0.5 && 0.4 < pT < 2.0
   if (TMath::Abs(particle->Rapidity()) < 0.5 && particle->Pt() > 0.4 && particle->Pt() < 2.0)
   {
      return true;
   }
   else
   {
      return false;
   }
}

bool centralityCut(TLorentzVector *particle)
{
   if (TMath::Abs(particle->Eta()) < 1 && particle->Pt() > 0.1)
   {
      return true;
   }
   else
   {
      return false;
   }
}

bool epRecoCut(TLorentzVector *particle)
{
   if (particle->Pt() > 0.2 && particle->Pt() < 2 && TMath::Abs(particle->Eta()) > 0.05 && TMath::Abs(particle->Eta()) < 1)
   {
      return true;
   }
   else
   {
      return false;
   }
}

void fillHistograms(TH2F **h2ImpNppPhi1Bin, TH2F *h2ImpNppPhi2Bin[][2], TH2F *h2ImpNppPhi3Bin[][3], TH2F *h2ImpNppPhi6Bin[][6], TH1F **h1EventPlane, bool onlyPrimary, double *resolution)
{
   TRandom *rndgen = new TRandom(time(0));
   // This function is for filling the histograms
   TTree *event = NULL;
   f14 reader(event);
   Long64_t nentries = reader.fChain->GetEntries();
   cout << "Total Entries:" << nentries << endl;
   Long64_t nbytes = 0, nb = 0;
   int eventCounter = 0;

   for (Long64_t jentry = 0; jentry < nentries; jentry++)
   {
      int primarieCounter = 0;
      if (jentry % 50000 == 0)
      {
         cout << "Progress %:" << jentry * 1.0 / nentries * 100 << endl;
      }
      // if (jentry == 50000)
      //    break;

      Long64_t ientry = reader.LoadTree(jentry);
      if (ientry < 0)
         break;

      double Psi[10];
      nb = reader.fChain->GetEntry(jentry);
      nbytes += nb;

      if (reader.fColHdr_timestep != 200 || reader.fColHdr_Ntrack == 0)
      {
         continue;
      }

      for (int i = 0; i < 10; i++)
      {
         // Get the Psi of this event according to the resolutions
         if (i == 0)
         {
            Psi[0] = 0.0;
            continue;
         }
         else
         {
            Psi[i] = rndgen->Gaus(0, resolution[i]);
         }
         h1EventPlane[i]->Fill(Psi[i]);
      }

      // Have known the Imp, fill the h2 with the phi of every particle
      int _nNetProtonPhi[10][6] = {{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
      for (Int_t j = 0; j < reader.fColHdr_Ntrack; j++)
      {
         if (onlyPrimary)
         {
            // Skip this particle if it is:
            // 1. Decay Secondaries;
            // 2. Spectators;
            // 3. Non-charged Paricles;
            if (reader.fTracksOut_pptype[j] == 20 || reader.fTracksOut_Nc[j] == 0 || reader.fTracksOut_chg[j] == 0)
            {
               continue;
            }
         }
         else
         {
            // Skip this particle if it is:
            // 2. Spectators;
            // 3. Non-charged Paricles;
            if (reader.fTracksOut_Nc[j] == 0 || reader.fTracksOut_chg[j] == 0)
            {
               continue;
            }
         }

         double pVector[4] = {reader.fTracksOut_px[j], reader.fTracksOut_py[j], reader.fTracksOut_pz[j], reader.fTracksOut_p0[j]};
         TLorentzVector *particle = new TLorentzVector(pVector);

         // Is it (anti-)proton?
         if (TMath::Abs(reader.fTracksOut_ityp[j]) == 1)
         {
            // (anti-)proton
            if (netProtonCut(particle))
            {
               for (int i = 0; i < 10; i++)
               {
                  int index = 0;
                  if (i != 0)
                  {
                     // Determine the index in (0,0.5pi)->(0,1,2,3,4,5)
                     index = (int)(atan(fabs((reader.fTracksOut_ry[j] * cos(Psi[i]) - reader.fTracksOut_rx[j] * sin(Psi[i])) / (reader.fTracksOut_rx[j] * cos(Psi[i]) + reader.fTracksOut_ry[j] * sin(Psi[i])))) / (TMath::Pi() * 1.0 / 12));
                  }
                  else
                  {
                     // if False, use Reaction Plane
                     index = (int)(atan(fabs(reader.fTracksOut_ry[j] / reader.fTracksOut_rx[j])) / (TMath::Pi() * 1.0 / 12));
                  }
                  // Net-Proton Multiplicity Counter
                  if (reader.fTracksOut_ityp[j] == 1)
                     _nNetProtonPhi[i][index]++;
                  else if (reader.fTracksOut_ityp[j] == -1)
                     _nNetProtonPhi[i][index]--;
               }
            }
         }
         // Delete this particle
         delete particle;
      }

      for (int i = 0; i < 10; i++)
      {
         // Fill the histograms
         for (int j = 0; j < 6; j++)
         {
            if (j < 3)
            {
               h2ImpNppPhi3Bin[i][j]->Fill(_nNetProtonPhi[i][2 * j] + _nNetProtonPhi[i][2 * j + 1], reader.fEvtHdr_imp);
               if (j < 2)
               {
                  h2ImpNppPhi2Bin[i][j]->Fill(_nNetProtonPhi[i][3 * j] + _nNetProtonPhi[i][3 * j + 1] + _nNetProtonPhi[i][3 * j + 2], reader.fEvtHdr_imp);
               }
            }
            h2ImpNppPhi6Bin[i][j]->Fill(_nNetProtonPhi[i][j], reader.fEvtHdr_imp);
         }
         h2ImpNppPhi1Bin[i]->Fill(_nNetProtonPhi[i][0] + _nNetProtonPhi[i][1] + _nNetProtonPhi[i][2] + _nNetProtonPhi[i][3] + _nNetProtonPhi[i][4] + _nNetProtonPhi[i][5], reader.fEvtHdr_imp);
         eventCounter++;
         // cout << _nNetProtonPhi[i][0] + _nNetProtonPhi[i][1] + _nNetProtonPhi[i][2] + _nNetProtonPhi[i][3] + _nNetProtonPhi[i][4] + _nNetProtonPhi[i][5]  << endl;
      }
   }
   cout << eventCounter << " events processed" << endl;
}

TH3I *drawEvent()
{
   TH3I *particles = new TH3I("particles", "particles", 100, -50, 50, 100, -50, 50, 100, -50, 50);
   TTree *event = NULL;

   f14 reader(event);

   Long64_t nentries = reader.fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry = 0; jentry < 50000; jentry++)
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

void getMoments(TH1D *netProtonDist, double *moments)
{
   // Loop on bins (possibly including underflows/overflows)
   Int_t binx;
   Double_t w;
   Double_t x;
   Int_t firstBinX = netProtonDist->GetXaxis()->GetFirst();
   Int_t lastBinX = netProtonDist->GetXaxis()->GetLast();

   double mean = netProtonDist->GetMean();

   int sum = 0;
   for (binx = firstBinX; binx <= lastBinX; binx++)
   {
      x = netProtonDist->GetXaxis()->GetBinCenter(binx);
      w = netProtonDist->GetBinContent(binx);
      sum += w;
      moments[2] += 1.0 * w * pow((x - mean), 2);
      moments[3] += 1.0 * w * pow((x - mean), 3);
      moments[4] += w * pow((x - mean), 4);
      moments[5] += w * pow((x - mean), 5);
      moments[6] += w * pow((x - mean), 6);
      moments[7] += w * pow((x - mean), 7);
      moments[8] += w * pow((x - mean), 8);
   }
   for (int i = 2; i < 9; i++)
   {
      moments[i] = moments[i] / sum;
   }
}

double varC2oC1(double *moments, double mean, int events)
{
   return (-(moments[2] * moments[2]) / (mean * mean) + moments[4] / (mean * mean) - (2 * moments[2] * moments[3]) / (mean * mean * mean) + (moments[2] * moments[2] * moments[2]) / (mean * mean * mean * mean)) / (1.0 * events);
}

double varC3oC2(double *moments, double mean, int events)
{
   return (9.0 * moments[2] - (6.0 * moments[4]) / moments[2] + (6.0 * moments[3] * moments[3]) / (moments[2] * moments[2]) + moments[6] / (moments[2] * moments[2]) - (2.0 * moments[3] * moments[5]) / (moments[2] * moments[2] * moments[2]) + moments[3] * moments[3] * moments[4] / (moments[2] * moments[2] * moments[2] * moments[2])) / (1.0 * events);
}

double varC4oC2(double *moments, double mean, int events)
{
   return (-9.0 * moments[2] * moments[2] + 9.0 * moments[4] + 40.0 * moments[3] * moments[3] / moments[2] - 6.0 * moments[6] / moments[2] - 8.0 * moments[3] * moments[5] / (moments[2] * moments[2]) + 6.0 * moments[4] * moments[4] / (moments[2] * moments[2]) + moments[8] / (moments[2] * moments[2]) + 8.0 * moments[3] * moments[3] * moments[4] / (moments[2] * moments[2] * moments[2]) - 2.0 * moments[4] * moments[6] / (moments[2] * moments[2] * moments[2]) + moments[4] * moments[4] * moments[4] / (moments[2] * moments[2] * moments[2] * moments[2])) / (1.0 * events);
}

void getCumulants(TH2F *h2ImpNpp, double *C2oC1, double *C3oC2, double *C4oC2, double *eC2oC1, double *eC3oC2, double *eC4oC2, double *centralityNchg)
{
   int binMax, binMin;
   int sum;

   // For centain phi, calculate Ci
   // CBWM included
   double C1[3] = {0, 0, 0}, C2[3] = {0, 0, 0}, C3[3] = {0, 0, 0}, C4[3] = {0, 0, 0};
   for (int l = 0; l < 3; l++)
   {
      sum = 0;
      if (l == 0)
      {
         // 0-10%, Set from centralityNchg to MinBin
         binMin = 0;
         binMax = h2ImpNpp->ProjectionY()->FindBin(centralityNchg[1]);
      }
      else if (l == 1)
      {
         // 10-60%, Set from centralityNchg to MinBin
         binMin = h2ImpNpp->ProjectionY()->FindBin(centralityNchg[1]);
         binMax = h2ImpNpp->ProjectionY()->FindBin(centralityNchg[6]);
      }
      else if (l == 2)
      {
         // 10-60%, Set from centralityNchg to MinBin
         binMin = h2ImpNpp->ProjectionY()->FindBin(centralityNchg[6]);
         binMax = h2ImpNpp->ProjectionY()->FindBin(centralityNchg[8]);
      }
      // Loop on every centrality bin and do the modification
      for (int j = binMin; j < binMax; j++)
      {
         TH1D *_netProtonDist = h2ImpNpp->ProjectionX("", j, j);
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
            continue;
         }
         _kappa = _netProtonDist->GetKurtosis();
         if (isnan(_kappa))
         {
            _kappa = 0;
            cout << "omitEvent: " << nEvents << endl;
            continue;
         }

         C1[l] += nEvents * _mean;
         C2[l] += nEvents * _sigma * _sigma;
         C3[l] += nEvents * _fS * _sigma * _sigma * _sigma;
         C4[l] += nEvents * _kappa * _sigma * _sigma * _sigma * _sigma;

         double _moments[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
         getMoments(_netProtonDist, _moments);
         eC2oC1[l] += nEvents * nEvents * varC2oC1(_moments, _mean, nEvents);
         eC3oC2[l] += nEvents * nEvents * varC3oC2(_moments, _mean, nEvents);
         eC4oC2[l] += nEvents * nEvents * varC4oC2(_moments, _mean, nEvents);

         sum += nEvents;
      }

      C1[l] = C1[l] / sum;
      C2[l] = C2[l] / sum;
      C3[l] = C3[l] / sum;
      C4[l] = C4[l] / sum;

      // Weighted Ratios
      C2oC1[l] = C2[l] / C1[l];
      C3oC2[l] = C3[l] / C2[l];
      C4oC2[l] = C4[l] / C2[l];

      eC2oC1[l] = TMath::Sqrt(eC2oC1[l]) / sum;
      eC3oC2[l] = TMath::Sqrt(eC3oC2[l]) / sum;
      eC4oC2[l] = TMath::Sqrt(eC4oC2[l]) / sum;
   }
}

void drawFigures(TH2F *h2ImpNppPhi1Bin, TH2F **h2ImpNppPhi2Bin, TH2F **h2ImpNppPhi3Bin, TH2F **h2ImpNppPhi6Bin, double *centralityNchg, TString outputplotsfolder)
{
   TCanvas *c1 = new TCanvas();
   double x1Bin = 0.5 * TMath::Pi() / 2;
   double x1err = 0.5 * TMath::Pi() / 2;
   double x2Bin[2] = {0.125 * TMath::Pi(), 0.375 * TMath::Pi()};
   double x2err[2] = {0.125 * TMath::Pi(), 0.125 * TMath::Pi()};
   double x3Bin[3] = {1.0 * TMath::Pi() / 12, 3.0 * TMath::Pi() / 12, 5.0 * TMath::Pi() / 12};
   double x3err[3] = {1.0 * TMath::Pi() / 12, 1.0 * TMath::Pi() / 12, 1.0 * TMath::Pi() / 12};
   double x6Bin[6] = {1.0 * TMath::Pi() / 24, 3.0 * TMath::Pi() / 24, 5.0 * TMath::Pi() / 24, 7.0 * TMath::Pi() / 24, 9.0 * TMath::Pi() / 24, 11.0 * TMath::Pi() / 24};
   double x6err[6] = {1.0 * TMath::Pi() / 24, 1.0 * TMath::Pi() / 24, 1.0 * TMath::Pi() / 24, 1.0 * TMath::Pi() / 24, 1.0 * TMath::Pi() / 24, 1.0 * TMath::Pi() / 24};

   double C2oC1VsPhi6Bin[3][6];
   double C3oC2VsPhi6Bin[3][6];
   double C4oC2VsPhi6Bin[3][6];

   double C2oC1VsPhi3Bin[3][3];
   double C3oC2VsPhi3Bin[3][3];
   double C4oC2VsPhi3Bin[3][3];

   double C2oC1VsPhi2Bin[3][2];
   double C3oC2VsPhi2Bin[3][2];
   double C4oC2VsPhi2Bin[3][2];

   double C2oC1VsPhi1Bin[3];
   double C3oC2VsPhi1Bin[3];
   double C4oC2VsPhi1Bin[3];

   double eC2oC1VsPhi6Bin[3][6];
   double eC3oC2VsPhi6Bin[3][6];
   double eC4oC2VsPhi6Bin[3][6];

   double eC2oC1VsPhi3Bin[3][3];
   double eC3oC2VsPhi3Bin[3][3];
   double eC4oC2VsPhi3Bin[3][3];

   double eC2oC1VsPhi2Bin[3][2];
   double eC3oC2VsPhi2Bin[3][2];
   double eC4oC2VsPhi2Bin[3][2];

   double eC2oC1VsPhi1Bin[3];
   double eC3oC2VsPhi1Bin[3];
   double eC4oC2VsPhi1Bin[3];

   // 6 Bins
   for (int i = 0; i < 6; i++)
   {
      double C2oC1[3] = {0}, C3oC2[3] = {0}, C4oC2[3] = {0};
      double eC2oC1[3] = {0}, eC3oC2[3] = {0}, eC4oC2[3] = {0};
      getCumulants(h2ImpNppPhi6Bin[i], C2oC1, C3oC2, C4oC2, eC2oC1, eC3oC2, eC4oC2, centralityNchg);

      for (int j = 0; j < 3; j++)
      {
         C2oC1VsPhi6Bin[j][i] = C2oC1[j];
         C3oC2VsPhi6Bin[j][i] = C3oC2[j];
         C4oC2VsPhi6Bin[j][i] = C4oC2[j];

         eC2oC1VsPhi6Bin[j][i] = eC2oC1[j];
         eC3oC2VsPhi6Bin[j][i] = eC3oC2[j];
         eC4oC2VsPhi6Bin[j][i] = eC4oC2[j];
      }
   }

   // 3 Bins
   for (int i = 0; i < 3; i++)
   {
      double C2oC1[3] = {0}, C3oC2[3] = {0}, C4oC2[3] = {0};
      double eC2oC1[3] = {0}, eC3oC2[3] = {0}, eC4oC2[3] = {0};
      getCumulants(h2ImpNppPhi3Bin[i], C2oC1, C3oC2, C4oC2, eC2oC1, eC3oC2, eC4oC2, centralityNchg);

      for (int j = 0; j < 3; j++)
      {
         C2oC1VsPhi3Bin[j][i] = C2oC1[j];
         C3oC2VsPhi3Bin[j][i] = C3oC2[j];
         C4oC2VsPhi3Bin[j][i] = C4oC2[j];

         eC2oC1VsPhi3Bin[j][i] = eC2oC1[j];
         eC3oC2VsPhi3Bin[j][i] = eC3oC2[j];
         eC4oC2VsPhi3Bin[j][i] = eC4oC2[j];
      }
   }

   // 2 Bins
   for (int i = 0; i < 2; i++)
   {
      double C2oC1[3] = {0}, C3oC2[3] = {0}, C4oC2[3] = {0};
      double eC2oC1[3] = {0}, eC3oC2[3] = {0}, eC4oC2[3] = {0};
      getCumulants(h2ImpNppPhi2Bin[i], C2oC1, C3oC2, C4oC2, eC2oC1, eC3oC2, eC4oC2, centralityNchg);

      for (int j = 0; j < 3; j++)
      {
         C2oC1VsPhi2Bin[j][i] = C2oC1[j];
         C3oC2VsPhi2Bin[j][i] = C3oC2[j];
         C4oC2VsPhi2Bin[j][i] = C4oC2[j];

         eC2oC1VsPhi2Bin[j][i] = eC2oC1[j];
         eC3oC2VsPhi2Bin[j][i] = eC3oC2[j];
         eC4oC2VsPhi2Bin[j][i] = eC4oC2[j];
      }
   }

   double C2oC1[3] = {0}, C3oC2[3] = {0}, C4oC2[3] = {0};
   double eC2oC1[3] = {0}, eC3oC2[3] = {0}, eC4oC2[3] = {0};
   getCumulants(h2ImpNppPhi1Bin, C2oC1, C3oC2, C4oC2, eC2oC1, eC3oC2, eC4oC2, centralityNchg);

   for (int j = 0; j < 3; j++)
   {
      C2oC1VsPhi1Bin[j] = C2oC1[j];
      C3oC2VsPhi1Bin[j] = C3oC2[j];
      C4oC2VsPhi1Bin[j] = C4oC2[j];

      eC2oC1VsPhi1Bin[j] = eC2oC1[j];
      eC3oC2VsPhi1Bin[j] = eC3oC2[j];
      eC4oC2VsPhi1Bin[j] = eC4oC2[j];
   }

   for (int i = 0; i < 3; i++)
   {
      c1->Clear();
      TGraphErrors *h1C2oC1VsPhi1Bin = new TGraphErrors(1, &x1Bin, &C2oC1VsPhi1Bin[i], &x1err, eC2oC1VsPhi1Bin);
      h1C2oC1VsPhi1Bin->SetLineColor(2);
      h1C2oC1VsPhi1Bin->SetMarkerSize(1);
      h1C2oC1VsPhi1Bin->SetMarkerStyle(20);
      TGraphErrors *h1C2oC1VsPhi2Bin = new TGraphErrors(2, x2Bin, C2oC1VsPhi2Bin[i], x2err, eC2oC1VsPhi2Bin[i]);
      h1C2oC1VsPhi2Bin->SetLineColor(3);
      h1C2oC1VsPhi2Bin->SetMarkerSize(1);
      h1C2oC1VsPhi2Bin->SetMarkerStyle(21);
      TGraphErrors *h1C2oC1VsPhi3Bin = new TGraphErrors(3, x3Bin, C2oC1VsPhi3Bin[i], x3err, eC2oC1VsPhi3Bin[i]);
      h1C2oC1VsPhi3Bin->SetLineColor(4);
      h1C2oC1VsPhi3Bin->SetMarkerSize(1);
      h1C2oC1VsPhi3Bin->SetMarkerStyle(22);
      TGraphErrors *h1C2oC1VsPhi6Bin = new TGraphErrors(6, x6Bin, C2oC1VsPhi6Bin[i], x6err, eC2oC1VsPhi6Bin[i]);
      h1C2oC1VsPhi6Bin->SetLineColor(7);
      h1C2oC1VsPhi6Bin->SetMarkerSize(1);
      h1C2oC1VsPhi6Bin->SetMarkerStyle(23);
      TMultiGraph *mg = new TMultiGraph();
      mg->Add(h1C2oC1VsPhi1Bin);
      mg->Add(h1C2oC1VsPhi2Bin);
      mg->Add(h1C2oC1VsPhi3Bin);
      mg->Add(h1C2oC1VsPhi6Bin);
      mg->SetTitle("C_{2}/C_{1}");
      mg->GetXaxis()->SetTitle("#phi-#psi_{EP}");

      mg->Draw("zpA");
      c1->SaveAs(Form(outputplotsfolder + "Centrality%d/h1C2oC1VsPhiCentral.png", i));
      c1->SetTitle(Form("h1C2oC1VsPhiCentral%d", i));

      TGraphErrors *h1C3oC2VsPhi1Bin = new TGraphErrors(1, &x1Bin, &C3oC2VsPhi1Bin[i], &x1err, eC3oC2VsPhi1Bin);
      h1C3oC2VsPhi1Bin->SetLineColor(2);
      h1C3oC2VsPhi1Bin->SetMarkerSize(1);
      h1C2oC1VsPhi1Bin->SetMarkerStyle(20);
      TGraphErrors *h1C3oC2VsPhi2Bin = new TGraphErrors(2, x2Bin, C3oC2VsPhi2Bin[i], x2err, eC3oC2VsPhi2Bin[i]);
      h1C3oC2VsPhi2Bin->SetLineColor(3);
      h1C3oC2VsPhi2Bin->SetMarkerSize(1);
      h1C3oC2VsPhi2Bin->SetMarkerStyle(21);
      TGraphErrors *h1C3oC2VsPhi3Bin = new TGraphErrors(3, x3Bin, C3oC2VsPhi3Bin[i], x3err, eC3oC2VsPhi3Bin[i]);
      h1C3oC2VsPhi3Bin->SetLineColor(4);
      h1C3oC2VsPhi3Bin->SetMarkerSize(1);
      h1C3oC2VsPhi3Bin->SetMarkerStyle(22);
      TGraphErrors *h1C3oC2VsPhi6Bin = new TGraphErrors(6, x6Bin, C3oC2VsPhi6Bin[i], x6err, eC3oC2VsPhi6Bin[i]);
      h1C3oC2VsPhi6Bin->SetLineColor(7);
      h1C3oC2VsPhi6Bin->SetMarkerSize(1);
      h1C3oC2VsPhi6Bin->SetMarkerStyle(23);
      mg = new TMultiGraph();
      mg->Add(h1C3oC2VsPhi1Bin);
      mg->Add(h1C3oC2VsPhi2Bin);
      mg->Add(h1C3oC2VsPhi3Bin);
      mg->Add(h1C3oC2VsPhi6Bin);
      mg->SetTitle("C_{3}/C_{2}");
      mg->GetXaxis()->SetTitle("#phi-#psi_{RP}");

      mg->Draw("zpA");
      c1->Size(0, 0);
      c1->SaveAs(Form(outputplotsfolder + "Centrality%d/h1C3oC2VsPhiCentral.png", i));
      c1->SetTitle(Form("h1vVsPhiCentral%d", i));

      c1->Clear();
      TGraphErrors *h1C4oC2VsPhi1Bin = new TGraphErrors(1, &x1Bin, &C4oC2VsPhi1Bin[i], &x1err, eC4oC2VsPhi1Bin);

      h1C4oC2VsPhi1Bin->SetLineColor(2);
      h1C4oC2VsPhi1Bin->SetMarkerSize(1);
      h1C4oC2VsPhi1Bin->SetMarkerStyle(20);
      TGraphErrors *h1C4oC2VsPhi2Bin = new TGraphErrors(2, x2Bin, C4oC2VsPhi2Bin[i], x2err, eC4oC2VsPhi2Bin[i]);
      h1C4oC2VsPhi2Bin->SetLineColor(3);
      h1C4oC2VsPhi2Bin->SetMarkerSize(1);
      h1C4oC2VsPhi2Bin->SetMarkerStyle(21);
      TGraphErrors *h1C4oC2VsPhi3Bin = new TGraphErrors(3, x3Bin, C4oC2VsPhi3Bin[i], x3err, eC4oC2VsPhi3Bin[i]);
      h1C4oC2VsPhi3Bin->SetLineColor(4);
      h1C4oC2VsPhi3Bin->SetMarkerSize(1);
      h1C4oC2VsPhi3Bin->SetMarkerStyle(22);
      TGraphErrors *h1C4oC2VsPhi6Bin = new TGraphErrors(6, x6Bin, C4oC2VsPhi6Bin[i], x6err, eC4oC2VsPhi6Bin[i]);
      h1C4oC2VsPhi6Bin->SetLineColor(7);
      h1C4oC2VsPhi6Bin->SetMarkerSize(1);
      h1C3oC2VsPhi6Bin->SetMarkerStyle(23);
      mg = new TMultiGraph();
      mg->Add(h1C4oC2VsPhi1Bin);
      mg->Add(h1C4oC2VsPhi2Bin);
      mg->Add(h1C4oC2VsPhi3Bin);
      mg->Add(h1C4oC2VsPhi6Bin);
      mg->SetTitle("C_{4}/C_{2}");
      mg->GetXaxis()->SetTitle("#phi-#psi_{RP}");

      mg->Draw("zpA");
      c1->Size(0, 0);
      c1->SaveAs(Form(outputplotsfolder + "Centrality%d/h1C4oC2VsPhiCentral.png", i));
      c1->SetTitle(Form("h1vVsPhiCentral%d", i));
   }

   TGraphErrors *h1C4oC2VsPhiNBin[3];
   auto legend = new TLegend(0.6,0.6,0.9,0.9);
   legend->SetHeader("The Legend Title","C");
   for (int i = 0; i < 3; i++)
   {
      double nBins[4] = {1, 2, 3, 6};
      double enBins[4] = {0.75, 0.75, 0.75, 0.75};
      double yNBins[4] = {C4oC2VsPhi1Bin[i], 0.5 * (C4oC2VsPhi2Bin[i][0] + C4oC2VsPhi2Bin[i][1]), 1.0 / 3.0 * (C4oC2VsPhi3Bin[i][0] + C4oC2VsPhi3Bin[i][1] + C4oC2VsPhi3Bin[i][2]), 1.0 / 6.0 * (C4oC2VsPhi6Bin[i][0] + C4oC2VsPhi6Bin[i][1] + C4oC2VsPhi6Bin[i][2] + C4oC2VsPhi6Bin[i][3] + C4oC2VsPhi6Bin[i][4] + C4oC2VsPhi6Bin[i][5])};
      double eyNBins[4] = {eC4oC2VsPhi1Bin[i], 0.5 * TMath::Sqrt(eC4oC2VsPhi2Bin[i][0] * eC4oC2VsPhi2Bin[i][0] + eC4oC2VsPhi2Bin[i][1] * eC4oC2VsPhi2Bin[i][1]), 1.0 / 3.0 * TMath::Sqrt(eC4oC2VsPhi3Bin[i][0] * eC4oC2VsPhi3Bin[i][0] + eC4oC2VsPhi3Bin[i][1] * eC4oC2VsPhi3Bin[i][1] + eC4oC2VsPhi3Bin[i][2] * eC4oC2VsPhi3Bin[i][2]), 1.0 / 6.0 * TMath::Sqrt(eC4oC2VsPhi6Bin[i][0] * eC4oC2VsPhi6Bin[i][0] + eC4oC2VsPhi6Bin[i][1] * eC4oC2VsPhi6Bin[i][1] + eC4oC2VsPhi6Bin[i][2] * eC4oC2VsPhi6Bin[i][2] + eC4oC2VsPhi6Bin[i][3] * eC4oC2VsPhi6Bin[i][3] + eC4oC2VsPhi6Bin[i][4] * eC4oC2VsPhi6Bin[i][4] + eC4oC2VsPhi6Bin[i][5] * eC4oC2VsPhi6Bin[i][5])};
      h1C4oC2VsPhiNBin[i] = new TGraphErrors(4, nBins, yNBins, enBins, eyNBins);
      if (i == 0)
      {
         h1C4oC2VsPhiNBin[i]->SetLineWidth(2);
         h1C4oC2VsPhiNBin[i]->SetLineColor(2);
         h1C4oC2VsPhiNBin[i]->SetMarkerSize(1);
         h1C4oC2VsPhiNBin[i]->SetMarkerStyle(20);
         legend->AddEntry(h1C4oC2VsPhiNBin[i],"0-10%","ep");
      }
      else if (i == 1)
      {
         h1C4oC2VsPhiNBin[i]->SetLineWidth(2);
         h1C4oC2VsPhiNBin[i]->SetLineColor(3);
         h1C4oC2VsPhiNBin[i]->SetMarkerSize(1);
         h1C4oC2VsPhiNBin[i]->SetMarkerStyle(21);
         legend->AddEntry(h1C4oC2VsPhiNBin[i],"10-60%","ep");
      }
      else
      {
         h1C4oC2VsPhiNBin[i]->SetLineWidth(2);
         h1C4oC2VsPhiNBin[i]->SetLineColor(7);
         h1C4oC2VsPhiNBin[i]->SetMarkerSize(1);
         h1C4oC2VsPhiNBin[i]->SetMarkerStyle(23);
         legend->AddEntry(h1C4oC2VsPhiNBin[i],"60-80%","ep");
      }
   }
   TMultiGraph *mg = new TMultiGraph();
   mg->Add(h1C4oC2VsPhiNBin[0]);
   mg->Add(h1C4oC2VsPhiNBin[1]);
   mg->Add(h1C4oC2VsPhiNBin[2]);
   mg->GetYaxis()->SetTitle("C_{4}/C_{2} Vs NBins");
   mg->GetXaxis()->SetTitle("Number of Bins");
   mg->SetTitle("");
   mg->Draw("zpA");
   legend->Draw();
   c1->Size(0, 0);
   c1->SaveAs(outputplotsfolder + "h1C4oC2VsNBins.png");

   c1->SetLogy(1);
   h2ImpNppPhi1Bin->ProjectionX()->DrawNormalized("E");
   h2ImpNppPhi1Bin->SetTitle("Net-proton production vs. N_{chg}");
   c1->SaveAs(outputplotsfolder + "netPro.pdf");
   c1->SaveAs(outputplotsfolder + "netPro.png");

   h2ImpNppPhi1Bin->ProjectionY()->DrawNormalized("E");
   c1->SaveAs(outputplotsfolder + "refmult3.pdf");
   c1->SaveAs(outputplotsfolder + "refmult3.png");

   TH1 *cumul = h2ImpNppPhi1Bin->ProjectionY()->GetCumulative(false);
   cumul->SetTitle("Cumulative Charged Particle Multiplicity");
   cumul->Scale(1.0 / h2ImpNppPhi1Bin->ProjectionY()->Integral());
   cumul->Draw("E");
   c1->SaveAs(outputplotsfolder + "refmult3Cumu.pdf");
   c1->SaveAs(outputplotsfolder + "refmult3Cumu.png");

   c1->SetLogz(1);
   c1->SetLogy(0);
   h2ImpNppPhi1Bin->Draw("COLZ");
   c1->SaveAs(outputplotsfolder + "h2ImpNppPhi1Bin.pdf");
   c1->SaveAs(outputplotsfolder + "h2ImpNppPhi1Bin.png");
   c1->Close();
}

void netProton(TString inputFile = "primaryOnly.root", TString outputFile = "primaryOutput.root", TString outputplotsfolder = "primaryOnly/", bool onlyPrimary = false)
{
   system("rm -rf " + outputplotsfolder);
   system("mkdir -p " + outputplotsfolder);
   for (int i = 0; i < 10; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         system(Form("mkdir -p " + outputplotsfolder + "Resolution%d/Centrality%d", i, j));
      }
   }

   double resolution[10];
   for (int i = 0; i < 10; i++)
   {
      resolution[i] = 1.0 / 3.0 * 0.5 * TMath::Pi() * i * 1.0 / 10.0;
   }

   TCanvas *c1 = new TCanvas();
   c1->cd();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0101);

   TFile f(inputFile);

   TH2F *h2ImpNppPhi6Bin[10][6];
   TH2F *h2ImpNppPhi3Bin[10][3];
   TH2F *h2ImpNppPhi2Bin[10][2];
   TH2F *h2ImpNppPhi1Bin[10];
   TH1F *h1EventPlane[10];

   for (int j = 0; j < 10; j++)
   {
      f.GetObject(Form("h2ImpNppPhi1Bin%d", j), h2ImpNppPhi1Bin[j]);
      for (int i = 0; i < 6; i++)
      {
         if (i < 3)
         {
            if (i < 2)
            {
               f.GetObject(Form("h2ImpNppPhi2Bin%d%d", i, j), h2ImpNppPhi2Bin[j][i]);
            }
            f.GetObject(Form("h2ImpNppPhi3Bin%d%d", i, j), h2ImpNppPhi3Bin[j][i]);
         }
         f.GetObject(Form("h2ImpNppPhi6Bin%d%d", i, j), h2ImpNppPhi6Bin[j][i]);
      }
      f.GetObject(Form("h1EventPlane%d", j), h1EventPlane[j]);
   }

   if (!h1EventPlane[0])
   {
      for (int j = 0; j < 10; j++)
      {
         h2ImpNppPhi1Bin[j] = new TH2F(Form("h2ImpNppPhi1Bin%d", j), Form("h2ImpNppPhi1Bin%d", j), 100, -5.5, 94.5, 130, 0, 13.0);
         h2ImpNppPhi1Bin[j]->GetYaxis()->SetTitle("Impact Parameter");
         h2ImpNppPhi1Bin[j]->GetXaxis()->SetTitle("Net-Proton");
         for (int i = 0; i < 6; i++)
         {
            if (i < 3)
            {
               if (i < 2)
               {
                  h2ImpNppPhi2Bin[j][i] = new TH2F(Form("h2ImpNppPhi2Bin%d%d", i, j), Form("h2ImpNppPhi2Bin%d%d", i, j), 100, -5.5, 94.5, 130, 0, 13.0);
                  h2ImpNppPhi2Bin[j][i]->GetYaxis()->SetTitle("Impact Parameter");
                  h2ImpNppPhi2Bin[j][i]->GetXaxis()->SetTitle("Net-Proton");
                  h2ImpNppPhi2Bin[j][i]->GetZaxis()->SetTitle("#phi");
               }
               h2ImpNppPhi3Bin[j][i] = new TH2F(Form("h2ImpNppPhi3Bin%d%d", i, j), Form("h2ImpNppPhi3Bin%d%d", i, j), 100, -5.5, 94.5, 130, 0, 13.0);
               h2ImpNppPhi3Bin[j][i]->GetYaxis()->SetTitle("Impact Parameter");
               h2ImpNppPhi3Bin[j][i]->GetXaxis()->SetTitle("Net-Proton");
               h2ImpNppPhi3Bin[j][i]->GetZaxis()->SetTitle("#phi");
            }
            h2ImpNppPhi6Bin[j][i] = new TH2F(Form("h2ImpNppPhi6Bin%d%d", i, j), Form("h2ImpNppPhi6Bin%d%d", i, j), 100, -5.5, 94.5, 130, 0, 13.0);
            h2ImpNppPhi6Bin[j][i]->GetYaxis()->SetTitle("Impact Parameter");
            h2ImpNppPhi6Bin[j][i]->GetXaxis()->SetTitle("Net-Proton");
            h2ImpNppPhi6Bin[j][i]->GetZaxis()->SetTitle("#phi");
         }
         h1EventPlane[j] = new TH1F(Form("h1EventPlane%d", j), Form("h1EventPlane%d", j), 100, -TMath::Pi(), TMath::Pi());
      }

      fillHistograms(h2ImpNppPhi1Bin, h2ImpNppPhi2Bin, h2ImpNppPhi3Bin, h2ImpNppPhi6Bin, h1EventPlane, onlyPrimary, resolution);
   }

   double centralityNchg[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
   centralityBound(centralityNchg, h2ImpNppPhi1Bin[0]->ProjectionY());
   ofstream logOut;
   logOut.open(outputplotsfolder + "centralityImp.log", ios_base::app);
   logOut << "Centrality Impact Parameter" << endl;

   for (int l = 0; l < 9; l++)
   {
      if (l == 0)
      {
         logOut << 0 << "%-5%: " << centralityNchg[l] << endl;
      }
      else
      {
         logOut << l * 10 << "%: " << centralityNchg[l] << endl;
      }
   }

   TFile *hfile = TFile::Open(outputFile, "RECREATE");

   for (int i = 0; i < 10; i++)
   {
      h1EventPlane[i]->Write();
      h1EventPlane[i]->GetXaxis()->SetTitle("cos 2(#psi_{2#eta+}-#psi_{2#eta-})");
      h1EventPlane[i]->Draw("COLZ");
      c1->SaveAs(Form(outputplotsfolder + "Resolution%d/h1EventPlane.png", i));

      h2ImpNppPhi1Bin[i]->Write();
      h2ImpNppPhi1Bin[i]->Draw("COLZ");
      c1->SaveAs(Form(outputplotsfolder + "Resolution%d/h2ImpNppPhi1Bin.png", i));
      for (int j = 0; j < 6; j++)
      {
         if (j < 3)
         {
            if (j < 2)
            {
               h2ImpNppPhi2Bin[i][j]->Write();
               h2ImpNppPhi2Bin[i][j]->Draw("COLZ");
               c1->SaveAs(Form(outputplotsfolder + "Resolution%d/h2ImpNppPhi2Bin%d.png", i, j));
            }
            h2ImpNppPhi3Bin[i][j]->Write();
            h2ImpNppPhi3Bin[i][j]->Draw("COLZ");
            c1->SaveAs(Form(outputplotsfolder + "Resolution%d/h2ImpNppPhi3Bin%d.png", i, j));
         }
         h2ImpNppPhi6Bin[i][j]->Write();
         h2ImpNppPhi6Bin[i][j]->Draw("COLZ");
         c1->SaveAs(Form(outputplotsfolder + "Resolution%d/h2ImpNppPhi6Bin%d.png", i, j));
      }
   }
   // Loop on every angle
   hfile->Close();
   c1->Close();
   for (int i = 0; i < 10; i++)
   {
      drawFigures(h2ImpNppPhi1Bin[i], h2ImpNppPhi2Bin[i], h2ImpNppPhi3Bin[i], h2ImpNppPhi6Bin[i], centralityNchg, Form(outputplotsfolder + "Resolution%d/", i));
   }
}
