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

void centralityBound(double *centralityNchg, TH1D *refmult3)
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
               centralityNchg[0] = cumul->GetBinLowEdge(k);
            }
         }
         else
         {
            if (cumul->GetBinContent(k) > l * 0.1)
            {
               centralityNchg[l] = cumul->GetBinLowEdge(k);
            }
         }
      }
   }
}

void fillHistograms(TH2F *h2NcgNppPhi1Bin, TH2F **h2NcgNppPhi2Bin, TH2F **h2NcgNppPhi3Bin, TH2F **h2NcgNppPhi6Bin,TH2F *h2EventPlane, bool reconstruction)
{
   TTree *event = NULL;

   f14 reader(event);

   Long64_t nentries = reader.fChain->GetEntries();
   cout << nentries << endl;
   Long64_t nbytes = 0, nb = 0;
   int eventCounter=0;
   for (Long64_t jentry = 0; jentry < nentries; jentry++)
   {
      if (jentry % 50000 == 0)
      {
         cout << "Progress %:" << jentry * 1.0 / nentries * 100 << endl;
      }
      Long64_t ientry = reader.LoadTree(jentry);
      if (ientry < 0)
         break;
      int _nNetProton = 0;
      int _nCharged = 0;
      int _nTotalCharged = 0;
      double sumQxP = 0, sumQyP = 0;
      double sumQxM = 0, sumQyM = 0;
      nb = reader.fChain->GetEntry(jentry);
      nbytes += nb;

      // if (Cut(ientry) < 0) continue;
      if (reader.fColHdr_timestep != 200 || reader.fColHdr_Ntrack == 0)
      {
         continue;
      }
      else
      {
         eventCounter++;
      }

      // Read a Event, Determine the Centrality and Event plane of this event
      for (Int_t j = 0; j < reader.fColHdr_Ntrack; j++)
      {
         // Charged Particle
         if ((reader.fTracksOut_chg[j] != 0))
         {
            float pT = TMath::Sqrt(reader.fTracksOut_px[j] * reader.fTracksOut_px[j] + reader.fTracksOut_py[j] * reader.fTracksOut_py[j]);
            float p = TMath::Sqrt(reader.fTracksOut_px[j] * reader.fTracksOut_px[j] + reader.fTracksOut_py[j] * reader.fTracksOut_py[j] + reader.fTracksOut_pz[j] * reader.fTracksOut_pz[j]);
            float eta = 0.5 * log((p + reader.fTracksOut_pz[j]) / (p - reader.fTracksOut_pz[j]));
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
               if (TMath::Abs(eta) < 1 && pT > 0.1)
               {
                  _nCharged++;
               }
            }
            if (fabs(eta)>0.05 && pT>0.2 && pT<2 && fabs(eta)<1)
            {
               if (eta>0)
               {
                  sumQxP += pT * cos(2 * atan2(reader.fTracksOut_ry[j], reader.fTracksOut_rx[j]));
                  sumQyP += pT * sin(2 * atan2(reader.fTracksOut_ry[j], reader.fTracksOut_rx[j]));
               }
               else
               {  
                  sumQxM += pT * cos(2 * atan2(reader.fTracksOut_ry[j], reader.fTracksOut_rx[j]));
                  sumQyM += pT * sin(2 * atan2(reader.fTracksOut_ry[j], reader.fTracksOut_rx[j]));
               }
            }
            _nTotalCharged++;
         }
      }
      // Event Plane Psi2
      double PsiP;
      double PsiM;
      // Make Psi in 0, pi
      if (atan2(sumQyP, sumQxP) > 0)
         PsiP = 0.5 * atan2(sumQyP, sumQxP);
      else
         PsiP = 0.5 * (atan2(sumQyP, sumQxP) + 2 * TMath::Pi());

      if (atan2(sumQyM, sumQxM) > 0)
         PsiM = 0.5 * atan2(sumQyM, sumQxM);
      else
         PsiM = 0.5 * (atan2(sumQyM, sumQxM) + 2 * TMath::Pi());
         
      // Fill the distribution
      h2EventPlane->Fill(cos(2*(PsiP-PsiM)),_nCharged);
      h2NcgNppPhi1Bin->Fill(_nNetProton, _nCharged);
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
               float p = TMath::Sqrt(reader.fTracksOut_px[j] * reader.fTracksOut_px[j] + reader.fTracksOut_py[j] * reader.fTracksOut_py[j] + reader.fTracksOut_pz[j] * reader.fTracksOut_pz[j]);
               float y = 0.5 * log((reader.fTracksOut_p0[j] + reader.fTracksOut_pz[j]) / (reader.fTracksOut_p0[j] - reader.fTracksOut_pz[j]));
               float eta = 0.5 * log((p + reader.fTracksOut_pz[j]) / (p - reader.fTracksOut_pz[j]));
               if (TMath::Abs(y) < 0.5 && pT < 2.0 && pT > 0.4)
               {
                  int index = 0;
                  double Psi=0;
                  if (reconstruction)
                  {
                     if (eta>0)
                     {
                        Psi=PsiM;
                     }
                     else
                     {
                        Psi=PsiP;
                     }
                     index = (int)(atan(fabs((reader.fTracksOut_ry[j] * cos(Psi) - reader.fTracksOut_rx[j] * sin(Psi)) / (reader.fTracksOut_rx[j] * cos(Psi) + reader.fTracksOut_ry[j] * sin(Psi)))) / (TMath::Pi() * 1.0 / 12));
                  }
                  else
                  {
                     index = (int)(atan(fabs(reader.fTracksOut_ry[j] / reader.fTracksOut_rx[j])) / (TMath::Pi() * 1.0 / 12));
                  }
                  if (reader.fTracksOut_ityp[j] == 1)
                     _nNetProtonPhi[index]++;
                  else if (reader.fTracksOut_ityp[j] == -1)
                     _nNetProtonPhi[index]--;
               }
            }
            else
            {
               // Other Particle
            }
         }
      }
      for (int j = 0; j < 6; j++)
      {
         if (j<3)
         {
            h2NcgNppPhi3Bin[j]->Fill(_nNetProtonPhi[2*j]+_nNetProtonPhi[2*j+1], _nCharged);
            if(j<2)
               {
                  h2NcgNppPhi2Bin[j]->Fill(_nNetProtonPhi[3*j]+_nNetProtonPhi[3*j+1]+_nNetProtonPhi[3*j+2], _nCharged);
               }
         }
         h2NcgNppPhi6Bin[j]->Fill(_nNetProtonPhi[j], _nCharged);
      }
   }
   cout << eventCounter << " events processed" <<endl;
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

void getCumulants(TH2F *h2NcgNpp, double *C1, double *C2, double *C3, double *C4, double *centralityNchg)
{
   int binMax, binMin;
   int sum;
   double mean[3] = {0}, sigma[3] = {0}, fS[3] = {0}, kappa[3] = {0};
   // For centain phi, calculate Ci
   // CBWM included
   for (int l = 0; l < 3; l++)
   {
      sum = 0;
      if (l == 0)
      {
         // 0-10%, Set from centralityNchg to MinBin
         binMin = h2NcgNpp->ProjectionY()->GetBin(centralityNchg[1]);
         binMax = h2NcgNpp->ProjectionY()->GetMinimumBin();
      }
      else if (l == 1)
      {
         // 10-60%, Set from centralityNchg to MinBin
         binMin = h2NcgNpp->ProjectionY()->GetBin(centralityNchg[6]);
         binMax = h2NcgNpp->ProjectionY()->GetBin(centralityNchg[1]);
      }
      else if (l == 2)
      {
         // 10-60%, Set from centralityNchg to MinBin
         binMin = h2NcgNpp->ProjectionY()->GetBin(centralityNchg[8]);
         binMax = h2NcgNpp->ProjectionY()->GetBin(centralityNchg[6]);
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

void calcResolution(TH2F *h2EventPlane, double *centralityNchg,double *Resolution)
{
   int binMax, binMin;
   int sum;
   double mean[3] = {0,0,0};
   // For centain phi, calculate Ci
   // CBWM include
   for (int l = 0; l < 3; l++)
   {
      sum = 0;
      if (l == 0)
      {
         // 0-10%, Set from centralityNchg to MinBin
         binMin = h2EventPlane->ProjectionY()->GetBin(centralityNchg[1]);
         binMax = h2EventPlane->ProjectionY()->GetMinimumBin();
      }
      else if (l == 1)
      {
         // 10-60%, Set from centralityNchg to MinBin
         binMin = h2EventPlane->ProjectionY()->GetBin(centralityNchg[6]);
         binMax = h2EventPlane->ProjectionY()->GetBin(centralityNchg[1]);
      }
      else if (l == 2)
      {
         // 10-60%, Set from centralityNchg to MinBin
         binMin = h2EventPlane->ProjectionY()->GetBin(centralityNchg[8]);
         binMax = h2EventPlane->ProjectionY()->GetBin(centralityNchg[6]);
      }
      // Loop on every centrality bin and do the modification
      for (int j = binMin; j < binMax; j++)
      {
         TH1D *_eventPlaneDist = h2EventPlane->ProjectionX("", j, j);

         int nEvents = _eventPlaneDist->Integral();

         if (nEvents == 0 or nEvents == 1)
            continue;

         float _mean = 0;
         if(!(_eventPlaneDist->GetMean()<0))
            _mean = sqrt(_eventPlaneDist->GetMean());
         else
            {
               cout << "ERROR"<< endl;
               continue;
            }
         mean[l] += nEvents * _mean;
         sum += nEvents;
      }
      Resolution[l] = mean[l] / sum;
   }
}


void netProton(TString outputplotsfolder = "outputFolder/", bool reconstructionRP = true)
{
   system("rm -rf " + outputplotsfolder);
   system("mkdir -p " + outputplotsfolder);
   for (int i = 0; i < 3; i++)
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

   TH2F *h2NcgNppPhi6Bin[6];
   TH2F *h2NcgNppPhi3Bin[3];
   TH2F *h2NcgNppPhi2Bin[2];
   TH2F *h2NcgNppPhi1Bin = NULL;
   TH2F *h2EventPlane=NULL;

   f.GetObject("h2NcgNppPhi1Bin", h2NcgNppPhi1Bin);
   f.GetObject("h2EventPlane", h2EventPlane);
   for (int i = 0; i < 6; i++)
   {
      if (i < 3)
      {
         if (i < 2)
         {
            f.GetObject(Form("h2NcgNppPhi2Bin%d", i), h2NcgNppPhi2Bin[i]);
         }
         f.GetObject(Form("h2NcgNppPhi3Bin%d", i), h2NcgNppPhi3Bin[i]);
      }
      f.GetObject(Form("h2NcgNppPhi6Bin%d", i), h2NcgNppPhi6Bin[i]);
   }

   if (!h2NcgNppPhi1Bin)
   {
      h2NcgNppPhi1Bin = new TH2F("h2NcgNppPhi1Bin", "h2NcgNppPhi1Bin", 100, -5.5, 94.5, 450, -0.5, 449.5);
      h2NcgNppPhi1Bin->GetYaxis()->SetTitle("Charged Multiplicity");
      h2NcgNppPhi1Bin->GetXaxis()->SetTitle("Net-Proton");
      for (int i = 0; i < 6; i++)
      {
         if (i < 3)
         {
            if (i < 2)
            {
               h2NcgNppPhi2Bin[i] = new TH2F(Form("h2NcgNppPhi2Bin%d", i), Form("h2NcgNppPhi2Bin%d", i), 100, -5.5, 94.5, 450, -0.5, 449.5);
               h2NcgNppPhi2Bin[i]->GetYaxis()->SetTitle("Charged Multiplicity");
               h2NcgNppPhi2Bin[i]->GetXaxis()->SetTitle("Net-Proton");
               h2NcgNppPhi2Bin[i]->GetZaxis()->SetTitle("#phi");
            }
            h2NcgNppPhi3Bin[i] = new TH2F(Form("h2NcgNppPhi3Bin%d", i), Form("h2NcgNppPhi3Bin%d", i), 100, -5.5, 94.5, 450, -0.5, 449.5);
            h2NcgNppPhi3Bin[i]->GetYaxis()->SetTitle("Charged Multiplicity");
            h2NcgNppPhi3Bin[i]->GetXaxis()->SetTitle("Net-Proton");
            h2NcgNppPhi3Bin[i]->GetZaxis()->SetTitle("#phi");
         }
         h2NcgNppPhi6Bin[i] = new TH2F(Form("h2NcgNppPhi6Bin%d", i), Form("h2NcgNppPhi6Bin%d", i), 100, -5.5, 94.5, 450, -0.5, 449.5);
         h2NcgNppPhi6Bin[i]->GetYaxis()->SetTitle("Charged Multiplicity");
         h2NcgNppPhi6Bin[i]->GetXaxis()->SetTitle("Net-Proton");
         h2NcgNppPhi6Bin[i]->GetZaxis()->SetTitle("#phi");
      }
      h2EventPlane = new TH2F("h2EventPlane", "h2EventPlane", 1000, -1, 1, 450, -0.5, 449.5);
      fillHistograms(h2NcgNppPhi1Bin, h2NcgNppPhi2Bin, h2NcgNppPhi3Bin, h2NcgNppPhi6Bin, h2EventPlane, reconstructionRP);
      h2EventPlane->GetXaxis()->SetTitle("cos 2(#psi_{2#eta+}-#psi_{2#eta-})");
      h2EventPlane->Draw();
      c1->SaveAs("h2EventPlane.png");
   }

   double centralityNchg[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
   double resolution[3] = {0, 0, 0};
   centralityBound(centralityNchg, h2NcgNppPhi1Bin->ProjectionY());
   TFile *hfile = TFile::Open("outputFile.root", "RECREATE");
   h2EventPlane->Write();
   calcResolution(h2EventPlane,centralityNchg,resolution);
   h2NcgNppPhi1Bin->Write();
   h2NcgNppPhi1Bin->Draw("COLZ");
   c1->SaveAs("h2NcgNppPhi1Bin.png");
   for (int i = 0; i < 6; i++)
      {
         if (i < 3)
         {
            if (i < 2)
            {
               h2NcgNppPhi2Bin[i]->Write();
               h2NcgNppPhi2Bin[i]->Draw("COLZ");
               c1->SaveAs(Form(outputplotsfolder + "h2NcgNppPhi2Bin%d.png", i));
            }
            h2NcgNppPhi3Bin[i]->Write();
            h2NcgNppPhi3Bin[i]->Draw("COLZ");
            c1->SaveAs(Form(outputplotsfolder + "h2NcgNppPhi3Bin%d.png", i));
         }
         h2NcgNppPhi6Bin[i]->Write();
         h2NcgNppPhi6Bin[i]->Draw("COLZ");
         c1->SaveAs(Form(outputplotsfolder + "h2NcgNppPhi6Bin%d.png", i));
      }

   ofstream logOut;
   logOut.open(outputplotsfolder + "centralityNchg.log", ios_base::app);
   logOut << "centralityNchg.log" <<endl;
   for (int l = 0; l < 9; l++)
   {
      if (l == 0)
      {
         logOut << 5 << "%-10%:" << centralityNchg[l] << endl;
      }
      else
      {
         logOut << l * 10 << "%" << centralityNchg[l] << endl;
      }
   }

      logOut << "resolution.log" <<endl;
   // for (int l = 0; l < 9; l++)
   // {
   //    if (l == 0)
   //    {
   //       logOut << 5 << "%-10%:" << resolution[l] << endl;
   //    }
   //    else
   //    {
   //       logOut << l * 10 << "%" << resolution[l] << endl;
   //    }
   // }

      for (int l = 0; l < 3; l++)
   {
      if (l == 0)
      {
         logOut << "0-10%:" << resolution[l] << endl;
      }
      else if (l == 1)
      {
         logOut << "10-60%" << resolution[l] << endl;
      }
      else if (l == 2)
      {
         logOut << "60-80%" << resolution[l] << endl;
      }
   }

   // Loop on every angle

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

   // 6 Bins
   for (int i = 0; i < 6; i++)
   {
      double C1[3] = {0}, C2[3] = {0}, C3[3] = {0}, C4[3] = {0};
      double C2oC1[3] = {0}, C3oC2[3] = {0}, C4oC2[3] = {0};
      getCumulants(h2NcgNppPhi6Bin[i], C1, C2, C3, C4, centralityNchg);
      for (int l = 0; l < 3; l++)
      {
         C2oC1[l] = C2[l] / C1[l];
         C3oC2[l] = C3[l] / C2[l];
         C4oC2[l] = C4[l] / C2[l];
      }
      for (int j = 0; j < 3; j++)
      {
         C2oC1VsPhi6Bin[j][i] = C2oC1[j];
         C3oC2VsPhi6Bin[j][i] = C3oC2[j];
         C4oC2VsPhi6Bin[j][i] = C4oC2[j];
      }
   }

   // 3 Bins
   for (int i = 0; i < 3;i++)
   {
      double C1[3] = {0}, C2[3] = {0}, C3[3] = {0}, C4[3] = {0};
      double C2oC1[3] = {0}, C3oC2[3] = {0}, C4oC2[3] = {0};
      getCumulants(h2NcgNppPhi3Bin[i], C1, C2, C3, C4, centralityNchg);
      for (int l = 0; l < 3; l++)
      {
         C2oC1[l] = C2[l] / C1[l];
         C3oC2[l] = C3[l] / C2[l];
         C4oC2[l] = C4[l] / C2[l];
      }
      for (int j = 0; j < 3; j++)
      {
         C2oC1VsPhi3Bin[j][i] = C2oC1[j];
         C3oC2VsPhi3Bin[j][i] = C3oC2[j];
         C4oC2VsPhi3Bin[j][i] = C4oC2[j];
      }
   }

   // 2 Bins
   for (int i = 0; i < 2;i++)
   {
      double C1[3] = {0}, C2[3] = {0}, C3[3] = {0}, C4[3] = {0};
      double C2oC1[3] = {0}, C3oC2[3] = {0}, C4oC2[3] = {0};
      getCumulants(h2NcgNppPhi2Bin[i], C1, C2, C3, C4, centralityNchg);
      for (int l = 0; l < 3; l++)
      {
         C2oC1[l] = C2[l] / C1[l];
         C3oC2[l] = C3[l] / C2[l];
         C4oC2[l] = C4[l] / C2[l];
      }
      for (int j = 0; j < 3; j++)
      {
         C2oC1VsPhi2Bin[j][i] = C2oC1[j];
         C3oC2VsPhi2Bin[j][i] = C3oC2[j];
         C4oC2VsPhi2Bin[j][i] = C4oC2[j];
      }
   }

   double C1[3] = {0}, C2[3] = {0}, C3[3] = {0}, C4[3] = {0};
   double C2oC1[3] = {0}, C3oC2[3] = {0}, C4oC2[3] = {0};
   getCumulants(h2NcgNppPhi1Bin, C1, C2, C3, C4, centralityNchg);
   for (int l = 0; l < 3; l++)
   {
      C2oC1[l] = C2[l] / C1[l];
      C3oC2[l] = C3[l] / C2[l];
      C4oC2[l] = C4[l] / C2[l];
   }
   for (int j = 0; j < 3; j++)
   {
      C2oC1VsPhi1Bin[j] = C2oC1[j];
      C3oC2VsPhi1Bin[j] = C3oC2[j];
      C4oC2VsPhi1Bin[j] = C4oC2[j];
   }

   for (int i = 0; i < 3; i++)
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
      mg->SetTitle("C_{2}/C_{1}");
      if(reconstructionRP)
      {
         mg->GetXaxis()->SetTitle("#phi-#psi_{EP}");
      }
      {
         mg->GetXaxis()->SetTitle("#phi-#psi_{RP}");
      }
      
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
      mg->SetTitle("C_{3}/C_{2}");
      if(reconstructionRP)
      {
         mg->GetXaxis()->SetTitle("#phi-#psi_{EP}");
      }
      {
         mg->GetXaxis()->SetTitle("#phi-#psi_{RP}");
      }
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
      mg->SetTitle("C_{4}/C_{2}");
      if(reconstructionRP)
      {
         mg->GetXaxis()->SetTitle("#phi-#psi_{EP}");
      }
      {
         mg->GetXaxis()->SetTitle("#phi-#psi_{RP}");
      }
      mg->Draw("zpA");
      c1->Size(0, 0);
      c1->SaveAs(Form(outputplotsfolder + "%d/h1C4oC2VsPhiCentral.png", i));
      c1->SetTitle(Form("h1vVsPhiCentral%d", i));
   }

   c1->SetLogy(1);
   h2NcgNppPhi1Bin->ProjectionX()->DrawNormalized("E");
   h2NcgNppPhi1Bin->SetTitle("Net-proton production vs. N_{chg}");
   c1->SaveAs(outputplotsfolder + "netPro.pdf");
   c1->SaveAs(outputplotsfolder + "netPro.png");

   h2NcgNppPhi1Bin->ProjectionY()->DrawNormalized("E");
   c1->SaveAs(outputplotsfolder + "refmult3.pdf");
   c1->SaveAs(outputplotsfolder + "refmult3.png");

   TH1 *cumul = h2NcgNppPhi1Bin->ProjectionY()->GetCumulative(false);
   cumul->SetTitle("Cumulative Charged Particle Multiplicity");
   cumul->Scale(1.0 / h2NcgNppPhi1Bin->ProjectionY()->Integral());
   cumul->Draw("E");
   c1->SaveAs(outputplotsfolder + "refmult3Cumu.pdf");
   c1->SaveAs(outputplotsfolder + "refmult3Cumu.png");

   c1->SetLogz(1);
   c1->SetLogy(0);
   h2NcgNppPhi1Bin->Draw("COLZ");
   c1->SaveAs(outputplotsfolder + "h2NcgNppPhi1Bin.pdf");
   c1->SaveAs(outputplotsfolder + "h2NcgNppPhi1Bin.png");

   int xbins = h2NcgNppPhi1Bin->GetXaxis()->GetNbins();
   int ybins = h2NcgNppPhi1Bin->GetYaxis()->GetNbins();

   // drawEvent()->Write();
   // hfile->Write();
   hfile->Close();
}
