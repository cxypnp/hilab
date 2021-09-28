#include <cmath>
#include <fstream>
#include <iostream>

#include "TH2I.h"
#include "TH1I.h"
#include "TH3F.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TMath.h"
#include "f14.h"

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

void fillHistograms(TH2I *h2NcgNpp, TH3F *h3NcgNppPhi)
{
   TTree *event = NULL;

   f14 reader(event);

   Long64_t nentries = reader.fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry = 0; jentry < nentries; jentry++)
   {
      if (jentry % 10000 == 0)
      {
         cout << jentry * 1.0 / nentries * 100 << endl;
      }
      Long64_t ientry = reader.LoadTree(jentry);
      if (ientry < 0)
         break;
      int _nNetProton = 0;
      int _nCharged = 0;
      nb = reader.fChain->GetEntry(jentry);
      nbytes += nb;

      // if (Cut(ientry) < 0) continue;
      if (reader.fColHdr_timestep != 200 || reader.fColHdr_Ntrack == 0)
      {
         continue;
      }

      // Read a Event
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
                  if (reader.fTracksOut_ityp[j] == 1)
                     _nNetProton++;
                  else if (reader.fTracksOut_ityp[j] == -1)
                     _nNetProton--;
               }
            }
            else
            {
               // Other Particle
               float pT = TMath::Sqrt(reader.fTracksOut_px[j] * reader.fTracksOut_px[j] + reader.fTracksOut_py[j] * reader.fTracksOut_py[j]);
               float p = TMath::Sqrt(reader.fTracksOut_px[j] * reader.fTracksOut_px[j] + reader.fTracksOut_py[j] * reader.fTracksOut_py[j] + reader.fTracksOut_pz[j] * reader.fTracksOut_pz[j]);
               float eta = 0.5 * log((p + reader.fTracksOut_pz[j]) / (p - reader.fTracksOut_pz[j]));
               if (TMath::Abs(eta) < 1 && pT > 0.1)
               {
                  _nCharged++;
               }
            }
         }
      }
      h2NcgNpp->Fill(_nNetProton, _nCharged);
      // Have known the Ncg and Npp, fill the h3 with the phi of every particle

      // _nNetProtonPhi[6]={0,0,0,0,0,0};
      // for (Int_t j = 0; j < reader.fColHdr_Ntrack; j++)
      // {
      //    // Charged Particle
      //    if ((reader.fTracksOut_chg[j] != 0))
      //    {
      //       // Is it (anti-)proton?
      //       if ((reader.fTracksOut_ityp[j] == 1) || (reader.fTracksOut_ityp[j] == -1))
      //       {
      //          // (anti-)proton
      //          float pT = TMath::Sqrt(reader.fTracksOut_px[j] * reader.fTracksOut_px[j] + reader.fTracksOut_py[j] * reader.fTracksOut_py[j]);
      //          float y = 0.5 * log((reader.fTracksOut_p0[j] + reader.fTracksOut_pz[j]) / (reader.fTracksOut_p0[j] - reader.fTracksOut_pz[j]));
      //          if (TMath::Abs(y) < 0.5 && pT < 2.0 && pT > 0.4)
      //          {
      //             int index = getPhi(reader.fTracksOut_rx,reader.fTracksOut_ry) / deltaphi;
      //             cout << index << endl;
      //             if (reader.fTracksOut_ityp[j] == 1)
      //                _nNetProtonPhi[index]++;
      //             else if (reader.fTracksOut_ityp[j] == -1)
      //                _nNetProtonPhi[index]--;
      //          }
      //       }
      //       else
      //       {
      //          // // Other Particle
      //       }
      //    }
      // }
      // for (int i = 0; i < 6; i++)
      // {
      //    h3NcgNppPhi->Fill(_nNetProton, _nCharged);
      // }
   }
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

void netProton(TString outputplotsfolder = "outputFolder/")
{
   system("rm -rf " + outputplotsfolder);
   system("mkdir -p " + outputplotsfolder);
   ofstream logOut;

   TCanvas *c1 = new TCanvas();
   c1->cd();
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0101);

   logOut.open(outputplotsfolder + "run.log", ios_base::app);

   TFile f("outputFile.root");
   TH2I *h2NcgNpp = NULL;
   TH3F *h3NcgNppPhi = NULL;
   f.GetObject("h2NcgNpp", h2NcgNpp);

   if (!h2NcgNpp)
   {
      h2NcgNpp = new TH2I("h2NcgNpp", "h2NcgNpp", 100, -5.5, 94.5, 450, -0.5, 449.5);
   h2NcgNpp->GetYaxis()->SetTitle("Charged Multiplicity");
   h2NcgNpp->GetXaxis()->SetTitle("Net-Proton");

   h3NcgNppPhi = new TH3F("h3NcgNppPhi", "h3NcgNppPhi", 100, -5.5, 94.5, 450, -0.5, 449.5, 100, 0, TMath::Pi()*0.5);
   h3NcgNppPhi->GetYaxis()->SetTitle("Charged Multiplicity");
   h3NcgNppPhi->GetXaxis()->SetTitle("Net-Proton");
   h3NcgNppPhi->GetZaxis()->SetTitle("#Phi");
      fillHistograms(h2NcgNpp,h3NcgNppPhi);
   }

   double centraclityNchg[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
   centralityBound(centraclityNchg, h2NcgNpp->ProjectionY());

   for (int l = 0; l < 9; l++)
   {
      if (l == 0)
      {
         logOut << 5 << "%-10%:" << centraclityNchg[l] << endl;
      }
      else
      {
         logOut << l * 10 << "%" << centraclityNchg[l] << endl;
      }
   }

   // Loop on every centrality bin and do the modification
   int binMax, binMin;
   double mean[9] = {0}, sigma[9] = {0}, fS[9] = {0}, kappa[9] = {0};
   int sum;
   double C1[9] = {0}, C2[9] = {0}, C3[9] = {0}, C4[9] = {0};
   double C2oC1[9] = {0}, C3oC2[9] = {0}, C4oC2[9] = {0};
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

      // {
      //    cout << "This" << l << endl;
      //    cout << " " << sum << endl;
      //    cout << "end" << endl;
      // }

      mean[l] = mean[l] / sum;
      sigma[l] = sigma[l] / sum;
      fS[l] = fS[l] / sum;
      kappa[l] = kappa[l] / sum;

      // cout << "This" << l << endl;
      // cout << mean[l]  << endl;
      // cout << sigma[l]  << endl;
      // cout << fS[l]  << endl;
      // cout << kappa[l]  << endl;

      C1[l] = mean[l];
      C2[l] = sigma[l] * sigma[l];
      C3[l] = fS[l] * C2[l] * sigma[l];
      C4[l] = kappa[l] * sigma[l] * sigma[l] * sigma[l] * sigma[l];
      C2oC1[l] = C2[l] / C1[l];
      C3oC2[l] = C3[l] / C2[l];
      C4oC2[l] = C4[l] / C2[l];
   }

   TGraph *gMean = new TGraph(9, centraclityNchg, mean);
   gMean->Draw("AC*");
   gMean->SetTitle("gMean");
   c1->SaveAs(outputplotsfolder + "gMean.png");
   TGraph *gSigma = new TGraph(9, centraclityNchg, sigma);
   gSigma->Draw("AC*");
   gSigma->SetTitle("gSigma");
   c1->SaveAs(outputplotsfolder + "gSigma.png");
   TGraph *gS = new TGraph(9, centraclityNchg, fS);
   gS->Draw("AC*");
   gS->SetTitle("gS");
   c1->SaveAs(outputplotsfolder + "gS.png");
   TGraph *gKappa = new TGraph(9, centraclityNchg, kappa);
   gKappa->Draw("AC*");
   gKappa->SetTitle("gKappa");
   c1->SaveAs(outputplotsfolder + "gKappa.png");

   TGraph *gC1 = new TGraph(9, centraclityNchg, C1);
   gC1->Draw("AC*");
   gC1->SetTitle("C1");
   c1->SaveAs(outputplotsfolder + "gC1.png");
   TGraph *gC2 = new TGraph(9, centraclityNchg, C2);
   gC2->Draw("AC*");
   gC2->SetTitle("C2");
   c1->SaveAs(outputplotsfolder + "gC2.png");
   TGraph *gC3 = new TGraph(9, centraclityNchg, C3);
   gC3->Draw("AC*");
   gC3->SetTitle("C3");
   c1->SaveAs(outputplotsfolder + "gC3.png");
   TGraph *gC4 = new TGraph(9, centraclityNchg, C4);
   gC4->Draw("AC*");
   gC4->SetTitle("C4");
   c1->SaveAs(outputplotsfolder + "gC4.png");

   TGraph *gC2oC1 = new TGraph(9, centraclityNchg, C2oC1);
   gC2oC1->Draw("AC*");
   gC2oC1->SetTitle("C2oC1");
   c1->SaveAs(outputplotsfolder + "gC2oC1.png");
   TGraph *gC3oC2 = new TGraph(9, centraclityNchg, C3oC2);
   gC3oC2->Draw("AC*");
   gC3oC2->SetTitle("C3oC2");
   c1->SaveAs(outputplotsfolder + "gC3oC2.png");
   TGraph *gC4oC2 = new TGraph(9, centraclityNchg, C4oC2);
   gC4oC2->Draw("AC*");
   gC4oC2->SetTitle("C4oC2");
   c1->SaveAs(outputplotsfolder + "gC4oC2.png");

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

   TFile *hfile = TFile::Open("outputFile.root", "RECREATE");
   h2NcgNpp->Write();
   h3NcgNppPhi->Write();
   // drawEvent()->Write();
   hfile->Write();
   hfile->Close();
}
