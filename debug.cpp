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

void centralityBoundImp(double *centralityNchg,double *centralityNchgCount, TH1D *refmult3)
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
    for (int i = 0; i < 9; i++)
    {
        if (i==0)
        {
            centralityNchgCount[0]=refmult3->Integral(centralityNchg[0],centralityNchg[0]);
        }
        else
        {
            centralityNchgCount[i]=refmult3->Integral(centralityNchg[i-1],centralityNchg[i]);
        }
    }
}

void centralityBoundNcg(double *centralityNchg,double *centralityNchgCount, TH1D *refmult3)
{
    // This function is for determine the centrality bins
    TH1 *cumul = refmult3->GetCumulative(false);
    cumul->Scale(1.0 / refmult3->Integral());
    for (int k = 0; k < refmult3->GetXaxis()->GetNbins(); k++)
    {
        for (int l = 0; l < 9; l++)
        {
            if (l == 0)
            {
                // 0-5% bin
                if (cumul->GetBinContent(k) > 0.05)
                {
                    centralityNchg[0] = cumul->GetBinLowEdge(k);
                }
            }
            else
            {
                // l*10% bins
                if (cumul->GetBinContent(k) > l * 0.1)
                {
                    centralityNchg[l] = cumul->GetBinLowEdge(k);
                }
            }
        }
    }
    for (int i = 0; i < 9; i++)
    {
        if (i==0)
        {
            centralityNchgCount[0]=refmult3->Integral(centralityNchg[0],refmult3->GetMinimumBin());
        }
        else
        {
            centralityNchgCount[i]=refmult3->Integral(centralityNchg[i],centralityNchg[i-1]);
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

void fillHistograms(TH2F *h2ImpNppPhi1Bin, TH2F *h2NcgNppPhi1Bin, TH2F *ImpVsNcg)
{
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
        if (jentry == 50000)
           break;

        Long64_t ientry = reader.LoadTree(jentry);
        if (ientry < 0)
            break;

        int _nCharged = 0;
        double sumQxP = 0, sumQyP = 0;
        double sumQxM = 0, sumQyM = 0;
        double PsiP;
        double PsiM;
        nb = reader.fChain->GetEntry(jentry);
        nbytes += nb;

        if (reader.fColHdr_timestep != 200 || reader.fColHdr_Ntrack == 0)
        {
            continue;
        }

        // Have known the Imp and Npp, fill the h2 with the phi of every particle
        int _nNetProtonPhi[6] = {0, 0, 0, 0, 0, 0};
        int nNetProton = 0;
        for (Int_t j = 0; j < reader.fColHdr_Ntrack; j++)
        {
            // Skip this particle if it is:
            // 1. Decay Secondaries;
            // 2. Spectators;
            // 3. Non-charged Paricles;
            if (reader.fTracksOut_pptype[j] == 20 || reader.fTracksOut_Nc[j] == 0 || reader.fTracksOut_chg[j] == 0)
            {
                continue;
            }
            double pVector[4] = {reader.fTracksOut_px[j], reader.fTracksOut_py[j], reader.fTracksOut_pz[j], reader.fTracksOut_p0[j]};
            TLorentzVector *particle = new TLorentzVector(pVector);

            // Is it (anti-)proton?
            if (TMath::Abs(reader.fTracksOut_ityp[j]) == 1)
            {
                // (anti-)proton
                if (netProtonCut(particle))
                {
                    // Net-Proton Multiplicity Counter
                    if (reader.fTracksOut_ityp[j] == 1)
                        nNetProton++;
                    else if (reader.fTracksOut_ityp[j] == -1)
                        nNetProton--;
                }
            }
            else if (centralityCut(particle))
            {
                // Not a proton, get the Nchg
                _nCharged++;
            }
            // Delete this particle
            delete particle;
        }

        // Fill the histograms
        h2ImpNppPhi1Bin->Fill(nNetProton, reader.fEvtHdr_imp);
        h2NcgNppPhi1Bin->Fill(nNetProton, _nCharged);
        ImpVsNcg->Fill(reader.fEvtHdr_imp,_nCharged);
        eventCounter++;
    }

    cout << eventCounter << " events processed" << endl;
}

void debug()
{
    TCanvas *c1=new TCanvas();
    TFile *hfile = TFile::Open("debugFile.root", "RECREATE");
    TH2F* h2ImpNppPhi1Bin = new TH2F("h2ImpNppPhi1Bin", "h2ImpNppPhi1Bin", 100, -5.5, 94.5, 130, 0, 13.0);
    h2ImpNppPhi1Bin->GetYaxis()->SetTitle("Impact Parameter");
    h2ImpNppPhi1Bin->GetXaxis()->SetTitle("Net-Proton");

    TH2F* h2NcgNppPhi1Bin = new TH2F("h2NcgNppPhi1Bin", "h2NcgNppPhi1Bin", 100, -5.5, 94.5, 100, -0.5, 99.5);
    h2NcgNppPhi1Bin->GetYaxis()->SetTitle("Charged Multiplicity");
    h2NcgNppPhi1Bin->GetXaxis()->SetTitle("Net-Proton");

    TH2F *ImpVsNcg = new TH2F("h2NcgNppPhi1Bin", "h2NcgNppPhi1Bin", 130, 0, 13.0, 100, -0.5, 99.5);
    h2NcgNppPhi1Bin->GetYaxis()->SetTitle("Charged Multiplicity");
    h2NcgNppPhi1Bin->GetXaxis()->SetTitle("Impact Parameter");

    fillHistograms(h2ImpNppPhi1Bin,h2NcgNppPhi1Bin,ImpVsNcg);
    double centralityNchg[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    double centralityNchgCount[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    double centralityImp[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    double centralityImpCount[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    centralityBoundImp(centralityImp,centralityImpCount, h2ImpNppPhi1Bin->ProjectionY());
    centralityBoundNcg(centralityNchg,centralityNchgCount, h2NcgNppPhi1Bin->ProjectionY());

    ofstream logOut;
    logOut.open("centralityNchg.log", ios_base::app);
   logOut << "centralityBoundImp.log" << endl;
   for (int l = 0; l < 9; l++)
   {
      if (l == 0)
      {
         logOut << 0 << "%-5%:" << centralityImp[l] << endl;
      }
      else
      {
         logOut << l * 10 << "%" << centralityImp[l] << endl;
      }
   }
   logOut << "centralityBoundNcg.log" << endl;
   for (int l = 0; l < 9; l++)
   {
      if (l == 0)
      {
         logOut << 0 << "%-5%:" << centralityNchg[l] << endl;
      }
      else
      {
         logOut << l * 10 << "%" << centralityNchg[l] << endl;
      }
   }

   logOut << "centralityBoundImpCount.log" << endl;
   for (int l = 0; l < 9; l++)
   {
      if (l == 0)
      {
         logOut << 0 << "%-5%:" << centralityImpCount[l] << endl;
      }
      else
      {
         logOut << l * 10 << "%" << centralityImpCount[l] << endl;
      }
   }
   logOut << "centralityBoundNcgCount.log" << endl;
   for (int l = 0; l < 9; l++)
   {
      if (l == 0)
      {
         logOut << 0 << "%-5%:" << centralityNchgCount[l] << endl;
      }
      else
      {
         logOut << l * 10 << "%" << centralityNchgCount[l] << endl;
      }
   }


    h2ImpNppPhi1Bin->Draw("COLZ");
    c1->SaveAs("h2ImpNppPhi1Bin.png");
    h2NcgNppPhi1Bin->Draw("COLZ");
    c1->SaveAs("h2NcgNppPhi1Bin.png");
    ImpVsNcg->Draw("COLZ");
    c1->SaveAs("ImpVsNcg.png");
    hfile->Write();
    hfile->Close();
}