#include <cmath>
#include <fstream>
#include <iostream>

#include "TF1.h"
#include "TH1D.h"
#include "TH3F.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TLegend.h"

using namespace std;

void fillHistograms(TH1D *h1NppPhi1Bin, TH1D **h1NppPhi2Bin, TH1D **h1NppPhi3Bin, TH1D **h1NppPhi6Bin)
{
    TRandom *rndgen = new TRandom(time(0));
    Long64_t nentries = 100000;
    cout << "Total Entries:" << nentries << endl;
    int eventCounter = 0;

    for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        int primarieCounter = 0;
        if (jentry % 50000 == 0)
        {
            cout << "Progress %:" << jentry * 1.0 / nentries * 100 << endl;
        }

        // Have known the Imp, fill the h2 with the phi of every particle
        int _nNetProtonPhi[6] = {0, 0, 0, 0, 0, 0};
        //   int NetProtonNum from distribution
        //   The phi of every proton also from distribution;
        int nProton = rndgen->Gaus(15,4);
        int nAntiProton = rndgen->Gaus(3,1);

        TF1 *phiDist=new TF1("phi","1+[0]*TMath::Cos(2*x)",0,2*TMath::Pi());
        phiDist->SetParameter(0,0);
        for (Int_t j = 0; j < nProton; j++)
        {
            //   get the phi of this proton
            int index = 0;
            index = (int)((phiDist->GetRandom()) / (TMath::Pi() * 1.0 / 12)) % 6;
            _nNetProtonPhi[index]++;
        }

        for (Int_t j = 0; j < nAntiProton; j++)
        {
            //   get the phi of this proton
            int index = 0;
            index = (int)((phiDist->GetRandom()) / (TMath::Pi() * 1.0 / 12)) % 6;
            _nNetProtonPhi[index]--;
        }

        // Fill the histograms
        for (int j = 0; j < 6; j++)
        {
            if (j < 3)
            {
                h1NppPhi3Bin[j]->Fill(_nNetProtonPhi[2 * j] + _nNetProtonPhi[2 * j + 1]);
                if (j < 2)
                {
                    h1NppPhi2Bin[j]->Fill(_nNetProtonPhi[3 * j] + _nNetProtonPhi[3 * j + 1] + _nNetProtonPhi[3 * j + 2]);
                }
            }
            h1NppPhi6Bin[j]->Fill(_nNetProtonPhi[j]);
        }
        h1NppPhi1Bin->Fill(_nNetProtonPhi[0] + _nNetProtonPhi[1] + _nNetProtonPhi[2] + _nNetProtonPhi[3] + _nNetProtonPhi[4] + _nNetProtonPhi[5]);
        eventCounter++;
        // cout << _nNetProtonPhi[i][0] + _nNetProtonPhi[i][1] + _nNetProtonPhi[i][2] + _nNetProtonPhi[i][3] + _nNetProtonPhi[i][4] + _nNetProtonPhi[i][5]  << endl;
    }
    cout << eventCounter << " events processed" << endl;
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

void getCumulants(TH1D *h1Npp, double *C2oC1, double *C3oC2, double *C4oC2, double *eC2oC1, double *eC3oC2, double *eC4oC2)
{
    double C1, C2, C3, C4;
    int nEvents = h1Npp->Integral();

    float _mean = 0, _sigma = 0, _fS = 0, _kappa = 0;
    _mean = h1Npp->GetMean();
    _sigma = h1Npp->GetStdDev();
    _fS = h1Npp->GetSkewness();

    if (isnan(_fS))
    {
        _fS = 0;
        cout << "omitEvent: " << nEvents << endl;
        // return;
    }
    _kappa = h1Npp->GetKurtosis();
    if (isnan(_kappa))
    {
        _kappa = 0;
        cout << "omitEvent: " << nEvents << endl;
        // return;
    }

    cout << "_mean" << _mean << endl;
    cout << "_sigma" << _sigma << endl;
    cout << "_fS" << _fS << endl;
    cout << "_kappa" << _kappa << endl;

    C1 = _mean;
    C2 = _sigma * _sigma;
    C3 = _fS * _sigma * _sigma * _sigma;
    C4 = _kappa * _sigma * _sigma * _sigma * _sigma;

    cout << "C1" << C1 << endl;
    cout << "C2" << C2 << endl;
    cout << "C3" << C3 << endl;
    cout << "C4" << C4 << endl;

    double _moments[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    getMoments(h1Npp, _moments);

    *C2oC1 = C2 / C1;
    *C3oC2 = C3 / C2;
    *C4oC2 = C4 / C2;
    *eC2oC1 = TMath::Sqrt(varC2oC1(_moments, _mean, nEvents));
    *eC3oC2 = TMath::Sqrt(varC3oC2(_moments, _mean, nEvents));
    *eC4oC2 = TMath::Sqrt(varC4oC2(_moments, _mean, nEvents));
}

void drawFigures(TH1D *h1NppPhi1Bin, TH1D **h1NppPhi2Bin, TH1D **h1NppPhi3Bin, TH1D **h1NppPhi6Bin, TString outputplotsfolder)
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

    double C2oC1VsPhi6Bin[6];
    double C3oC2VsPhi6Bin[6];
    double C4oC2VsPhi6Bin[6];

    double C2oC1VsPhi3Bin[3];
    double C3oC2VsPhi3Bin[3];
    double C4oC2VsPhi3Bin[3];

    double C2oC1VsPhi2Bin[2];
    double C3oC2VsPhi2Bin[2];
    double C4oC2VsPhi2Bin[2];

    double C2oC1VsPhi1Bin;
    double C3oC2VsPhi1Bin;
    double C4oC2VsPhi1Bin;

    double eC2oC1VsPhi6Bin[6];
    double eC3oC2VsPhi6Bin[6];
    double eC4oC2VsPhi6Bin[6];

    double eC2oC1VsPhi3Bin[3];
    double eC3oC2VsPhi3Bin[3];
    double eC4oC2VsPhi3Bin[3];

    double eC2oC1VsPhi2Bin[2];
    double eC3oC2VsPhi2Bin[2];
    double eC4oC2VsPhi2Bin[2];

    double eC2oC1VsPhi1Bin;
    double eC3oC2VsPhi1Bin;
    double eC4oC2VsPhi1Bin;

    // 6 Bins
    for (int i = 0; i < 6; i++)
    {
        double C2oC1, C3oC2, C4oC2;
        double eC2oC1, eC3oC2, eC4oC2;
        getCumulants(h1NppPhi6Bin[i], &C2oC1, &C3oC2, &C4oC2, &eC2oC1, &eC3oC2, &eC4oC2);

        C2oC1VsPhi6Bin[i] = C2oC1;
        C3oC2VsPhi6Bin[i] = C3oC2;
        C4oC2VsPhi6Bin[i] = C4oC2;
        eC2oC1VsPhi6Bin[i] = eC2oC1;
        eC3oC2VsPhi6Bin[i] = eC3oC2;
        eC4oC2VsPhi6Bin[i] = eC4oC2;
        cout << C2oC1 << " " << C3oC2 << " " << C4oC2 << endl;
        cout <<"error" << eC2oC1 << " " << eC3oC2 << " " << eC4oC2 << endl;
    }

    // 3 Bins
    for (int i = 0; i < 3; i++)
    {
        double C2oC1, C3oC2, C4oC2;
        double eC2oC1, eC3oC2, eC4oC2;
        getCumulants(h1NppPhi3Bin[i], &C2oC1, &C3oC2, &C4oC2, &eC2oC1, &eC3oC2, &eC4oC2);

        C2oC1VsPhi3Bin[i] = C2oC1;
        C3oC2VsPhi3Bin[i] = C3oC2;
        C4oC2VsPhi3Bin[i] = C4oC2;
        eC2oC1VsPhi3Bin[i] = eC2oC1;
        eC3oC2VsPhi3Bin[i] = eC3oC2;
        eC4oC2VsPhi3Bin[i] = eC4oC2;
                cout << C2oC1 << " " << C3oC2 << " " << C4oC2 << endl;
        cout <<"error" << eC2oC1 << " " << eC3oC2 << " " << eC4oC2 << endl;
    }

    // 2 Bins
    for (int i = 0; i < 2; i++)
    {
        double C2oC1, C3oC2, C4oC2;
        double eC2oC1, eC3oC2, eC4oC2;
        getCumulants(h1NppPhi2Bin[i], &C2oC1, &C3oC2, &C4oC2, &eC2oC1, &eC3oC2, &eC4oC2);

        C2oC1VsPhi2Bin[i] = C2oC1;
        C3oC2VsPhi2Bin[i] = C3oC2;
        C4oC2VsPhi2Bin[i] = C4oC2;
        eC2oC1VsPhi2Bin[i] = eC2oC1;
        eC3oC2VsPhi2Bin[i] = eC3oC2;
        eC4oC2VsPhi2Bin[i] = eC4oC2;
        cout << C2oC1 << " " << C3oC2 << " " << C4oC2 << endl;
        cout <<"error" << eC2oC1 << " " << eC3oC2 << " " << eC4oC2 << endl;
    }

    double C2oC1, C3oC2, C4oC2;
    double eC2oC1, eC3oC2, eC4oC2;
    getCumulants(h1NppPhi1Bin, &C2oC1, &C3oC2, &C4oC2, &eC2oC1, &eC3oC2, &eC4oC2);
    C2oC1VsPhi1Bin = C2oC1;
    C3oC2VsPhi1Bin = C3oC2;
    C4oC2VsPhi1Bin = C4oC2;
    eC2oC1VsPhi1Bin = eC2oC1;
    eC3oC2VsPhi1Bin = eC3oC2;
    eC4oC2VsPhi1Bin = eC4oC2;
            cout << C2oC1 << " " << C3oC2 << " " << C4oC2 << endl;
        cout <<"error" << eC2oC1 << " " << eC3oC2 << " " << eC4oC2 << endl;

    c1->Clear();
    TGraphErrors *h1C2oC1VsPhi1Bin = new TGraphErrors(1, &x1Bin, &C2oC1VsPhi1Bin, &x1err, &eC2oC1VsPhi1Bin);
    h1C2oC1VsPhi1Bin->SetLineColor(2);
    h1C2oC1VsPhi1Bin->SetMarkerSize(1);
    h1C2oC1VsPhi1Bin->SetMarkerStyle(20);
    TGraphErrors *h1C2oC1VsPhi2Bin = new TGraphErrors(2, x2Bin, C2oC1VsPhi2Bin, x2err, eC2oC1VsPhi2Bin);
    h1C2oC1VsPhi2Bin->SetLineColor(3);
    h1C2oC1VsPhi2Bin->SetMarkerSize(1);
    h1C2oC1VsPhi2Bin->SetMarkerStyle(21);
    TGraphErrors *h1C2oC1VsPhi3Bin = new TGraphErrors(3, x3Bin, C2oC1VsPhi3Bin, x3err, eC2oC1VsPhi3Bin);
    h1C2oC1VsPhi3Bin->SetLineColor(4);
    h1C2oC1VsPhi3Bin->SetMarkerSize(1);
    h1C2oC1VsPhi3Bin->SetMarkerStyle(22);
    TGraphErrors *h1C2oC1VsPhi6Bin = new TGraphErrors(6, x6Bin, C2oC1VsPhi6Bin, x6err, eC2oC1VsPhi6Bin);
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
    c1->SaveAs(outputplotsfolder + "h1C2oC1VsPhiCentral.png");
    c1->SetTitle("h1C2oC1VsPhiCentral");

    TGraphErrors *h1C3oC2VsPhi1Bin = new TGraphErrors(1, &x1Bin, &C3oC2VsPhi1Bin, &x1err, &eC3oC2VsPhi1Bin);
    h1C3oC2VsPhi1Bin->SetLineColor(2);
    h1C3oC2VsPhi1Bin->SetMarkerSize(1);
    h1C2oC1VsPhi1Bin->SetMarkerStyle(20);
    TGraphErrors *h1C3oC2VsPhi2Bin = new TGraphErrors(2, x2Bin, C3oC2VsPhi2Bin, x2err, eC3oC2VsPhi2Bin);
    h1C3oC2VsPhi2Bin->SetLineColor(3);
    h1C3oC2VsPhi2Bin->SetMarkerSize(1);
    h1C3oC2VsPhi2Bin->SetMarkerStyle(21);
    TGraphErrors *h1C3oC2VsPhi3Bin = new TGraphErrors(3, x3Bin, C3oC2VsPhi3Bin, x3err, eC3oC2VsPhi3Bin);
    h1C3oC2VsPhi3Bin->SetLineColor(4);
    h1C3oC2VsPhi3Bin->SetMarkerSize(1);
    h1C3oC2VsPhi3Bin->SetMarkerStyle(22);
    TGraphErrors *h1C3oC2VsPhi6Bin = new TGraphErrors(6, x6Bin, C3oC2VsPhi6Bin, x6err, eC3oC2VsPhi6Bin);
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
    c1->SaveAs(outputplotsfolder + "h1C3oC2VsPhiCentral.png");
    c1->SetTitle("h1vVsPhiCentral");

    c1->Clear();
    TGraphErrors *h1C4oC2VsPhi1Bin = new TGraphErrors(1, &x1Bin, &C4oC2VsPhi1Bin, &x1err, &eC4oC2VsPhi1Bin);

    h1C4oC2VsPhi1Bin->SetLineColor(2);
    h1C4oC2VsPhi1Bin->SetMarkerSize(1);
    h1C4oC2VsPhi1Bin->SetMarkerStyle(20);
    TGraphErrors *h1C4oC2VsPhi2Bin = new TGraphErrors(2, x2Bin, C4oC2VsPhi2Bin, x2err, eC4oC2VsPhi2Bin);
    h1C4oC2VsPhi2Bin->SetLineColor(3);
    h1C4oC2VsPhi2Bin->SetMarkerSize(1);
    h1C4oC2VsPhi2Bin->SetMarkerStyle(21);
    TGraphErrors *h1C4oC2VsPhi3Bin = new TGraphErrors(3, x3Bin, C4oC2VsPhi3Bin, x3err, eC4oC2VsPhi3Bin);
    h1C4oC2VsPhi3Bin->SetLineColor(4);
    h1C4oC2VsPhi3Bin->SetMarkerSize(1);
    h1C4oC2VsPhi3Bin->SetMarkerStyle(22);
    TGraphErrors *h1C4oC2VsPhi6Bin = new TGraphErrors(6, x6Bin, C4oC2VsPhi6Bin, x6err, eC4oC2VsPhi6Bin);
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
    c1->SaveAs(outputplotsfolder + "h1C4oC2VsPhiCentral.png");
    c1->SetTitle("h1vVsPhiCentral");

    // TGraphErrors *h1C4oC2VsPhiNBin;
    // auto legend = new TLegend(0.6, 0.6, 0.9, 0.9);
    // legend->SetHeader("The Legend Title", "C");
    // double nBins[4] = {1, 2, 3, 6};
    // double enBins[4] = {0.75, 0.75, 0.75, 0.75};
    // double yNBins[4] = {C4oC2VsPhi1Bin, 0.5 * (C4oC2VsPhi2Bin[0] + C4oC2VsPhi2Bin[1]), 1.0 / 3.0 * (C4oC2VsPhi3Bin[0] + C4oC2VsPhi3Bin[1] + C4oC2VsPhi3Bin[2]), 1.0 / 6.0 * (C4oC2VsPhi6Bin[0] + C4oC2VsPhi6Bin[1] + C4oC2VsPhi6Bin[2] + C4oC2VsPhi6Bin[3] + C4oC2VsPhi6Bin[4] + C4oC2VsPhi6Bin[5])};
    // double eyNBins[4] = {eC4oC2VsPhi1Bin, 0.5 * TMath::Sqrt(eC4oC2VsPhi2Bin[0] * eC4oC2VsPhi2Bin[0] + eC4oC2VsPhi2Bin[1] * eC4oC2VsPhi2Bin[1]), 1.0 / 3.0 * TMath::Sqrt(eC4oC2VsPhi3Bin[0] * eC4oC2VsPhi3Bin[0] + eC4oC2VsPhi3Bin[1] * eC4oC2VsPhi3Bin[1] + eC4oC2VsPhi3Bin[2] * eC4oC2VsPhi3Bin[2]), 1.0 / 6.0 * TMath::Sqrt(eC4oC2VsPhi6Bin[0] * eC4oC2VsPhi6Bin[0] + eC4oC2VsPhi6Bin[1] * eC4oC2VsPhi6Bin[1] + eC4oC2VsPhi6Bin[2] * eC4oC2VsPhi6Bin[2] + eC4oC2VsPhi6Bin[3] * eC4oC2VsPhi6Bin[3] + eC4oC2VsPhi6Bin[4] * eC4oC2VsPhi6Bin[4] + eC4oC2VsPhi6Bin[5] * eC4oC2VsPhi6Bin[5])};
    // h1C4oC2VsPhiNBin = new TGraphErrors(4, nBins, yNBins, enBins, eyNBins);
    // h1C4oC2VsPhiNBin->SetLineWidth(2);
    // h1C4oC2VsPhiNBin->SetLineColor(2);
    // h1C4oC2VsPhiNBin->SetMarkerSize(1);
    // h1C4oC2VsPhiNBin->SetMarkerStyle(20);
    // legend->AddEntry(h1C4oC2VsPhiNBin, "Toy MC", "ep");
    // mg = new TMultiGraph();
    // mg->Add(h1C4oC2VsPhiNBin);
    // mg->GetYaxis()->SetTitle("C_{4}/C_{2} Vs NBins");
    // mg->GetXaxis()->SetTitle("Number of Bins");
    // mg->SetTitle("");
    // mg->Draw("zpA");
    // legend->Draw();
    // c1->Size(0, 0);
    // c1->SaveAs(outputplotsfolder + "h1C4oC2VsNBins.png");

    c1->SetLogy(1);
    h1NppPhi1Bin->DrawNormalized("E");
    h1NppPhi1Bin->SetTitle("Net-proton production vs. N_{chg}");
    c1->SaveAs(outputplotsfolder + "netPro.pdf");
    c1->SaveAs(outputplotsfolder + "netPro.png");

    c1->SetLogz(1);
    c1->SetLogy(0);
    h1NppPhi1Bin->Draw("COLZ");
    c1->SaveAs(outputplotsfolder + "h1NppPhi1Bin.pdf");
    c1->SaveAs(outputplotsfolder + "h1NppPhi1Bin.png");
    c1->Close();
}

void toyMC(TString inputFile = "primaryOnly.root", TString outputFile = "primaryOutput.root", TString outputplotsfolder = "primaryOnly/")
{
    system("rm -rf " + outputplotsfolder);
    system("mkdir -p " + outputplotsfolder);

    TCanvas *c1 = new TCanvas();
    c1->cd();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0101);

    TFile f(inputFile);

    TH1D *h1NppPhi6Bin[6];
    TH1D *h1NppPhi3Bin[3];
    TH1D *h1NppPhi2Bin[2];
    TH1D *h1NppPhi1Bin;

    f.GetObject("h1NppPhi1Bin", h1NppPhi1Bin);
    for (int i = 0; i < 6; i++)
    {
        if (i < 3)
        {
            if (i < 2)
            {
                f.GetObject(Form("h1NppPhi2Bin%d", i), h1NppPhi2Bin[i]);
            }
            f.GetObject(Form("h1NppPhi3Bin%d", i), h1NppPhi3Bin[i]);
        }
        f.GetObject(Form("h1NppPhi6Bin%d", i), h1NppPhi6Bin[i]);
    }

    if (!h1NppPhi1Bin)
    {
        h1NppPhi1Bin = new TH1D("h1NppPhi1Bin", "h1NppPhi1Bin", 100, -5.5, 94.5);
        h1NppPhi1Bin->GetYaxis()->SetTitle("Impact Parameter");
        h1NppPhi1Bin->GetXaxis()->SetTitle("Net-Proton");
        for (int i = 0; i < 6; i++)
        {
            if (i < 3)
            {
                if (i < 2)
                {
                    h1NppPhi2Bin[i] = new TH1D(Form("h1NppPhi2Bin%d", i), Form("h1NppPhi2Bin%d", i), 100, -5.5, 94.5);
                    h1NppPhi2Bin[i]->GetYaxis()->SetTitle("Impact Parameter");
                    h1NppPhi2Bin[i]->GetXaxis()->SetTitle("Net-Proton");
                    h1NppPhi2Bin[i]->GetZaxis()->SetTitle("#phi");
                }
                h1NppPhi3Bin[i] = new TH1D(Form("h1NppPhi3Bin%d", i), Form("h1NppPhi3Bin%d", i), 100, -5.5, 94.5);
                h1NppPhi3Bin[i]->GetYaxis()->SetTitle("Impact Parameter");
                h1NppPhi3Bin[i]->GetXaxis()->SetTitle("Net-Proton");
                h1NppPhi3Bin[i]->GetZaxis()->SetTitle("#phi");
            }
            h1NppPhi6Bin[i] = new TH1D(Form("h1NppPhi6Bin%d", i), Form("h1NppPhi6Bin%d", i), 100, -5.5, 94.5);
            h1NppPhi6Bin[i]->GetYaxis()->SetTitle("Impact Parameter");
            h1NppPhi6Bin[i]->GetXaxis()->SetTitle("Net-Proton");
            h1NppPhi6Bin[i]->GetZaxis()->SetTitle("#phi");
        }
        fillHistograms(h1NppPhi1Bin, h1NppPhi2Bin, h1NppPhi3Bin, h1NppPhi6Bin);
    }
    TFile *hfile = TFile::Open(outputFile, "RECREATE");

    h1NppPhi1Bin->Write();
    h1NppPhi1Bin->Draw("COLZ");
    c1->SaveAs(outputplotsfolder + "h1NppPhi1Bin.png");
    for (int j = 0; j < 6; j++)
    {
        if (j < 3)
        {
            if (j < 2)
            {
                h1NppPhi2Bin[j]->Write();
                h1NppPhi2Bin[j]->Draw("COLZ");
                c1->SaveAs(Form(outputplotsfolder + "h1NppPhi2Bin%d.png", j));
            }
            h1NppPhi3Bin[j]->Write();
            h1NppPhi3Bin[j]->Draw("COLZ");
            c1->SaveAs(Form(outputplotsfolder + "h1NppPhi3Bin%d.png", j));
        }
        h1NppPhi6Bin[j]->Write();
        h1NppPhi6Bin[j]->Draw("COLZ");
        c1->SaveAs(Form(outputplotsfolder + "h1NppPhi6Bin%d.png", j));
    }
    // Loop on every angle
    hfile->Close();
    c1->Close();
    drawFigures(h1NppPhi1Bin, h1NppPhi2Bin, h1NppPhi3Bin, h1NppPhi6Bin, outputplotsfolder);
}
