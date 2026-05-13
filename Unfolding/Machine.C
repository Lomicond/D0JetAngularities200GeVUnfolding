#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TSystem.h"
#include "TH1F.h"
#include "TChain.h"
#include "TObject.h"
#include "TClonesArray.h"

#include "TParticle.h"
#include "TDatabasePDG.h"
#include <TLorentzVector.h>

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

using namespace std;

#ifndef __CINT__
#include "TROOT.h"
#include "TFile.h"
#include "TError.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "Riostream.h"
#include <cstdlib>
#include "TH3F.h"
#include "TH2F.h"
#include "THn.h"
#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <RooUnfold.h>
#include "Riostream.h"
#include "TGraph.h"
#include "TStopwatch.h"
#include "TPaveText.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TLatex.h"
#include <vector>
#include "THnSparse.h"
#include "TGraphErrors.h"
#include "TPaletteAxis.h"
#include "../config.h"
#include "../config_hist.h"

#endif

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"


template<typename T>
void DeleteArray(T *array[], int size) {
    for (int i = 0; i < size; ++i) {
        delete array[i];
    }
}

double SetMinD0Pt = 0;

map <Int_t, Int_t> centralityMap = {
        {8, 0}, // 8->   0-10%
        {7, 1}, // 7->  10-20%
        {6, 1}, // 6->  20-30%
        {5, 1}, // 5->  30-40%
        {4, 2}, // 4->  40-50%
        {3, 2}, // 3->  50-60%
        {2, 2}, // 2->  60-70%
        {1, 2}  // 1->  70-80%
};

vector <Double_t> getAxisVector(TAxis *axis) {
    vector <Double_t> axisVector;
    for (Int_t i = 1; i <= axis->GetNbins(); i++) {
        axisVector.push_back(axis->GetBinLowEdge(i));
    }
    axisVector.push_back(axis->GetBinUpEdge(axis->GetNbins()));
    return axisVector;
}

vector <Double_t> getEqualBining(TH1D *hist, const Int_t &nBins) {
    Int_t nBinsHist = hist->GetNbinsX();
    Double_t integral = hist->Integral();
    Double_t sum = 0;
    vector <Double_t> binEdges;
    binEdges.push_back(hist->GetBinLowEdge(1));
    for (Int_t i = 1; i <= nBinsHist; i++) {
        sum += hist->GetBinContent(i);
        if (sum > integral / nBins) {
            binEdges.push_back(hist->GetBinLowEdge(i + 1));
            sum = 0;
        }
    }
    binEdges.push_back(hist->GetBinLowEdge(nBinsHist + 1));
    return binEdges;
}

template <typename T>
TString VecToTString(const std::vector<T>& v, int precision = 3)
{
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(precision);

    for (size_t i = 0; i < v.size(); ++i) {
        ss << v[i];
        if (i + 1 < v.size()) ss << ", ";
    }

    return TString(ss.str());
}

#include "TCanvas.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"

#include <vector>
#include <sstream>
#include <iomanip>

// pomocná funkce: rozdělí binning do více řádků
template <typename T>
std::vector<TString> VecToTStringLines(const std::vector<T>& v,
                                       int precision = 2,
                                       int nPerLine = 6)
{
    std::vector<TString> lines;

    std::ostringstream ss;
    ss << std::fixed << std::setprecision(precision);

    int counter = 0;

    for (size_t i = 0; i < v.size(); ++i) {
        ss.str("");
        ss.clear();
        ss << v[i];

        TString val = ss.str();

        if (counter == 0) lines.push_back(val);
        else              lines.back() += ", " + val;

        counter++;
        if (counter == nPerLine) counter = 0;
    }

    return lines;
}
#include "TCanvas.h"
#include "TLatex.h"
#include "TString.h"
#include "TStyle.h"

#include <vector>
#include <sstream>
#include <iomanip>

static const char* kAngNames[6] = {
    "#lambda_{1}^{1}",
    "#lambda_{1.5}^{1}",
    "#lambda_{2}^{1}",
    "#lambda_{3}^{1}",
    "#lambda_{0.5}^{1}",
    "p_{T}^{D}"
};

static const char* kCentNames[3] = {
    "0-10%",
    "10-40%",
    "40-80%"
};

template <typename T>
std::vector<TString> FormatVectorLines(const std::vector<T>& v, int maxPerLine = 8, int precision = 2)
{
    std::vector<TString> out;
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(precision);

    for (size_t i = 0; i < v.size(); ++i) {
        if (i > 0 && i % maxPerLine == 0) {
            out.push_back(TString(ss.str()));
            ss.str("");
            ss.clear();
            ss << std::fixed << std::setprecision(precision);
        }

        ss << v[i];
        if (i + 1 < v.size() && (i + 1) % maxPerLine != 0) ss << ", ";
    }

    if (!ss.str().empty()) out.push_back(TString(ss.str()));
    return out;
}

template <typename T>
void DrawVectorBlock(TLatex& latex, double x, double& y,
                     const TString& label, const std::vector<T>& v,
                     int maxPerLine = 8, int precision = 2,
                     double lineStep = 0.030, double extraGap = 0.010)
{
    std::vector<TString> lines = FormatVectorLines(v, maxPerLine, precision);

    if (lines.empty()) {
        latex.DrawLatex(x, y, label);
        y -= (lineStep + extraGap);
        return;
    }

    latex.DrawLatex(x, y, Form("%s %s", label.Data(), lines[0].Data()));
    y -= lineStep;

    for (size_t i = 1; i < lines.size(); ++i) {
        latex.DrawLatex(x + 0.04, y, lines[i]);
        y -= lineStep;
    }

    y -= extraGap;
}

bool NearlyEqual(double a, double b, double eps = 1e-9)
{
    return fabs(a - b) < eps;
}
void CheckCacheCompatibility1D(const char* what,
                               int nCache, double cacheMin, double cacheMax,
                               const vector<Double_t>& finalBins)
{
    if (nCache <= 0) {
        cerr << "ERROR: " << what << ": nCache <= 0" << endl;
        gSystem->Exit(1);
    }

    if (finalBins.size() < 2) {
        cerr << "ERROR: " << what << ": final binning has less than 2 edges" << endl;
        gSystem->Exit(1);
    }

    const double step = (cacheMax - cacheMin) / nCache;

    if (!(cacheMin <= finalBins.front() + 1e-12 &&
          cacheMax >= finalBins.back()  - 1e-12)) {
        cerr << "ERROR: " << what << ": cache range ["
             << cacheMin << ", " << cacheMax
             << "] does not cover final range ["
             << finalBins.front() << ", " << finalBins.back() << "]" << endl;
        gSystem->Exit(1);
    }

    for (size_t i = 0; i < finalBins.size(); i++) {
        double x = finalBins[i];
        double pos = (x - cacheMin) / step;
        double nearest = TMath::Nint(pos);

        if (fabs(pos - nearest) > 1e-8) {
            cerr << "ERROR: " << what
                 << ": final edge " << x
                 << " is not aligned with cache binning." << endl;
            cerr << "       cache: n=" << nCache
                 << " range=[" << cacheMin << ", " << cacheMax << "]"
                 << " step=" << step << endl;
            gSystem->Exit(1);
        }
    }
}
void CheckAllCacheCompatibility()
{


    for (int iCent = 0; iCent < nCentralityBins; iCent++) {

        CheckCacheCompatibility1D(Form("pT reco cent %d", iCent),
                                  nPtRecoCache, ptRecoMin, ptRecoMax,
                                  ptRecoBinsVec[iCent]);

        CheckCacheCompatibility1D(Form("pT true cent %d", iCent),
                                  nPtTrueCache, ptTrueMin, ptTrueMax,
                                  ptMcBinsVecCustom[iCent]);

        CheckCacheCompatibility1D(Form("z reco cent %d", iCent),
                                  nZRecoCache, zRecoMin, zRecoMax,
                                  zRecoBinsVec[iCent]);

        CheckCacheCompatibility1D(Form("z true cent %d", iCent),
                                  nZTrueCache, zTrueMin, zTrueMax,
                                  zMcBinsVecCustom[iCent]);

        for (int iAng = 0; iAng < nAngularities; iAng++) {


            CheckCacheCompatibility1D(Form("ang reco cent %d ang %d", iCent, iAng),
                                      nAngRecoCache, angRecoMin, angRecoMax,
                                      angRecoBinsVec[iCent][iAng]);

            CheckCacheCompatibility1D(Form("ang true cent %d ang %d", iCent, iAng),
                                      nAngTrueCache, angTrueMin, angTrueMax,
                                      angMcBinsVecCustom[iCent][iAng]);
        }
    }
}
void DrawAllBinning()
{
    gStyle->SetOptStat(0);

 
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);

    double xL = 0.03;
    double xR = 0.53;
    double yL = 0.97;
    double yR = 0.97;

    // ---------------- LEFT COLUMN ----------------
    latex.SetTextSize(0.030);
    latex.DrawLatex(xL, yL, "TRUE level");
    yL -= 0.045;

    latex.SetTextSize(0.024);
    DrawVectorBlock(latex, xL, yL, "p_{T,Jet}:", ptMcBinsVecCustom[0], 12, 1);
    DrawVectorBlock(latex, xL, yL, "z:", zMcBinsVecCustom[0], 12, 2);

    for (int ia = 0; ia < 6; ++ia) {
        DrawVectorBlock(latex, xL, yL, Form("%s:", kAngNames[ia]), angMcBinsVecCustom[0][ia], 12, 2);
    }

    yL -= 0.02;
    latex.SetTextSize(0.030);
    latex.DrawLatex(xL, yL, "RECO level");
    yL -= 0.045;

    latex.SetTextSize(0.022);
    for (int ic = 0; ic < 3; ++ic) {
        latex.DrawLatex(xL, yL, Form("%s", kCentNames[ic]));
        yL -= 0.035;

        DrawVectorBlock(latex, xL + 0.02, yL, "p_{T,Jet}:", ptRecoBinsVec[ic], 14, 1);
        DrawVectorBlock(latex, xL + 0.02, yL, "z:", zRecoBinsVec[ic], 12, 2);

        yL -= 0.01;
    }

    // ---------------- RIGHT COLUMN ----------------
    latex.SetTextSize(0.030);
    latex.DrawLatex(xR, yR, "RECO angularities");
    yR -= 0.045;

    latex.SetTextSize(0.022);
    for (int ic = 0; ic < 3; ++ic) {
        latex.DrawLatex(xR, yR, Form("%s", kCentNames[ic]));
        yR -= 0.032;

        for (int ia = 0; ia < 6; ++ia) {
            DrawVectorBlock(latex, xR + 0.02, yR, Form("%s:", kAngNames[ia]), angRecoBinsVec[ic][ia], 14, 2, 0.026, 0.006);
        }

        yR -= 0.015;
    }

}
void DrawLineAround(int iCent) {
    // Vykreslení čáry mezi biny
    TLine *lineA = new TLine(ptRecoBinsVec[iCent][0], ptMcBinsVecCustom[iCent][0],
                             ptRecoBinsVec[iCent][ptRecoBinsVec[iCent].size() - 1], ptMcBinsVecCustom[iCent][0]);
    lineA->SetLineColor(kRed);
    lineA->SetLineWidth(1);
    lineA->Draw();
    TLine *lineB = new TLine(ptRecoBinsVec[iCent][0], ptMcBinsVecCustom[iCent][ptMcBinsVecCustom[iCent].size() - 1],
                             ptRecoBinsVec[iCent][ptRecoBinsVec[iCent].size() - 1],
                             ptMcBinsVecCustom[iCent][ptMcBinsVecCustom[iCent].size() - 1]);
    lineB->SetLineColor(kRed);
    lineB->SetLineWidth(1);
    lineB->Draw();
    TLine *lineC = new TLine(ptRecoBinsVec[iCent][0], ptMcBinsVecCustom[iCent][0], ptRecoBinsVec[iCent][0],
                             ptMcBinsVecCustom[iCent][ptMcBinsVecCustom[iCent].size() - 1]);
    lineC->SetLineColor(kRed);
    lineC->SetLineWidth(1);
    lineC->Draw();
    TLine *lineD = new TLine(ptRecoBinsVec[iCent][ptRecoBinsVec[iCent].size() - 1], ptMcBinsVecCustom[iCent][0],
                             ptRecoBinsVec[iCent][ptRecoBinsVec[iCent].size() - 1],
                             ptMcBinsVecCustom[iCent][ptMcBinsVecCustom[iCent].size() - 1]);
    lineD->SetLineColor(kRed);
    lineD->SetLineWidth(1);
    lineD->Draw();
}

void DrawLineAround2(int iCent) {
    // Vykreslení čáry mezi biny
    TLine *lineA = new TLine(ptMcBinsVecCustom[iCent][0], ptMcBinsVecCustom[iCent][0],
                             ptMcBinsVecCustom[iCent][ptMcBinsVecCustom[iCent].size() - 1],
                             ptMcBinsVecCustom[iCent][0]);
    lineA->SetLineColor(kRed);
    lineA->SetLineWidth(1);
    lineA->Draw();
    TLine *lineB = new TLine(ptMcBinsVecCustom[iCent][0], ptMcBinsVecCustom[iCent][ptMcBinsVecCustom[iCent].size() - 1],
                             ptMcBinsVecCustom[iCent][ptMcBinsVecCustom[iCent].size() - 1],
                             ptMcBinsVecCustom[iCent][ptMcBinsVecCustom[iCent].size() - 1]);
    lineB->SetLineColor(kRed);
    lineB->SetLineWidth(1);
    lineB->Draw();
    TLine *lineC = new TLine(ptMcBinsVecCustom[iCent][0], ptMcBinsVecCustom[iCent][0], ptMcBinsVecCustom[iCent][0],
                             ptMcBinsVecCustom[iCent][ptMcBinsVecCustom[iCent].size() - 1]);
    lineC->SetLineColor(kRed);
    lineC->SetLineWidth(1);
    lineC->Draw();
    TLine *lineD = new TLine(ptMcBinsVecCustom[iCent][ptMcBinsVecCustom[iCent].size() - 1], ptMcBinsVecCustom[iCent][0],
                             ptMcBinsVecCustom[iCent][ptMcBinsVecCustom[iCent].size() - 1],
                             ptMcBinsVecCustom[iCent][ptMcBinsVecCustom[iCent].size() - 1]);
    lineD->SetLineColor(kRed);
    lineD->SetLineWidth(1);
    lineD->Draw();
}

TH2D *invertAxis(TH2D *hist) {

    vector <Double_t> xAxisVector = getAxisVector(hist->GetXaxis());
    vector <Double_t> yAxisVector = getAxisVector(hist->GetYaxis());

    Int_t nBinsX = xAxisVector.size() - 1;
    Int_t nBinsY = yAxisVector.size() - 1;

    TH2D *histInv = new TH2D(hist->GetName(), hist->GetTitle(), nBinsY, &yAxisVector[0], nBinsX, &xAxisVector[0]);
    for (Int_t iBinX = 1; iBinX <= nBinsX; iBinX++) {
        for (Int_t iBinY = 1; iBinY <= nBinsY; iBinY++) {
            histInv->SetBinContent(iBinY, iBinX, hist->GetBinContent(iBinX, iBinY));
            histInv->SetBinError(iBinY, iBinX, hist->GetBinError(iBinX, iBinY));
        }
    }
    return histInv;
}

void NormalizeByBinWidth(TH1D *hist, const Int_t color) {

    for (int i = 1; i <= hist->GetNbinsX(); i++) {
        if (hist->GetBinWidth(i) != 0) {
            hist->SetBinContent(i, hist->GetBinContent(i) / hist->GetBinWidth(i));
            hist->SetBinError(i, hist->GetBinError(i) / hist->GetBinWidth(i));
        }
    }

    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(20);
}

void RecolorHistogram(TH1D *hist, const Int_t color) {

    for (int i = 1; i <= hist->GetNbinsX(); i++) {
        if (hist->GetBinWidth(i) != 0) {
            hist->SetBinContent(i, hist->GetBinContent(i));
            hist->SetBinError(i, hist->GetBinError(i));
        }
    }
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(20);
}
Double_t findMin(TH1D hist, double Amin) {

        for (int j = 1; j <= hist.GetNbinsX(); j++) {
            if (hist.GetBinContent(j) < Amin && Amin > 0.1 && hist.GetBinContent(j) > 0) {
                Amin = hist.GetBinContent(j);
            }
    }
    return Amin;
}
Double_t findMin(TH1D *hist[], double Amin, int size) {

    for (int i = 0; i < size; i++) {
        for (int j = 1; j <= hist[i]->GetNbinsX(); j++) {
            if (hist[i]->GetBinContent(j) < Amin && Amin > 0.1 && hist[i]->GetBinContent(j) > 0) {
                Amin = hist[i]->GetBinContent(j);
            }
        }
    }
    return Amin;
}

Double_t findMin(TH1D hist[], double Amin, int size) {

    for (int i = 0; i < size; i++) {
        for (int j = 1; j <= hist[i].GetNbinsX(); j++) {
            if (hist[i].GetBinContent(j) < Amin && Amin > 0.1 && hist[i].GetBinContent(j) > 0) {
                Amin = hist[i].GetBinContent(j);
            }
        }
    }
    return Amin;
}

void DrawTextAbove(int siter = 0, int size = 20) {
    TLatex *texMin = new TLatex();
    texMin->SetNDC();
    texMin->SetTextFont(43);
    texMin->SetTextSize(size);
    texMin->SetTextAlign(22);
    TString BinningReco[3];
    TString BinningMc[3];
    //Vypíšu všechny okraje z ptRecoBinsVec[0] do BinningReco[0]
    BinningReco[0] = "p_{T}^{Reco}(0-10%) = {";
    BinningReco[1] = "p_{T}^{Reco}(10-40%) = {";
    BinningReco[2] = "p_{T}^{Reco}(40-80%) = {";
    BinningMc[0] = "p_{T}^{True}(0-10%) = {";
    BinningMc[1] = "p_{T}^{True}(10-40%) = {";
    BinningMc[2] = "p_{T}^{True}(40-80%) = {";

    for (int i = 0; i < ptRecoBinsVec[0].size() - 1; i++) BinningReco[0] += Form("%.0f, ", ptRecoBinsVec[0][i]);
    for (int i = 0; i < ptRecoBinsVec[1].size() - 1; i++) BinningReco[1] += Form("%.0f, ", ptRecoBinsVec[1][i]);
    for (int i = 0; i < ptRecoBinsVec[2].size() - 1; i++) BinningReco[2] += Form("%.0f, ", ptRecoBinsVec[2][i]);

    for (int i = 0; i < ptMcBinsVecCustom[0].size() - 1; i++) BinningMc[0] += Form("%.0f, ", ptMcBinsVecCustom[0][i]);
    for (int i = 0; i < ptMcBinsVecCustom[1].size() - 1; i++) BinningMc[1] += Form("%.0f, ", ptMcBinsVecCustom[1][i]);
    for (int i = 0; i < ptMcBinsVecCustom[2].size() - 1; i++) BinningMc[2] += Form("%.0f, ", ptMcBinsVecCustom[2][i]);

    BinningReco[0] += Form("%.0f}", ptRecoBinsVec[0][ptRecoBinsVec[0].size() - 1]);
    BinningReco[1] += Form("%.0f}", ptRecoBinsVec[1][ptRecoBinsVec[1].size() - 1]);
    BinningReco[2] += Form("%.0f}", ptRecoBinsVec[2][ptRecoBinsVec[2].size() - 1]);
    BinningMc[0] += Form("%.0f}", ptMcBinsVecCustom[0][ptMcBinsVecCustom[0].size() - 1]);
    BinningMc[1] += Form("%.0f}", ptMcBinsVecCustom[1][ptMcBinsVecCustom[1].size() - 1]);
    BinningMc[2] += Form("%.0f}", ptMcBinsVecCustom[2][ptMcBinsVecCustom[2].size() - 1]);

    const char *Word;
    if (ClosureTest) Word = "Closure test:";
    else Word = "Real data:";

//    TString textAbove = Form("#color[4]{Closure test:} p_{T}(D^{0}) > #color[2]{%.2f GeV/c}; Truth Bins: #color[2]{%d} (%.1f - %.1f), ecb: #color[2]{%s}; Reco Bins: #color[2]{%d}, ecb: #color[2]{%s}, w. prior: #color[2]{%s}", SetMinD0Pt, nMcBins, TruthJetPtMin, TruthJetPtMax ,useCustomPtMcBins?"false":"true", nRecoBins, useCustomPtRecoBins?"false":"true", WeightedPrior?"true":"false");
    if (siter != -1)
        texMin->DrawLatex(0.5, 0.98, Form("#color[4]{%s} p_{T}(D^{0}) > #color[2]{%.2f GeV/c}", Word, SetMinD0Pt));
    else
        texMin->DrawLatex(0.5, 0.98, Form("#color[4]{%s} p_{T}(D^{0}) > #color[2]{%.2f GeV/c}", Word, SetMinD0Pt));
    texMin->DrawLatex(0.28, 0.95, BinningReco[0]);
    texMin->DrawLatex(0.28, 0.92, BinningReco[1]);
    texMin->DrawLatex(0.28, 0.89, BinningReco[2]);
    texMin->DrawLatex(0.73, 0.95, BinningMc[0]);
    texMin->DrawLatex(0.73, 0.92, BinningMc[1]);
    texMin->DrawLatex(0.73, 0.89, BinningMc[2]);
}

Double_t findMax(TH1D *hist[], double Amax, int size) {

    for (int i = 0; i < size; i++) {
        for (int j = 1; j <= hist[i]->GetNbinsX(); j++) {
            if (hist[i]->GetBinContent(j) > Amax) {
                Amax = hist[i]->GetBinContent(j);
            }
        }
    }
    return Amax;
}

Double_t findMax(TH1D hist[], double Amax, int size) {

    for (int i = 0; i < size; i++) {
        for (int j = 1; j <= hist[i].GetNbinsX(); j++) {
            if (hist[i].GetBinContent(j) > Amax) {
                Amax = hist[i].GetBinContent(j);
            }
        }
    }
    return Amax;
}

//double centBins[nCentralityBins + 1] = {0, 10, 40, 80}; // in icreasing order
//double centBins[nCentralityBins + 1] = {0, 10, 40, 80}; // in icreasing order









//----2D----------
void plotComparison(TCanvas *can,
                    RooUnfoldResponse *hResponse,
                    TH2D *hUnfolded[],
                    TH2D *hRealData,
                    TH2D *hMcMeasured,
                    TH2D *hMcA,
                    const Int_t &iCent,
                    TString var,
                    const char *OutputFile,
                    TH2D *hMc,
                    TH2D *hMcMes) {

    int variable = 0;
    if (var == "z") variable = 1;
    if (var == "#lambda^{1}_{1}") variable = 2;
    if (var == "#lambda^{1}_{1.5}") variable = 3;
    if (var == "#lambda^{1}_{2}") variable = 4;
    if (var == "#lambda^{1}_{3}") variable = 5;
    if (var == "#lambda^{1}_{0.5}") variable = 6;
    if (var == "p_{T}^{D}") variable = 7;


    //vykreslím hResponse
    can->Clear();
    //gPad->SetTopMargin(0.15);

    can->Divide(2, 2);


    TLegend *leg1 = new TLegend(0.28, 0.70, 0.43, 0.92);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);

    TLegend *leg2 = new TLegend(0.75, 0.8, 0.95, 0.9);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.055);

    //x-axis
    TH1D *hRealDataProjX = (TH1D *) hRealData->ProjectionX(
            TString("hRealDataProjX_") + var + Form("_%i", iCent))->Clone(
            TString("Clone_hRealDataProjX_") + var + Form("_%i", iCent));
    //check if hMc is empty

    TH1D *hMcProjX = (TH1D *) hMc->ProjectionX(TString("hMcProjX_") + var + Form("_%i", iCent))->Clone(
            TString("ChMcProjX_") + var + Form("_%i", iCent));
    TH1D *hMcMeasuredProjX = (TH1D *) hMcMes->ProjectionX(
            TString("hMcMeasuredProjX_") + var + Form("_%i", iCent))->Clone(
            TString("Clone_hMcMeasuredProjX_") + var + Form("_%i", iCent));

    NormalizeByBinWidth(hRealDataProjX, 2001);
    NormalizeByBinWidth(hMcProjX, 2002);
    NormalizeByBinWidth(hMcMeasuredProjX, 2003);

    leg2->AddEntry(hRealDataProjX, "Real Data", "lp");
    if(ClosureTest) leg1->AddEntry(hMcProjX, "Mc (scaled)", "lp");
    else leg1->AddEntry((TH1D *) 0, "Real unfolded:", "");
    leg2->AddEntry(hMcMeasuredProjX, "Mc Reco (scaled)", "lp");

    can->cd(1);
/*
    TPad *myPad1 = new TPad("myPad_11_X", "myPad_11_Y", 0, 0.4, 1.0, 0.92);
    myPad1->SetBottomMargin(0.0);
    myPad1->SetBorderMode(0);
    TPad *myPad2 = new TPad("myPad_22_X", "myPad_22_Y", 0, 0.09, 1.0, 0.4);
    myPad2->SetTopMargin(0.0);
    myPad2->SetBorderMode(0);
    myPad1->Draw();
    myPad2->Draw();
    myPad1->cd();
*/
    TPad *Ratio_A = new TPad("Ratio_A", "Ratio_A", 0, 0.08, 1.0, 0.5);
    Ratio_A->SetTopMargin(0.0);

    TPad *Ratio_B = new TPad("Ratio_B", "Ratio_B", 0, 0.5, 1.0, 0.92);
    //bottom pad
    Ratio_B->SetBottomMargin(0.0);
    Ratio_B->SetBorderMode(0);
    TPad *Ratio_A2 = new TPad("Ratio_A2", "Ratio_A2", 0, 0.08, 1.0, 0.5);
    Ratio_A2->SetTopMargin(0.0);

    TPad *Ratio_B2 = new TPad("Ratio_B2", "Ratio_B2", 0, 0.5, 1.0, 0.92);
    Ratio_B2->SetBottomMargin(0.0);
    Ratio_B2->SetBorderMode(0);
    gPad->SetLogy();
    if(ClosureTest){
        hMcProjX->SetMarkerStyle(21);
        hMcProjX->Draw();
    } else hRealDataProjX->Draw();
    //od druhého binu do předposledního
    //hMcProjX->GetXaxis()->SetRange(2, hMcProjX->GetNbinsX()-1);





    //pause();
    // cout << "A" << endl;
    can->cd(2);
    gPad->SetBottomMargin(0.2);


    //y-axis
    //Udělat Clone projekce TH1D *hRealDataProjY = (TH1D *)hRealData->ProjectionY(TString("hRealDataProjY_")+var+Form("_%i", iCent));
    TH1D *hRealDataProjY = (TH1D *) hRealData->ProjectionY(
            TString("hRealDataProjY_") + var + Form("_%i", iCent))->Clone(
            TString("Clone_hRealDataProjY_") + var + Form("_%i", iCent));
    TH1D *hMcProjY = (TH1D *) hMc->ProjectionY(TString("hMcProjY_") + var + Form("_%i", iCent))->Clone(
            TString("Clone_hMcProjY_") + var + Form("_%i", iCent));
    TH1D *hMcMeasuredProjY = (TH1D *) hMcMes->ProjectionY(
            TString("hMcMeasuredProjY_") + var + Form("_%i", iCent))->Clone(
            TString("Clone_hMcMeasuredProjY_") + var + Form("_%i", iCent));

    Ratio_A2->Draw();
    Ratio_B2->Draw();
    NormalizeByBinWidth(hRealDataProjY, 2001);
    NormalizeByBinWidth(hMcProjY, 2002);
    NormalizeByBinWidth(hMcMeasuredProjY, 2003);
    ////can->SaveAs("./OutputPdf/"+TString(OutputFile)+"Rcp.pdf");


    can->cd(4);
    Ratio_A->Draw();
    Ratio_B->Draw();
    Ratio_B2->cd();
    can->cd(3);

    // cout << "B" << endl;
    gPad->SetLogy();
    hMcProjY->SetMarkerStyle(21);

    if (ClosureTest) {
        hMcProjY->Draw();
        //set x axis title
        hMcProjY->GetXaxis()->SetTitle(var);

    }
    else {
        hRealDataProjY->Draw();
        //set x axis title
        hRealDataProjY->GetXaxis()->SetTitle(var);

    }



    //hMcProjY->GetXaxis()->SetRange(1, hMcProjY->GetNbinsX());


    // first draw and compare y projections of hRealData and hMcMeasured
    can->cd(4);

    can->cd();
    leg1->Draw("same");
    //leg2->Draw("same");
    tex->DrawLatex(0.2, 0.85, "p_{T}");

    tex->DrawLatex(0.1, 0.35, centralityTitles[iCent] + "  " + var);

    TH1D *hUnfoldedProjX[nIter];
    TH1D *hUnfoldedProjY[nIter];

    PrintCheckNumbers[iCent][variable][0] = hRealData->Integral();

    TH1D *hBayRatiosX[nIter];
    TH1D *hBayRatiosX2[nIter];

    TH1D *hBayRatiosY[nIter];
    TH1D *hBayRatiosY2[nIter];



    for (int iter = 0; iter < nIter; iter++) {
        RooUnfoldBayes unfolding(hResponse, hRealData, PlotIterations[iter]);
        hUnfolded[iter] = (TH2D *) unfolding.Hreco();
        ////hUnfolded[iter] = (TH2D *) unfolding.Hunfold();

        PrintCheckNumbers[iCent][variable][iter + 1] = hUnfolded[iter]->Integral();


        //x-axis
        hUnfoldedProjX[iter] = (TH1D * )(
                hUnfolded[iter]->ProjectionX(TString("hUnfoldedProjX_") + var + Form("_%i_%i", iter, iCent)))->Clone(
                TString("Clone_hUnfoldedProjX_") + var + Form("_%i_%i", iter, iCent));


        NormalizeByBinWidth(hUnfoldedProjX[iter], 2005 + iter);
        hUnfoldedProjX[iter]->SetMarkerStyle(27);

        leg1->AddEntry(hUnfoldedProjX[iter], Form("Iter%i", PlotIterations[iter]), "lep");


        // first draw and compare x projections of hUnfolded and hMc
        can->cd(1);
       // myPad1->cd();

        hUnfoldedProjX[iter]->GetYaxis()->SetTitle("dN/dp_{T}");
        if (iter == 4) {
            hUnfoldedProjX[iter]->SetMarkerStyle(20);
            // hUnfoldedProjX[iter]->SetMarkerSize(1);
        }
        hUnfoldedProjX[iter]->Draw("same");

        //y-axis
        hUnfoldedProjY[iter] = (TH1D *) hUnfolded[iter]->ProjectionY(
                TString("hUnfoldedProjY_") + var + Form("_%i_%i", iter, iCent))->Clone(
                TString("Clone_hUnfoldedProjY_") + var + Form("_%i_%i", iter, iCent));

        NormalizeByBinWidth(hUnfoldedProjY[iter], 2005 + iter);
        hUnfoldedProjY[iter]->SetMarkerStyle(27);

        // first draw and compare y projections of hUnfolded and hMc
        can->cd(3);
        //gPad->SetLogx();

        //gPad->SetLogx();
        //zobrazi od druhého binu
        hUnfoldedProjY[iter]->GetYaxis()->SetTitle("dN/d" + var);
        //set x axis
        hUnfoldedProjY[iter]->GetXaxis()->SetTitle(var);

        if (iter == 4) {
            hUnfoldedProjY[iter]->SetMarkerStyle(20);
            // hUnfoldedProjY[iter]->SetMarkerSize(1);
        }
        hUnfoldedProjY[iter]->Draw("same");
        //hUnfoldedProjY[iter]->GetXaxis()->SetRange(1, hUnfoldedProjY[iter]->GetNbinsX());
        //hUnfoldedProjY[iter]->GetXaxis()->SetRangeUser(1, 2);


        hBayRatiosX[iter] = (TH1D *) hUnfoldedProjX[iter]->Clone(
                TString("hBayRatiosX_") + var + Form("_%i_%i", iter, iCent));
        hBayRatiosX[iter]->GetYaxis()->SetTitle("Unfolded i-th/Unfolded 1st");

        hBayRatiosX2[iter] = (TH1D *) hUnfoldedProjX[iter]->Clone(
                TString("hBayRatiosX2_") + var + Form("_%i_%i", iter, iCent));
        hBayRatiosX2[iter]->GetYaxis()->SetTitle("Unfolded i-th/Unfolded 1st");

        hBayRatiosY[iter] = (TH1D *) hUnfoldedProjY[iter]->Clone(
                TString("hBayRatiosY_") + var + Form("_%i_%i", iter, iCent));
        hBayRatiosY[iter]->GetYaxis()->SetTitle("Unfolded i-th/Unfolded 1st");

        hBayRatiosY2[iter] = (TH1D *) hUnfoldedProjY[iter]->Clone(
                TString("hBayRatiosY2_") + var + Form("_%i_%i", iter, iCent));

        hBayRatiosY2[iter]->GetYaxis()->SetTitle("Unfolded i-th/Unfolded 1st");


    }

   // myPad2->cd();
    Ratio_A2->cd();
    //bottom margin


    for (int iter = (UseRelativeP ? 1 : 0); iter < nIter; iter++) {

        hBayRatiosX[iter]->Draw(iter == (UseRelativeP ? 1 : 0) ? "" : "same");
        hBayRatiosX[iter]->Draw(iter == 0 ? "" : "same");
        hBayRatiosX[iter]->GetYaxis()->SetRangeUser(0.77, 1.23);
        //title size
        hBayRatiosX[iter]->GetYaxis()->SetTitleSize(0.09);
        //center
        hBayRatiosX[iter]->GetYaxis()->CenterTitle();
        //label size
        hBayRatiosX[iter]->GetYaxis()->SetLabelSize(0.08);
        //Vynechat každé druhé číslo na y-ové ose
        hBayRatiosX[iter]->GetYaxis()->SetNdivisions(505);
        //ylabel offset
        hBayRatiosX[iter]->GetYaxis()->SetTitleOffset(0.5);
        //x-axis
        hBayRatiosX[iter]->GetXaxis()->SetTitleSize(0.09);
        //label
        hBayRatiosX[iter]->GetXaxis()->SetLabelSize(0.08);
        //decreasing line width
        hBayRatiosX[(UseRelativeP ? 1 : 0)]->GetYaxis()->SetRangeUser(0.77, 1.23);
        if (iter == 0) continue;
        if (UseRelativeP) hBayRatiosX[iter]->Divide(hUnfoldedProjX[iter - 1]);
        else hBayRatiosX[iter]->Divide(hUnfoldedProjX[0]);

        //set range
        if (UseRelativeP) hBayRatiosX[iter]->GetYaxis()->SetTitle("i-th/(i-1)-th iteration");
        else hBayRatiosX[iter]->GetYaxis()->SetTitle("Unfolded/(1st iter)");

    }
    for (int iter = nIter - 1; iter >= (UseRelativeP ? 1 : 0); iter--) {
        hBayRatiosX[iter]->Draw("same");
    }
    DrawLineOne2(hBayRatiosX[0]->GetBinLowEdge(1), hBayRatiosX[0]->GetBinLowEdge(hBayRatiosX[0]->GetNbinsX() + 1));

    Ratio_B2->cd();
    for (int iter = 0; iter < nIter; iter++) {
        hBayRatiosX2[iter]->Draw(iter == 0 ? "" : "same");
        hBayRatiosX2[iter]->GetYaxis()->SetRangeUser(0.50, 1.5);
        hBayRatiosX2[iter]->Divide(hUnfoldedProjX[0]);
        hBayRatiosX2[iter]->GetYaxis()->SetRangeUser(0.77, 1.23);
        //title size
        hBayRatiosX2[iter]->GetYaxis()->SetTitleSize(0.09);
        //center
        hBayRatiosX2[iter]->GetYaxis()->CenterTitle();
        //label size
        hBayRatiosX2[iter]->GetYaxis()->SetLabelSize(0.08);
        //Vynechat každé druhé číslo na y-ové ose
        hBayRatiosX2[iter]->GetYaxis()->SetNdivisions(505);
        //ylabel offset
        hBayRatiosX2[iter]->GetYaxis()->SetTitleOffset(0.5);
        //x-axis
        hBayRatiosX2[iter]->GetXaxis()->SetTitleSize(0.09);
        //label
        hBayRatiosX2[iter]->GetXaxis()->SetLabelSize(0.08);
        //set range
        hBayRatiosX2[iter]->GetYaxis()->SetTitle("Unfolded/(1st iter)");
    }

    double xMax2 = findMax(hBayRatiosX2, 0, nIter);
    double xMin2 = findMin(hBayRatiosX2, 1, nIter);
    hBayRatiosX2[0]->GetYaxis()->SetRangeUser(0.77, 1.23);


    Ratio_A->cd();
    //myPad2Y->cd();
    //gPad->SetLogx();

    for (int iter = (UseRelativeP ? 1 : 0); iter < nIter; iter++) {
        hBayRatiosY[iter]->Draw(iter == 0 ? "" : "same");
        hBayRatiosY[iter]->GetYaxis()->SetRangeUser(0.77, 1.23);
        //title size
        hBayRatiosY[iter]->GetYaxis()->SetTitleSize(0.09);
        //center
        hBayRatiosY[iter]->GetYaxis()->CenterTitle();
        //label size
        hBayRatiosY[iter]->GetYaxis()->SetLabelSize(0.08);
        //Vynechat každé druhé číslo na y-ové ose
        hBayRatiosY[iter]->GetYaxis()->SetNdivisions(505);
        //ylabel offset
        hBayRatiosY[iter]->GetYaxis()->SetTitleOffset(0.5);
        //x-axis
        hBayRatiosY[iter]->GetXaxis()->SetTitleSize(0.09);
        //label
        hBayRatiosY[iter]->GetXaxis()->SetLabelSize(0.08);

        hBayRatiosY[(UseRelativeP ? 1 : 0)]->GetYaxis()->SetRangeUser(0.77, 1.23);
        if (iter == 0) continue;

        if (UseRelativeP) hBayRatiosY[iter]->Divide(hUnfoldedProjY[iter - 1]);
        else hBayRatiosY[iter]->Divide(hUnfoldedProjY[0]);
        //set range
        if (UseRelativeP) hBayRatiosY[iter]->GetYaxis()->SetTitle("i-th/(i-1)-th iteration");
        else hBayRatiosY[iter]->GetYaxis()->SetTitle("Unfolded/(1st iter)");
    }

    for (int iter = nIter - 1; iter >= (UseRelativeP ? 1 : 0); iter--) {
        hBayRatiosY[iter]->Draw("same");
    }
    DrawLineOne2(hBayRatiosY[0]->GetBinLowEdge(1), hBayRatiosY[0]->GetBinLowEdge(hBayRatiosY[0]->GetNbinsX() + 1));

    hMcProjX->Scale(hUnfoldedProjX[0]->Integral() / hMcProjX->Integral());
    hMcProjY->Scale(hUnfoldedProjY[0]->Integral() / hMcProjY->Integral());
    can->cd();
    //Nakreslím vertiáklní čáru
    TLine *line = new TLine(16, -8, 16, 12);
    line->SetLineColor(kBlack);
    line->SetLineColor(kBlack);
    line->SetLineStyle(1);
    line->SetLineWidth(2);
    //line->Draw();



    Double_t Xmax = findMax(hUnfoldedProjX, hMcProjX->GetMaximum(), nIter);
    Double_t Xmin = findMin(hUnfoldedProjX, hMcProjX->GetMinimum(), nIter);

    if (!ClosureTest){
        Xmax = max(Xmax, hRealDataProjX->GetMaximum());
        Xmin = min(Xmin, hRealDataProjX->GetMinimum());
    }

    Double_t Ymax = findMax(hUnfoldedProjY, hMcProjY->GetMaximum(), nIter);
    Double_t Ymin = findMin(hUnfoldedProjY, hMcProjY->GetMaximum(), nIter);

    if (!ClosureTest){
        Ymax = max(Ymax, hRealDataProjY->GetMaximum());
        Ymin = min(Ymin, hRealDataProjY->GetMinimum());
    }




    Ratio_B->cd();

for (int iter = 0; iter < nIter; iter++) {
    hBayRatiosY2[iter]->Draw(iter == 0 ? "" : "same");
    hBayRatiosY2[iter]->GetYaxis()->SetRangeUser(0.77, 1.23);
    hBayRatiosY2[iter]->Divide(hUnfoldedProjY[0]);
    }
    double yMax2 = findMax(hBayRatiosY2, 0, nIter);
    double yMin2 = findMin(hBayRatiosY2, 1, nIter);
    hBayRatiosY2[0]->GetYaxis()->SetRangeUser(0.77, 1.23);
    //title size
    hBayRatiosY2[0]->GetYaxis()->SetTitleSize(0.08);
    //center
    hBayRatiosY2[0]->GetYaxis()->CenterTitle();
    //label size
    hBayRatiosY2[0]->GetYaxis()->SetLabelSize(0.08);
    //Vynechat každé druhé číslo na y-ové ose
    hBayRatiosY2[0]->GetYaxis()->SetNdivisions(505);
    //offset
    hBayRatiosY2[0]->GetYaxis()->SetTitleOffset(0.5);
    if (ClosureTest) hMcProjX->GetYaxis()->SetRangeUser(Xmin>0?Xmin*0.01:0.001, Xmax*2);
    else hRealDataProjX->GetYaxis()->SetRangeUser(Xmin>0?Xmin*0.01:0.001, Xmax*2);
    if (ClosureTest) hMcProjY->GetYaxis()->SetRangeUser(Ymin>0?Ymin*0.5:0.5, Ymax*2);
    else hRealDataProjY->GetYaxis()->SetRangeUser(Ymin>0?Ymin*0.5:0.5, Ymax*2);

    hMcProjY->GetYaxis()->SetRangeUser(Ymin > 0 ? Ymin * 0.5 : 500000, Ymax * 2);
    //hMcProjY->GetYaxis()->SetRangeUser(500000, Ymax*2);

    can->cd();

    //DrawTextAbove();
    can->SaveAs("./OutputPdf/" + TString(OutputFile) + "Rcp.pdf");
    can->Clear();

    //Divid into two
    can->Divide(2, 1);
    can->cd(1);
    //left margin
    gPad->SetLeftMargin(0.18);
    gPad->SetLogy();

    //show underflow and overflow
    //set stats
    //hRealDataProjX->SetStats(kTRUE); // Povolit statistický box pro histogram
    //gStyle->SetOptStat("uo");
    hRealDataProjX->GetYaxis()->SetTitle("dN/dp_{T,Jet}");
    //offset
    hRealDataProjX->GetYaxis()->SetTitleOffset(1.5);
    hRealDataProjX->GetXaxis()->SetTitle("p_{T,Jet}");
    hRealDataProjX->Draw();
    hMcMeasuredProjX->Scale(hRealDataProjX->Integral("width") / hMcMeasuredProjX->Integral("width"));

    hMcMeasuredProjX->Draw("same");
    //hMcMeasuredProjX->SetStats(kTRUE); // Povolit statistický box pro histogram
    //hMcMeasuredProjX->GetXaxis()->SetRange(1, hMcMeasuredProjX->GetNbinsX());

    //Vypíšu hodnoty binů obou histogramů:
    for (int i = 1; i <= hRealDataProjX->GetNbinsX(); i++) {
        cout << "2DBin " << i << ": Real Data: " << hRealDataProjX->GetBinContent(i) << ", Mc Measured: "
             << hMcMeasuredProjX->GetBinContent(i) << endl;
    }


    TLegend leg_1(0.2, 0.2, 0.45, 0.28);
    leg_1.SetBorderSize(0);
    leg_1.SetFillStyle(0);
    //text size
    leg_1.SetTextSize(0.04);
    leg_1.AddEntry(hRealDataProjX, "Real Data", "lp");
    leg_1.AddEntry(hMcMeasuredProjX, "RM reco (scaled)", "lp");

    //draw
    leg_1.Draw("same");



    can->cd(2);
    //left margin
    gPad->SetLeftMargin(0.18);
    gPad->SetLogy();
    hRealDataProjY->GetYaxis()->SetTitle("dN/d" + var);
    //offset
    hRealDataProjY->GetYaxis()->SetTitleOffset(1.5);
    hRealDataProjY->GetXaxis()->SetTitle(var);
    hRealDataProjY->Draw();
    hMcMeasuredProjY->Scale(hRealDataProjY->Integral("width") / hMcMeasuredProjY->Integral("width"));
    hMcMeasuredProjY->Draw("same");
    //hMcMeasuredProjY->GetXaxis()->SetRange(1, hMcMeasuredProjY->GetNbinsX());

    can->SaveAs("./OutputPdf/" + TString(OutputFile) + "Rcp.pdf");
    can->Clear();

    //----------------------------------------------------------------------------
    if (true) {
        can->cd();
        can->Clear();
        can->Divide(2, 2);
        can->cd(1);
        gPad->SetTopMargin(0.25);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.15);

        //Tlatex text
        TLatex *texA = new TLatex();
        texA->SetNDC();
        texA->SetTextFont(42);
        texA->SetTextSize(0.055);
        texA->SetTextAlign(22);


        hMc->GetYaxis()->SetTitle(var);
        Stejn(hMc, TString("hUnfoldedPraaaojY_") + var + Form("_%i_%i", 20, iCent));
        //hMc->Draw("colz");
        //set title of histogram
        //hMc->SetTitle("MC reco");

        //setlogz
        gPad->SetLogz();

        texA->DrawLatex(0.5, 0.5, "MC truth");
        /////////////////////////////////////////////
        can->cd(3);
        gPad->SetTopMargin(0.25);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.15);
        hMcMes->GetYaxis()->SetTitle(var);

        Stejn(hMcMes, TString("hUnfoldedPraaaasdojY_") + var + Form("_%i_%i", 20, iCent));

        //hMcMes->SetTitle("MC truth");

        //hMcMes->Draw("colz");

        gPad->SetLogz();
        texA->DrawLatex(0.5, 0.5, "MC reco");

        /////////////////////////////////////////////
        can->cd(4);
        gPad->SetTopMargin(0.25);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.15);
        hRealData->GetYaxis()->SetTitle(var);
        //draw stats
        hRealData->SetStats(1);

        Stejn(hRealData, TString("hUnfolasdojY_") + var + Form("_%i_%i", 20, iCent));

        gPad->SetLogz();
        texA->DrawLatex(0.5, 0.5, "Real reco");
        ////////////////////////////////////////////////
        can->cd(2);
        gPad->SetTopMargin(0.25);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.15);
        hUnfolded[0]->GetYaxis()->SetTitle(var);
        Stejn(hUnfolded[0], TString("hUnfoldedPsraaaojY_") + var + Form("_%i_%i", 20, iCent));
        gPad->SetLogz();
        texA->DrawLatex(0.5, 0.5, "Unfolded");


        can->cd();
        texA->DrawLatex(0.5, 0.5, "(p_{T}, " + var + ") " + centralityTitles[iCent]);

        can->SaveAs("./OutputPdf/" + TString(OutputFile) + "Rcp.pdf");
        can->Clear();

        int histN;
        if (var == "z") {
            histN = 1;
        } else if (var == "#lambda^{1}_{1}") {
            histN = 2;
        } else if (var == "#lambda^{1}_{1.5}") {
            histN = 3;
        } else if (var == "#lambda^{1}_{2}") {
            histN = 4;
        } else if (var == "#lambda^{1}_{3}") {
            histN = 5;
        } else if (var == "#lambda^{1}_{0.5}") {
            histN = 6;
        } else if (var == "p_{T}^{D}") {
            histN = 7;
        }

        can->Divide(2, 2);
        can->cd(1);
        //leftmargin
        gPad->SetLeftMargin(0.15);
        //rightmargin
        gPad->SetRightMargin(0.15);
        Stejn(hRespZ[variable][iCent], TString("hresponzea") + var + Form("_%i_%i", variable, iCent));
        gPad->SetLogz();


        can->cd(3);
        //leftmargin
        gPad->SetLeftMargin(0.15);
        //rightmargin
        gPad->SetRightMargin(0.15);
        Stejn(hRespZ[0][iCent], TString("hresponzeapT") + var + Form("_%i_%i", variable, iCent));
        gPad->SetLogz();
        can->cd(2);
        //leftmargin
        gPad->SetLeftMargin(0.15);
        //rightmargin
        gPad->SetRightMargin(0.15);
        hRespZHighRes[variable][iCent]->Draw("colz");
        gPad->SetLogz();
        can->cd(4);
        //leftmargin
        gPad->SetLeftMargin(0.15);
        //rightmargin
        gPad->SetRightMargin(0.15);
        hRespZHighRes[0][iCent]->Draw("colz");
        gPad->SetLogz();
        can->cd();
        texA->DrawLatex(0.5, 0.5, "2D Response matrix (p_{T}, " + var + ") " + centralityTitles[iCent]);
        can->SaveAs("./OutputPdf/" + TString(OutputFile) + "Rcp.pdf");
        can->Clear();
        can->cd();
        //hResponse->Hresponse()->Draw();
        //uložím to do nového histrogramu
        TH2D *hResponseClone = (TH2D *) hResponse->Hresponse()->Clone(Form("hResponseClone_%i_%s", iCent, var.Data()));
        hResponseClone->Draw("colz");
        //set x axis label
        hResponseClone->GetXaxis()->SetTitle("N_{" + var + "}^{reco} (N_{p_{T}}^{reco})");
        hResponseClone->GetYaxis()->SetTitle("N_{" + var + "}^{true} (N_{p_{T}}^{true})");

        gPad->SetLogz();

        int nxbinsZ = hMcMes->GetYaxis()->GetNbins();
        int nybinsZ = hUnfolded[0]->GetYaxis()->GetNbins();
        int nxbinsPt = hMcMes->GetXaxis()->GetNbins();
        int nybinsPt = hUnfolded[0]->GetXaxis()->GetNbins();

        //draw vertical lines
        for (int i = 0; i < nybinsZ + 1; i++) {
            TLine *line = new TLine(0, nybinsPt * i, nxbinsZ * nxbinsPt, nybinsPt * i);

            line->SetLineColor(kBlack);
            line->SetLineStyle(1);
            line->SetLineWidth(1);
            line->Draw(i == 0 ? "" : "same");


        }

        for (int i = 0; i < nxbinsZ + 1; i++) {

            TLine *line2 = new TLine(nxbinsPt * i, 0, nxbinsPt * i, nybinsZ * nybinsPt);

            line2->SetLineColor(kBlack);
            line2->SetLineStyle(1);
            line2->SetLineWidth(1);
            line2->Draw(i == 0 ? "" : "same");
        }
        texA->DrawLatex(0.5, 0.93, "4D Response matrix (p_{T}, " + var + ") " + centralityTitles[iCent]);


        can->SaveAs("./OutputPdf/" + TString(OutputFile) + "Rcp.pdf");
        can->Clear();
    }
//
    if (ClosureTest) {

        TH1D *hClosureMCpT2 = (TH1D *) hMc->ProjectionX(TString("__hClosureMCpT2") + Form("_%i", iCent))->Clone(
                TString("__hClosureMCpTA2") + Form("_%i", iCent));
        TH1D *hClosureMCpT2Y = (TH1D *) hMc->ProjectionY(TString("__hClosureMCpT2Y") + Form("_%i", iCent))->Clone(
                TString("__hClosureMCpTA2Y") + Form("_%i", iCent));

        TH1D *hUnfoldedPt2[nIter];
        TH1D *hUnfoldedPt2Y[nIter];

        for (int iter = 0; iter < nIter; iter++) {
            hUnfoldedPt2[iter] = (TH1D *) hUnfoldedProjX[iter]->Clone(
                    TString("_hUnfoldedPt2X") + Form("_%i_%i", iter, iCent));
            hUnfoldedPt2Y[iter] = (TH1D *) hUnfoldedProjY[iter]->Clone(
                    TString("_hUnfoldedPt2Y") + Form("_%i_%i", iter, iCent));

        }

        can->cd();
        can->Clear();
        can->Divide(2, 2);
        //------------------------------------
        can->cd(1);
        gPad->SetTopMargin(0.25);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.15);
        hClosureMCpT2->GetYaxis()->SetTitle("dN/dp_{T, Jet} [(GeV/c)^{-1}]");
        hClosureMCpT2->GetXaxis()->SetTitle("p_{T, Jet} [GeV/c]");
        //set logy
        gPad->SetLogy();
        NormalizeByBinWidth(hClosureMCpT2, 4);

        hClosureMCpT2->Draw();
        //legenda
        TLegend *leg = new TLegend(0.55, 0.35, 0.95, 0.62);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(hClosureMCpT2, "MC (scaled)", "lep");


        for (int iter = 0; iter < nIter; iter++) {
            //NormalizeByBinWidth(hUnfoldedPt2[iter], 2000 + iter);
            //scale by maximum
            hUnfoldedPt2[iter]->Scale(hClosureMCpT2->Integral() / hUnfoldedPt2[iter]->Integral());

            hUnfoldedPt2[iter]->Draw("same");
            leg->AddEntry(hUnfoldedPt2[iter], Form("Iter%i", PlotIterations[iter]), "lep");
        }

        leg->Draw("same");
        //---------------------------------------
        can->cd(3);
        gPad->SetTopMargin(0.15);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.15);
        hClosureMCpT2Y->GetYaxis()->SetTitle("dN/dp_{T, Jet} [(GeV/c)^{-1}]");
        hClosureMCpT2Y->GetXaxis()->SetTitle("p_{T, Jet} [GeV/c]");
        //set logy
        gPad->SetLogy();
        //gPad->SetLogx();
        NormalizeByBinWidth(hClosureMCpT2Y, 4);

        hClosureMCpT2Y->Draw();
        //legenda
        TLegend *legY = new TLegend(0.55, 0.55, 0.95, 0.82);
        legY->SetBorderSize(0);
        legY->SetFillStyle(0);
        legY->AddEntry(hClosureMCpT2Y, "MC (scaled)", "lep");

        hClosureMCpT2Y->Scale(1. * hUnfoldedPt2Y[0]->Integral("width") / hClosureMCpT2Y->Integral("width"));

        for (int iter = 0; iter < nIter; iter++) {
            //NormalizeByBinWidth(hUnfoldedPt2[iter], 2000 + iter);
            //scale by maximum
            //hUnfoldedPt2Y[iter]->Scale(hClosureMCpT2Y->GetMaximum() / hUnfoldedPt2Y[iter]->GetMaximum());

            hUnfoldedPt2Y[iter]->Draw("same");
            legY->AddEntry(hUnfoldedPt2Y[iter], Form("Iter%i", PlotIterations[iter]), "lep");
        }

        legY->Draw("same");
        //----------------------------------------
        can->cd(2);
        gPad->SetTopMargin(0.25);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.15);

        //legenda2
        TLegend *leg2 = new TLegend(0.55, 0.47, 0.95, 0.74);
        leg2->SetBorderSize(0);
        leg2->SetFillStyle(0);

        TH1D *hRatiosMc[nIter];
        for (int iter = 0; iter < nIter; iter++) {
            hRatiosMc[iter] = (TH1D *) hUnfoldedPt2[iter]->Clone(
                    TString("hRyatiosMc") + Form("_%i_%i_%i", iter, iCent, 0));
            hRatiosMc[iter]->Divide(hClosureMCpT2);
            hRatiosMc[iter]->GetYaxis()->SetTitle("Unfolded/MC");
            hRatiosMc[iter]->GetYaxis()->SetRangeUser(0.5, 1.5);
            hRatiosMc[iter]->Draw(iter == 0 ? "" : "same");
            leg2->AddEntry(hRatiosMc[iter], Form("Iter%i", PlotIterations[iter]), "lep");
        }
        leg2->Draw("same");
        DrawLineOne();
//-----------------------------------------------------------
        can->cd(4);
        gPad->SetTopMargin(0.15);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.15);

        //legenda2
        TLegend *leg2Y = new TLegend(0.55, 0.55, 0.95, 0.82);
        leg2Y->SetBorderSize(0);
        leg2Y->SetFillStyle(0);

        TH1D *hRatiosMcY[nIter];
        for (int iter = 0; iter < nIter; iter++) {
            hRatiosMcY[iter] = (TH1D *) hUnfoldedPt2Y[iter]->Clone(
                    TString("hRyatiosMcY") + Form("_%i_%i_%i", iter, iCent, 0));
            hRatiosMcY[iter]->Divide(hClosureMCpT2Y);
            hRatiosMcY[iter]->GetYaxis()->SetTitle("Unfolded/MC");
            hRatiosMcY[iter]->GetYaxis()->SetRangeUser(0.5, 1.5);
            hRatiosMcY[iter]->Draw(iter == 0 ? "" : "same");
            leg2Y->AddEntry(hRatiosMcY[iter], Form("Iter%i", PlotIterations[iter]), "lep");
        }
        leg2Y->Draw("same");

        //left bin edge of hRatiosMcY
        DrawLineOne2(hRatiosMcY[0]->GetBinLowEdge(2), hRatiosMcY[0]->GetBinLowEdge(hRatiosMcY[0]->GetNbinsX() + 1));
        //------------------------------------------------------
        can->cd();
        DrawTextAbove(0);
        tex->DrawLatex(0.5, 0.05, centralityTitles[iCent] + "  " + var);
//
        //save
        can->SaveAs("./OutputPdf/" + TString(OutputFile) + "Rcp.pdf");

        //---------------------------------
        can->Clear();
        can->cd();
        can->Divide(2, 1);

        TLegend *leg3 = new TLegend(0.55, 0.55, 0.95, 0.82);
        leg3->SetBorderSize(0);
        leg3->SetFillStyle(0);

        can->cd(1);
        //upper margin
        gPad->SetTopMargin(0.15);
        gPad->SetLeftMargin(0.15);

        //vypnu logscale
        gPad->SetLogy(0);
        TH1D *hRelativeRatios[nIter];
        for (int iter = 1; iter < nIter; iter++) {
            hRelativeRatios[iter] = (TH1D *) hUnfoldedPt2[iter]->Clone(
                    TString("hyRelativeRatios") + Form("_%i_%i_%i", iter, iCent, 0));
            hRelativeRatios[iter]->Divide(hUnfoldedPt2[iter - 1]);
            hRelativeRatios[iter]->GetYaxis()->SetTitle("Unfolded n-th/Unfolded (n-1)-th");
            hRelativeRatios[iter]->GetYaxis()->SetRangeUser(0.5, 1.5);
            hRelativeRatios[iter]->GetYaxis()->SetTitleOffset(1.1);
            hRelativeRatios[iter]->Draw(iter == 1 ? "" : "same");
            leg3->AddEntry(hRelativeRatios[iter], Form("Iter%i", PlotIterations[iter]), "lep");
        }
        leg3->Draw("same");
        DrawLineOne();

        can->cd(2);
        //upper margin
        gPad->SetTopMargin(0.15);
        gPad->SetLeftMargin(0.15);

        TH1D *hRelativeUncertRatiosVs3[nIter];
        TH1D *hRelativeUncertRatios;
        TH1D *hRelativeUncertRatiosA[nIter];

        TLegend *leg4 = new TLegend(0.55, 0.55, 0.95, 0.82);
        leg4->SetBorderSize(0);
        leg4->SetFillStyle(0);

        for (int iter = 0; iter < nIter; iter++) {
            hRelativeUncertRatiosA[iter] = (TH1D *) hUnfoldedPt2[iter]->Clone(
                    TString("ahRelativeUncertRatiosVs3") + Form("_%i_%i_%i", iter, iCent, 0));
            //Vynuluji všechny hodnoty a chyby;


            for (int i = 0; i < hRelativeUncertRatiosA[iter]->GetNbinsX(); i++) {
                hRelativeUncertRatiosA[iter]->SetBinContent(i + 1, 0);
                hRelativeUncertRatiosA[iter]->SetBinError(i + 1, 0.00001);
                double value = hUnfoldedPt2[iter]->GetBinError(i + 1) / hUnfoldedPt2[2]->GetBinError(i + 1);
                hRelativeUncertRatiosA[iter]->SetBinContent(i + 1, value);
            }


            hRelativeUncertRatiosA[iter]->GetYaxis()->SetTitle("Unfolded n-th/Unfolded 3-rd (uncertainty)");
            //left offset
            hRelativeUncertRatiosA[iter]->GetYaxis()->SetTitleOffset(1.1);
            hRelativeUncertRatiosA[iter]->GetYaxis()->SetRangeUser(0, 3);
            hRelativeUncertRatiosA[iter]->Draw(iter == 0 ? "" : "same");
            leg4->AddEntry(hRelativeUncertRatiosA[iter], Form("Iter%i", PlotIterations[iter]), "lep");
        }

        leg4->Draw("same");


        can->cd();
        DrawTextAbove(0);
        tex->DrawLatex(0.5, 0.05, centralityTitles[iCent] + "  " + var);

        can->SaveAs("./OutputPdf/" + TString(OutputFile) + "Rcp.pdf");

        //---------------------------------
        DeleteArray(hRatiosMc, nIter);
        DeleteArray(hRatiosMcY, nIter);

        DeleteArray(hUnfoldedPt2, nIter);
        DeleteArray(hUnfoldedPt2Y, nIter);

        //delete hClosureMCpT2;
        delete hClosureMCpT2;
        delete hClosureMCpT2Y;


    }





    //---------------------------------------------------------------------------/
    delete leg1;
    delete leg2;
    delete tex;
    delete hRealDataProjX;
    delete hMcProjX;
    delete hMcMeasuredProjX;
    delete hRealDataProjY;
    delete hMcProjY;
    delete hMcMeasuredProjY;
    DeleteArray(hUnfoldedProjX, nIter);
    DeleteArray(hUnfoldedProjY, nIter);
}

void plotComparisonSVD(TCanvas *can, RooUnfoldResponse *hResponse, TH1D *hUnfoldedSVD[], TH1D *hRealData, TH1D *hMc,
                       TH1D *hMcMeasured, const Int_t &iCent, TString var, const char *OutputFile) {


    can->Clear();
    can->Divide(2, 1);
    TLegend *leg1 = new TLegend(0.29, 0.65, 0.45, 0.83);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);

    TLegend *leg2 = new TLegend(0.75, 0.75, 0.95, 0.82);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);


    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.055);

    TH1D *hRealDataProjX = (TH1D *) hRealData->Clone(TString("SVDhRealDataProjXPt") + Form("_%i", iCent));
    TH1D *hMcProjX = (TH1D *) hMc->Clone(TString("SVDhMcProjXPt") + Form("_%i__%i", 1, iCent));
    TH1D *hMcMeasuredProjX = (TH1D *) hMcMeasured->Clone(TString("SVDhMcMeasuredProjXPt") + Form("_%i", iCent));

    NormalizeByBinWidth(hRealDataProjX, 2000 + 1);
    NormalizeByBinWidth(hMcProjX, 2000 + 2);
    NormalizeByBinWidth(hMcMeasuredProjX, 2000 + 3);

    leg2->AddEntry(hRealDataProjX, "Real Data", "lep");
    leg1->AddEntry(hMcProjX, "Mc", "lep");
    leg1->AddEntry((TH1D *) 0, "Real unfolded:", "");
    leg2->AddEntry(hMcMeasuredProjX, "Mc Reco (scaled)", "lp");

    can->cd(1);
    gPad->SetTopMargin(0.15);

    TPad *myPad1 = new TPad("myPad_1", "myPad_1", 0, 0.4, 1.0, 0.92);
    myPad1->SetBottomMargin(0.0);
    myPad1->SetBorderMode(0);
    TPad *myPad2 = new TPad("myPad_2", "myPad_2", 0, 0.09, 1.0, 0.4);
    myPad2->SetTopMargin(0.0);
    myPad2->SetBorderMode(0);
    myPad1->Draw();
    myPad2->Draw();
    myPad1->cd();
    //Nastavím marker na hvězdičku
    hMcProjX->SetMarkerStyle(21);

    //hMcProjX->SetMarkerSize(2);
    gPad->SetLogy();
    hMcProjX->Draw();
    //hMcProjX->GetYaxis()->SetRangeUser(100, 10000);

    //hRealMeasured->Draw("same");
    // first draw and compare x projections of hRealData and hMcMeasured
    can->cd(2);
    gPad->SetTopMargin(0.15);
    gPad->SetLogy();

    hRealDataProjX->GetYaxis()->SetTitle("dN/dp_{T, Jet} [(GeV/c)^{-1}]");
    hRealDataProjX->GetXaxis()->SetTitle("p_{T, Jet} [GeV/c]");
    hRealDataProjX->Draw();
    hMcMeasuredProjX->Scale(hRealDataProjX->Integral() / hMcMeasuredProjX->Integral());
    hMcMeasuredProjX->Draw("same");

    can->cd();
    leg1->Draw("same");
    leg2->Draw("same");
    //tex->DrawLatex(0.2, 0.68, "p_{T}");

    tex->DrawLatex(0.15, 0.05, "SVD " + centralityTitles[iCent] + "  " + var);


    TH1D *hUnfoldedProjXSVD[nKterm];
    TH1D *hSvdRatios[nKterm];
    for (int iKterm = 0; iKterm < nKterm; iKterm++) {
        //RooUnfoldSvd unfoldingSVD(hResponse, hRealData, kterm[iKterm]);
        RooUnfoldSvd *unfoldingSVD = new RooUnfoldSvd(hResponse, hRealData, kterm[iKterm]);

        hUnfoldedSVD[iKterm] = (TH1D *) unfoldingSVD->Hreco();
        ////hUnfoldedSVD[iKterm] = (TH1D *) unfoldingSVD->Hunfold();

        hUnfoldedSVD[iKterm]->SetName(TString("hSVDUnfoldedPt") + Form("_%i_%i", iKterm, iCent));
        hUnfoldedProjXSVD[iKterm] = (TH1D *) hUnfoldedSVD[iKterm]->Clone(
                TString("hUnfoldedSVDProjXPt") + Form("_%i_%i", iKterm, iCent));
        NormalizeByBinWidth(hUnfoldedProjXSVD[iKterm], 2005 + iKterm);
        hUnfoldedProjXSVD[iKterm]->SetMarkerStyle(27);
        leg1->AddEntry(hUnfoldedProjXSVD[iKterm], Form("kterm%i", kterm[iKterm]), "lep");
        myPad1->cd();
        hUnfoldedProjXSVD[iKterm]->GetYaxis()->SetTitle("dN/dp_{T, Jet} [(GeV/c)^{-1}]");
        if (iKterm == 2) {
            hUnfoldedProjXSVD[iKterm]->SetMarkerStyle(20);
            // hUnfoldedProjXSVD[iKterm]->SetMarkerSize(1);
        }
        hUnfoldedProjXSVD[iKterm]->Draw("same");
        myPad2->cd();
        hSvdRatios[iKterm] = (TH1D *) hUnfoldedProjXSVD[iKterm]->Clone(
                TString("hSvdRatios") + Form("_%i_%i", iKterm, iCent));
        hSvdRatios[iKterm]->Divide(hMcProjX);
        hSvdRatios[iKterm]->Draw(iKterm == 0 ? "" : "same");
        hSvdRatios[iKterm]->GetYaxis()->SetTitle("Unfolded/Mc");
    }
    double MinA = findMin(hSvdRatios, 2000, nKterm);
    double MaxA = findMax(hSvdRatios, 0, nKterm);
    hSvdRatios[0]->GetYaxis()->SetRangeUser(MinA - 0.1 * MinA, MaxA + 0.1 * MaxA);

    //Najdu který prvek v kterm[nKterm] je roven 2

    //hMcProjX->Scale(hUnfoldedProjXSVD[ChosenKterm]->GetMaximum() / hMcProjX->GetMaximum());
    //hRealMeasured->Scale(hMcProjX->GetMaximum() / hRealMeasured->GetMaximum());

    Double_t Xmax = findMax(hUnfoldedProjXSVD, hMcProjX->GetMaximum(), nKterm);
    Double_t Xmin = findMin(hUnfoldedProjXSVD, hMcProjX->GetMinimum(), nKterm);
    hMcProjX->GetYaxis()->SetRangeUser(Xmin > 0 ? Xmin * 0.5 : 10, Xmax * 2);
    can->cd();
    //Nakreslím vertiáklní čáru
    TLine *line = new TLine(16, -8, 16, 12);
    line->SetLineColor(kBlack);
    line->SetLineColor(kBlack);
    line->SetLineStyle(1);
    line->SetLineWidth(2);
    //line->Draw();
    can->cd();
    DrawTextAbove();
    can->SaveAs("./OutputPdf/" + TString(OutputFile) + "Rcp.pdf");

}
void OpenStabilityFile(const char* scanDir)
{
    if (!scanDir || scanDir[0] == '\0') scanDir = ".";

    // ofstream složku nevytvoří, proto ji vytvoříme ručně
    if (gSystem->AccessPathName(scanDir)) {
        int status = gSystem->mkdir(scanDir, kTRUE);  // kTRUE = recursive mkdir
        if (status != 0 && gSystem->AccessPathName(scanDir)) {
            cout << "[stability] ERROR: cannot create directory: " << scanDir << endl;
            gSystem->Exit(1);
        }
    }

    if (fout.is_open()) fout.close();

    TString stabilityFile = TString(scanDir) + "/stability.tsv";
    fout.open(stabilityFile.Data(), std::ios::out | std::ios::app);

    if (!fout.is_open()) {
        cout << "[stability] ERROR: cannot open file: " << stabilityFile << endl;
        gSystem->Exit(1);
    }

    cout << "[stability] writing to: " << stabilityFile << endl;
}

TH2D *getPearsonCoeffs1D(const TMatrixD &covMatrix) {
    Int_t nrows = covMatrix.GetNrows();
    Int_t ncols = covMatrix.GetNcols();
    //cout << "nrows: " << nrows << " ncols: " << ncols << endl;
 //   TString histName = Form("PearsonCoeffs_1D_%s", covMatrix.GetName()); // Oprava názvu histogramu
    TString histName = Form("PearsonCoeffs_1D_%s_%d", covMatrix.GetName(), (int)time(0));
    TH2D *PearsonCoeffs = new TH2D(histName, "Pearson Coefficients 1D;Coefficient Index;Pearson Coefficient", ncols, 0,
                                   ncols, nrows, 0, nrows);

    for (Int_t row = 0; row < nrows; row++) {
        for (Int_t col = 0; col < ncols; col++) {
            Double_t pearson = 0.;
            if (covMatrix(row, row) != 0. && covMatrix(col, col) != 0.)
                pearson = covMatrix(row, col) / TMath::Sqrt(covMatrix(row, row) * covMatrix(col, col));

            // Calculate the linear index for the 1D histogram
            //Int_t binIndex = (row * ncols) + col + 1; // Oprava výpočtu lineárního indexu

            PearsonCoeffs->SetBinContent(col + 1, row + 1, pearson);
        }
    }

    //PearsonCoeffs->GetYaxis()->SetRangeUser(-1, 1);
    return PearsonCoeffs;
}

TH2D getPearsonCoeffs1D_(const TMatrixD &covMatrix) {
    Int_t nrows = covMatrix.GetNrows();
    Int_t ncols = covMatrix.GetNcols();
    //cout << "nrows: " << nrows << " ncols: " << ncols << endl;
    //   TString histName = Form("PearsonCoeffs_1D_%s", covMatrix.GetName()); // Oprava názvu histogramu
    TString histName = Form("PearsonCoeffs_1D_%s_%d", covMatrix.GetName(), (int)time(0));
    TH2D PearsonCoeffs(histName, "Pearson Coefficients 1D;Coefficient Index;Pearson Coefficient", ncols, 0,
                                  ncols, nrows, 0, nrows);

    for (Int_t row = 0; row < nrows; row++) {
        for (Int_t col = 0; col < ncols; col++) {
            Double_t pearson = 0.;
            if (covMatrix(row, row) != 0. && covMatrix(col, col) != 0.)
                pearson = covMatrix(row, col) / TMath::Sqrt(covMatrix(row, row) * covMatrix(col, col));

            // Calculate the linear index for the 1D histogram
            //Int_t binIndex = (row * ncols) + col + 1; // Oprava výpočtu lineárního indexu

            PearsonCoeffs.SetBinContent(col + 1, row + 1, pearson);
        }
    }

    //PearsonCoeffs->GetYaxis()->SetRangeUser(-1, 1);
    return PearsonCoeffs;
}
/*
void superIterEvol(TCanvas* can, TH1D* RatioWeight[], TString OutputFile){
    can->Clear();
    //Vykreslím všechny RatioWeight
    for (int i = 0; i < nIter; i++){
        RatioWeight[i]->SetMarkerStyle(20);
        RatioWeight[i]->SetMarkerSize(1);
        RatioWeight[i]->SetMarkerColor(i+2000);
        RatioWeight[i]->SetLineColor(i+2000);
        RatioWeight[i]->Draw(i==0?"":"same");
    }
    can->SaveAs("./OutputPdf/"+TString(OutputFile)+"RatioWeight.pdf");

}
*/
//---1D---------------
//plotComparison( can, responsePt[iCent], hUnfoldedPt[iCent], hMeasuredPt[iCent], iCent, "p_{T}",OutputFile);

void plotComparison1D(TCanvas *can, const Int_t &iCent){

    //TO DELETE LATER
    TString var = "p_{T,Jet}";

    can->Clear();
    can->Divide(2, 1);

    TLegend leg1(0.30, 0.63, 0.47, 0.81);
    leg1.SetBorderSize(0);
    leg1.SetFillStyle(0);
    leg1.SetTextSize(0.035);

    TLegend leg2(0.70, 0.75, 0.79, 0.83);
    leg2.SetBorderSize(0);
    leg2.SetFillStyle(0);
    leg2.SetTextSize(0.035);

    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    tex.SetTextSize(0.050);

    TH1D hUnfoldedPtCopy[nIter];
    RooUnfoldBayes rubUnfoldingPt[nIter];
    TH1D hBayRatiosNth[nIter];
    TH1D hBayRatiosStep[nIter];
    TH1D hBackFoldBayRatios[nIter];
    TH2D hPearsonCoeffs[nIter];

    TH1D hRealDataPtCopy = *((TH1D*)hRealData[iCent].Clone(Form("hRealDataPtCopy%i", iCent)));
    TH2D hResponse1D = *((TH2D *)rurResponse[iCent].Hresponse()->Clone(Form("hResponse1D_%i", iCent)));
    TH1D hResponseReco = *((TH1D *)hResponse1D.ProjectionX()->Clone(Form("hResponseReco_%i", iCent)));
    TH1D hResponseTruth = *((TH1D *)hResponse1D.ProjectionY()->Clone(Form("hResponseTrue_%i", iCent)));
    TH1D hResponseTruthNoFakes = *((TH1D*) rurResponse[iCent].Htruth()->Clone(Form("hTrueResponseTrue_%i", iCent)));

    NormalizeByBinWidth(&hRealDataPtCopy, 2001);
    NormalizeByBinWidth(&hResponseReco, 2002);
    NormalizeByBinWidth(&hResponseTruth, 2003);
    NormalizeByBinWidth(&hResponseTruthNoFakes, 2003);

    PrintCheckNumbers[iCent][0][0] = hRealData[iCent].Integral();

    if (ClosureTest) leg1.AddEntry(&hResponseTruth, "MC true", "lep");
    else {leg1.AddEntry(&hRealDataPtCopy, "Data reco", "lep");}

    leg2.AddEntry(&hRealDataPtCopy, "Data reco", "lep");
    leg2.AddEntry(&hResponseReco, "MC reco (scaled)", "lep");

    //upper pad
    TPad *padRecoComp = new TPad("RecoComparisons", "RecoComparisons", 0.00, 0.5, 1.0, 0.9);
    padRecoComp->SetBottomMargin(0.0);
    padRecoComp->SetTopMargin(0.0);
    padRecoComp->SetBorderMode(0);
    padRecoComp->SetLeftMargin(0.15);

    //bottom pad
    TPad *padRelToNth = new TPad("NextStepRatios", "NextStepRatios", 0.00, 0.29, 1.0, 0.5);
    padRelToNth->SetBottomMargin(0.0);
    padRelToNth->SetTopMargin(0.0);
    padRelToNth->SetBorderMode(0);
    padRelToNth->SetLeftMargin(0.15);

    //middle pad
    TPad *padNextStep = new TPad("RelativeToNthRatios", "RelativeToNthRatios", 0.00, 0.02, 1.0, 0.29);
    padNextStep->SetTopMargin(0.0);
    padNextStep->SetBorderMode(0);
    padNextStep->SetBottomMargin(0.21);
    padNextStep->SetLeftMargin(0.15);


    // procentuální fake podíl: fakes / measured * 100

        TH1D* hFakes = (TH1D*) rurResponse[iCent].Hfakes()->Clone(Form("hFakes_%d", iCent));
        TH1D* hFakesPct = (TH1D*) hFakes->Clone(Form("hFakesPct_%d", iCent));
        hFakesPct->SetTitle(Form("Fake fraction, cent %d", iCent));
        TH1D* hMeasured = (TH1D*) rurResponse[iCent].Hmeasured()->Clone(Form("hMeasured_%d", iCent));
        hFakesPct->Divide(hMeasured);
        TH1D* hOneMinusFake = (TH1D*) hFakesPct->Clone(Form("hOneMinusFake_%d", iCent));
        for (int iBin = 1; iBin <= hOneMinusFake->GetNbinsX(); ++iBin) {
            hOneMinusFake->SetBinContent(iBin, 1.0 - hFakesPct->GetBinContent(iBin));
            hOneMinusFake->SetBinError(iBin, hFakesPct->GetBinError(iBin));
        }

    //left side
    can->cd(1);

        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.16);
        gPad->SetBottomMargin(0.08);
        gPad->SetLeftMargin(0.15);
        gPad->SetLogy();

        if (ClosureTest) {
            hResponseTruthNoFakes.GetYaxis()->SetTitle("dN/dp_{T, Jet} [(GeV/c)^{-1}]");
            hResponseTruthNoFakes.GetXaxis()->SetTitle("p_{T, Jet} [GeV/c]");
            hResponseTruthNoFakes.GetYaxis()->SetTitleOffset(1.5);
            hResponseTruthNoFakes.GetYaxis()->SetTitleSize(0.045);
            hResponseTruthNoFakes.GetXaxis()->SetTitleOffset(0.7);
            hResponseTruthNoFakes.GetXaxis()->SetTitleSize(0.045);
            hResponseTruthNoFakes.SetMarkerStyle(21);
            hResponseTruthNoFakes.Draw();
            hResponseTruthNoFakes.GetXaxis()->SetLabelOffset(-0.015);
        } else {

            hRealDataPtCopy.Draw();

        }

    //right side
    can->cd(2);

        padRecoComp->Draw();
        padNextStep->Draw();
        padRelToNth->Draw();

            //Upper pad
            padRecoComp->cd();

            Double_t padHeight = padRelToNth->GetHNDC(); //high of the pad in NDC coordinates
            Double_t desiredLabelSize = 0.015; //aimed label size in NDC coordinates

            gPad->SetTopMargin(0.15);
            gPad->SetLogy();

            //Real data
            hRealDataPtCopy.GetYaxis()->SetTitle("dN/dp_{T,Jet} [(GeV/c)^{-1}]");
            hRealDataPtCopy.GetXaxis()->SetTitle("p_{T,Jet} [GeV/c]");
            hRealDataPtCopy.Draw();
            hRealDataPtCopy.GetYaxis()->SetLabelSize(desiredLabelSize / padHeight);
            hRealDataPtCopy.GetYaxis()->SetTitleSize(desiredLabelSize / padHeight); // Nastavení velikosti názvu osy X
            hRealDataPtCopy.GetYaxis()->SetTitleOffset(0.80);

            //MC reco
            hResponseReco.Scale(1.*hRealDataPtCopy.Integral("width") / hResponseReco.Integral("width"));
            hResponseReco.Draw("same");

    can->cd();

        leg1.Draw("same");
        leg2.Draw("same");
        tex.SetTextAngle(90);
        tex.DrawLatex(0.978, 0.4, "1D " + centralityTitles[iCent] + " p_{T,Jet}");

    can->cd(1);




            //Unfolding
            for (int iter = 0; iter < nIter; iter++) {

                rubUnfoldingPt[iter] = RooUnfoldBayes(&rurResponse[iCent], &hRealData[iCent], PlotIterations[iter]);
                hUnfoldedPt[iCent][iter] = *(TH1D *) rubUnfoldingPt[iter].Hreco();
                ////hUnfoldedPt[iCent][iter] = (TH1D *) rubUnfoldingPt[iter]->Hunfold();
                hUnfoldedPt[iCent][iter].SetName(Form("hUnfoldedPt_%i_%i", iter, iCent));

                hUnfoldedPtCopy[iter] = *((TH1D *) hUnfoldedPt[iCent][iter].Clone(Form("hUnfoldedPtCopy_%i_%i", iter, iCent)));
                hUnfoldedPtCopy[iter].GetYaxis()->SetTitle("dN/dp_{T,Jet} [(GeV/c)^{-1}]");
                NormalizeByBinWidth(&hUnfoldedPtCopy[iter], 2005 + iter);
                hUnfoldedPtCopy[iter].SetMarkerStyle(27);
                hUnfoldedPtCopy[iter].Draw("same");

                //Set automatic range for unfolded histograms
                if (!ClosureTest){
                    double minAA = min(findMin(hUnfoldedPtCopy, 2000, nIter), hRealDataPtCopy.GetMinimum());
                    double maxAA = max(findMax(hUnfoldedPtCopy, 0, nIter), hRealDataPtCopy.GetMaximum());
                    hRealDataPtCopy.GetYaxis()->SetRangeUser(minAA>0?0.1*minAA:1, 10 * maxAA);
                    if (iter == 0) hResponseTruthNoFakes.Scale(1.*hUnfoldedPtCopy[iter].Integral("width") / hResponseTruthNoFakes.Integral("width"));

                }

                leg1.AddEntry(&hUnfoldedPtCopy[iter], Form("%i iter", PlotIterations[iter]), "lep");

                PrintCheckNumbers[iCent][0][iter + 1] = hUnfoldedPt[iCent][iter].Integral();

                //Backfolding
                hBackfoldedPt[iCent][iter] = *(TH1D*)  rurResponse[iCent].ApplyToTruth(&hUnfoldedPt[iCent][iter]);
                hBackfoldedPt[iCent][iter].Divide(hOneMinusFake);

                NormalizeByBinWidth(&hBackfoldedPt[iCent][iter], 2000 + iter);
                hBackFoldBayRatios[iter] = *((TH1D *) hBackfoldedPt[iCent][iter].Clone(Form("hBackFoldBayRatios_%i_%i", iter, iCent)));
                
                hBackFoldBayRatios[iter].Divide(&hRealDataPtCopy);

            //Pearson coefficients
            hPearsonCoeffs[iter] = getPearsonCoeffs1D_(rubUnfoldingPt[iter].Ereco(RooUnfold::kCovariance));
            hPearsonCoeffs[iter].SetName(Form("hPearsonCoeffs_%i%i", PlotIterations[iter], iCent));

            //Ratios Step by step
            hBayRatiosStep[iter] = *((TH1D *) hUnfoldedPtCopy[iter].Clone(Form("hBayRatiosStep_%i_%i", iter, iCent)));
            hBayRatiosStep[iter].GetYaxis()->SetTitle("Unfolded i-th/Unfolded 1st");
            hBayRatiosStep[iter].GetYaxis()->SetTitleOffset(0.80);

            //Ratio to the first iteration
            hBayRatiosNth[iter] = *(TH1D *) hUnfoldedPtCopy[iter].Clone(Form("hBayRatiosNth_%i_%i", iter, iCent));
            hBayRatiosNth[iter].GetYaxis()->SetTitle("Unfolded i-th/Unfolded 1st");

        }

    if(ClosureTest){
        Double_t Xmax = findMax(hUnfoldedPtCopy, hResponseTruth.GetMaximum(), nIter);
        Double_t Xmin = findMin(hUnfoldedPtCopy, hResponseTruth.GetMinimum(), nIter);
        hResponseTruth.GetYaxis()->SetRangeUser(Xmin > 0 ? Xmin * 0.5 : 0.001, Xmax * 2);
    }

    //Middle pad
    padRelToNth->cd();

        gPad->SetLogy(0);

        //Unfolded ratios
        for (int iter = 0; iter < nIter; iter++) {

            hBayRatiosNth[iter].Draw(iter == 0 ? "" : "same");
            hBayRatiosNth[iter].Divide(&hUnfoldedPtCopy[0]);
            hBayRatiosNth[iter].GetYaxis()->SetTitle("Unfolded/(1st iter)");
            double MinA = findMin(hBayRatiosNth, 2000, nIter);
            double MaxA = findMax(hBayRatiosNth, 0, nIter);
            Double_t padHeight = padRelToNth->GetHNDC(); // výška padu v NDC (0-1)

            Double_t desiredLabelSize = 0.025; // cílená velikost pro plnou výšku 1.0
            hBayRatiosNth[iter].GetYaxis()->SetLabelSize(desiredLabelSize / padHeight);
            hBayRatiosNth[iter].GetXaxis()->SetLabelSize(desiredLabelSize / padHeight);
            hBayRatiosNth[iter].GetXaxis()->SetTitleSize(desiredLabelSize / padHeight*1.1);
            hBayRatiosNth[iter].GetYaxis()->SetTitleSize(desiredLabelSize / padHeight*1.1);
            hBayRatiosNth[iter].GetYaxis()->SetTitleOffset(0.45);
            hBayRatiosNth[iter].GetYaxis()->CenterTitle(); // Center the title
            hBayRatiosNth[iter].GetXaxis()->SetTitleOffset(0.85);
            hBayRatiosNth[0].GetYaxis()->SetRangeUser(0.74, 1.26);

        }

        DrawLineOne();

    //Lower  pad
    padNextStep->cd();

        gPad->SetLogy(0);

        double WorstBinMetric    = 0.0;
        double MeanAbsoluteDrift = 0.0;
        double RMSDrift          = 0.0;
        double WeightedDrift     = 0.0;

        int refIter = nIter - 1;   // nebo GivenIter-1, pokud chces konkretni iteraci

        for (int iter = 1; iter < nIter; iter++) {
            hBayRatiosStep[iter].Draw(iter == 1 ? "" : "same");
            hBayRatiosStep[iter].GetYaxis()->SetRangeUser(0.74, 1.26);
            hBayRatiosStep[iter].Divide(&hUnfoldedPtCopy[iter - 1]);

            hBayRatiosStep[iter].GetYaxis()->SetTitle("i-th/(i-1)-th iteration");
            hBayRatiosStep[iter].GetXaxis()->SetLabelSize(0.05);
            hBayRatiosStep[iter].GetYaxis()->SetLabelSize(0.05);
            hBayRatiosStep[iter].GetXaxis()->SetTitleSize(0.06);
            hBayRatiosStep[iter].GetYaxis()->SetTitleSize(0.06);

            if (iter == refIter) {

                double sumAbsDrift      = 0.0;
                double sumSquaredDrift  = 0.0;
                double sumWeights       = 0.0;
                double sumWeightedDrift = 0.0;
                int countBins           = 0;

                for (int iBin = 1; iBin <= hBayRatiosStep[iter].GetNbinsX(); iBin++) {
                    double x = hBayRatiosStep[iter].GetBinCenter(iBin);
                    if (x <= 5 || x >= 20) continue;
                    

                    double ratio = hBayRatiosStep[iter].GetBinContent(iBin);
                    if (!std::isfinite(ratio) || ratio <= 0) continue;

                    double drift = std::abs(1.0 - ratio);

                    if (drift > WorstBinMetric) WorstBinMetric = drift;

                    sumAbsDrift     += drift;
                    sumSquaredDrift += drift * drift;
                    countBins++;

                    // na zacatek radsi vazit obsahem predchozi iterace
                    double w = hUnfoldedPtCopy[iter - 1].GetBinContent(iBin);
                    if (w > 0) {
                        sumWeights       += w;
                        sumWeightedDrift += drift * w;
                    }
                }

                if (countBins > 0) {
                    MeanAbsoluteDrift = sumAbsDrift / countBins;
                    RMSDrift          = std::sqrt(sumSquaredDrift / countBins);
                }

                if (sumWeights > 0) {
                    WeightedDrift = sumWeightedDrift / sumWeights;
                }
            }
        }



        

        DrawLineOne();

    can->cd();

    DrawTextAbove(0);
    can->SaveAs(outPdf);

    double UnfoldedToMc = 0;
    if (ClosureTest||true){
        can->cd();
        can->Clear();
        can->SetCanvasSize(1200, 1200);

        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.15);
        gPad->SetTopMargin(0.29);

        TH1D hUnfoldeToMc[nIter];

        TLegend leg3(0.20, 0.50, 0.45, 0.70);
        leg3.SetBorderSize(0);
        leg3.SetFillStyle(0);
        leg3.SetTextSize(0.035);
        leg3.SetTextFont(42);

        for (int iter = 0; iter < nIter; iter++) {
            hUnfoldeToMc[iter] = *((TH1D *) hUnfoldedPtCopy[iter].Clone(Form("hUnfoldeToMc_%i_%i", iter, iCent)));
            hUnfoldeToMc[iter].Scale(1. * hResponseTruthNoFakes.Integral("width") / hUnfoldeToMc[iter].Integral("width"));
            hUnfoldeToMc[iter].Divide(&hResponseTruthNoFakes);
            hUnfoldeToMc[iter].SetName(Form("hUnfoldeToMc_%i_%i", iter, iCent));
            hUnfoldeToMc[iter].GetYaxis()->SetTitle("Unfolded/MC true");
            hUnfoldeToMc[iter].GetXaxis()->SetTitle("p_{T,Jet} [GeV/c]");
            hUnfoldeToMc[iter].GetYaxis()->SetTitleOffset(1.3);
            hUnfoldeToMc[iter].GetYaxis()->SetTitleSize(0.04);
            hUnfoldeToMc[iter].GetXaxis()->SetTitleOffset(1.2);
            hUnfoldeToMc[iter].GetXaxis()->SetTitleSize(0.04);
            hUnfoldeToMc[iter].Draw(iter == 0 ? "" : "same");
            hUnfoldeToMc[iter].GetYaxis()->SetRangeUser(0.5, 1.5);
            leg3.AddEntry(&hUnfoldeToMc[iter], Form("%i iter", PlotIterations[iter]), "lep");

            UnfoldedToMc = ComputeTVShapeDistance(hUnfoldedPtCopy[iter], hResponseTruthNoFakes, 5, 20);
        }

        leg3.Draw("same");

        DrawTextAbove(0);
        DrawLineOne();
        can->SaveAs(outPdf);

    }

            fout << runId << "\t"
            << iCent << "\t"
            << "1D" << "\t"
            << "" << "\t"
            << "PT" << "\t"
            << std::fixed << std::setprecision(4)
            << WorstBinMetric * 100    << "\t"
            << MeanAbsoluteDrift * 100 << "\t"
            << RMSDrift * 100          << "\t"
            << WeightedDrift * 100 << "\t"
            << UnfoldedToMc * 100
            << "\n";

    //Response matrix
    if (ResponseMatrix) {

        can->Clear();
        can->SetCanvasSize(1200, 1000);
        can->Divide(2, 2);

        can->cd(1);

        gPad->SetLogz(1);
        gPad->SetRightMargin(0.15);
        gPad->SetTopMargin(0.29);

        TH2D hUnfoldingMatrix = (TH2D)(rubUnfoldingPt[2].UnfoldingMatrix());
        hUnfoldingMatrix.SetTitleOffset(0.5, "Y");
        hUnfoldingMatrix.SetName(Form("UnfoldingMatrix_%i", iCent));
        hUnfoldingMatrix.GetXaxis()->SetTitle("N^{(" + var + ", reco)}_{bin}");
        hUnfoldingMatrix.GetYaxis()->SetTitle("N^{(" + var + ", unfold)}_{bin}");
        hUnfoldingMatrix.GetZaxis()->SetTitleOffset(0.7);
        hUnfoldingMatrix.Draw("colz");
        gPad->Update();

        tex.SetTextAngle(0);
        tex.DrawLatex(0.3, 0.35, centralityTitles[iCent] + "  " + var);
        tex.DrawLatex(0.25, 0.25, "Bin Migration probability (iter 3)");

        can->cd(2);

        gPad->SetRightMargin(0.15);
        gPad->SetTopMargin(0.29);
        gPad->SetLogz(1);

        //Show even overflow and underflow bins
        vector <Double_t> ptRecoBinsVecExtended[nCentralityBins];
        vector <Double_t> ptMcBinsVecExtended[nCentralityBins];

        int underflowVal = -10;
        int overflowVal = 1000;
        for (int j = 0; j < 3; j++) {

            if (UseOverflow) {
                ptRecoBinsVecExtended[j].push_back(underflowVal);
                ptMcBinsVecExtended[j].push_back(underflowVal);
            }

            for (int i = 0; i < ptRecoBinsVec[j].size(); i++) ptRecoBinsVecExtended[j].push_back(ptRecoBinsVec[j][i]);

            for (int i = 0; i < ptMcBinsVecCustom[j].size(); i++)
                ptMcBinsVecExtended[j].push_back(ptMcBinsVecCustom[j][i]);

            if (UseOverflow) {
                ptRecoBinsVecExtended[j].push_back(overflowVal);
                ptMcBinsVecExtended[j].push_back(overflowVal);
            }

        }


        TH2D hResponseRebinned(Form("hResponseRebinned_%d", iCent),
                               "Response Matrix;reco;true", ptRecoBinsVecExtended[iCent].size() - 1,
                               &ptRecoBinsVecExtended[iCent][0], ptMcBinsVecExtended[iCent].size() - 1,
                               &ptMcBinsVecExtended[iCent][0]);

        for (int i = 0; i < (ptRecoBinsVecExtended[iCent].size() - 1); i++) {
            for (int j = 0; j < (ptMcBinsVecExtended[iCent].size() - 1); j++) {
                hResponseRebinned.SetBinContent(i + 1, j + 1,
                                                hResponse1D.GetBinContent(i + !UseOverflow, j + !UseOverflow));
                hResponseRebinned.SetBinError(i + 1, j + 1, hResponse1D.GetBinError(i + !UseOverflow, j + !UseOverflow));
            }
        }

        hResponseRebinned.SetTitle("");
        hResponseRebinned.GetXaxis()->SetTitle("p_{T}^{reco} [GeV/c]");
        hResponseRebinned.GetYaxis()->SetTitle("p_{T}^{true} [GeV/c]");
        hResponseRebinned.Draw("colz");

        tex.DrawLatex(0.3, 0.35, centralityTitles[iCent] + "  " + var);
        tex.DrawLatex(0.25, 0.25, "Response matrix");

        //Draw a rectangle around the center of the matrix
        if (UseOverflow) DrawLineAround(iCent);

        //X projection of RM
        can->cd(3);

        gPad->SetLogy(1);

        TH1D hResponseXCopy = *(TH1D *) hResponseRebinned.ProjectionX("hResponseX")->Clone("hResponseXCopy");
        hResponseXCopy.Sumw2();
        NormalizeByBinWidth(&hResponseXCopy, 3);
        hResponseXCopy.Draw("E");
        hResponseXCopy.GetYaxis()->SetTitle("dN/dp_{T,Jet}");

        TLegend legend(0.6, 0.6, 0.8, 0.8);
        legend.AddEntry(&hResponseXCopy, "reco", "lep");
        legend.SetBorderSize(0);
        legend.SetFillStyle(0);
        legend.SetTextSize(0.08);
        legend.Draw();

        tex.DrawLatex(0.3, 0.35, centralityTitles[iCent] + "  " + var);

        //Y projection of RM
        can->cd(4);

        gPad->SetLogy(1);

        TH1D hResponseYCopy = *(TH1D *) hResponseRebinned.ProjectionY("hResponseY")->Clone("hResponseYCopy");
        hResponseYCopy.Sumw2();
        NormalizeByBinWidth(&hResponseYCopy, 4);
        hResponseYCopy.Draw("E");
        hResponseYCopy.GetYaxis()->SetTitle("dN/dp_{T,Jet}");

        TLegend legend2(0.6, 0.6, 0.8, 0.8);
        legend2.AddEntry(&hResponseYCopy, "true", "lep");
        legend2.SetBorderSize(0);
        legend2.SetFillStyle(0);
        legend2.SetTextSize(0.08);
        legend2.Draw();

        can->cd();

        DrawTextAbove(0);

        can->SaveAs(outPdf);
        can->Clear();

    }

    //Pearson coefficients
    if (PearsonCoeff) {
        can->SetCanvasSize(1200, 1000);
        can->Clear();
        can->cd();

        gPad->SetTopMargin(0.15);
        gPad->SetRightMargin(0.2);
        gPad->SetLeftMargin(0.2);
        can->SetLogz(0);

        TH2D hPearsRebinned[nIter];

        //Show even overflow and underflow bins
        vector <Double_t> ptMcBinsVecExtended[nCentralityBins];

        int underflowVal = -10;
        int overflowVal = 1000;
        for (int j = 0; j < 3; j++) {

            if (UseOverflow) {
                ptMcBinsVecExtended[j].push_back(underflowVal);
            }

            for (int i = 0; i < ptMcBinsVecCustom[j].size(); i++)
                ptMcBinsVecExtended[j].push_back(ptMcBinsVecCustom[j][i]);

            if (UseOverflow) {
                ptMcBinsVecExtended[j].push_back(overflowVal);
            }

        }

        can->SetCanvasSize(3 * 400, (TMath::Ceil(nIter / 3.) + 1) * 400);
        can->Divide(3, TMath::Ceil(nIter / 3.) + 1);

        for (int iter = 0; iter < nIter; iter++) {
            can->cd(iter + 1 + 3);
            gPad->SetRightMargin(0.15);
            gPad->SetLeftMargin(0.15);
            hPearsRebinned[iter] = TH2D(Form("hPearsRebinned_%.d_%.d", iCent, iter),
                                        "Response Matrix;reco;true", ptMcBinsVecExtended[iCent].size() - 1,
                                        &ptMcBinsVecExtended[iCent][0], ptMcBinsVecExtended[iCent].size() - 1,
                                        &ptMcBinsVecExtended[iCent][0]);
            //přendám biny
            for (int i = 0; i < (ptMcBinsVecExtended[iCent].size() - 1); i++) {
                for (int j = 0; j < (ptMcBinsVecExtended[iCent].size() - 1); j++) {
                    hPearsRebinned[iter].SetBinContent(i + 1, j + 1, hPearsonCoeffs[iter].GetBinContent(i + 1, j + 1));
                }
            }

            hPearsRebinned[iter].SetTitle((TString) ";p_{T,Jet};p_{T,Jet}");
            hPearsRebinned[iter].Draw("colz");
            hPearsRebinned[iter].GetZaxis()->SetRangeUser(-1, 1);
            gPad->Update();
            hPearsRebinned[iter].GetZaxis()->SetTitleOffset(0.7);
            tex.DrawLatex(0.25, 0.65, centralityTitles[iCent] + "  " + var);
            tex.DrawLatex(0.25, 0.35, Form("Pearson Coefficients Iter %i", PlotIterations[iter]));
            if (UseOverflow) DrawLineAround2(iCent);
        }

        can->cd();
        DrawTextAbove(0);
        can->SaveAs(outPdf);

    }

    //Backfold
    if (true) {

        can->Clear();
        can->SetCanvasSize(1200, 1000);

        can->Divide(2, 1);

        can->cd(1);

        gPad->SetTopMargin(0.15);

        TPad *padBackfoldDistr = new TPad("padBackfoldDistr", "padBackfoldDistr", 0, 0.4, 1.0, 0.92);
        padBackfoldDistr->SetBottomMargin(0.0);
        padBackfoldDistr->SetBorderMode(0);
        TPad *padBackfoldRatios = new TPad("padBackfoldRatios", "padBackfoldRatios", 0, 0.09, 1.0, 0.4);
        padBackfoldRatios->SetTopMargin(0.0);
        padBackfoldRatios->SetBorderMode(0);

        padBackfoldDistr->Draw();
        padBackfoldRatios->Draw();

        padBackfoldDistr->cd();

        hRealDataPtCopy.GetYaxis()->SetTitle("dN/dp_{T,Jet}");
        hRealDataPtCopy.GetXaxis()->SetTitle("p_{T,Jet} [GeV/c]");

        gPad->SetLogy();

        hRealDataPtCopy.Draw();
        hRealDataPtCopy.SetMarkerStyle(20);
        hRealDataPtCopy.SetMarkerSize(1);
        hRealDataPtCopy.SetMarkerColor(kBlack);
        hRealDataPtCopy.SetLineColor(kBlack);

        TLegend legBF(0.65, 0.55, 1.0, 0.82);
        legBF.SetBorderSize(0);
        legBF.SetFillStyle(0);
        legBF.AddEntry(&hRealDataPtCopy, "Real Data", "lep");

        for (int iter = 0; iter < nIter; iter++) {
            hBackfoldedPt[iCent][iter].Draw("same");
            legBF.AddEntry(&hBackfoldedPt[iCent][iter], Form("Iter%i", PlotIterations[iter]), "lep");
        }

        legBF.Draw("same");

        padBackfoldRatios->cd();

        for (int iter = 0; iter < nIter; iter++) {
            hBackFoldBayRatios[iter].Draw(iter == 0 ? "" : "same");
            hBackFoldBayRatios[iter].GetYaxis()->SetRangeUser(0.5, 1.5);
            hBackFoldBayRatios[iter].GetYaxis()->SetTitleOffset(0.20);
        }

        DrawLineOne2(ptRecoBinsVec[iCent][0], ptRecoBinsVec[iCent][ptRecoBinsVec[iCent].size() - 1]);

        can->cd(2);
        gPad->SetLogy(0);

        gPad->SetLeftMargin(0.25);
        TH1D hChiSquaredNdf("hChiSquaredNdf", "hChiSquaredNdf", nIter, 0, nIter);

        for (int iter = 0; iter < nIter; iter++) {

            double chi2ndf = 0;
            for (int i = 1; i <= hRealDataPtCopy.GetNbinsX(); i++) {
                if (hRealDataPtCopy.GetBinError(i) == 0) continue;
                chi2ndf += (hRealDataPtCopy.GetBinContent(i) - hBackfoldedPt[iCent][iter].GetBinContent(i)) *
                           (hRealDataPtCopy.GetBinContent(i) - hBackfoldedPt[iCent][iter].GetBinContent(i)) /
                           (hRealDataPtCopy.GetBinError(i) * hRealDataPtCopy.GetBinError(i));
            }

            chi2ndf /= hRealDataPtCopy.GetNbinsX();
            hChiSquaredNdf.SetBinContent(iter + 1, chi2ndf);
            hChiSquaredNdf.GetXaxis()->SetBinLabel(iter + 1, Form("Iter %i", PlotIterations[iter]));

        }

        hChiSquaredNdf.Draw("ph");
        hChiSquaredNdf.GetYaxis()->SetTitle("#chi^{2}/ndf = #frac{1}{N} #sum_{i}^{N} #frac{(Data_{i} - Backfolded_{i})^{2}}{Error_{i}^{2}}");
        hChiSquaredNdf.GetYaxis()->SetTitleOffset(1.6);

        can->SaveAs(outPdf);
        can->Clear();
    }


    if (true) {

        TH1D* hTruth = (TH1D*) rurResponse[iCent].Htruth()->Clone(
            Form("hTruth_%d", iCent));

        TH2D* hResp = (TH2D*) rurResponse[iCent].Hresponse();

        // matched truth = projekce response na truth osu
        TH1D* hMatchedTruth = (TH1D*) hResp->ProjectionY(
            Form("hMatchedTruth_%d", iCent));

        // misses = all truth - matched truth
        TH1D* hMisses = (TH1D*) hTruth->Clone(Form("hMisses_%d", iCent));
        hMisses->Add(hMatchedTruth, -1.0);


        // procentuální miss podíl: misses / truth * 100
        TH1D* hMissesPct = (TH1D*) hMisses->Clone(Form("hMissesPct_%d", iCent));
        hMissesPct->SetTitle(Form("Miss fraction, cent %d", iCent));
        hMissesPct->Divide(hTruth);
        //hMissesPct->Scale(100.0);

        can->Divide(2,1);

        // Fakes %
        can->cd(1);
        gPad->SetLogy(0);
        hFakesPct->GetXaxis()->SetTitle("p_{T,Jet}^{reco}");
        hFakesPct->GetYaxis()->SetTitle("Fake fraction");
        hFakesPct->SetLineColor(kRed+1);
        hFakesPct->SetMinimum(0.0);
        hFakesPct->SetMaximum(1.0);
        hFakesPct->Draw("hist");

        // Misses %
        can->cd(2);
        gPad->SetLogy(0);
        hMissesPct->GetXaxis()->SetTitle("p_{T,Jet}^{true}");
        hMissesPct->GetYaxis()->SetTitle("Miss fraction");
        hMissesPct->SetLineColor(kBlue+1);
        hMissesPct->SetMinimum(0.0);
        hMissesPct->SetMaximum(1.0);
        hMissesPct->Draw("hist");

        can->SaveAs(outPdf);
        can->Clear();

        delete hTruth;
        delete hMatchedTruth;
        delete hMisses;
        delete hMissesPct;

    }
        delete hFakes;
        delete hMeasured;
        delete hFakesPct;
        delete hOneMinusFake;

}

void plotComparison2D(TCanvas *can, const Int_t &iCent, TString var) {

    int variable = getVariable(var) - 1;
    can->Clear();

    can->Divide(2, 2);

    TLegend leg1(0.28, 0.70, 0.43, 0.92);
    leg1.SetBorderSize(0);
    leg1.SetFillStyle(0);

    TLegend leg2(0.50, 0.80, 0.85, 0.87);
    leg2.SetBorderSize(0);
    leg2.SetFillStyle(0);
    leg2.SetTextSize(0.035);


    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    tex.SetTextSize(0.055);

    TH1D hUnfolded2DCopy_X[nIter];
    TH1D hUnfolded2DCopy_Y[nIter];
    RooUnfoldBayes rubUnfolding2D[nIter];
    TH1D hBayRatiosNth_X[nIter];
    TH1D hBayRatiosNth_Y[nIter];
    TH1D hBayRatiosStep_X[nIter];
    TH1D hBayRatiosStep_Y[nIter];
    TH1D hBackFoldBayRatios_X[nIter];
    TH1D hBackFoldBayRatios_Y[nIter];

    hRealData2D[iCent][variable].GetXaxis()->SetTitle("p_{T,Jet} [GeV/c]");
    hRealData2D[iCent][variable].GetYaxis()->SetTitle(var);

    TH1D hRealDataCopyX = *(TH1D *) hRealData2D[iCent][variable].ProjectionX()->Clone(
            Form("hRealDataCopyX_%i_%i", iCent, variable));
    TH1D hRealDataCopyY = *(TH1D *) hRealData2D[iCent][variable].ProjectionY()->Clone(
            Form("hRealDataCopyY_%i_%i", iCent, variable));

    TH2D hResponseReco2D = *(TH2D *) rurResponse2D[iCent][variable].Hmeasured()->Clone(
            Form("hResponseReco2D_%i_%i", iCent, variable));
    hResponseReco2D.GetXaxis()->SetTitle("p_{T,Jet}^{reco} [GeV/c]");
    hResponseReco2D.GetYaxis()->SetTitle(var+" ^{,reco}");
    TH2D hResponseTruth2D = *(TH2D *) rurResponse2D[iCent][variable].Htruth()->Clone(
            Form("hResponseTruth2D_%i_%i", iCent, variable));
    hResponseTruth2D.GetXaxis()->SetTitle("p_{T,Jet}^{true} [GeV/c]");
    hResponseTruth2D.GetYaxis()->SetTitle(var+" ^{,true}");

    TH1D hResponseReco2D_X = *(TH1D *) hResponseReco2D.ProjectionX()->Clone(
            Form("hResponseReco2D_X_%i_%i", iCent, variable));
    TH1D hResponseReco2D_Y = *(TH1D *) hResponseReco2D.ProjectionY()->Clone(
            Form("hResponseReco2D_Y_%i_%i", iCent, variable));

    TH1D hResponseTruth2D_X = *(TH1D *) hResponseTruth2D.ProjectionX()->Clone(
            Form("hResponseTruth2D_X%i_%i", iCent, variable));
    TH1D hResponseTruth2D_Y = *(TH1D *) hResponseTruth2D.ProjectionY()->Clone(
            Form("hResponseTruth2D_Y%i_%i", iCent, variable));

    NormalizeByBinWidth(&hRealDataCopyX, 2001);
    NormalizeByBinWidth(&hResponseTruth2D_X, 2002);
    NormalizeByBinWidth(&hResponseReco2D_X, 2003);

    NormalizeByBinWidth(&hRealDataCopyY, 2001);
    NormalizeByBinWidth(&hResponseTruth2D_Y, 2002);
    NormalizeByBinWidth(&hResponseReco2D_Y, 2003);

    leg2.AddEntry(&hRealDataCopyX, "Real Data", "lp");
    if (ClosureTest) leg1.AddEntry(&hResponseTruth2D_X, "Mc (scaled)", "lp");
    else leg1.AddEntry((TH1D *) 0, "Real unfolded:", "");
    leg2.AddEntry(&hResponseReco2D_X, "Mc Reco (scaled)", "lp");

    //Left Up
    can->cd(1);
    gPad->SetLeftMargin(0.15);

    TPad *padPtRatioStep = new TPad("padPtRatioStep", "padPtRatioStep", 0, 0.08, 1.0, 0.5);
    padPtRatioStep->SetTopMargin(0.0);

    TPad *padPtRatioNth = new TPad("padPtRatioNth", "padPtRatioNth", 0, 0.5, 1.0, 0.92);
    padPtRatioNth->SetBottomMargin(0.0);
    padPtRatioNth->SetBorderMode(0);

    TPad *padVarRatioStep = new TPad("padVarRatioStep", "padVarRatioStep", 0, 0.08, 1.0, 0.5);
    padVarRatioStep->SetTopMargin(0.0);

    TPad *padVarRatioNth = new TPad("padVarRatioNth", "padVarRatioNth", 0, 0.5, 1.0, 0.92);
    padVarRatioNth->SetBottomMargin(0.0);
    padVarRatioNth->SetBorderMode(0);

    gPad->SetLogy();

    if (ClosureTest) {
        hResponseTruth2D_X.SetMarkerStyle(21);
        hResponseTruth2D_X.Scale(1. * hRealDataCopyX.Integral("width") / hResponseTruth2D_X.Integral("width"));
        hResponseTruth2D_X.Draw();
        hResponseTruth2D_X.GetYaxis()->SetTitle("dN/dp_{T}^{true}");
        hResponseTruth2D_X.GetYaxis()->SetTitleOffset(1.3);

    } else {
        hRealDataCopyX.Draw();
        hRealDataCopyX.GetYaxis()->SetTitle("dN/dp_{T}^{true}");

    }

    //Right Up
    can->cd(2);

    padPtRatioStep->Draw();
    padPtRatioNth->Draw();

    //Left Down
    can->cd(3);
    gPad->SetLeftMargin(0.15);

    gPad->SetLogy();
    hResponseTruth2D_Y.SetMarkerStyle(21);

    if (ClosureTest) {
        hResponseTruth2D_Y.Draw();
        hResponseTruth2D_Y.GetXaxis()->SetTitle(var+"^{, true}");
        hResponseTruth2D_Y.GetYaxis()->SetTitle("dN/d"+var+"^{, true}");
        hResponseTruth2D_Y.GetYaxis()->SetTitleOffset(1.3);
        hResponseTruth2D_Y.Scale(1. * hRealDataCopyY.Integral("width") / hResponseTruth2D_Y.Integral("width"));
    } else {
        hRealDataCopyY.Draw();
        //set x axis title
        hRealDataCopyY.GetXaxis()->SetTitle(var);
        hRealDataCopyY.GetXaxis()->SetTitle(var+"^{, true}");
        hRealDataCopyY.GetYaxis()->SetTitle("dN/d"+var+"^{, true}");
        hRealDataCopyY.GetYaxis()->SetTitleOffset(1.3);


    }

    //Right Down
    can->cd(4);

    padVarRatioStep->Draw();
    padVarRatioNth->Draw();

    //The main one
    can->cd();

    PrintCheckNumbers[iCent][variable][0] = hRealData2D[iCent][variable].Integral();

    leg1.Draw("same");

    for (int iter = 0; iter < nIter; iter++) {

        rubUnfolding2D[iter] = RooUnfoldBayes(&rurResponse2D[iCent][variable], &hRealData2D[iCent][variable],
                                              PlotIterations[iter]);
        ////hUnfolded2D[iCent][variable][iter] = *(TH2D *) rubUnfolding2D[iCent].Hreco();
        ////hUnfolded2D[iter] = (TH2D *) rubUnfolding2D.Hunfold();
        TH1 *hrecoBase = rubUnfolding2D[iter].Hreco();
        TH2D *hrecoTemp = dynamic_cast<TH2D *>(hrecoBase);
        hUnfolded2D[iCent][variable][iter] = *(TH2D *) hrecoTemp->Clone(
                Form("hUnfolded2D_%i_%i_%i", iCent, variable, iter));
        hUnfolded2D[iCent][variable][iter].GetXaxis()->SetTitle("p_{T,Jet} [GeV/c]");
        hUnfolded2D[iCent][variable][iter].GetYaxis()->SetTitle(var);

        PrintCheckNumbers[iCent][variable][iter + 1] = hUnfolded2D[iCent][variable][iter].Integral();
        hUnfolded2DCopy_X[iter] = *((TH1D *) hUnfolded2D[iCent][variable][iter].ProjectionX()->Clone(
                Form("hUnfoldedPtCopy_X_%i_%i", iter, iCent)));
        hUnfolded2DCopy_Y[iter] = *((TH1D *) hUnfolded2D[iCent][variable][iter].ProjectionY()->Clone(
                Form("hUnfoldedPtCopy_X_%i_%i", iter, iCent)));

        NormalizeByBinWidth(&hUnfolded2DCopy_X[iter], 2005 + iter);
        NormalizeByBinWidth(&hUnfolded2DCopy_Y[iter], 2005 + iter);

        hUnfolded2DCopy_X[iter].SetMarkerStyle(27);
        hUnfolded2DCopy_X[iter].GetYaxis()->SetTitle("dN/dp_{T,Jet}");
        hUnfolded2DCopy_Y[iter].SetMarkerStyle(27);
        hUnfolded2DCopy_Y[iter].GetYaxis()->SetTitle("dN/d" + var);
        hUnfolded2DCopy_Y[iter].GetXaxis()->SetTitle(var);

        leg1.AddEntry(&hUnfolded2DCopy_X[iter], Form("Iter%i", PlotIterations[iter]), "lep");

        can->cd(1);

        hUnfolded2DCopy_X[iter].Draw("same");

        can->cd(3);

 
        hUnfolded2DCopy_Y[iter].Draw(iter==0?"":"same");
        if(iter == 0)        hRealDataCopyY.Draw("same");


        hBayRatiosNth_X[iter] = *(TH1D *) hUnfolded2DCopy_X[iter].Clone(
                Form("hBayRatiosNth_X_%i_%i_", iter, iCent) + var);
        hBayRatiosNth_X[iter].GetYaxis()->SetTitle("Unfolded i-th/Unfolded 1st");

        hBayRatiosStep_X[iter] = *(TH1D *) hUnfolded2DCopy_X[iter].Clone(
                Form("hBayRatiosStep_X_%i_%i_", iter, iCent) + var);
        hBayRatiosStep_X[iter].GetYaxis()->SetTitle("Unfolded i-th/Unfolded (i-1)-th");

        hBayRatiosNth_Y[iter] = *(TH1D *) hUnfolded2DCopy_Y[iter].Clone(
                Form("hBayRatiosNth_Y_%i_%i_", iter, iCent) + var);
        hBayRatiosNth_Y[iter].GetYaxis()->SetTitle("Unfolded i-th/Unfolded 1st");

        hBayRatiosStep_Y[iter] = *(TH1D *) hUnfolded2DCopy_Y[iter].Clone(
                Form("hBayRatiosStep_Y_%i_%i_", iter, iCent) + var);
        hBayRatiosStep_Y[iter].GetYaxis()->SetTitle("Unfolded i-th/Unfolded (i-1)-th");


    }

    if (!ClosureTest) {
        Double_t Xmax = findMax(hUnfolded2DCopy_X, hRealDataCopyX.GetMaximum(), nIter);
        Double_t Xmin = findMin(hUnfolded2DCopy_X, hRealDataCopyX.GetMinimum()==0?1000:hRealDataCopyX.GetMinimum(), nIter);
        hRealDataCopyX.GetYaxis()->SetRangeUser(Xmin > 0 ? Xmin * 0.5 : 1, Xmax * 2);
        Double_t Xmax2 = findMax(hUnfolded2DCopy_Y, hRealDataCopyY.GetMaximum(), nIter);
        Double_t Xmin3 = findMin(hRealDataCopyY,hRealDataCopyY.GetMaximum());
        Double_t Xmin2 = findMin(hUnfolded2DCopy_Y, Xmin3, nIter);
        hRealDataCopyY.GetYaxis()->SetRangeUser(Xmin2 > 0 ? Xmin2 * 0.5 : 1, Xmax2 * 2);


    }


    padPtRatioStep->cd();

double WorstBinMetric_X    = 0.0;
double MeanAbsoluteDrift_X = 0.0;
double RMSDrift_X          = 0.0;
double WeightedDrift_X     = 0.0;

int refIterX = nIter - 1;   // nebo GivenIter-1

for (int iter = nIter - 1; iter >= 1; iter--) {

    hBayRatiosStep_X[iter].Draw(iter == (nIter - 1) ? "" : "same");
    hBayRatiosStep_X[iter].GetYaxis()->SetRangeUser(0.77, 1.23);
    hBayRatiosStep_X[iter].Divide(&hUnfolded2DCopy_X[iter - 1]);

    hBayRatiosStep_X[iter].GetYaxis()->SetTitleSize(0.09);
    hBayRatiosStep_X[iter].GetYaxis()->CenterTitle();
    hBayRatiosStep_X[iter].GetYaxis()->SetLabelSize(0.08);
    hBayRatiosStep_X[iter].GetYaxis()->SetNdivisions(505);
    hBayRatiosStep_X[iter].GetYaxis()->SetTitleOffset(0.5);
    hBayRatiosStep_X[iter].GetXaxis()->SetTitleSize(0.09);
    hBayRatiosStep_X[iter].GetXaxis()->SetLabelSize(0.08);

    if (iter == refIterX) {

        double sumAbsDrift_X      = 0.0;
        double sumSquaredDrift_X  = 0.0;
        double sumWeights_X       = 0.0;
        double sumWeightedDrift_X = 0.0;
        int countBins_X           = 0;

        for (int iBin = 1; iBin <= hBayRatiosStep_X[iter].GetNbinsX(); iBin++) {
            double x = hBayRatiosStep_X[iter].GetBinCenter(iBin);
            if (x <= 5 || x >= 20) continue;

            double ratio = hBayRatiosStep_X[iter].GetBinContent(iBin);
            if (!std::isfinite(ratio) || ratio <= 0.0) continue;

            double drift = std::abs(1.0 - ratio);

            if (drift > WorstBinMetric_X) WorstBinMetric_X = drift;

            sumAbsDrift_X     += drift;
            sumSquaredDrift_X += drift * drift;
            countBins_X++;

            double w = hUnfolded2DCopy_X[iter - 1].GetBinContent(iBin);
            if (w > 0.0) {
                sumWeights_X       += w;
                sumWeightedDrift_X += drift * w;
            }
        }

        if (countBins_X > 0) {
            MeanAbsoluteDrift_X = sumAbsDrift_X / countBins_X;
            RMSDrift_X          = std::sqrt(sumSquaredDrift_X / countBins_X);
        }

        if (sumWeights_X > 0.0) {
            WeightedDrift_X = sumWeightedDrift_X / sumWeights_X;
        }
    }
}



    DrawLineOne2(hBayRatiosStep_X[0].GetBinLowEdge(1),
                 hBayRatiosStep_X[0].GetBinLowEdge(hBayRatiosStep_X[0].GetNbinsX() + 1));

    padPtRatioNth->cd();

    for (int iter = nIter - 1; iter >= 1; iter--) {

        hBayRatiosNth_X[iter].Draw(iter == (nIter - 1) ? "" : "same");
        hBayRatiosNth_X[iter].GetYaxis()->SetRangeUser(0.77, 1.23);
        hBayRatiosNth_X[iter].Divide(&hUnfolded2DCopy_X[0]);

        hBayRatiosNth_X[iter].GetYaxis()->SetTitleSize(0.09);
        hBayRatiosNth_X[iter].GetYaxis()->CenterTitle();
        hBayRatiosNth_X[iter].GetYaxis()->SetLabelSize(0.08);
        hBayRatiosNth_X[iter].GetYaxis()->SetNdivisions(505);
        hBayRatiosNth_X[iter].GetYaxis()->SetTitleOffset(0.5);
        hBayRatiosNth_X[iter].GetXaxis()->SetTitleSize(0.09);
        hBayRatiosNth_X[iter].GetXaxis()->SetLabelSize(0.08);

    }

    DrawLineOne2(hBayRatiosNth_X[0].GetBinLowEdge(1),
                 hBayRatiosNth_X[0].GetBinLowEdge(hBayRatiosNth_X[0].GetNbinsX() + 1));

    padVarRatioStep->cd();

double WorstBinMetric_Y    = 0.0;
double MeanAbsoluteDrift_Y = 0.0;
double RMSDrift_Y          = 0.0;
double WeightedDrift_Y     = 0.0;

int refIterY = nIter - 1;   // nebo GivenIter-1

for (int iter = nIter - 1; iter >= 1; iter--) {

    hBayRatiosStep_Y[iter].Draw(iter == (nIter - 1) ? "" : "same");
    hBayRatiosStep_Y[iter].GetYaxis()->SetRangeUser(0.77, 1.23);
    hBayRatiosStep_Y[iter].Divide(&hUnfolded2DCopy_Y[iter - 1]);

    hBayRatiosStep_Y[iter].GetYaxis()->SetTitleSize(0.09);
    hBayRatiosStep_Y[iter].GetYaxis()->CenterTitle();
    hBayRatiosStep_Y[iter].GetYaxis()->SetLabelSize(0.08);
    hBayRatiosStep_Y[iter].GetYaxis()->SetNdivisions(505);
    hBayRatiosStep_Y[iter].GetYaxis()->SetTitleOffset(0.5);
    hBayRatiosStep_Y[iter].GetXaxis()->SetTitleSize(0.09);
    hBayRatiosStep_Y[iter].GetXaxis()->SetLabelSize(0.08);

    if (iter == refIterY) {

        double sumAbsDrift_Y      = 0.0;
        double sumSquaredDrift_Y  = 0.0;
        double sumWeights_Y       = 0.0;
        double sumWeightedDrift_Y = 0.0;
        int countBins_Y           = 0;

        for (int iBin = 1; iBin <= hBayRatiosStep_Y[iter].GetNbinsX(); iBin++) {

            double ratio = hBayRatiosStep_Y[iter].GetBinContent(iBin);
            if (!std::isfinite(ratio) || ratio <= 0.0) continue;

             //pro z a momdips uvažujeme pouze biny nad 0.2

             if (var == "z" || var == "p_{T}^{D}") {
                 double x = hBayRatiosStep_Y[iter].GetBinCenter(iBin);
                 if (x <= 0.2) continue;
             } else //skip the last bin
                if (iBin == hBayRatiosStep_Y[iter].GetNbinsX()) continue;


            double drift = std::abs(1.0 - ratio);

            if (drift > WorstBinMetric_Y) WorstBinMetric_Y = drift;

            sumAbsDrift_Y     += drift;
            sumSquaredDrift_Y += drift * drift;
            countBins_Y++;

            double w = hUnfolded2DCopy_Y[iter - 1].GetBinContent(iBin);
            if (w > 0.0) {
                sumWeights_Y       += w;
                sumWeightedDrift_Y += drift * w;
            }
        }

        if (countBins_Y > 0) {
            MeanAbsoluteDrift_Y = sumAbsDrift_Y / countBins_Y;
            RMSDrift_Y          = std::sqrt(sumSquaredDrift_Y / countBins_Y);
        }

        if (sumWeights_Y > 0.0) {
            WeightedDrift_Y = sumWeightedDrift_Y / sumWeights_Y;
        }
    }
}



    DrawLineOne2(hBayRatiosStep_Y[0].GetBinLowEdge(1),
                 hBayRatiosStep_Y[0].GetBinLowEdge(hBayRatiosStep_Y[0].GetNbinsX() + 1));

    padVarRatioNth->cd();

    for (int iter = nIter - 1; iter >= 1; iter--) {

        hBayRatiosNth_Y[iter].Draw(iter == (nIter - 1) ? "" : "same");
        hBayRatiosNth_Y[iter].GetYaxis()->SetRangeUser(0.77, 1.23);
        hBayRatiosNth_Y[iter].Divide(&hUnfolded2DCopy_Y[0]);

        hBayRatiosNth_Y[iter].GetYaxis()->SetTitleSize(0.09);
        hBayRatiosNth_Y[iter].GetYaxis()->CenterTitle();
        hBayRatiosNth_Y[iter].GetYaxis()->SetLabelSize(0.08);
        hBayRatiosNth_Y[iter].GetYaxis()->SetNdivisions(505);
        hBayRatiosNth_Y[iter].GetYaxis()->SetTitleOffset(0.5);
        hBayRatiosNth_Y[iter].GetXaxis()->SetTitleSize(0.09);
        hBayRatiosNth_Y[iter].GetXaxis()->SetLabelSize(0.08);

    }

    DrawLineOne2(hBayRatiosNth_Y[0].GetBinLowEdge(1),
                 hBayRatiosNth_Y[0].GetBinLowEdge(hBayRatiosNth_Y[0].GetNbinsX() + 1));


    can->cd();

    tex.DrawLatex(0.38, 0.48, "(p_{T}, " + var + ") " + centralityTitles[iCent]);

    can->SaveAs(outPdf);
    can->Clear();
    gPad->SetLogz(0);


    //2D iteration ratios

    TH2D hBayRatios2D[nIter];
    for (int iter = 1; iter < nIter; iter++) {

        hBayRatios2D[iter] = *(TH2D *) hUnfolded2D[iCent][variable][iter].Clone(
                Form("hBayRatios2D_%i_%i_%i", iCent, variable, iter));
        hBayRatios2D[iter].Divide(&hUnfolded2D[iCent][variable][iter-1]);
        hBayRatios2D[iter].GetXaxis()->SetTitle("p_{T,Jet} [GeV/c]");
        hBayRatios2D[iter].GetYaxis()->SetTitle(var);
    }
    gStyle->SetPaintTextFormat("4.2f");  // šířka 4, 2 desetinná místa

    //rozdělím canvas na 5 oken
    can->SetCanvasSize(3 * 400, (TMath::Ceil(nIter / 3.) + 1) * 400);
    can->Divide(3, TMath::Ceil(nIter / 3.) + 1);
    for (int iter = 1; iter < nIter; iter++) {
        can->cd(iter + 1 + 3);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.15);
        hBayRatios2D[iter].Draw("colztext");
        hBayRatios2D[iter].GetZaxis()->SetRangeUser(0, 2);

        gPad->Update();
        hBayRatios2D[iter].GetZaxis()->SetTitleOffset(0.7);
        tex.DrawLatex(0.25, 0.05, centralityTitles[iCent] + "  " + var);
        tex.DrawLatex(0.25, 0.15, Form("Bayes Ratio Iter %i/%i", PlotIterations[iter], PlotIterations[iter-1]));
    }


    can->SaveAs(outPdf);
    can->Clear();

    //Distributions Real vs Mc
    if(!ClosureTest && false) {
        can->SetCanvasSize(800, 400);
        can->Divide(2, 1);
    } else {
        can->SetCanvasSize(1200, 1200);
        can->Divide(2, 2);
    }


    can->cd(1);

    gPad->SetLeftMargin(0.18);
    gPad->SetLogy();

    hRealDataCopyX.GetYaxis()->SetTitle("dN/dp_{T,Jet}");
    hRealDataCopyX.GetYaxis()->SetTitleOffset(1.5);
    hRealDataCopyX.GetXaxis()->SetTitle("p_{T,Jet} [GeV/c]");
    hRealDataCopyX.Draw();
    hResponseReco2D_X.Scale(hRealDataCopyX.Integral("width") / hResponseReco2D_X.Integral("width"));
    hResponseReco2D_X.Draw("same");
    leg2.Draw("same");

    can->cd(2);

    gPad->SetLeftMargin(0.18);
    gPad->SetLogy();

    hRealDataCopyY.GetYaxis()->SetTitle("dN/d" + var + "");
    hRealDataCopyY.GetYaxis()->SetTitleOffset(1.5);
    hRealDataCopyY.GetXaxis()->SetTitleOffset(1.1);
    hRealDataCopyY.GetXaxis()->SetTitle(var);
    hRealDataCopyY.Draw();
    hResponseReco2D_Y.Scale(hRealDataCopyY.Integral("width") / hResponseReco2D_Y.Integral("width"));
    hResponseReco2D_Y.Draw("same");
    leg2.Draw("same");

    Double_t Xmin = TMath::Min(hRealDataCopyY.GetMinimum(), hResponseReco2D_Y.GetMinimum());
    Double_t Xmax = TMath::Max(hRealDataCopyY.GetMaximum(), hResponseReco2D_Y.GetMaximum());
    hRealDataCopyY.GetYaxis()->SetRangeUser(Xmin > 0 ? Xmin * 0.5 : 1, Xmax * 2);

    TH1D hUnfoldeToMc_X[nIter];
    TH1D hUnfoldeToMc_Y[nIter];

    Double_t UnfoldedToMc_x=0.0;
    Double_t UnfoldedToMc_y=0.0;
    if (ClosureTest || true){

        can->cd(3);
        gPad->SetLeftMargin(0.18);
        gPad->SetLogy(0);

        for (Int_t iter = 0; iter < nIter; iter++) {

            hUnfoldeToMc_X[iter] = *(TH1D *) hUnfolded2DCopy_X[iter].Clone(
                    Form("hUnfoldeToMc_X_%i_%i", iter, iCent) + var);
            hUnfoldeToMc_X[iter].Divide(&hResponseTruth2D_X);
            hUnfoldeToMc_X[iter].Scale(1./hUnfolded2DCopy_X[iter].Integral()*hResponseTruth2D_X.Integral());
            hUnfoldeToMc_X[iter].GetYaxis()->SetTitle("Unfolded/MC True");
            hUnfoldeToMc_X[iter].GetXaxis()->SetTitle("p_{T,Jet} [GeV/c]");
            hUnfoldeToMc_X[iter].GetYaxis()->SetTitleOffset(1.5);
            hUnfoldeToMc_X[iter].Draw(iter == 0 ? "" : "same");
            hUnfoldeToMc_X[iter].GetYaxis()->SetRangeUser(0.5,1.5);

            if(iter == 5) UnfoldedToMc_x = ComputeTVShapeDistance(hUnfolded2DCopy_X[5], hResponseTruth2D_X, 5.0, 20.0);
        }


        DrawLineOne();

        can->cd(4);
        gPad->SetLeftMargin(0.18);
        gPad->SetLogy(0);
        for (Int_t iter = 0; iter < nIter; iter++) {

            hUnfoldeToMc_Y[iter] = *(TH1D *) hUnfolded2DCopy_Y[iter].Clone(
                    Form("hUnfoldeToMc_Y_%i_%i", iter, iCent) + var);
            hUnfoldeToMc_Y[iter].Divide(&hResponseTruth2D_Y);
            hUnfoldeToMc_Y[iter].Scale(1./hUnfolded2DCopy_Y[iter].Integral()*hResponseTruth2D_Y.Integral());
            hUnfoldeToMc_Y[iter].GetYaxis()->SetTitle("Unfolded/Mc true");
            hUnfoldeToMc_Y[iter].GetXaxis()->SetTitle(var);
            hUnfoldeToMc_Y[iter].GetYaxis()->SetTitleOffset(1.5);
            hUnfoldeToMc_Y[iter].GetXaxis()->SetTitleOffset(1.1);
            hUnfoldeToMc_Y[iter].Draw(iter == 0 ? "" : "same");
            hUnfoldeToMc_Y[iter].GetYaxis()->SetRangeUser(0.5,1.5);

            if(iter == 5) UnfoldedToMc_y = ComputeTVShapeDistance(hUnfolded2DCopy_Y[5], hResponseTruth2D_Y);
        }
        DrawLineOne2(hUnfoldeToMc_Y[0].GetXaxis()->GetBinLowEdge(1),
                     hUnfoldeToMc_Y[0].GetXaxis()->GetBinLowEdge(hUnfoldeToMc_Y[0].GetNbinsX() + 1));

        can->cd();
        tex.DrawLatex(0.35, 0.48, "(p_{T}, " + var + ") " + centralityTitles[iCent]);

    }

    fout << runId << "\t"
     << iCent << "\t"
     << "2D" << "\t"
     << var.Data() << "\t"
     << "PT" << "\t"
     << std::fixed << std::setprecision(4)
     << WorstBinMetric_X * 100    << "\t"
     << MeanAbsoluteDrift_X * 100 << "\t"
     << RMSDrift_X * 100          << "\t"
     << WeightedDrift_X * 100 << "\t"
     << UnfoldedToMc_x * 100 
     << "\n";

    fout << runId << "\t"
     << iCent << "\t"
     << "2D" << "\t"
     << var.Data() << "\t"
     << "var" << "\t"
     << std::fixed << std::setprecision(4)
     << WorstBinMetric_Y * 100    << "\t"
     << MeanAbsoluteDrift_Y * 100 << "\t"
     << RMSDrift_Y * 100          << "\t"
     << WeightedDrift_Y * 100 << "\t"
     << UnfoldedToMc_y * 100
     << "\n";

    can->SaveAs(outPdf);

    can->Clear();
    can->SetCanvasSize(1200, 1000);

    if (true) {

        can->Divide(2, 2);

        can->cd(1);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);

        TLatex texA;
        texA.SetNDC();
        texA.SetTextFont(42);
        texA.SetTextSize(0.055);
        texA.SetTextAlign(22);

        Stejn(&hResponseTruth2D, TString("hResponseTruth2D_plot_") + var + Form("_%i", iCent));
        gPad->SetLogz();
        texA.DrawLatex(0.5, 0.5, "Mc true");

        can->cd(3);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);

        Stejn(&hResponseReco2D, TString("hResponseReco2D_plot_") + var + Form("_%i", iCent));
        gPad->SetLogz();
        texA.DrawLatex(0.5, 0.5, "MC reco");

        can->cd(2);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);

        int showIter = GivenIter-1; //3rd it.
        Stejn(&hUnfolded2D[iCent][variable][showIter], TString("hUnfolded2D_") + var + Form("_%i", iCent));
        gPad->SetLogz();
        texA.DrawLatex(0.5, 0.5, Form("Real Unfolded (it. %i)", showIter));

        can->cd(4);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);

        Stejn(&hRealData2D[iCent][variable], TString("hRealData2D_") + var + Form("_%i", iCent));
        gPad->SetLogz();
        texA.DrawLatex(0.5, 0.5, "Real reco");

        can->cd();
        texA.DrawLatex(0.5, 0.49, "(p_{T}, " + var + ") " + centralityTitles[iCent]);

        can->SaveAs(outPdf);
        can->Clear();

        can->Divide(2, 2);

        can->cd(1);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);
        TH2D Response2DRecoTrueVar = *ProjectRecoTrue2D(&rurResponse2D[iCent][variable], "reco", "y", "true", "y");
        Response2DRecoTrueVar.GetXaxis()->SetTitle(var + " ^{,reco}");
        Response2DRecoTrueVar.GetYaxis()->SetTitle(var + " ^{,true}");
        Stejn(&Response2DRecoTrueVar, TString("Response2DRecoTrueVar") + var + Form("_%i", iCent));
        texA.DrawLatex(0.5, 0.5, "RM projection");
        gPad->SetLogz();

        can->cd(3);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);
        TH2D Response2DRecoTruePt = *ProjectRecoTrue2D(&rurResponse2D[iCent][variable], "reco", "x", "true", "x");
        Response2DRecoTruePt.GetXaxis()->SetTitle("p_{T,Jet}^{reco} [GeV/c]");
        Response2DRecoTruePt.GetYaxis()->SetTitle("p_{T,Jet}^{true} [GeV/c]");
        Stejn(&Response2DRecoTruePt, TString("Response2DRecoTruePt") + var + Form("_%i", iCent));
        texA.DrawLatex(0.5, 0.5, "RM projection");
        gPad->SetLogz();

        can->cd(2);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);
        hRespZHighRes[variable + 1][iCent]->Draw("colz");
        hRespZHighRes[variable + 1][iCent]->GetYaxis()->SetTitleOffset(1.3);
        texA.DrawLatex(0.5, 0.5, "RM High Res");
        gPad->SetLogz();

        can->cd(4);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);
        hRespZHighRes[0][iCent]->Draw("colz");
        hRespZHighRes[0][iCent]->GetYaxis()->SetTitleOffset(1.3);
        texA.DrawLatex(0.5, 0.5, "RM High Res");
        gPad->SetLogz();

        can->cd();
        texA.DrawLatex(0.5, 0.49, "(p_{T}, " + var + ") " + centralityTitles[iCent]);

        can->SaveAs(outPdf);

    }

    if(ClosureTest){

        can->cd();
        can->Clear();
        can->SetCanvasSize(800, 400);
        can->Divide(2, 1);






    }



    if(true){
        can->Clear();
        can->SetCanvasSize(1200, 1000);

        can->cd();
        TH2D hResponseClone = *(TH2D *) rurResponse2D[iCent][variable].Hresponse()->Clone(Form("hResponseClone_%i_%s", iCent, var.Data()));
        //set z axis range
        hResponseClone.GetZaxis()->SetRangeUser(0.000001,10e5);
        hResponseClone.Draw("colz");
        hResponseClone.GetXaxis()->SetLabelSize(0);
        hResponseClone.GetYaxis()->SetLabelSize(0);
        hResponseClone.GetXaxis()->SetTickLength(0);
        hResponseClone.GetYaxis()->SetTickLength(0);
        gPad->SetLogz();


        int nybinsZ = hResponseTruth2D_Y.GetXaxis()->GetNbins();
        int nxbinsZ = hResponseReco2D_Y.GetXaxis()->GetNbins();
        int nybinsPt = hResponseTruth2D_X.GetXaxis()->GetNbins();
        int nxbinsPt = hResponseReco2D_X.GetXaxis()->GetNbins();

        TLine line[nybinsZ+1];
        TLine line2[nxbinsZ+1];

        for (int i = 0; i < nybinsZ + 1; i++) {
            line[i] = TLine(0, nybinsPt * i, nxbinsZ * nxbinsPt, nybinsPt * i);
            line[i].SetLineColor(kBlack);
            line[i].SetLineStyle(1);
            line[i].SetLineWidth(1);
            line[i].Draw(i == 0 ? "" : "same");
        }

        for (int i = 0; i < nxbinsZ + 1; i++) {

            line2[i] = TLine(nxbinsPt * i, 0, nxbinsPt * i, nybinsZ * nybinsPt);
            line2[i].SetLineColor(kBlack);
            line2[i].SetLineStyle(1);
            line2[i].SetLineWidth(1);
            line2[i].Draw(i == 0 ? "" : "same");
        }

        //smal axes (X)
        for (int iz = 0; iz < nxbinsZ; ++iz) {
            int firstBinIndex = iz * nxbinsPt + 1;
            int lastBinIndex  = (iz + 1) * nxbinsPt;

            float x_start = hResponseClone.GetXaxis()->GetBinLowEdge(firstBinIndex);
            float x_end   = hResponseClone.GetXaxis()->GetBinLowEdge(lastBinIndex + 1); // +1 for upper edge

            // Fyzikální rozsah datové osy z jednotlivé části (vždy stejný!)
            float pt_low  = hResponseReco2D_X.GetBinLowEdge(1)-1;
            float pt_high = hResponseReco2D_X.GetBinLowEdge(nxbinsPt + 1);

            TGaxis* ptAxis = new TGaxis(x_start, 0, x_end, 0,
                                        pt_low, pt_high,
                                        nxbinsPt, "-N");
             ptAxis->SetLabelSize(0.005);
            ptAxis->SetLabelOffset(-0.012);
            ptAxis->SetLabelFont(42);
            ptAxis->SetTextAngle(90);
            ptAxis->Draw("same");

            for (int i = 1; i <= nxbinsPt + 1; ++i) {
                if (i == 1 || i == nxbinsPt + 1) {
                    ptAxis->ChangeLabel(i, -1, 0); // skryj label
                } else {
                    double val = hResponseReco2D_X.GetBinLowEdge(i);
                    ptAxis->ChangeLabel(i, -1, -1, 12, -1, -1, Form("%.2g", val));
                }
            }

        }

        float x_start =  hResponseClone.GetXaxis()->GetBinLowEdge(1);
        float x_end   = hResponseClone.GetXaxis()->GetBinLowEdge(nxbinsZ * nxbinsPt  + 1);

        TGaxis* zAxis =  new TGaxis(x_start, 0, x_end, 0,
                                      0, 10,
                                    hResponseReco2D_Y.GetXaxis()->GetNbins(), "-N");
        zAxis->SetLabelSize(0.015);
        zAxis->SetLabelOffset(-0.04);
        zAxis->SetLabelFont(42);
        zAxis->Draw("same");

        for (int i = 1; i <= nxbinsZ + 1; ++i) {
                double val = hResponseReco2D_Y.GetBinLowEdge(i);
                zAxis->ChangeLabel(i, -1, -1, 23, -1, -1, Form("%.3g", val));
        }

        zAxis->SetTitle(var+" ^{,reco}");
        zAxis->SetTitleOffset(-1.1);
        zAxis->SetTitleSize(0.03);

        for (int iz = 0; iz < nybinsZ; ++iz) {
            int firstBinIndex = iz * nybinsPt + 1;  // ROOT bin index starts at 1
            int lastBinIndex  = (iz + 1) * nybinsPt;

            float y_start = hResponseClone.GetYaxis()->GetBinLowEdge(firstBinIndex);
            float y_end   = hResponseClone.GetYaxis()->GetBinLowEdge(lastBinIndex + 1);

            float pt_low  = hResponseTruth2D_X.GetBinLowEdge(1) - 1;
            float pt_high = hResponseTruth2D_X.GetBinLowEdge(nybinsPt + 1);

            TGaxis* ptAxisY = new TGaxis(0, y_start, 0, y_end,
                                         pt_low, pt_high,
                                         nybinsPt, "-N");

            ptAxisY->SetLabelSize(0.006);
            ptAxisY->SetLabelOffset(0.008);
            ptAxisY->SetLabelFont(42);
            ptAxisY->SetTextAngle(0);
            ptAxisY->Draw("same");

            for (int i = 1; i <= nybinsPt + 1; ++i) {
                if (i == 1 || i == nybinsPt + 1) {
                    ptAxisY->ChangeLabel(i, -1, 0);  // hide labels at edges
                } else {
                    double val = hResponseTruth2D_X.GetBinLowEdge(i);
                    ptAxisY->ChangeLabel(i, -1, -1, 12, -1, -1, Form("%.2g", val));
                }
            }
        }

        float y_start = hResponseClone.GetYaxis()->GetBinLowEdge(1);
        float y_end   = hResponseClone.GetYaxis()->GetBinLowEdge(nybinsZ * nybinsPt + 1);

        TGaxis* zAxisY = new TGaxis(0, y_start, 0, y_end,
                                    0, 10, hResponseTruth2D_Y.GetXaxis()->GetNbins(), "+N");

        zAxisY->SetLabelSize(0.015);
        zAxisY->SetLabelOffset(-0.05);
        zAxisY->SetLabelFont(42);
        zAxisY->Draw("same");

        for (int i = 1; i <= nybinsZ + 1; ++i) {
            double val = hResponseTruth2D_Y.GetBinLowEdge(i);
            zAxisY->ChangeLabel(i, -1, -1, 12, -1, -1, Form("%.3g", val));
        }

        zAxisY->SetTitle(var + " ^{,true}");
        zAxisY->SetTitleOffset(-1.3);
        zAxisY->SetTitleSize(0.03);

        //přidám přes TLatex popisek k ose x
        TLatex texX;
        texX.SetNDC();
        texX.SetTextFont(42);
        texX.SetTextAngle(90);
        texX.SetTextSize(0.01);
        texX.DrawLatex(0.19, 0.08, "p_{T,Jet}^{reco} [GeV/c]");
        texX.SetTextAngle(0);
        texX.DrawLatex(0.14, 0.135, "p_{T,Jet}^{true} [GeV/c]");

        tex.DrawLatex(0.15, 0.9, "4D Response matrix " + centralityTitles[iCent] + " (p_{T,Jet}, " + var + ")");

        can->SaveAs(outPdf);






    }

    if (true){

        TH2D* hFakes = (TH2D*) rurResponse2D[iCent][variable].Hfakes()->Clone(
            Form("hFakes_%d_%d", iCent, variable));

        TH2D* hTruth = (TH2D*) rurResponse2D[iCent][variable].Htruth()->Clone(
            Form("hTruth_%d_%d", iCent, variable));

        TH2D* hMeasured = (TH2D*) rurResponse2D[iCent][variable].Hmeasured()->Clone(
            Form("hMeasured_%d_%d", iCent, variable));

        TH2D* hResp = (TH2D*) rurResponse2D[iCent][variable].Hresponse();

        // matched truth ve flattenutém prostoru
        TH1D* hMatchedTruthFlat = (TH1D*) hResp->ProjectionY(
            Form("hMatchedTruthFlat_%d_%d", iCent, variable));

        // rozbalení zpět do 2D
        TH2D* hMatchedTruth = (TH2D*) hTruth->Clone(Form("hMatchedTruth_%d_%d", iCent, variable));
        hMatchedTruth->Reset("ICES");

        int nXBins = hTruth->GetNbinsX();
        int nYBins = hTruth->GetNbinsY();

        for (int iY = 1; iY <= nYBins; ++iY) {
            for (int iX = 1; iX <= nXBins; ++iX) {

                int flatBin = iX + nXBins * (iY - 1);

                hMatchedTruth->SetBinContent(iX, iY, hMatchedTruthFlat->GetBinContent(flatBin));
                hMatchedTruth->SetBinError  (iX, iY, hMatchedTruthFlat->GetBinError(flatBin));
            }
        }

        // misses = all truth - matched truth
        TH2D* hMisses = (TH2D*) hTruth->Clone(Form("hMisses_%d_%d", iCent, variable));
        hMisses->Add(hMatchedTruth, -1.0);

        // fake fraction [%]
        TH2D* hFakesPct = (TH2D*) hFakes->Clone(Form("hFakesPct_%d_%d", iCent, variable));
        hFakesPct->SetTitle(Form("Fake fraction, cent %d", iCent));
        hFakesPct->Divide(hMeasured);

        // miss fraction [%]
        TH2D* hMissesPct = (TH2D*) hMisses->Clone(Form("hMissesPct_%d_%d", iCent, variable));
        hMissesPct->SetTitle(Form("Miss fraction, cent %d", iCent));
        hMissesPct->Divide(hTruth);
        can->Clear();
        can->Divide(2,1);

        // Fakes %
        can->cd(1);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);
        gPad->SetTopMargin(0.15);
        gPad->SetLogz(0);
        hFakesPct->SetTitle("Fakes");
        hFakesPct->GetXaxis()->SetTitle("p_{T,Jet}^{reco}");
        hFakesPct->GetYaxis()->SetTitle(var+"^{reco}");
        hFakesPct->SetMinimum(0.0);
        hFakesPct->SetMaximum(1.0);
        hFakesPct->Draw("COLZ TEXT");

        // Misses %
        can->cd(2);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);
        gPad->SetTopMargin(0.15);
        gPad->SetLogz(0);
        hMissesPct->SetTitle("Misses");
        hMissesPct->GetXaxis()->SetTitle("p_{T,Jet}^{true}");
        hMissesPct->GetYaxis()->SetTitle(var+"^{true}");
        hMissesPct->SetMinimum(0.0);
        hMissesPct->SetMaximum(1.0);
        hMissesPct->Draw("COLZ TEXT");

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.05);
        latex.SetTextAlign(22); // center

        // Fakes
        can->cd(1);
        latex.DrawLatex(0.5, 0.92, "Fakes");

        // Misses
        can->cd(2);
        latex.DrawLatex(0.5, 0.92, "Misses");

        can->cd();
        latex.DrawLatex(0.5,0.92,var + " " + centralityTitles[iCent]);

        can->SaveAs(outPdf);
        can->Clear();

        delete hFakes;
        delete hTruth;
        delete hMeasured;
        delete hMatchedTruth;
        delete hMisses;
        delete hFakesPct;
        delete hMissesPct;


    }


}


































void plotComparison(TCanvas *can, RooUnfoldResponse *hResponse, RooUnfoldResponse *hResponseOLD, TH1D *hUnfoldedPt[],
                    TH1D *hBackfoldedPt[], TH1D *hRealData, const Int_t &iCent, TString var, const char *OutputFile,
                    Int_t iSuper, TH1D *Prior, TH1D *hResponseTruthProj, TH2D *Hist1DPTResponse, TH1D *hChiSquared,
                    TH1D *hClosureMCpT) {


    can->Clear();
    can->Divide(2, 1);

    TLegend *leg1 = new TLegend(0.30, 0.63, 0.47, 0.81);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.035);


    TLegend *leg2 = new TLegend(0.70, 0.75, 0.79, 0.83);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    //text size
    leg2->SetTextSize(0.035);

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.050);


    vector <Double_t> ptRecoBinsVec2[nCentralityBins];
    vector <Double_t> ptMcBinsVecCustom2[nCentralityBins];
    for (int j = 0; j < 3; j++) {
        if (UseOverflow) {
            ptRecoBinsVec2[j].push_back(underflowplot);
            ptMcBinsVecCustom2[j].push_back(underflowplot);
        }

        for (int i = 0; i < ptRecoBinsVec[j].size(); i++) ptRecoBinsVec2[j].push_back(ptRecoBinsVec[j][i]);

        for (int i = 0; i < ptMcBinsVecCustom[j].size(); i++) ptMcBinsVecCustom2[j].push_back(ptMcBinsVecCustom[j][i]);

        if (UseOverflow) {
            ptRecoBinsVec2[j].push_back(overflowplot);
            ptMcBinsVecCustom2[j].push_back(overflowplot);
        }
    }

    TH1D *hRealDataProjXPt = (TH1D *) hRealData->Clone(TString("hRealDataProjXPt") + Form("_%i", iCent));
    TH2D *hhResponse = (TH2D *) hResponse->Hresponse();
    TH2D *hhResponseRebinned = new TH2D(Form("hhResponseRebinned_%.d_%i", iCent, iSuper),
                                        "Response Matrix (iter 4);reco;true", ptRecoBinsVec2[iCent].size() - 1,
                                        &ptRecoBinsVec2[iCent][0], ptMcBinsVecCustom2[iCent].size() - 1,
                                        &ptMcBinsVecCustom2[iCent][0]);


    for (int i = 0; i < (ptRecoBinsVec2[iCent].size() - 1); i++) {
        for (int j = 0; j < (ptMcBinsVecCustom2[iCent].size() - 1); j++) {
            hhResponseRebinned->SetBinContent(i + 1, j + 1,
                                              hhResponse->GetBinContent(i + !UseOverflow, j + !UseOverflow));
            hhResponseRebinned->SetBinError(i + 1, j + 1, hhResponse->GetBinError(i + !UseOverflow, j + !UseOverflow));

        }
    }
    TH1D *hhResponseReco = (TH1D *) hhResponseRebinned->ProjectionX(
            TString("hhResponseReco_") + Form("_%i_%i", iCent, iSuper))->Clone(
            TString("Clone_hhResponseReco_") + Form("_%i_%i", iCent, iSuper));
    TH1D *hhResponseTruth = (TH1D *) hhResponseRebinned->ProjectionY(
            TString("hhResponseTruth_") + Form("_%i_%i", iCent, iSuper))->Clone(
            TString("Clone_hhResponseTruth_") + Form("_%i_%i", iCent, iSuper));
    //přepíšu biny z hhResponseTruth do hResponseTruthProj, pokud je histogram prázdný
    //pokud je první bin prázdný
    if (hResponseTruthProj->GetBinContent(1) == 0) {
        for (int i = 1; i <= hResponseTruthProj->GetNbinsX(); i++) {
            hResponseTruthProj->SetBinContent(i, hhResponseTruth->GetBinContent(i));
            hResponseTruthProj->SetBinError(i, hhResponseTruth->GetBinError(i));
        }
    }


    NormalizeByBinWidth(hRealDataProjXPt, 2000 + 1);
    NormalizeByBinWidth(hhResponseReco, 2000 + 2);
    NormalizeByBinWidth(hhResponseTruth, 2000 + 2);


    leg2->AddEntry(hRealDataProjXPt, "Data reco", "lep");
    leg2->AddEntry(hhResponseReco, "RM reco (scaled)", "lep");


    if (ClosureTest) leg1->AddEntry(hhResponseTruth, "RM true", "lep");
    else {leg1->AddEntry(hRealDataProjXPt, "Real data (reco)", "lep");}

    can->cd(1);
    //right margin
    gPad->SetRightMargin(0.05);
    //top
    gPad->SetTopMargin(0.16);
    //bottom
    gPad->SetBottomMargin(0.08);
    //left
    gPad->SetLeftMargin(0.15);

    TPad *Ratio_1 = new TPad("Ratio_1", "Ratio_1", 0.00, 0.5, 1.0, 0.9); //padRecoComp
    Ratio_1->SetBottomMargin(0.0);
    Ratio_1->SetTopMargin(0.0);
    Ratio_1->SetBorderMode(0);
    //left margin
    Ratio_1->SetLeftMargin(0.15);

    TPad *Ratio_2 = new TPad("Ratio_2", "Ratio_2", 0.00, 0.29, 1.0, 0.5); //padRelToNth
    Ratio_2->SetBottomMargin(0.0);
    Ratio_2->SetTopMargin(0.0);
    Ratio_2->SetBorderMode(0);
    Ratio_2->SetLeftMargin(0.15);

    TPad *Ratio_3 = new TPad("Ratio_3", "Ratio_3", 0.00, 0.02, 1.0, 0.29); //padNextStep
    Ratio_3->SetTopMargin(0.0);
    Ratio_3->SetBorderMode(0);
    Ratio_3->SetBottomMargin(0.21);
    Ratio_3->SetLeftMargin(0.15);

/*

    TPad *myPad1 = new TPad("myPad_11", "myPad_11", 0, 0.00, 1.0, 0.90);
    myPad1->SetBottomMargin(0.10);
   // myPad1->SetBorderMode(0);
    //TPad *myPad2 = new TPad("myPad_22", "myPad_22", 0, 0.09, 1.0, 0.4);
    myPad1->SetTopMargin(0.070);
    //myPad2->SetBorderMode(0);
    myPad1->SetLeftMargin(0.15);
    myPad1->Draw();
    //myPad2->Draw();
    myPad1->cd();
*/

    can->cd(2);
    Ratio_1->Draw();
    Ratio_2->Draw();
    Ratio_3->Draw();

    Ratio_1->cd();



    gPad->SetTopMargin(0.15);
    gPad->SetLogy();

    hRealDataProjXPt->GetYaxis()->SetTitle("dN/dp_{T, Jet} [(GeV/c)^{-1}]");//
    hRealDataProjXPt->GetXaxis()->SetTitle("p_{T, Jet} [GeV/c]");
    hRealDataProjXPt->Draw();

    Double_t padHeight = Ratio_2->GetHNDC(); // výška padu v NDC (0-1)

    Double_t desiredLabelSize = 0.015; // cílená velikost pro plnou výšku 1.0
    hRealDataProjXPt->GetYaxis()->SetLabelSize(desiredLabelSize / padHeight);
    hRealDataProjXPt->GetYaxis()->SetTitleSize(desiredLabelSize / padHeight); // Nastavení velikosti názvu osy X
    hRealDataProjXPt->GetYaxis()->SetTitleOffset(0.80);

    hhResponseReco->Scale(hRealDataProjXPt->Integral("width") / hhResponseReco->Integral("width"));
    hhResponseReco->Draw("same");

    ////Vypíšu hodnoty binů obou histogramů:
    for (int i = 1; i <= hhResponseReco->GetNbinsX(); i++) {
        cout << "1DBin " << i << ": Real Data: " << hRealDataProjXPt->GetBinContent(i)
             << " RM reco: " << hhResponseReco->GetBinContent(i) << endl;
    }


    can->cd(1);
    gPad->SetLogy();
    if (ClosureTest) {
        hhResponseTruth->GetYaxis()->SetTitle("dN/dp_{T, Jet} [(GeV/c)^{-1}]");
        hhResponseTruth->GetXaxis()->SetTitle("p_{T, Jet} [GeV/c]");
        hhResponseTruth->GetYaxis()->SetTitleOffset(1.5);
        hhResponseTruth->GetYaxis()->SetTitleSize(0.045);
        hhResponseTruth->GetXaxis()->SetTitleOffset(0.7);
        hhResponseTruth->GetXaxis()->SetTitleSize(0.045);
        hhResponseTruth->SetMarkerStyle(21);
        hhResponseTruth->Draw();
        //offset čísel na ose x
        hhResponseTruth->GetXaxis()->SetLabelOffset(-0.015);
        hhResponseTruth->GetXaxis()->SetRange(1 + UseOverflow, hhResponseTruth->GetNbinsX() - UseOverflow);
        //hhResponseTruth->Scale(hRealDataProjXPt->Integral("width") / hhResponseTruth->Integral("width"));/////TEST
    } else
        hRealDataProjXPt->Draw();

    //myPad2->cd();
    //gPad->SetLogy();

    can->cd();
    leg1->Draw("same");
    leg2->Draw("same");
    //tex->DrawLatex(0.2, 0.68, "p_{T}");
    tex->SetTextAngle(90);
    //center text
    tex->DrawLatex(0.978, 0.4, "1D " + centralityTitles[iCent] + "  " + var);


    TH1D *hUnfoldedProjXPt[nIter];
    RooUnfoldBayes *unfoldingPt[nIter];
    TH1D *hBayRatios[nIter];
    TH1D *hBayRatios3[nIter];

    TH1D *hBackFoldBayRatios[nIter];
    TH2D *fPearsonCoeffs[nIter];

    PrintCheckNumbers[iCent][0][0] = hRealData->Integral();
    for (int iter = 0; iter < nIter; iter++) {


        unfoldingPt[iter] = new RooUnfoldBayes(hResponse, hRealData, PlotIterations[iter]);

        hUnfoldedPt[iter] = (TH1D *) unfoldingPt[iter]->Hreco();
        ////hUnfoldedPt[iter] = (TH1D *) unfoldingPt[iter]->Hunfold();


        hBackfoldedPt[iter] = (TH1D *) hResponse->ApplyToTruth(hUnfoldedPt[iter]);
        NormalizeByBinWidth(hBackfoldedPt[iter], 2000 + iter);

        //-------------------------------------
        PrintCheckNumbers[iCent][0][iter + 1] = hUnfoldedPt[iter]->Integral();

        //-------------------------------------

        //nastavím jméno hUnfoldedPt
        hUnfoldedPt[iter]->SetName(TString("hUnfoldedPt") + Form("_%i_%i", iter, iCent));

        fPearsonCoeffs[iter] = getPearsonCoeffs1D(unfoldingPt[iter]->Ereco(RooUnfold::kCovariance));
        /////fPearsonCoeffs[iter] = getPearsonCoeffs1D(unfoldingPt[iter]->Eunfold(RooUnfold::kCovariance));

        //Nastavím jiné jméno
        fPearsonCoeffs[iter]->SetName(Form("PearsonCoeffsIter%i%d%i", PlotIterations[iter], iCent, iSuper));
        //copy of TH1D
        hUnfoldedProjXPt[iter] = (TH1D *) hUnfoldedPt[iter]->Clone(
                TString("hUnfoldedProjXPt") + Form("_%i_%i", iter, iCent));

        NormalizeByBinWidth(hUnfoldedProjXPt[iter], 2005 + iter);
        hUnfoldedProjXPt[iter]->SetMarkerStyle(27);

        leg1->AddEntry(hUnfoldedProjXPt[iter], Form("%i iter", PlotIterations[iter]), "lep");


        // first draw and compare x projections of hUnfolded and hMc
        can->cd(1);
        gPad->SetLogy();
        hUnfoldedProjXPt[iter]->GetYaxis()->SetTitle("dN/dp_{T}");

        if (iter == 4) {
            hUnfoldedProjXPt[iter]->SetMarkerStyle(20);
            // hUnfoldedProjX[iter]->SetMarkerSize(1);
        }
        hUnfoldedProjXPt[iter]->Draw("same");
        //set title size



       //// myPad2->cd();
        Ratio_3->cd();
        gPad->SetLogy(0);
        hBayRatios[iter] = (TH1D *) hUnfoldedProjXPt[iter]->Clone(TString("hBayRatios") + Form("_%i_%i", iter, iCent));
        hBayRatios[iter]->GetYaxis()->SetTitle("Unfolded i-th/Unfolded 1st");
        //offset
        hBayRatios[iter]->GetYaxis()->SetTitleOffset(0.80);
        hBackFoldBayRatios[iter] = (TH1D *) hBackfoldedPt[iter]->Clone(
                TString("hBackFoldBayRatios") + Form("_%i_%i", iter, iCent));
        hBackFoldBayRatios[iter]->Divide(hRealDataProjXPt);
        hBackFoldBayRatios[iter]->GetYaxis()->SetTitle("Backfolded/Real Data");

        hBayRatios3[iter] = (TH1D *) hUnfoldedProjXPt[iter]->Clone(TString("hBayRatios3") + Form("_%i_%i", iter, iCent));
        hBayRatios3[iter]->GetYaxis()->SetTitle("Unfolded i-th/Unfolded 1st");
        //hBayRatios3[iter]->Divide(hUnfoldedProjXPt[0]);

    }

    //Nastavím rozsah hRealDataProjXPt, který je vykreslen spolu s hUnfoldedProjXPt[iter], tak aby byl vidět celý graf

    if (!ClosureTest){
    double minAA = min(findMin(hUnfoldedProjXPt, 2000, nIter), hRealDataProjXPt->GetMinimum());
    double maxAA = max(findMax(hUnfoldedProjXPt, 0, nIter), hRealDataProjXPt->GetMaximum());
    hRealDataProjXPt->GetYaxis()->SetRangeUser(minAA>0?0.1*minAA:0.001, 10 * maxAA);
    }
    for (int iter = (UseRelativeP ? 1 : 0); iter < nIter; iter++) {
        hBayRatios[iter]->Draw(iter == (UseRelativeP ? 1 : 0) ? "" : "same");
        hBayRatios[(UseRelativeP ? 1 : 0)]->GetYaxis()->SetRangeUser(0.74, 1.26);
        if (iter == 0) continue;
        if (UseRelativeP) hBayRatios[iter]->Divide(hUnfoldedProjXPt[iter - 1]);
        else hBayRatios[iter]->Divide(hUnfoldedProjXPt[0]);

        //set range
        if (UseRelativeP) hBayRatios[iter]->GetYaxis()->SetTitle("i-th/(i-1)-th iteration");
        else hBayRatios[iter]->GetYaxis()->SetTitle("Unfolded/(1st iter)");
        hBayRatios[iter]->GetXaxis()->SetLabelSize(0.05); // velikost čísel na ose X
        hBayRatios[iter]->GetYaxis()->SetLabelSize(0.05); // velikost čísel na ose Y
        hBayRatios[iter]->GetXaxis()->SetTitleSize(0.06); // velikost titulku osy X
        hBayRatios[iter]->GetYaxis()->SetTitleSize(0.06); // velikost titulku osy Y
        //offset

    }


    //    for (int iter = (UseRelativeP?1:0); iter < nIter; iter++) {
    //po zpátku
    for (int iter = nIter - 1; iter >= (UseRelativeP ? 1 : 0); iter--) {
        Double_t padHeight = Ratio_3->GetHNDC(); // výška padu v NDC (0-1)

        Double_t desiredLabelSize = 0.025; // cílená velikost pro plnou výšku 1.0
        hBayRatios[iter]->GetYaxis()->SetLabelSize(desiredLabelSize / padHeight);
        hBayRatios[iter]->GetXaxis()->SetLabelSize(desiredLabelSize / padHeight);
        hBayRatios[iter]->GetXaxis()->SetTitleSize(desiredLabelSize / padHeight);
        hBayRatios[iter]->GetYaxis()->SetTitleSize(desiredLabelSize / padHeight);
        hBayRatios[iter]->GetYaxis()->SetTitleOffset(0.40);
        hBayRatios[iter]->GetYaxis()->CenterTitle(); // Center the title


        hBayRatios[iter]->Draw("same");
    }


    hBayRatios[0]->Divide(hBayRatios[0]);


    //hBayRatios[0]->GetYaxis()->SetRangeUser(MinA-0.1*MinA,MaxA+0.1*MaxA);

    DrawLineOne();
    //hMcProjXPt->Scale(hUnfoldedProjXPt[3]->GetMaximum() / hMcProjXPt->GetMaximum());
    //hRealMeasured->Scale(hMcProjX->GetMaximum() / hRealMeasured->GetMaximum());


    Ratio_2->cd();

    for (int iter = 0; iter < nIter; iter++) {


        hBayRatios3[iter]->Draw(iter == 0 ? "" : "same");
        hBayRatios3[iter]->Divide(hUnfoldedProjXPt[0]);
        hBayRatios3[iter]->GetYaxis()->SetTitle("Unfolded/(1st iter)");
        DrawLineOne();
        double MinA = findMin(hBayRatios3, 2000, nIter);
        double MaxA = findMax(hBayRatios3, 0, nIter);
        Double_t padHeight = Ratio_2->GetHNDC(); // výška padu v NDC (0-1)

        Double_t desiredLabelSize = 0.025; // cílená velikost pro plnou výšku 1.0
        hBayRatios3[iter]->GetYaxis()->SetLabelSize(desiredLabelSize / padHeight);
        hBayRatios3[iter]->GetXaxis()->SetLabelSize(desiredLabelSize / padHeight);
        hBayRatios3[iter]->GetXaxis()->SetTitleSize(desiredLabelSize / padHeight*1.1);
        hBayRatios3[iter]->GetYaxis()->SetTitleSize(desiredLabelSize / padHeight*1.1);
        hBayRatios3[iter]->GetYaxis()->SetTitleOffset(0.45);
        hBayRatios3[iter]->GetYaxis()->CenterTitle(); // Center the title
        hBayRatios3[iter]->GetXaxis()->SetTitleOffset(0.85);

        //y raaange
        //hBayRatios3[0]->GetYaxis()->SetRangeUser(MinA - 0.1 * MinA, MaxA + 0.1 * MaxA);
        hBayRatios3[0]->GetYaxis()->SetRangeUser(0.74, 1.26);


    }



    Double_t Xmax = findMax(hUnfoldedProjXPt, hhResponseTruth->GetMaximum(), nIter);
    Double_t Xmin = findMin(hUnfoldedProjXPt, hhResponseTruth->GetMinimum(), nIter);

    hhResponseTruth->GetYaxis()->SetRangeUser(Xmin > 0 ? Xmin * 0.5 : 0.001, Xmax * 2);

    can->cd();
    //Nakreslím vertiáklní čáru
    TLine *line = new TLine(16, -8, 16, 12);
    line->SetLineColor(kBlack);
    line->SetLineColor(kBlack);
    line->SetLineStyle(1);
    line->SetLineWidth(2);
    //line->Draw();
    can->cd();
    DrawTextAbove(iSuper);
    can->SaveAs(outPdf);
    can->Clear();
    //-------------------------------------
    if (ResponseMatrix) {

        can->Divide(2, 2);

        //RooUnfoldBayes* lastUnfolding = new RooUnfoldBayes(hResponse, hRealData,1);

        TH2D *hUnfoldingMatrix = new TH2D(unfoldingPt[2]->UnfoldingMatrix());

        hUnfoldingMatrix->SetTitleOffset(0.5, "Y");
        hUnfoldingMatrix->SetName(Form("UnfoldingMatrix_%i", iCent));
        hUnfoldingMatrix->GetXaxis()->SetTitle("N^{(" + var + ", reco)}_{bin}");
        hUnfoldingMatrix->GetYaxis()->SetTitle("N^{(" + var + ", unfold)}_{bin}");
        can->cd(1);

        //Pravý okraj
        gPad->SetLogz(1);
        gPad->SetRightMargin(0.15);
        gPad->SetTopMargin(0.29);
        hUnfoldingMatrix->GetZaxis()->SetTitleOffset(0.7);

        hUnfoldingMatrix->Draw("colz");
        gPad->Update();
        //TPaletteAxis *palette = (TPaletteAxis *)hUnfoldingMatrix->GetListOfFunctions()->FindObject("palette");
        //palette->SetX1NDC(0.90);
        //palette->SetX2NDC(0.92);

        tex->DrawLatex(0.3, 0.35, centralityTitles[iCent] + "  " + var);
        //tex->DrawLatex(0.5, 0.25, "Unfolding Matrix");
        tex->DrawLatex(0.25, 0.25, "Bin Migration probability (iter 3)");
        can->cd(2);
        gPad->SetRightMargin(0.15);
        gPad->SetTopMargin(0.29);
        gPad->SetLogz(1);

        //Vykreslím response jako clon a převod na th2D
        TH2D *hhResponse = (TH2D *) hResponse->Hresponse();
        hhResponse->Sumw2();

//Vytvořím vektor vector<Double_t> ptRecoBinsVec[nCentralityBins+2]

        vector <Double_t> ptRecoBinsVec2[nCentralityBins];
        vector <Double_t> ptMcBinsVecCustom2[nCentralityBins];
        for (int j = 0; j < 3; j++) {
            if (UseOverflow) {
                ptRecoBinsVec2[j].push_back(underflowplot);
                ptMcBinsVecCustom2[j].push_back(underflowplot);
            }

            for (int i = 0; i < ptRecoBinsVec[j].size(); i++) ptRecoBinsVec2[j].push_back(ptRecoBinsVec[j][i]);

            for (int i = 0; i < ptMcBinsVecCustom[j].size(); i++)
                ptMcBinsVecCustom2[j].push_back(ptMcBinsVecCustom[j][i]);

            if (UseOverflow) {
                ptRecoBinsVec2[j].push_back(overflowplot);
                ptMcBinsVecCustom2[j].push_back(overflowplot);
            }
        }

        //Binování z ptRecoBinsVec[nCentralityBins] a ptMcBinsVecCustom[nCentralityBins]
        TH2D *hResponseRebinned = new TH2D(Form("hResponseRebinned_%.d_%i", iCent, iSuper),
                                           "Response Matrix (iter 4);reco;true", ptRecoBinsVec2[iCent].size() - 1,
                                           &ptRecoBinsVec2[iCent][0], ptMcBinsVecCustom2[iCent].size() - 1,
                                           &ptMcBinsVecCustom2[iCent][0]);
        //přepíšu biny
        for (int i = 0; i < (ptRecoBinsVec2[iCent].size() - 1); i++) {
            for (int j = 0; j < (ptMcBinsVecCustom2[iCent].size() - 1); j++) {
                hResponseRebinned->SetBinContent(i + 1, j + 1,
                                                 hhResponse->GetBinContent(i + !UseOverflow, j + !UseOverflow));
                //Nastavím chyby
                hResponseRebinned->SetBinError(i + 1, j + 1,
                                               hhResponse->GetBinError(i + !UseOverflow, j + !UseOverflow));
            }
        }


        hResponseRebinned->SetTitle("");
        hResponseRebinned->GetXaxis()->SetTitle("p_{T}^{reco} [GeV/c]");
        hResponseRebinned->GetYaxis()->SetTitle("p_{T}^{true} [GeV/c]");
        hResponseRebinned->Draw("colz");


        //Normalize hResponseNwRebinned to one
        //hResponseNwRebinned->Scale(1/hResponseNwRebinned->Integral());
        if (UseOverflow) DrawLineAround(iCent);
        tex->DrawLatex(0.3, 0.35, centralityTitles[iCent] + "  " + var);
        //tex->DrawLatex(0.5, 0.25, "Unfolding Matrix");
        tex->DrawLatex(0.25, 0.25, "Response matrix");
        can->cd(3);
        gPad->SetLogy(1);
        //Vytvořím x-ovou projekci hResponse
        TH1D *hResponseX_ = (TH1D *) hResponseRebinned->ProjectionX("hResponseX_")->Clone("Clone_hResponseX_");
        //Nastavím jinou barvu
        //hResponseX_->SetTitle("Response Matrix (iter 4) X projection");
        tex->DrawLatex(0.3, 0.35, centralityTitles[iCent] + "  " + var);
        //zapnu chyby:
        hResponseX_->Sumw2();
        //tex->DrawLatex(0.25, 0.25, "Response matrix (iter 4)");
        //normalizuji
        NormalizeByBinWidth(hResponseX_, 3);
        hResponseX_->Draw("E");
        //Nastavím popisek y-ovou osu
        hResponseX_->GetYaxis()->SetTitle("dN/dp_{T}");
        //dám legendu
        TLegend *legend = new TLegend(0.6, 0.6, 0.8, 0.8);
        legend->AddEntry(hResponseX_, "reco", "l");
        //žádný rámeček
        legend->SetBorderSize(0);
        //žádné pozadí
        legend->SetFillStyle(0);
        legend->SetTextSize(0.08);
        legend->Draw();

        can->cd(4);
        gPad->SetLogy(1);
        //Vytvořím y-ovou projekci hResponse
        TH1D *hResponseY_ = (TH1D *) hResponseRebinned->ProjectionY("hResponseY_")->Clone("Clone_hResponseY_");
        //Nastavím jinou barvu
        //hResponseY_->SetTitle("Response Matrix (iter 4) Y projection");
        NormalizeByBinWidth(hResponseY_, 4);
        hResponseY_->Sumw2();
        hResponseY_->Draw("E");
        //Přepíšu biny do Prior
        if (Prior) {
            for (int i = 0; i < hResponseY_->GetNbinsX(); i++) {
                Prior->SetBinContent(i + 1, hResponseY_->GetBinContent(i + 1 + UseOverflow));
                Prior->SetBinError(i + 1, hResponseY_->GetBinError(i + 1 + UseOverflow));
            }
        }
        //Nastavím popisek y-ovou osu
        hResponseY_->GetYaxis()->SetTitle("dN/dp_{T}");
        //dám legendu
        TLegend *legend2 = new TLegend(0.6, 0.6, 0.8, 0.8);
        //Nastavím velikost
        legend2->SetTextSize(0.08);
        legend2->AddEntry(hResponseY_, "true", "l");
        //žádný rámeček
        legend2->SetBorderSize(0);
        //žádné pozadí
        legend2->SetFillStyle(0);
        legend2->Draw();

        can->cd();

        DrawTextAbove(iSuper);
        can->SaveAs(outPdf);
        //delete hUnfoldingMatrixB;
        delete hUnfoldingMatrix;
        //delete hEfficiency;
        //delete hPurity;
        //delete ResolutionMatrix;
        can->Clear();
    }
    if (PearsonCoeff) {
        can->cd();
        //topmargin
        gPad->SetTopMargin(0.15);
        //right margin
        gPad->SetRightMargin(0.2);
        gPad->SetLeftMargin(0.2);
        can->SetLogz(0);
        TH2D *hPearsRebinned[nIter];

        //vector<Double_t> ptRecoBinsVec2[nCentralityBins];
        vector <Double_t> ptMcBinsVecCustom2[nCentralityBins];

        for (int j = 0; j < 3; j++) {
            if (UseOverflow) {
                ptMcBinsVecCustom2[j].push_back(underflowplot);
            }

            for (int i = 0; i < ptMcBinsVecCustom[j].size(); i++)
                ptMcBinsVecCustom2[j].push_back(ptMcBinsVecCustom[j][i]);
            if (UseOverflow) {
                ptMcBinsVecCustom2[j].push_back(overflowplot);
            }
        }
        //Rozdělím canvas na nIter částí
        can->SetCanvasSize(3 * 400, (TMath::Ceil(nIter / 3.) + 1) * 400);

        can->Divide(3, TMath::Ceil(nIter / 3.) + 1);

        for (int iter = 0; iter < nIter; iter++) {
            can->cd(iter + 1 + 3);
            gPad->SetRightMargin(0.15);
            gPad->SetLeftMargin(0.15);
            hPearsRebinned[iter] = new TH2D(Form("hPearsRebinned_%.d_%.d", iCent, iter),
                                            "Response Matrix (iter 4);reco;true", ptMcBinsVecCustom2[iCent].size() - 1,
                                            &ptMcBinsVecCustom2[iCent][0], ptMcBinsVecCustom2[iCent].size() - 1,
                                            &ptMcBinsVecCustom2[iCent][0]);
            //přendám biny
            for (int i = 0; i < (ptMcBinsVecCustom2[iCent].size() - 1); i++) {
                for (int j = 0; j < (ptMcBinsVecCustom2[iCent].size() - 1); j++) {
                    hPearsRebinned[iter]->SetBinContent(i + 1, j + 1,
                                                        fPearsonCoeffs[iter]->GetBinContent(i + 1, j + 1));
                }
            }
            //Nastavím na x-ové ose, aby se zobrazovaly čísla po 5
            //hPearsRebinned[iter]->GetXaxis()->SetNdivisions(-605);
            hPearsRebinned[iter]->SetTitle((TString) ";p_{T};p_{T}");
            hPearsRebinned[iter]->Draw("colz");
            //Nastavím range z-ové osy na -1 až 1
            hPearsRebinned[iter]->GetZaxis()->SetRangeUser(-1, 1);
            gPad->Update();
            hPearsRebinned[iter]->GetZaxis()->SetTitleOffset(0.7);
            tex->DrawLatex(0.25, 0.65, centralityTitles[iCent] + "  " + var);
            tex->DrawLatex(0.25, 0.35, Form("Pearson Coefficients Iter %i", PlotIterations[iter]));
            if (UseOverflow) DrawLineAround2(iCent);
            // TPaletteAxis *palette = (TPaletteAxis *)fPearsonCoeffs[iter]->GetListOfFunctions()->FindObject("palette");
            //palette->SetX1NDC(0.91);
            // palette->SetX2NDC(0.93);

        }
        can->cd();
        DrawTextAbove(iSuper);
        can->SaveAs(outPdf);
        can->SetCanvasSize(1200, 1000);
        for (int iter = 0; iter < nIter; iter++) {
            ////delete fPearsonCoeffs[iter];
            delete hPearsRebinned[iter];
        }

        can->Clear();
        can->Divide(2, 1);
        can->cd(1);
        /************************************************************/
        gPad->SetTopMargin(0.15);

        TPad *myPad1BF = new TPad("myPad_11BF", "myPad_11BF", 0, 0.4, 1.0, 0.92);
        myPad1BF->SetBottomMargin(0.0);
        myPad1BF->SetBorderMode(0);
        TPad *myPad2BF = new TPad("myPad_22BF", "myPad_22BF", 0, 0.09, 1.0, 0.4);
        myPad2BF->SetTopMargin(0.0);
        myPad2BF->SetBorderMode(0);
        myPad1BF->Draw();
        myPad2BF->Draw();
        myPad1BF->cd();

        //there will be backfolded
        hRealDataProjXPt->GetYaxis()->SetTitle("dN/dp_{T}");
        hRealDataProjXPt->GetXaxis()->SetTitle("p_{T, Jet} [GeV/c]");
        //set logy
        gPad->SetLogy();
        //Draw all backfolded
        hRealDataProjXPt->Draw();
        hRealDataProjXPt->SetMarkerStyle(20);
        hRealDataProjXPt->SetMarkerSize(1);
        hRealDataProjXPt->SetMarkerColor(kBlack);
        hRealDataProjXPt->SetLineColor(kBlack);
        //new
        TLegend *legBF = new TLegend(0.65, 0.55, 1.0, 0.82);
        //no boder
        legBF->SetBorderSize(0);
        //no fill
        legBF->SetFillStyle(0);

        legBF->AddEntry(hRealDataProjXPt, "Real Data", "lep");


        for (int iter = 0; iter < nIter; iter++) {
            hBackfoldedPt[iter]->Draw("same");
            legBF->AddEntry(hBackfoldedPt[iter], Form("Iter%i", PlotIterations[iter]), "lep");
        }
        legBF->Draw("same");

        myPad2BF->cd();
        //hBackFoldBayRatios all

        for (int iter = 0; iter < nIter; iter++) {
            hBackFoldBayRatios[iter]->Draw(iter == 0 ? "" : "same");
            hBackFoldBayRatios[iter]->GetYaxis()->SetRangeUser(0.5, 1.5);
        }

        //draw one line
        DrawLineOne2(ptRecoBinsVec[iCent][0], ptRecoBinsVec[iCent][ptRecoBinsVec[iCent].size() - 1]);

        can->cd(2);
        //left margin
        gPad->SetLeftMargin(0.25);
        TH1D *chisquaredndf = new TH1D("chisquaredndf", "chisquaredndf", nIter, 0, nIter);
        for (int iter = 0; iter < nIter; iter++) {
            //chi2/ndf = 1/N * sum((hRealDataProjXPt->GetBinContent(i) - hBackfoldedPt[iter]->GetBinContent(i))^2 / hRealDataProjXPt->GetBinError(i)^2)
            double chi2ndf = 0;
            for (int i = 1; i <= hRealDataProjXPt->GetNbinsX(); i++) {
                chi2ndf += (hRealDataProjXPt->GetBinContent(i) - hBackfoldedPt[iter]->GetBinContent(i)) *
                           (hRealDataProjXPt->GetBinContent(i) - hBackfoldedPt[iter]->GetBinContent(i)) /
                           (hRealDataProjXPt->GetBinError(i) * hRealDataProjXPt->GetBinError(i));
            }
            chi2ndf /= hRealDataProjXPt->GetNbinsX();

            chisquaredndf->SetBinContent(iter + 1, chi2ndf);

            //change bin name to iter
            chisquaredndf->GetXaxis()->SetBinLabel(iter + 1, Form("Iter %i", PlotIterations[iter]));
        }
        chisquaredndf->Draw("ph");
        //set y label
        chisquaredndf->GetYaxis()->SetTitle(
                "#chi^{2}/ndf = #frac{1}{N} #sum_{i}^{N} #frac{(Data_{i} - Backfolded_{i})^{2}}{Error_{i}^{2}}");
        //offset
        chisquaredndf->GetYaxis()->SetTitleOffset(1.6);







        /*********************************************************/
        //  gPad->SetTopMargin(0.15);
        //gPad->SetRightMargin(0.15);
        //  gPad->SetLeftMargin(0.15);



        can->SaveAs(outPdf);

    }




    //---------------------------------
    delete leg1;
    delete leg2;
    delete tex;

    delete hRealDataProjXPt;


    if (ClosureTest) {
        //kopie hClosureMCpT
        TH1D *hClosureMCpT2 = (TH1D *) hClosureMCpT->Clone(TString("_hClosureMCpT2") + Form("_%i_%i", iCent, iSuper));
        TH1D *hUnfoldedPt2[nIter];
        for (int iter = 0; iter < nIter; iter++) {
            hUnfoldedPt2[iter] = (TH1D *) hUnfoldedPt[iter]->Clone(
                    TString("_hUnfoldedPt2") + Form("_%i_%i", iter, iCent));
        }

        can->cd();
        can->Clear();
        gPad->SetTopMargin(0.15);
        can->Divide(2, 1);
        can->cd(1);
        gPad->SetTopMargin(0.15);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.15);
        hClosureMCpT2->GetYaxis()->SetTitle("dN/dp_{T}");
        hClosureMCpT2->GetXaxis()->SetTitle("p_{T, Jet} [GeV/c]");
        //set logy
        gPad->SetLogy();
        NormalizeByBinWidth(hClosureMCpT2, 1);

        hClosureMCpT2->Draw();
        //legenda
        TLegend *leg = new TLegend(0.55, 0.55, 0.95, 0.82);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(hClosureMCpT2, "MC (scaled)", "lep");


        for (int iter = 0; iter < nIter; iter++) {
            NormalizeByBinWidth(hUnfoldedPt2[iter], 2000 + iter);
            //scale by maximum
            hUnfoldedPt2[iter]->Scale(hClosureMCpT2->Integral() / hUnfoldedPt2[iter]->Integral());

            hUnfoldedPt2[iter]->Draw("same");
            leg->AddEntry(hUnfoldedPt2[iter], Form("Iter%i", PlotIterations[iter]), "lep");
        }

        leg->Draw("same");

        can->cd(2);
        gPad->SetTopMargin(0.15);
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.15);

        //legenda2
        TLegend *leg2 = new TLegend(0.55, 0.55, 0.95, 0.82);
        leg2->SetBorderSize(0);
        leg2->SetFillStyle(0);

        TH1D *hRatiosMc[nIter];
        for (int iter = 0; iter < nIter; iter++) {
            hRatiosMc[iter] = (TH1D *) hUnfoldedPt2[iter]->Clone(
                    TString("hRatiosMc") + Form("_%i_%i_%i", iter, iCent, iSuper));
            hRatiosMc[iter]->Divide(hClosureMCpT2);
            hRatiosMc[iter]->GetYaxis()->SetTitle("Unfolded/MC");
            hRatiosMc[iter]->GetYaxis()->SetRangeUser(0.5, 1.5);
            hRatiosMc[iter]->Draw(iter == 0 ? "" : "same");
            leg2->AddEntry(hRatiosMc[iter], Form("Iter%i", PlotIterations[iter]), "lep");
        }
        leg2->Draw("same");
        DrawLineOne();
        can->cd();
        DrawTextAbove(iSuper);
        tex->DrawLatex(0.5, 0.05, centralityTitles[iCent] + "  " + var);

        //save
        can->SaveAs(outPdf);

        //---------------------------------
        can->Clear();
        can->cd();
        can->Divide(2, 1);

        TLegend *leg3 = new TLegend(0.55, 0.55, 0.95, 0.82);
        leg3->SetBorderSize(0);
        leg3->SetFillStyle(0);

        can->cd(1);
        //upper margin
        gPad->SetTopMargin(0.15);
        gPad->SetLeftMargin(0.15);

        //vypnu logscale
        gPad->SetLogy(0);
        TH1D *hRelativeRatios[nIter];
        for (int iter = 1; iter < nIter; iter++) {
            hRelativeRatios[iter] = (TH1D *) hUnfoldedPt2[iter]->Clone(
                    TString("hRelativeRatios") + Form("_%i_%i_%i", iter, iCent, iSuper));
            hRelativeRatios[iter]->Divide(hUnfoldedPt2[iter - 1]);
            hRelativeRatios[iter]->GetYaxis()->SetTitle("Unfolded n-th/Unfolded (n-1)-th");
            hRelativeRatios[iter]->GetYaxis()->SetRangeUser(0.5, 1.5);
            hRelativeRatios[iter]->GetYaxis()->SetTitleOffset(1.1);
            hRelativeRatios[iter]->Draw(iter == 1 ? "" : "same");
            leg3->AddEntry(hRelativeRatios[iter], Form("Iter%i", PlotIterations[iter]), "lep");
        }
        leg3->Draw("same");
        DrawLineOne();

        can->cd(2);
        //upper margin
        gPad->SetTopMargin(0.15);
        gPad->SetLeftMargin(0.15);

        TH1D *hRelativeUncertRatiosVs3[nIter];
        TH1D *hRelativeUncertRatios;
        TH1D *hRelativeUncertRatiosA[nIter];

        TLegend *leg4 = new TLegend(0.55, 0.55, 0.95, 0.82);
        leg4->SetBorderSize(0);
        leg4->SetFillStyle(0);

        for (int iter = 0; iter < nIter; iter++) {
            hRelativeUncertRatiosA[iter] = (TH1D *) hUnfoldedPt2[iter]->Clone(
                    TString("hRelativeUncertRatiosVs3") + Form("_%i_%i_%i", iter, iCent, iSuper));
            //Vynuluji všechny hodnoty a chyby;


            for (int i = 0; i < hRelativeUncertRatiosA[iter]->GetNbinsX(); i++) {
                hRelativeUncertRatiosA[iter]->SetBinContent(i + 1, 0);
                hRelativeUncertRatiosA[iter]->SetBinError(i + 1, 0.00001);
                double value = hUnfoldedPt2[iter]->GetBinError(i + 1) / hUnfoldedPt2[2]->GetBinError(i + 1);
                hRelativeUncertRatiosA[iter]->SetBinContent(i + 1, value);
            }


            hRelativeUncertRatiosA[iter]->GetYaxis()->SetTitle("Unfolded n-th/Unfolded 3-rd (uncertainty)");
            //left offset
            hRelativeUncertRatiosA[iter]->GetYaxis()->SetTitleOffset(1.1);
            hRelativeUncertRatiosA[iter]->GetYaxis()->SetRangeUser(0, 3);
            hRelativeUncertRatiosA[iter]->Draw(iter == 0 ? "" : "same");
            leg4->AddEntry(hRelativeUncertRatiosA[iter], Form("Iter%i", PlotIterations[iter]), "lep");
        }

        leg4->Draw("same");


        can->cd();
        DrawTextAbove(iSuper);
        tex->DrawLatex(0.5, 0.05, centralityTitles[iCent] + "  " + var);

        can->SaveAs(outPdf);

        //---------------------------------
        DeleteArray(hRatiosMc, nIter);
        DeleteArray(hUnfoldedPt2, nIter);
        //delete hClosureMCpT2;
        delete hClosureMCpT2;


    }

    DeleteArray(hUnfoldedProjXPt, nIter);
    DeleteArray(unfoldingPt, nIter);

}
void plotRcp1D(TCanvas *can){

    can->Clear();
    can->Divide(3, 1);
    can->SetCanvasSize(1200, 400);


    TLegend leg1 = TLegend(0.63, 0.55, 0.80, 0.82);
    leg1.SetBorderSize(0);
    leg1.SetFillStyle(0);

    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    tex.SetTextSize(0.055);
    tex.SetTextAlign(22);

    Int_t numOfComb = factorial(nCentralityBins) / (factorial(nCentralityBins - 2) * factorial(2));
    TH1D hRcpPt[numOfComb][nIter];

    TLegend *leg[numOfComb];
    for (Int_t iComb = 0; iComb < numOfComb; iComb++) {
        leg[iComb] = new TLegend(0.13, 0.60, 0.30, 0.87);
        leg[iComb]->SetBorderSize(0);
        leg[iComb]->SetFillStyle(0);
        leg[iComb]->SetTextSize(0.055);
        leg[iComb]->SetTextFont(42);
    }
    //takovou smyčku, aby vždy centralita v čitateli byla menší než ve jmenovateli
    for (int iCent = 0; iCent < nCentralityBins; iCent++) {
        for (int jCent = iCent + 1; jCent < nCentralityBins; jCent++) {
            for (int iter = 0; iter < nIter; iter++) {
                Int_t nthParam = GetCombIndex(iCent,jCent,nCentralityBins);
                hRcpPt[nthParam][iter] = *(TH1D*)hUnfoldedPt[iCent][iter].Clone(Form("hRcpPt_%i_%i_%i", iCent, jCent, iter));
                hRcpPt[nthParam][iter].Divide(&hUnfoldedPt[jCent][iter]);
                hRcpPt[nthParam][iter].Scale(1.*Ncoll[jCent] / Ncoll[iCent]);
                hRcpPt[nthParam][iter].Scale(1.*NumberOfWEvents[jCent] / NumberOfWEvents[iCent]);
                hRcpPt[nthParam][iter].GetYaxis()->SetRangeUser(0.0, 2);
                can->cd(nthParam+1);
                hRcpPt[nthParam][iter].Draw(iter == 0 ? "" : "same");
                hRcpPt[nthParam][iter].SetLineColor(2000 + iter);
                hRcpPt[nthParam][iter].SetMarkerColor(2000 + iter);
                leg[nthParam]->AddEntry(&hRcpPt[nthParam][iter], Form("Iter %i", PlotIterations[iter]), "lep");
                if (iter == nIter-1){
                    tex.DrawLatex(0.5, 0.95, TString("1D R_{cp} ") + RcpTitles2[iCent] + "/" + RcpTitles2[jCent]);
                    leg[nthParam]->Draw();
                    DrawLineOne();
                }
            }
        }
    }

    can->SaveAs(outPdf);

}

void plotRcp2D(TCanvas *can){



    TLegend leg1 = TLegend(0.63, 0.55, 0.80, 0.82);
    leg1.SetBorderSize(0);
    leg1.SetFillStyle(0);

    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    tex.SetTextSize(0.055);
    tex.SetTextAlign(22);

    Int_t numOfComb = factorial(nCentralityBins) / (factorial(nCentralityBins - 2) * factorial(2));
    Int_t nVar = 7;
    TH2D hRcpPtVar[numOfComb][nVar][nIter];
    TH2D hRcpPtVarDen[numOfComb][nVar][nIter];
    TH1D hRcpPtVar_X[numOfComb][nVar][nIter];
    TH1D hRcpPtVar_Y[numOfComb][nVar][nIter];
    TH1D hRcpPtVarDen_X[numOfComb][nVar][nIter];
    TH1D hRcpPtVarDen_Y[numOfComb][nVar][nIter];
    TH1D hRcpPtVarDenReduced_Y[numOfComb][nVar][nIter];
    //Parameters for reduced
    Double_t minJetPt = 5;
    Double_t maxJetPt = 20-0.001;

    TLegend *leg[numOfComb];
    for (Int_t iComb = 0; iComb < numOfComb; iComb++) {
        leg[iComb] = new TLegend(0.13, 0.60, 0.30, 0.87);
        leg[iComb]->SetBorderSize(0);
        leg[iComb]->SetFillStyle(0);
        leg[iComb]->SetTextSize(0.055);
        leg[iComb]->SetTextFont(42);
    }


    for (int iCent  = 0; iCent < nCentralityBins; iCent++) {
        for (int iVar = 0; iVar < nVar; iVar++) {
            for (int iter = 0; iter < nIter; iter++) {
                int binMin = hUnfolded2D[iCent][0][iter].GetXaxis()->FindBin(5);
                int binMax = hUnfolded2D[iCent][0][iter].GetXaxis()->FindBin(19.9);
                hVarReduced_Y[iCent][iVar][iter] = *(TH1D *) hUnfolded2D[iCent][iVar][iter].ProjectionY(
                        Form("hVarReduced_Y_%i_%i_%i_%i", iCent, iCent, iVar, iter), binMin, binMax)->Clone(Form("hVarReduced_Y_clone_%d_%d_%d", iCent, iVar, iter));
            }
        }
    }





    //načtu OutputFile.root
    TFile *fLoad = new TFile("./Data/OutputFile.root");
    //načtu TGraphErrors z OutputFile.root
    //Unfolded Wide Z Cent = 0 Default 1 < D0pT < 10 GeV gStat
    //Unfolded Wide Z Cent = 0 Default 1 < D0pT < 10 GeV gSys
    //načtu je a uložím do jednoho TGraphErrors, kde gStat je statistická chyba a gSys je systémová chyba
    TGraphErrors *Z_gStat[3];
    TGraphErrors *Z_gSys[3];
    TGraphErrors *pT_gStat[3];
    TGraphErrors *pT_gSys[3];
    for (int iCent = 0; iCent < 2; iCent++) {
        Z_gStat[iCent] = (TGraphErrors *) fLoad->Get(Form("Results/Unfolded Wide Z Cent = %i Default 1 < D0pT < 10 GeV RCP gStat", iCent));
        Z_gSys[iCent] = (TGraphErrors *) fLoad->Get(Form("Results/Unfolded Wide Z Cent = %i Default 1 < D0pT < 10 GeV RCP gSys", iCent));
        Z_gStat[iCent]->SetMarkerStyle(20);
        Z_gStat[iCent]->SetMarkerColor(kBlack);
        Z_gStat[iCent]->SetLineColor(kBlack);
        Z_gSys[iCent]->SetMarkerStyle(20);
        Z_gSys[iCent]->SetMarkerColor(kBlack);
        Z_gSys[iCent]->SetLineColor(kBlack);
        pT_gStat[iCent] = (TGraphErrors *) fLoad->Get(Form("Results/Unfolded Wide p_{T} Cent = %i Default 1 < D0pT < 10 GeV RCP gStat", iCent));
        pT_gSys[iCent] = (TGraphErrors *) fLoad->Get(Form("Results/Unfolded Wide p_{T} Cent = %i Default 1 < D0pT < 10 GeV RCP gSys", iCent));
        pT_gStat[iCent]->SetMarkerStyle(20);
        pT_gStat[iCent]->SetMarkerColor(kBlack);
        pT_gStat[iCent]->SetLineColor(kBlack);
        pT_gSys[iCent]->SetMarkerStyle(20);
        pT_gSys[iCent]->SetMarkerColor(kBlack);
        pT_gSys[iCent]->SetLineColor(kBlack);
    }
fLoad->Close();



    //takovou smyčku, aby vždy centralita v čitateli byla menší než ve jmenovateli
    for (int iCent = 0; iCent < nCentralityBins; iCent++) {
        for (int jCent = iCent + 1; jCent < nCentralityBins; jCent++) {
            Int_t nthParam = GetCombIndex(iCent,jCent,nCentralityBins);
            for (int iVar = 0; iVar < nVar; iVar++) {

                can->Clear();
                can->Divide(3, 1);
                can->SetCanvasSize(1200, 400);
                        for (int iter = 0; iter < nIter; iter++) {
                            Int_t nthParam = GetCombIndex(iCent, jCent, nCentralityBins);
                            hRcpPtVar[nthParam][iVar][iter] = *(TH2D *) hUnfolded2D[iCent][iVar][iter].Clone(
                                    Form("hRcpPtVar_%i_%i_%i_%i", iCent, jCent, iVar, iter));
                            hRcpPtVar_X[nthParam][iVar][iter] = *(TH1D *) hRcpPtVar[nthParam][iVar][iter].ProjectionX(
                                    Form("hRcpPtVar_X_%i_%i_%i_%i", iCent, jCent, iVar, iter));
                            hRcpPtVar_Y[nthParam][iVar][iter] = *(TH1D *) hRcpPtVar[nthParam][iVar][iter].ProjectionY(
                                    Form("hRcpPtVar_Y_%i_%i_%i_%i", iCent, jCent, iVar, iter));
                            hRcpPtVarDen[nthParam][iVar][iter] = *(TH2D *) hUnfolded2D[jCent][iVar][iter].Clone(
                                    Form("hRcpPtVarDen_%i_%i_%i_%i", iCent, jCent, iVar, iter));
                            hRcpPtVarDen_X[nthParam][iVar][iter] = *(TH1D *) hRcpPtVarDen[nthParam][iVar][iter].ProjectionX(
                                    Form("hRcpPtVarDen_X_%i_%i_%i_%i", iCent, jCent, iVar, iter));
                            hRcpPtVarDen_Y[nthParam][iVar][iter] = *(TH1D *) hRcpPtVarDen[nthParam][iVar][iter].ProjectionY(
                                    Form("hRcpPtVarDen_Y_%i_%i_%i_%i", iCent, jCent, iVar, iter));
                            int binMin = hRcpPtVar_X[nthParam][iVar][iter].GetXaxis()->FindBin(minJetPt);
                            int binMax = hRcpPtVar_X[nthParam][iVar][iter].GetXaxis()->FindBin(maxJetPt);
                            hRcpPtVarReduced_Y[nthParam][iVar][iter] = *(TH1D *) hRcpPtVar[nthParam][iVar][iter].ProjectionY(
                                    Form("hRcpPtVarReduced_Y_%i_%i_%i_%i", iCent, jCent, iVar, iter), binMin, binMax);
                            //hVarReduced_Y[nthParam][iVar][iter] = *(TH1D *) hUnfolded2D[iCent][iVar][iter].ProjectionY(
                            //        Form("hVarReduced_Y_%i_%i_%i_%i", iCent, jCent, iVar, iter),binMin, binMax)->Clone(Form("asdasdas_%d",iCent)); //Fix this
                            hRcpPtVarDenReduced_Y[nthParam][iVar][iter] = *(TH1D *) hRcpPtVarDen[nthParam][iVar][iter].ProjectionY(
                                    Form("hRcpPtVarDenReduced_Y_%i_%i_%i_%i", iCent, jCent, iVar, iter), binMin,
                                    binMax);
                            can->cd(1);
                            hRcpPtVar_X[nthParam][iVar][iter].Divide(&hRcpPtVarDen_X[nthParam][iVar][iter]);
                            hRcpPtVar_X[nthParam][iVar][iter].Scale(1. * Ncoll[jCent] / Ncoll[iCent]);
                            hRcpPtVar_X[nthParam][iVar][iter].Scale(1. * NumberOfWEvents[jCent] / NumberOfWEvents[iCent]);
                            hRcpPtVar_X[nthParam][iVar][iter].GetYaxis()->SetRangeUser(0.0, 2);
                            hRcpPtVar_X[nthParam][iVar][iter].Draw(iter == 0 ? "" : "same");
                            hRcpPtVar_X[nthParam][iVar][iter].SetLineColor(2000 + iter);
                            hRcpPtVar_X[nthParam][iVar][iter].SetMarkerColor(2000 + iter);
                           if(jCent == 2 && iter == 0){
                                pT_gStat[iCent]->Draw("P same");
                                pT_gSys[iCent]->Draw("E2 same");
                                hRcpPtVar_X[nthParam][iVar][iter].Draw("same");

                            }
                            if (iter == nIter-1){
                                tex.DrawLatex(0.5, 0.95, TString("2D R_{cp} ") + RcpTitles2[iCent] + "/" + RcpTitles2[jCent]);
                                DrawLineOne();
                            }
                            can->cd(2);
                            hRcpPtVar_Y[nthParam][iVar][iter].Divide(&hRcpPtVarDen_Y[nthParam][iVar][iter]);
                            hRcpPtVar_Y[nthParam][iVar][iter].Scale(1. * Ncoll[jCent] / Ncoll[iCent]);
                            hRcpPtVar_Y[nthParam][iVar][iter].Scale(1. * NumberOfWEvents[jCent] / NumberOfWEvents[iCent]);
                            hRcpPtVar_Y[nthParam][iVar][iter].GetYaxis()->SetRangeUser(0.0, 2);
                            hRcpPtVar_Y[nthParam][iVar][iter].Draw(iter == 0 ? "" : "same");
                            hRcpPtVar_Y[nthParam][iVar][iter].SetLineColor(2000 + iter);
                            hRcpPtVar_Y[nthParam][iVar][iter].SetMarkerColor(2000 + iter);
                            if (iter == nIter-1){
                                tex.DrawLatex(0.5, 0.95, TString("2D R_{cp} ") + RcpTitles2[iCent] + "/" + RcpTitles2[jCent]);
                                DrawLineOne2(hRcpPtVar_Y[nthParam][iVar][iter].GetXaxis()->GetXmin(),
                                             hRcpPtVar_Y[nthParam][iVar][iter].GetXaxis()->GetXmax());
                            }
                            can->cd(3);
                            hRcpPtVarReduced_Y[nthParam][iVar][iter].Divide(&hRcpPtVarDenReduced_Y[nthParam][iVar][iter]);
                            hRcpPtVarReduced_Y[nthParam][iVar][iter].Scale(1. * Ncoll[jCent] / Ncoll[iCent]);
                            hRcpPtVarReduced_Y[nthParam][iVar][iter].Scale(1. * NumberOfWEvents[jCent] / NumberOfWEvents[iCent]);
                            hRcpPtVarReduced_Y[nthParam][iVar][iter].GetYaxis()->SetRangeUser(0.0, 2);
                            hRcpPtVarReduced_Y[nthParam][iVar][iter].Draw(iter == 0 ? "" : "same");
                            hRcpPtVarReduced_Y[nthParam][iVar][iter].SetLineColor(2000 + iter);
                            hRcpPtVarReduced_Y[nthParam][iVar][iter].SetMarkerColor(2000 + iter);
                            if (iter == nIter-1){
                                tex.DrawLatex(0.5, 0.95, TString("2D R_{cp} ") + RcpTitles2[iCent] + "/" + RcpTitles2[jCent]+ ", 5 < p_{T,Jet} < 20 GeV/c");
                                DrawLineOne2(
                                        hRcpPtVarReduced_Y[nthParam][iVar][iter].GetXaxis()->GetXmin(),
                                        hRcpPtVarReduced_Y[nthParam][iVar][iter].GetXaxis()->GetXmax());
                            }
                            if(iVar == 0 && jCent == 2 && iter == 0){
                                Z_gStat[iCent]->Draw("P same");
                                Z_gSys[iCent]->Draw("E2 same");
                                hRcpPtVarReduced_Y[nthParam][iVar][iter].Draw("same");
                            }
                        }
                        can->SaveAs(outPdf);
        }
    }
}
}


void ProcessSpectra(TH1 *R){
    for(int i = 1;i<=R->GetNbinsX()+1;i++){
        double val = R->GetBinContent(i);
        double er = R->GetBinError(i);
        double width = R->GetBinWidth(i);
        double center = fabs(R->GetBinCenter(i));

        R->SetBinContent(i,val/width/2./1.2/TMath::Pi()/center/0.0395/2);
        R->SetBinError(i,er/width/2./1.2/TMath::Pi()/center/0.0395/2);
    }
}
void plotFinalComp(TCanvas *can, const char *ScanDir = "default") {

    can->Clear();
    can->SetCanvasSize(1200, 1200);

    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    tex.SetTextSize(0.055);
    tex.SetTextAlign(22);

    can->Clear();
    can->Divide(3, 3);
    TGraph *gr[3];
    TGraph *grZ[3];

    Double_t x1[6] = {6, 8, 10, 12, 14, 17.5};
    Double_t y1[6] = {0.00010675022383376, 1.15832328625471E-05, 1.80019333289155E-06, 3.6331530050798E-07,
                      1.01646446037012E-07, 1.74234740059355E-08};
    Double_t x2[6] = {6, 8, 10, 12, 14, 17.5};
    Double_t y2[6] = {5.28937374675731E-05, 6.3302146943594E-06, 1.05021107963667E-06, 2.26261118182508E-07,
                      6.75751835538622E-08, 1.82983254467949E-08};
    Double_t x3[6] = {6, 8, 10, 12, 14, 17.5};
    Double_t y3[6] = {7.82739086833751E-06, 0.000001, 1.71412525329131E-07, 3.24069603578717E-08,
                      7.45315967438229E-09, 1.36380341328775E-09};
    /*
    gr[0] = new TGraph(6, x1, y1);
    gr[1] = new TGraph(6, x2, y2);
    gr[2] = new TGraph(6, x3, y3);
*/
    Double_t x2z[7] = {0.1, 0.3, 0.5, 0.65, 0.75, 0.85,0.95};

    TFile *file4 = TFile::Open("./Data/4.root");
    TH1 *h4[3];
    TH1 *h4Z[3];

    double nevents[3] = {1.1248969e+08, 3.5062620e+08, 4.7563863e+08};

    double x[3][6];
    double y[3][6];
    double x_2[3][7];
    double y_2[3][7];
    for (int ic = 0; ic < 3; ic++) {
        TString name4 = Form("Unfolded Wide p_{T} Cent = %d SI = 0 Iter = 4", ic);
        h4[ic] = (TH1F*)file4->Get(name4);
        h4[ic]->SetDirectory(0);
        h4[ic]->Scale(1.0 / nevents[ic]);
        ProcessSpectra(h4[ic]);
        for (int u = 0; u < 6; u++) {
            x[ic][u] = x1[u];
            y[ic][u] = h4[ic]->GetBinContent(h4[ic]->FindBin(x1[u]));
        }


        gr[ic] = new TGraph(6, x[ic], y[ic]);

        TString name4z = Form("Unfolded Wide Z Cent = %d SI = 0 Iter = 4", ic);
        h4Z[ic] = (TH1F*)file4->Get(name4z);
        h4Z[ic]->SetDirectory(0);
        h4Z[ic]->Scale(1.0 / nevents[ic]);
        ProcessSpectra(h4Z[ic]);
        for (int u = 0; u < 7; u++) {
            x_2[ic][u] = x2z[u];
            y_2[ic][u] = h4Z[ic]->GetBinContent(h4Z[ic]->FindBin(x2z[u]));
        }

        grZ[ic] = new TGraph(7, x_2[ic], y_2[ic]);



    }

    file4->Close();


    //OutputFile.root
    /*
    TFile *file7 = TFile::Open("./Data/NEIL.root");
    TGraph *Neil[3];
    Neil[0] = (TGraph*)file7->Get("Results/Unfolded Wide p_{T} Cent = 0 Default 1 < D0pT < 10 GeV gStat");
    Neil[1] = (TGraph*)file7->Get("Results/Unfolded Wide p_{T} Cent = 1 Default 1 < D0pT < 10 GeV gStat");
    Neil[2] = (TGraph*)file7->Get("Results/Unfolded Wide p_{T} Cent = 2 Default 1 < D0pT < 10 GeV gStat");
    file7->Close();*/
    TFile *file7 = TFile::Open("./Data/Komp.root");
    TGraph *Neil[3];
    Neil[0] = (TGraph*)file7->Get("neil_0");
    Neil[1] = (TGraph*)file7->Get("neil_1");
    Neil[2] = (TGraph*)file7->Get("neil_2");
    TGraph *NeilZ[3];
    NeilZ[0] = (TGraph*)file7->Get("neilZ_0");
    NeilZ[1] = (TGraph*)file7->Get("neilZ_1");
    NeilZ[2] = (TGraph*)file7->Get("neilZ_2");
    file7->Close();
    //Komp.root

    TLegend *legend = new TLegend(0.02, 0.64, 0.98, 0.68);
    legend->SetNColumns(9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.015);


    TString Variables[8] = { "p_{T}", "(p_{T},z)", "(p_{T},#lambda_{1}^{1})",
                             "(p_{T},#lambda_{1.5}^{1})", "(p_{T},#lambda_{2}^{1})",
                             "(p_{T},#lambda_{3}^{1})", "(p_{T},#lambda_{0.5}^{1})", "(p_{T},p_{T}^{D})" };



    //create results file
    TFile *fOutput = new TFile("Output/Results.root", "RECREATE");
        fOutput->cd();
    
    for (int iCent = 0; iCent < 3; iCent++) {

        can->cd(iCent + 7);
        gPad->SetLeftMargin(0.17);

        can->cd(iCent + 4);
        gPad->SetLeftMargin(0.17);

        can->cd(iCent + 1);
        gPad->SetLeftMargin(0.17);
        gPad->SetLogy();

        //int iter = 2;// niter
        int iter = GivenIter-1;// niter
        ////iter = savedIter-1;

    

        Unfolded1D[iCent] = *(TH1D*)hUnfoldedPt[iCent][iter].Clone(Form("Unfolded1D_Copy_%i", iCent));
        Unfolded1D[iCent].SetMarkerStyle(20);
        Unfolded1D[iCent].GetXaxis()->SetTitle("p_{T, Jet} [GeV/c]");
        Unfolded1D[iCent].GetYaxis()->SetTitle("1/N_{ev}d^{2}N/(2#timesB.R.#times2#pip_{T,Jet}dp_{T,Jet}d#eta_{Jet})");
        Unfolded1D[iCent].GetYaxis()->SetTitleOffset(1.6);
        Unfolded1D[iCent].Scale(1. / NumberOfWEvents[iCent]);
        NormalizeFinalSpectraPt(&Unfolded1D[iCent], 2000);

        Unfolded1D[iCent].Draw("");
        Unfolded1D[iCent].GetYaxis()->SetRangeUser(1e-11, 0.1);
        Unfolded1D[iCent].SetMarkerColor(kRed);
        Unfolded1D[iCent].SetLineColor(kRed);
        Unfolded1D[iCent].Write();


        if (iCent == 0)legend->AddEntry(&Unfolded1D[iCent],  Form("%s (%i it.)", Variables[0].Data(), iter+1), "lep");


        //new file
        //TFile *f = new TFile(Form("UnfoldedICS", iCent, iter), "READ");

        Unfolded2D_Y[iCent] = *((TH2D *) hUnfolded2D[iCent][0][iter].Clone())->ProjectionY(Form("Unfolded2D_Y_%i", iCent));
        Unfolded2D_Y[iCent].SetMarkerStyle(20);
        Unfolded2D_Y[iCent].GetXaxis()->SetTitle("z");
        Unfolded2D_Y[iCent].GetYaxis()->SetTitle("1/N_{ev}d^{2}N/(2#timesB.R.#times2#pip_{T,Jet}dp_{T,Jet}d#eta_{Jet})");
        Unfolded2D_Y[iCent].GetYaxis()->SetTitleOffset(1.6);
        //Unfolded2D_Y[iCent].Draw("same");
        Unfolded2D_Y[iCent].Scale(1. / NumberOfWEvents[iCent]);
        NormalizeFinalSpectraPt(&Unfolded2D_Y[iCent], 2001);
        Unfolded2D_Y[iCent].GetYaxis()->SetRangeUser(1e-5, 1);
        //Unfolded2D_Y[iCent].Write();

        for (int iVar = 0; iVar < 7; iVar++) {
          //  if (iVar >0) continue; // TO CHANGE


            Unfolded2D_X[iVar][iCent] = *((TH2D *) hUnfolded2D[iCent][iVar][iter].Clone())->ProjectionX(Form("Comparison_%i_%i_%i", iCent, iVar, iter));
            Unfolded2D_X[iVar][iCent].SetMarkerStyle(20);
            Unfolded2D_X[iVar][iCent].GetXaxis()->SetTitle("p_{T, Jet} [GeV/c]");
            Unfolded2D_X[iVar][iCent].GetYaxis()->SetTitle("1/N_{ev}d^{2}N/(2#dot2#pip_{T,Jet}dp_{T,Jet}d#eta_{Jet})");
            Unfolded2D_X[iVar][iCent].GetYaxis()->SetTitleOffset(1.6);
            Unfolded2D_X[iVar][iCent].Draw("same");
            Unfolded2D_X[iVar][iCent].Scale(1. / NumberOfWEvents[iCent]);
            NormalizeFinalSpectraPt(&Unfolded2D_X[iVar][iCent], 2001 + iVar);
            Unfolded2D_X[iVar][iCent].GetYaxis()->SetRangeUser(1e-11, 0.1);
 
            Unfolded2D_X[iVar][iCent].Write()

;
            if (iCent == 0) legend->AddEntry(&Unfolded2D_X[iVar][iCent], Form("%s (%i it.)", Variables[iVar+1].Data(), iter+1), "lep");

            can->cd(iCent + 1);
        }

        Neil[iCent]->SetLineColor(kBlack);
        Neil[iCent]->SetMarkerColor(kBlack);
        //set line width
        Neil[iCent]->SetLineWidth(1);
        Neil[iCent]->Draw("same");




        can->cd(iCent + 4);

        Unfolded1DRatio[iCent] = *(TH1D*)Unfolded1D[iCent].Clone(Form("Unfolded1DRatio_%i", iCent));
        Unfolded1DRatio[iCent].GetYaxis()->SetTitle("(Unfolded-Neil's)/Neil's");
        for (int i = 0; i < Unfolded1DRatio[iCent].GetXaxis()->GetNbins(); i++) {
            double x = Unfolded1DRatio[iCent].GetBinCenter(i+1);
            double y = Unfolded1DRatio[iCent].GetBinContent(i + 1);
            double norm = Neil[iCent]->Eval(x, nullptr, "S");
            Unfolded1DRatio[iCent].SetBinContent(i + 1, (Unfolded1DRatio[iCent].GetBinContent(i + 1) -norm)/ norm);
        }
        for (int i = 0; i < Unfolded1DRatio[iCent].GetXaxis()->GetNbins(); i++) {
            if (Unfolded1DRatio[iCent].FindBin(Unfolded1DRatio[iCent].GetXaxis()->GetBinCenter(i + 1)) < Unfolded1DRatio[iCent].FindBin(1)) {
                Unfolded1DRatio[iCent].SetBinContent(i + 1, 0);
                Unfolded1DRatio[iCent].SetBinError(i + 1, 0);
            }
        }
        Unfolded1DRatio[iCent].SetMarkerStyle(20);
        Unfolded1DRatio[iCent].Draw("");
        Unfolded1DRatio[iCent].GetYaxis()->SetRangeUser(-1.0, 1.0);

        for (int iVar = 0; iVar < 7; iVar++) {
            //if (iVar >0) continue; // TO CHANGE

            can->cd(iCent + 4);

            Unfolded2DRatio[iVar][iCent] = *(TH1D*)Unfolded2D_X[iVar][iCent].Clone(Form("Unfolded2DRatio_%i_%i", iVar, iCent));
            Unfolded2DRatio[iVar][iCent].GetYaxis()->SetTitle("(Unfolded-Neil's)/Neil's");
            for (int i = 0; i < Unfolded2DRatio[iVar][iCent].GetXaxis()->GetNbins(); i++) {
                double x = Unfolded2DRatio[iVar][iCent].GetBinCenter(i+1);
                double y = Unfolded2DRatio[iVar][iCent].GetBinContent(i + 1);
                double norm = Neil[iCent]->Eval(x, nullptr, "S");
                Unfolded2DRatio[iVar][iCent].SetBinContent(i + 1, (Unfolded2DRatio[iVar][iCent].GetBinContent(i + 1) -norm) / norm);
              }
            for (int i = 0; i < Unfolded2DRatio[iVar][iCent].GetXaxis()->GetNbins(); i++) {
                if (Unfolded2DRatio[iVar][iCent].FindBin(Unfolded2DRatio[iVar][iCent].GetXaxis()->GetBinCenter(i + 1)) < Unfolded2DRatio[iVar][iCent].FindBin(1)) {
                    Unfolded2DRatio[iVar][iCent].SetBinContent(i + 1, 0);
                    Unfolded2DRatio[iVar][iCent].SetBinError(i + 1, 0);
                }
            }


            Unfolded2DRatio[iVar][iCent].SetMarkerStyle(20);
            Unfolded2DRatio[iVar][iCent].Draw("same");
            Unfolded2DRatio[iVar][iCent].GetYaxis()->SetRangeUser(0.0, 2.0);


            can->cd(iCent + 7);

            Unfolded2DRatioVar[iVar][iCent] = *(TH1D*)Unfolded2D_X[iVar][iCent].Clone(Form("Unfolded2DRatioVar_%i_%i", iVar, iCent));
            Unfolded2DRatioVar[iVar][iCent].Reset();

            for (int iBin = 1; iBin <= Unfolded2DRatioVar[iVar][iCent].GetNbinsX(); iBin++) {


                double x  = 0;
                double y  = 0;
                double ex = 0;
                double ey = 0;
                if (iVar == 0){
                x  = Unfolded1D[iCent].GetBinContent(iBin);
                y  = Unfolded2D_X[0][iCent].GetBinContent(iBin);
                ex = Unfolded1D[iCent].GetBinError(iBin);
                ey = Unfolded2D_X[0][iCent].GetBinError(iBin);

                Unfolded2DRatioVar[iVar][iCent].SetMarkerStyle(20);
                Unfolded2DRatioVar[iVar][iCent].SetMarkerColor(2000);
                Unfolded2DRatioVar[iVar][iCent].SetLineColor(2000);

                } else { // skip pT vs pT case
                x  = Unfolded2D_X[iVar][iCent].GetBinContent(iBin);
                y  = Unfolded2D_X[0][iCent].GetBinContent(iBin);
                ex = Unfolded2D_X[iVar][iCent].GetBinError(iBin);
                ey = Unfolded2D_X[0][iCent].GetBinError(iBin);
                }

                if (y != 0) {
                    double ratio = (x - y) / y;
                    double err   = std::sqrt( (ex*ex)/(y*y) + (x*x*ey*ey)/(y*y*y*y) ); // error propagation
                    Unfolded2DRatioVar[iVar][iCent].SetBinContent(iBin, ratio);
                    Unfolded2DRatioVar[iVar][iCent].SetBinError(iBin, err);
                } else {
                    Unfolded2DRatioVar[iVar][iCent].SetBinContent(iBin, 0);
                    Unfolded2DRatioVar[iVar][iCent].SetBinError(iBin, 0);
                }
            }

           // Unfolded2DRatioVar[iVar][iCent] = *(TH1D*)Unfolded2D_X[iVar][iCent].Clone(Form("Unfolded2DRatioVar_%i_%i", iVar, iCent));
            //Unfolded2DRatioVar[iVar][iCent].Divide(&Unfolded1D[iCent]);
            Unfolded2DRatioVar[iVar][iCent].GetYaxis()->SetTitle("(Unfolded-2D(pt,z))/2D(pt,z)");
            Unfolded2DRatioVar[iVar][iCent].GetYaxis()->SetTitleOffset(1.2);
            Unfolded2DRatioVar[iVar][iCent].Draw(iVar==0?"":"same");
            Unfolded2DRatioVar[iVar][iCent].GetYaxis()->SetRangeUser(-1.0, 1.0);


        }
        can->cd(iCent + 4);
        DrawLineZero();
        DrawLineVertical(5);
        DrawLineVertical(20);
        can->cd(iCent + 7);
        DrawLineZero();
        DrawLineVertical(5);
        DrawLineVertical(20);

}


fOutput->Close();

    can->cd();
    legend->AddEntry(gr[0], "Neil's", "l");
    legend->Draw("same");
    can->SaveAs(outPdf);

    /*
    //Update file OutputSpectra.root, if it doesnt exist, create it
    TFile *fMach = new TFile("./OutputPdf/pTCheck/OutputSpectra"+TString(runId)+".root", "RECREATE");
       Bool_t oldAddDir = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);   // NIC se nebude automaticky lepit do fMach
    if (fMach->IsZombie()) {
        fMach = new TFile("./OutputPdf/pTCheck/OutputSpectra"+TString(runId)+".root", "RECREATE");
    }
    fMach->cd();
    */

    TString fMachDir = TString(ScanDir) + "/Output";

    if (gSystem->AccessPathName(fMachDir.Data())) {
        int status = gSystem->mkdir(fMachDir.Data(), kTRUE);  // kTRUE = create parent dirs too
        if (status != 0 && gSystem->AccessPathName(fMachDir.Data())) {
            cout << "[fMach] ERROR: cannot create directory: " << fMachDir << endl;
            gSystem->Exit(1);
        }
    }

    TString fMachName = fMachDir + "/OutputSpectra" + TString(runId) + ".root";
    TFile *fMach = new TFile(fMachName.Data(), "RECREATE");

    Bool_t oldAddDir = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);   // NIC se nebude automaticky lepit do fMach

    if (!fMach || fMach->IsZombie()) {
        cout << "[fMach] ERROR: cannot open file: " << fMachName << endl;
        gSystem->Exit(1);
    }

    cout << "[fMach] writing to: " << fMachName << endl;
    fMach->cd();
    
    TH1D _Unfolded1D[3];
    TH1D _Unfolded2D_X[7][3];

 for (int iCent = 0; iCent < 3; iCent++) {

    for (int iTer = 2; iTer < 5; iTer++) {
        //clear Unfolded1D and Unfolded2D_X
        _Unfolded1D[iCent].Reset();

        _Unfolded1D[iCent] = *(TH1D*)hUnfoldedPt[iCent][iTer].Clone(Form("Unfolded1D_Copy_%i", iCent));
        _Unfolded1D[iCent].SetMarkerStyle(20);
        _Unfolded1D[iCent].GetXaxis()->SetTitle("p_{T, Jet} [GeV/c]");
        _Unfolded1D[iCent].GetYaxis()->SetTitle("1/N_{ev}d^{2}N/(2#timesB.R.#times2#pip_{T,Jet}dp_{T,Jet}d#eta_{Jet})");
        _Unfolded1D[iCent].GetYaxis()->SetTitleOffset(1.6);
        _Unfolded1D[iCent].Scale(1. / NumberOfWEvents[iCent]);
        NormalizeFinalSpectraPt(&_Unfolded1D[iCent], 2000);

        TString baseName = Form("it%d_%s",iTer, Method.Data());

        _Unfolded1D[iCent].SetDirectory(0); 
        _Unfolded1D[iCent].SetLineColor(kRed);
        _Unfolded1D[iCent].SetMarkerColor(kRed);
        TString fullTitle1D = Form("d0pt_%d_%s", iCent, baseName.Data());
        _Unfolded1D[iCent].SetTitle(fullTitle1D);
        _Unfolded1D[iCent].Write(fullTitle1D, TObject::kOverwrite);

        for (int iVar = 0; iVar < 7; iVar++) {
        _Unfolded2D_X[iVar][iCent].Reset();

        _Unfolded2D_X[iVar][iCent] = *((TH2D *) hUnfolded2D[iCent][iVar][iTer].Clone())->ProjectionX(Form("Comparison_%i_%i_%i", iCent, iVar, iTer));
        _Unfolded2D_X[iVar][iCent].SetMarkerStyle(20);
        _Unfolded2D_X[iVar][iCent].GetXaxis()->SetTitle("p_{T, Jet} [GeV/c]");
        _Unfolded2D_X[iVar][iCent].GetYaxis()->SetTitle("1/N_{ev}d^{2}N/(2#dot2#pip_{T,Jet}dp_{T,Jet}d#eta_{Jet})");
        _Unfolded2D_X[iVar][iCent].GetYaxis()->SetTitleOffset(1.6);
        _Unfolded2D_X[iVar][iCent].Draw("same");
        _Unfolded2D_X[iVar][iCent].Scale(1. / NumberOfWEvents[iCent]);
        NormalizeFinalSpectraPt(&_Unfolded2D_X[iVar][iCent], 2001 + iVar);
        _Unfolded2D_X[iVar][iCent].GetYaxis()->SetRangeUser(1e-11, 0.1);

        _Unfolded2D_X[iVar][iCent].SetDirectory(0);
        _Unfolded2D_X[iVar][iCent].SetLineColor(kGreen);
        _Unfolded2D_X[iVar][iCent].SetMarkerColor(kGreen);

        TString fullTitle2D_X = Form("d0ptLambda%d_%d_%s", iVar, iCent, baseName.Data());
        _Unfolded2D_X[iVar][iCent].SetTitle(fullTitle2D_X);
        _Unfolded2D_X[iVar][iCent].Write(fullTitle2D_X, TObject::kOverwrite);
        }

        }
    fMach->Write();
}
    can->Clear();
    can->SetCanvasSize(1200, 800);
    can->Divide(3, 2);

    //TGraphErrors *NeilZ[3];
    TH2D hResponseTruth2D[3];
    TH1D hResponseTruth2D_Y[3];


    for (int iCent = 0; iCent < 3; iCent++){
        int binMin = hUnfolded2D[iCent][0][0].GetXaxis()->FindBin(5);
        int binMax = hUnfolded2D[iCent][0][0].GetXaxis()->FindBin(19.9);
     hResponseTruth2D[iCent] = *(TH2D *) rurResponse2D[iCent][0].Htruth()->Clone(
             Form("hResponseTruth22D_%i_%i", iCent,0 ));
    hResponseTruth2D_Y[iCent] = *(TH1D *) hResponseTruth2D[iCent].ProjectionY(Form("Copy__%i",iCent),binMin,binMax)->Clone(
            Form("hResponseTruth22D_Y%i_%i", iCent, 0));
    }



    TLegend *legend2 = new TLegend(0.02, 0.44, 0.98, 0.58);
    legend2->SetNColumns(4);
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(0);
    legend2->SetTextSize(0.015);

    legend2->AddEntry(NeilZ[0], "Neil's", "l");


    TH1D hVarReduced_YCopy[3];
    TH1D hVarReduced_Ratio[3];
    int iter = GivenIter-1;

    for (int iCent = 0; iCent < 3; iCent++) {
        can->cd(iCent + 1);
        gPad->SetLeftMargin(0.17);
        gPad->SetLogy();

        hVarReduced_YCopy[iCent] = *(TH1D *) hVarReduced_Y[iCent][0][iter].Clone(
                    Form("hRcpPtVarReduced_YCopy_%i", iCent));
        hVarReduced_YCopy[iCent].Draw("same");
        hVarReduced_YCopy[iCent].Scale(1. / NumberOfWEvents[iCent]);

        NormalizeFinalSpectraPt(&hVarReduced_YCopy[iCent], 2000);
        NeilZ[iCent]->SetLineColor(kGreen + 4);
        NeilZ[iCent]->Draw("same");
        hVarReduced_YCopy[iCent].GetYaxis()->SetTitle("1/N_{ev}d^{2}N/(2#timesB.R.#times2#pip_{T,Jet}dzd#eta_{Jet})");
        hVarReduced_YCopy[iCent].GetYaxis()->SetRangeUser(1e-6, 1e-1);
        //offset
        hVarReduced_YCopy[iCent].GetYaxis()->SetTitleOffset(1.6);

        NormalizeFinalSpectraPt(&hResponseTruth2D_Y[iCent], 2005);
        hResponseTruth2D_Y[iCent].GetYaxis()->SetTitle("1/N_{ev}d^{2}N/(2#timesB.R.#timesdzd#eta_{Jet})");
        hResponseTruth2D_Y[iCent].SetLineStyle(3);
        hResponseTruth2D_Y[iCent].Scale(1.*hVarReduced_YCopy[iCent].Integral()/hResponseTruth2D_Y[iCent].Integral());

        jetZ[iCent]->Scale(1.*hVarReduced_YCopy[iCent].Integral()/jetZ[iCent]->Integral());
        jetZ[iCent]->SetLineColor(kGreen+2);
        jetZ[iCent]->SetMarkerColor(kGreen+2);
        jetZ[iCent]->SetLineStyle(3);

    can->cd(iCent + 4);
        gPad->SetLogy(0);

        gPad->SetLeftMargin(0.17);
    hVarReduced_Ratio[iCent] = *(TH1D *) hVarReduced_YCopy[iCent].Clone(
            Form("hVarReduced_Ratio_%i", iCent));
    hVarReduced_Ratio[iCent].GetYaxis()->SetTitle("(Unfolded-Neil's)/Neil's");
    for (int i = 0; i < hVarReduced_Ratio[iCent].GetXaxis()->GetNbins(); i++) {
        double x = hVarReduced_Ratio[iCent].GetBinCenter(i + 1);
        double y = hVarReduced_Ratio[iCent].GetBinContent(i + 1);
        double norm = NeilZ[iCent]->Eval(x, nullptr, "S");
        hVarReduced_Ratio[iCent].SetBinContent(i + 1, (hVarReduced_Ratio[iCent].GetBinContent(i + 1)-norm) / norm);
    }
    for (int i = 0; i < hVarReduced_Ratio[iCent].GetXaxis()->GetNbins(); i++) {
        if (hVarReduced_Ratio[iCent].FindBin(hVarReduced_Ratio[iCent].GetXaxis()->GetBinCenter(i + 1)) < hVarReduced_Ratio[iCent].FindBin(0.1)) {
            hVarReduced_Ratio[iCent].SetBinContent(i + 1, 0);
            hVarReduced_Ratio[iCent].SetBinError(i + 1, 0);
        }
    }
    hVarReduced_Ratio[iCent].SetMarkerStyle(20);
    hVarReduced_Ratio[iCent].GetYaxis()->SetRangeUser(-1, 1);
    hVarReduced_Ratio[iCent].SetLineColor(2000);
    hVarReduced_Ratio[iCent].SetMarkerColor(2000);
    hVarReduced_Ratio[iCent].Draw("");

    //přerušovaná čára od 0 do 1
        DrawLineZero2(0,1.01);
    }
    legend2->AddEntry(&hVarReduced_YCopy[0], "z unfolded (5<p_{T,Jet}<20)", "lp");
   //// legend2->AddEntry(&hResponseTruth2D_Y[0], "MC true (weighted)", "lp");
    ////legend2->AddEntry(jetZ[0], "MC true (raw)", "lp");


    can->cd();
    legend2->Draw("same");
    can->SaveAs(outPdf);

    TLegend *legLam[7][3];

    for(int iLamb = 0; iLamb < 7; iLamb++){
        can->Clear();
        can->SetCanvasSize(1200, 400);
        can->Divide(3, 1);

        for (int iCent = 0; iCent < 3; iCent++) {

            for (int iIter = 0; iIter < nIter; iIter++) {

            can->cd(iCent + 1);
            hVarReduced_Y_ang[iLamb][iCent][iIter] = *(TH1D *) hVarReduced_Y[iCent][iLamb][iIter].Clone(
                    Form("hRcpPtVarReduced_YCopy_%i_%i", iCent, iLamb));
            //logscale
            gPad->SetLogy();
            gPad->SetLeftMargin(0.17);

            hVarReduced_Y_ang[iLamb][iCent][iIter].Scale(1. / NumberOfWEvents[iCent]);
            NormalizeFinalSpectra(&hVarReduced_Y_ang[iLamb][iCent][iIter], 2000);
            TString xtitle = hVarReduced_Y_ang[iLamb][iCent][iIter].GetXaxis()->GetTitle();
            hVarReduced_Y_ang[iLamb][iCent][iIter].GetYaxis()->SetTitle(Form("1/N_{ev}d^{2}N/(2#timesB.R.#timesd%sd#eta_{Jet})", xtitle.Data()));
            hVarReduced_Y_ang[iLamb][iCent][iIter].GetYaxis()->SetRangeUser(1e-6, 1e+1);
            hVarReduced_Y_ang[iLamb][iCent][iIter].GetYaxis()->SetTitleOffset(1.6);

            if (iIter != iter) continue; // only plot the given iter
            hVarReduced_Y_ang[iLamb][iCent][iIter].Draw();

            //save to fMach
            hVarReduced_Y_ang[iLamb][iCent][iIter].SetDirectory(0); // být autonomní – pro jistotu
            //set color blue
            hVarReduced_Y_ang[iLamb][iCent][iIter].SetLineColor(kBlue);
            hVarReduced_Y_ang[iLamb][iCent][iIter].SetMarkerColor(kBlue);
            //set title
/*
            hVarReduced_Y_ang[iLamb][iCent][iIter].SetTitle(Form("Lambda%i_%d_fj_%d_CONeg_%d_minJetPtRecoCut_%.1f_it%d_%s",
                                                         iLamb, iCent, FONLLjet, CutOfNegative, (double)minJetPtRecoCut, savedIter,Method.Data()));
            fMach->cd();

            hVarReduced_Y_ang[iLamb][iCent][iIter].Write(Form("Lambda%i_%d_fj_%d_CONeg_%d_minJetPtRecoCut_%.1f_it%d_%s",
                                                       iLamb, iCent, FONLLjet, CutOfNegative, (double)minJetPtRecoCut, savedIter,Method.Data()),
                                              TObject::kOverwrite);
            fMach->Write();
*/
            legLam[iLamb][iCent] = new TLegend(0.50, 0.75, 0.65, 0.87);
            legLam[iLamb][iCent]->SetBorderSize(0);
            legLam[iLamb][iCent]->SetFillStyle(0);
            legLam[iLamb][iCent]->SetTextSize(0.03);
            legLam[iLamb][iCent]->SetTextFont(42);
           // legLam[iLamb-1][iCent]->SetHeader(Form("Angularity %i", iLamb));
            legLam[iLamb][iCent]->AddEntry((TH1D *) 0, Form( "%i-%i%%", centrality[iCent][0],centrality[iCent][1]), "");
            legLam[iLamb][iCent]->AddEntry(&hVarReduced_Y_ang[iLamb][iCent][iIter],"5 < p_{T,Jet} [GeV/c] < 20 ", "lp");
            legLam[iLamb][iCent]->AddEntry((TH1D *) 0, Form( "Iter %i", iter + 1), "");

            legLam[iLamb][iCent]->Draw("same");

            }
        }

        can->SaveAs(outPdf);

    }

TH1::AddDirectory(oldAddDir);

fMach->cd();
for(int iLamb = 0; iLamb < 7; iLamb++){
    for (int iCent = 0; iCent < 3; iCent++) {
        for (int iIter = 2; iIter < 5; iIter++) {


        TString baseName = Form("%d_it%d_%s",iCent,iIter, Method.Data());
        hVarReduced_Y_ang[iLamb][iCent][iIter].SetDirectory(0);
        hVarReduced_Y_ang[iLamb][iCent][iIter].SetLineColor(kBlue);
        hVarReduced_Y_ang[iLamb][iCent][iIter].SetMarkerColor(kBlue);
        TString xTitle = hVarReduced_Y_ang[iLamb][iCent][iIter].GetXaxis()->GetTitle();
        if (xTitle.IsNull()) xTitle = "x";  // fallback, kdyby bylo prázdné

        TString fullTitle = Form("Lambda%i_%s", iLamb, baseName.Data());

        hVarReduced_Y_ang[iLamb][iCent][iIter].SetTitle(fullTitle);
        hVarReduced_Y_ang[iLamb][iCent][iIter].Write(fullTitle, TObject::kOverwrite);
        }

}
}


// hRcpPtVarReduced_Y[nthParam][iVar][iter].Draw("same");Int_t nthParam = GetCombIndex(iCent, jCent, nCentralityBins);
//uložím hRcpPtVarReduced_Y[nthParam][iVar][iter] do FinalOutput.root, kde název Form("RCP_5_20_%i_%d_fj_%d_CONeg_%d_minJetPtRecoCut_%.1f_it%d_%d/%d_%s", iLamb, iCent, FONLLjet, CutOfNegative, (double)minJetPtRecoCut, savedIter , iCent, jCent,Method.Data()));
    for (int iCent = 0; iCent < nCentralityBins; iCent++) {
        for (int jCent = iCent + 1; jCent < nCentralityBins; jCent++) {
            Int_t nthParam = GetCombIndex(iCent, jCent, nCentralityBins);
            for (int iVar = 0; iVar < 7; iVar++) {
                for (int iter = 2; iter < 5; iter++) {
                                        //set all colours to black
        
            //TString baseName = Form("fj_%d_CONeg_%d_minJetPtRecoCut_%.1f_it%d_%s",
             //           FONLLjet, CutOfNegative, (double)minJetPtRecoCut, iter, Method.Data());
                        TString baseName = Form("%d_it%d_%s",nthParam,iter, Method.Data());
                                        //set all colours to black
                    hRcpPtVarReduced_Y[nthParam][iVar][iter].SetLineColor(kBlack);
                    hRcpPtVarReduced_Y[nthParam][iVar][iter].SetMarkerColor(kBlack);
                    TString fullTitle = Form("RCP_5_20_Lambda%i_%s", iVar, baseName.Data());
                    hRcpPtVarReduced_Y[nthParam][iVar][iter].Write(fullTitle, TObject::kOverwrite);

                }
            }
        }
    }


fMach->Close();
}

//-----1D-----------------
void plotRCP(int step, bool svd, TCanvas *can, TH1D *hRcpPt_[],
             TH1D *hUnfoldedNumerator[], TH1D *hUnfoldedDenominator[],
             const Double_t Ncoll[], Int_t iCent, Int_t jCent,
             TString RcpTitles2[], const char *OutputFile, Int_t iSuper) {

    // Počet iterací
    const int nIterations = svd ? nKterm : nIter;

    // Přednastavení legend a textů
    TLegend *leg1 = new TLegend(0.63, 0.55, 0.80, 0.82);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);

    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    tex.SetTextSize(0.055);

    TH1D *hzPt[nIterations];

    // Nastavení canvasu
    can->cd(step);
    gPad->SetLeftMargin(0.15);

    for (int iter = 0; iter < nIterations; iter++) {
        // Výpočet histogramů
        TH1D *numerator = hUnfoldedNumerator[iter];
        TH1D *denominator = hUnfoldedDenominator[iter];

        hRcpPt_[iter] = (TH1D *) numerator->Clone(Form("hRcpPt_%i", iter));
        hRcpPt_[iter]->Divide(denominator);
        hRcpPt_[iter]->Scale(Ncoll[jCent] / Ncoll[iCent]);
        hRcpPt_[iter]->Scale(NumberOfWEvents[jCent] / NumberOfWEvents[iCent]);

        hzPt[iter] = (TH1D *) hRcpPt_[iter]->Clone(Form("hzPt_%i", iter));
        hzPt[iter]->SetLineColor(2000 + iter);
        hzPt[iter]->SetMarkerColor(2000 + iter);
        hzPt[iter]->SetMarkerStyle(20);

        TString label = svd ? Form("kterm%i", kterm[iter]) : Form("Iter%i", PlotIterations[iter]);
        leg1->AddEntry(hzPt[iter], label, "lep");
    }

    // Určení maxima a minima
    double max = findMax(hzPt, 0, nIterations);
    double min = findMin(hzPt, 100, nIterations);

    // Kreslení histogramů
    for (int iter = 0; iter < nIterations; iter++) {
        hzPt[iter]->Draw(iter == 0 ? "" : "same");
        hzPt[iter]->GetXaxis()->SetTitleOffset(0.5);
        hzPt[iter]->GetXaxis()->SetLabelOffset(-0.023);
    }

    gPad->SetTopMargin(0.127);
    hzPt[0]->GetYaxis()->SetRangeUser(0, std::max(max + 0.6, 2.0));
    gPad->SetLogy(0);
    DrawLineOne();

    // Popisky
    tex.DrawLatex(0.19, 0.05, (svd ? "SVD " : "") + TString("R_{cp}") +
                              RcpTitles2[iCent] + "/" + RcpTitles2[jCent] +
                              Form(", sIter: %i", iSuper));

    leg1->Draw();
}

//------2D----------------
void
plotRCP(TCanvas *can, TH2D *hRcp[], TH2D *hUnfoldedNumerator[], TH2D *hUnfoldedDenominator[], const Double_t Ncoll[],
        Int_t iCent, Int_t jCent, TString RcpTitles2[], TString var, const char *OutputFile) {

    can->cd();
    can->Clear();
    can->Divide(2, 1);

    TLegend *leg = new TLegend(0.6, 0.6, 0.9, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.055);

    TLegend *leg1 = new TLegend(0.60, 0.55, 1.2, 0.9);
    TLegend *leg2 = new TLegend(0.60, 0.55, 1.2, 0.9);
    for (auto leg: {leg1, leg2}) {
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
    }

    TH1D *hzX[nIter];
    TH1D *hzY[nIter];


    gPad->SetLeftMargin(0.15);

    TH1D *hDenomX[nIter];
    TH1D *hDenomY[nIter];

    for (int iter = 0; iter < nIter; iter++) {

        hRcp[iter] = (TH2D *) hUnfoldedNumerator[iter]->Clone(
                TString("hUnfProjXPt") + Form("_%i_c%i_%i_%s", iter, iCent, jCent, var.Data()));

        hDenomX[iter] = (TH1D *) hUnfoldedDenominator[iter]->ProjectionX(
                TString("hDenomX") + Form("_%i_c%i_%i_%s", iter, iCent, jCent, var.Data()));
        hDenomY[iter] = (TH1D *) hUnfoldedDenominator[iter]->ProjectionY(
                TString("hDenomY") + Form("_%i_c%i_%i_%s", iter, iCent, jCent, var.Data()));

        //Projekce na x-ovou osu
        hzX[iter] = (TH1D *) hRcp[iter]->ProjectionX(
                TString("hzX") + Form("_%i_c%i_%i_%s", iter, iCent, jCent, var.Data()))->Clone(
                TString("C_hzX") + Form("_%i_c%i_%i_%s", iter, iCent, jCent, var.Data()));
        hzX[iter]->Divide(hDenomX[iter]);
        hzX[iter]->Scale(1. * Ncoll[jCent] / Ncoll[iCent]);
        hzX[iter]->Scale(1. * NumberOfWEvents[jCent] / NumberOfWEvents[iCent]);
        hzX[iter]->SetLineColor(2000 + iter);
        hzX[iter]->SetMarkerColor(2000 + iter);
        hzX[iter]->SetMarkerStyle(20);

        //Projekce na y-ovou osu
        hzY[iter] = (TH1D *) hRcp[iter]->ProjectionY(
                TString("hzY") + Form("_%i_c%i_%i_%s", iter, iCent, jCent, var.Data()))->Clone(
                TString("C_hzY") + Form("_%i_c%i_%i_%s", iter, iCent, jCent, var.Data()));
        hzY[iter]->Divide(hDenomY[iter]);
        hzY[iter]->Scale(1. * Ncoll[jCent] / Ncoll[iCent]);
        hzY[iter]->Scale(1. * NumberOfWEvents[jCent] / NumberOfWEvents[iCent]);
        hzY[iter]->SetLineColor(2000 + iter);
        hzY[iter]->SetMarkerColor(2000 + iter);
        hzY[iter]->SetMarkerStyle(20);

        leg1->AddEntry(hzX[iter], Form("Iter%i", PlotIterations[iter]), "lep");
        leg2->AddEntry(hzY[iter], Form("Iter%i", PlotIterations[iter]), "lep");

        can->cd(1);
        //gPad->SetLogy();
        hzX[iter]->Draw(iter == 0 ? "" : "same");

        can->cd(2);
        //gPad->SetLogx();

        hzY[iter]->Draw(iter == 0 ? "" : "same");
    }

    Double_t Xmax = findMax(hzX, 0, nIter);
    Double_t Xmin = findMin(hzX, 100, nIter);
    Double_t Ymax = findMax(hzY, 0, nIter);
    Double_t Ymin = findMin(hzY, 100, nIter);

    ////hzX[0]->GetYaxis()->SetRangeUser(Xmin-0.2, Xmax+0.6);
    ////hzY[0]->GetYaxis()->SetRangeUser(Ymin-0.2, Ymax+0.6);

    hzX[0]->GetYaxis()->SetRangeUser(0, 2);
    hzY[0]->GetYaxis()->SetRangeUser(0, 2);

    can->cd(1);
    if (Xmax > 10 && false) {
        can->cd(1)->Clear();
        can->cd();
        //Nový tex, kde bude číslo 10 o stejné velikosti, jako text na y-ové ose
        TLatex *tex2 = new TLatex();
        tex2->SetNDC();
        tex2->SetTextFont(42);
        tex2->SetTextSize(0.02);
        tex2->DrawLatex(0.0355, 0.4925, "10");

        //Rozdělím cd(1) na dva grafy pod sebou
        can->cd(1);


        TPad *myPad1 = new TPad("myPad1", "myPad1", 0, 0.5, 1.0, 0.945);
        myPad1->SetBottomMargin(0.0);
        myPad1->SetBorderMode(0);
        TPad *myPad2 = new TPad("myPad2", "myPad2", 0, 0.09, 1.0, 0.5);
        myPad2->SetTopMargin(0.0);
        myPad2->SetBorderMode(0);


        myPad1->Draw();
        myPad1->cd();


        //Nastavím logscale
        //gPad->SetLogy();
        //Nakreslím první graf
        TH1D *hzXCopy[nIter];
        TH1D *hzXCopy2[nIter];

        hzXCopy[0] = (TH1D *) hzX[0]->Clone(TString("hzXCopy") + Form("_%i_c%i_%i_%s", 0, iCent, jCent, var.Data()));

        hzXCopy[0]->Draw("E0");
        hzXCopy[0]->GetYaxis()->SetRangeUser(10, Xmax * 3);
        //Odstraním popisek na x-ové ose
        hzXCopy[0]->GetXaxis()->SetTitle("");
        //Odstraním x-ovou osu
        hzXCopy[0]->GetXaxis()->SetLabelSize(0);
        //Odstraním x-ový okraj grafu
        hzXCopy[0]->GetXaxis()->SetTickLength(0);

        //Posunu graf níže, aby se dotýkal grafu v cd(2)
        hzXCopy[0]->GetYaxis()->SetTitleOffset(1.5);
        //smyčka přes zbytek iterací
        for (int i = 1; i < nIter; i++) {
            hzXCopy[i] = (TH1D *) hzX[i]->Clone(TString("hzXCopy") + Form("_%i_c%i_%i", 0, iCent, jCent));
            hzXCopy[i]->Draw("sameE0");
        }


        //druhý pad v cd(1)
        can->cd(1);

        myPad2->Draw();
        myPad2->cd();
        hzXCopy2[0] = (TH1D *) hzX[0]->Clone(TString("hzXCopy2") + Form("_%i_c%i_%i_%s", 0, iCent, jCent, var.Data()));
        hzXCopy2[0]->Draw("E0");
        hzXCopy2[0]->GetYaxis()->SetRangeUser(-1, 10);
        //smyčka přes zbytek iterací
        for (int i = 1; i < nIter; i++) {
            hzXCopy2[i] = (TH1D *) hzX[i]->Clone(
                    TString("hzXCopy2") + Form("_%i_c%i_%i_%s", i, iCent, jCent, var.Data()));
            hzXCopy2[i]->Draw("sameE0");
        }

    }


    //čára kolem jedničky
    double LevyOkrajPrvnihoBinu = hzX[0]->GetXaxis()->GetBinLowEdge(1);
    double PravyOkrajPoslednihoBinu = hzX[0]->GetXaxis()->GetBinUpEdge(hzX[0]->GetNbinsX());

    TLine *line = new TLine(LevyOkrajPrvnihoBinu, 1, PravyOkrajPoslednihoBinu, 1);
    line->SetLineStyle(2);
    line->Draw("same");

    leg1->Draw("same");
    can->cd(2);


    //čára kolem jedničky
    double LevyOkrajPrvnihoBinu2 = hzY[0]->GetXaxis()->GetBinLowEdge(1);
    double PravyOkrajPoslednihoBinu2 = hzY[0]->GetXaxis()->GetBinUpEdge(hzY[0]->GetNbinsX());

    TLine *line2 = new TLine(LevyOkrajPrvnihoBinu2, 1, PravyOkrajPoslednihoBinu2, 1);
    line2->SetLineStyle(2);
    line2->Draw("same");

    leg2->Draw("same");

    can->cd();
    tex->DrawLatex(0.35, 0.93, TString(var + " R_{cp}") + RcpTitles2[iCent] + TString("/") + RcpTitles2[jCent]);

    can->cd();
    //DrawTextAbove();
    can->SaveAs(outPdf);
    can->Clear();
    delete leg;
    delete leg1;
    delete leg2;
    delete tex;
    DeleteArray(hzX, nIter);
    DeleteArray(hzY, nIter);

    //TFile *fFinal = new TFile("./OutputPdf/FinalOutput.root", "UPDATE");
    //save RCP plot



}

void Centr(std::vector <std::vector<int>> centralityRange, const char *InputFile) {
    //Otevřu soubor OUTPUT_2014_20012024.root
    TFile *file = TFile::Open(InputFile);
    //check
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << InputFile << std::endl;
        exit(1);
    }
    //Načtu histogram hcentr
    TH1D *hcentr = (TH1D *) file->Get("event/hCentralityW");
    //set directory
    hcentr->SetDirectory(0);
    if (!hcentr) {
        std::cerr << "Error: Histogram hCentralityW not found in file: " << InputFile << std::endl;
        file->Close();
        exit(1);
    }
    //Vykreslím histogram, nad každým binem bude počet událostí
    hcentr->Draw();

    // Vytvoření textových objektů pro každý bin
    for (int i = 1; i <= hcentr->GetNbinsX(); ++i) {
        TLatex *latex = new TLatex();
        latex->SetTextSize(0.03);
        latex->SetTextAlign(22); // Centrování textu
        double binContent = hcentr->GetBinContent(i);
        double binCenter = hcentr->GetBinCenter(i);
        latex->DrawLatex(binCenter, binContent + 3 * ((hcentr->GetMaximum() - hcentr->GetMinimum()) / 100),
                         Form("%.0f", binContent));
    }
    int NCent = centralityRange.size();
    int NEvents[NCent];

    for (int iCent = 0; iCent < NCent; iCent++) {
        NEvents[iCent] = 0;
        for (int i = 1; i <= hcentr->GetNbinsX(); i++) {
            if (hcentr->GetBinCenter(i) >= centralityRange[iCent][0] &&
                hcentr->GetBinCenter(i) <= centralityRange[iCent][1]) {
                NEvents[iCent] += hcentr->GetBinContent(i);
            }
        }
        cout << "NEvents[" << iCent << "] = " << NEvents[iCent] << endl;
        NumberOfWEvents[iCent] = NEvents[iCent];

    }

    //
    //Vykreslím vertikální čáry podle centralityRange
    for (int i = 0; i < NCent - 1; i++) {

        TLine *line = new TLine(centralityRange[i][0] - 0.5, hcentr->GetMinimum() * 0.8, centralityRange[i][0] - 0.5,
                                hcentr->GetMaximum() * 1.2);
        line->SetLineColor(kRed);
        //line->Draw();
        /*   TLine *line2 = new TLine(centralityRange[i][1], 0, centralityRange[i][1], hcentr->GetMaximum());
           line2->SetLineColor(kRed);
           line2->Draw("same");*/
    }

    hcentr->GetYaxis()->SetRangeUser(hcentr->GetMinimum() * 0.8, hcentr->GetMaximum() * 1.2);
    for (int i = 0; i < NCent; i++) {
        TLatex *text = new TLatex();
        text->SetTextSize(0.03);
        text->SetTextColor(kRed);
        text->SetTextAlign(22);
        text->DrawLatex((centralityRange[i][0] + centralityRange[i][1]) / 2., hcentr->GetMinimum() * 1.2,
                        Form("%d-%d% %", centrality[i][0], centrality[i][1]));

    }
    file->Close();
}

double OutFlow(double Value, vector <Double_t> BinEdgesVector) {
    double ModValue = 0;

    int VecSize = BinEdgesVector.size();

    if (Value > BinEdgesVector[BinEdgesVector.size() - 2]) {
        ModValue = 0.5 * (BinEdgesVector[BinEdgesVector.size() - 1] + BinEdgesVector[BinEdgesVector.size() - 2]);
    } else if (Value < BinEdgesVector[1]) {
        ModValue = 0.5 * (BinEdgesVector[0] + BinEdgesVector[1]);
    } else {
        ModValue = Value;
    }

    return ModValue;
}

void plotResolution(TCanvas *can){

    can->cd();
    can->Clear();
    can->SetCanvasSize(1500, 900);
    gPad->SetLogy();

    TString pTTitle[5] = {
        "p_{T,Jet}^{True} #in (0-5) GeV/c",
        "p_{T,Jet}^{True} #in (5-10) GeV/c",
        "p_{T,Jet}^{True} #in (10-15) GeV/c",
        "p_{T,Jet}^{True} #in (15-20) GeV/c",
        "p_{T,Jet}^{True} #in (20-30) GeV/c"
    };

    TString varTitle[8] = {
        "p_{T,Jet}",
        "z",
        "#lambda_{1}^{1}",
        "#lambda_{1.5}^{1}",
        "#lambda_{2}^{1}",
        "#lambda_{3}^{1}",
        "#lambda_{0.5}^{1}",
        "p_{T}^{D}",
    };

    // loop over all variables jet p_T, z, lambda_..., momDispersion
    for (int iVar = 0; iVar < 8; iVar++) {

        // ---- Text objects (create once per iVar) ----
        TLatex *tex = new TLatex();
        tex->SetNDC();
        tex->SetTextFont(42);
        tex->SetTextSize(0.04);

        TLatex *tex2 = new TLatex();
        tex2->SetNDC();
        tex2->SetTextFont(42);
        tex2->SetTextSize(0.05);
        tex2->SetTextAngle(-90);

        // ---- 2D maps: mean and sigma (RMS) ----
        TH2D *hMeanMap  = new TH2D(Form("hMeanMap_var%d",  iVar),
                                  "Mean of resolution; p_{T,Jet}^{True} bin; Centrality",
                                  5, 0.5, 5.5,   // x bins = 1..5
                                  3, 0.5, 3.5);  // y bins = 1..3

        TH2D *hSigmaMap = new TH2D(Form("hSigmaMap_var%d", iVar),
                                  "Sigma (RMS) of resolution; p_{T,Jet}^{True} bin; Centrality",
                                  5, 0.5, 5.5,
                                  3, 0.5, 3.5);

        // axis bin labels
        for (int iCut = 0; iCut < 5; iCut++) {
            hMeanMap ->GetXaxis()->SetBinLabel(iCut+1, pTTitle[iCut]);
            hSigmaMap->GetXaxis()->SetBinLabel(iCut+1, pTTitle[iCut]);
        }
        for (int iCent = 0; iCent < 3; iCent++) {
            TString centLabel = Form("%d-%d%%", centrality[iCent][0], centrality[iCent][1]);
            hMeanMap ->GetYaxis()->SetBinLabel(iCent+1, centLabel);
            hSigmaMap->GetYaxis()->SetBinLabel(iCent+1, centLabel);
        }

        // ---- Page 1: your 5x3 panel of resolution histograms ----
        can->cd();
        can->Clear();
        gPad->SetLeftMargin(0.15);
        can->Divide(5, 3);

        for (int iCent = 0; iCent < 3; iCent++) {

            for (int iCut = 0; iCut < 5; iCut++) {

                can->cd(5*iCent + iCut + 1);

                TH1D *h = hResVar[iCent][iVar][iCut];
                if (!h) continue;
                // 1) vykresli axis + error bary (volitelně)
                h->SetMarkerStyle(20);
                h->SetMarkerSize(0.5);
                h->SetLineWidth(2);

                // 2) z histogramu udělej graph
                TGraphErrors *g = new TGraphErrors(h);
                g->SetMarkerStyle(20);
                g->SetMarkerSize(1.0);
                g->SetLineWidth(1);

                // 3) draw
                g->Draw("APLE");   // A=axis, P=points, L=line, E=errors

                // fill mean/sigma maps
                const double mean  = h->GetMean();
                const double sigma = h->GetRMS();

                hMeanMap ->SetBinContent(iCut+1, iCent+1, mean);
                hSigmaMap->SetBinContent(iCut+1, iCent+1, sigma);

                tex2->DrawLatex(0.96, 0.75, pTTitle[iCut]);
            }

            can->cd();
            tex->DrawLatex(0.42, 0.97-iCent*0.33,
                           Form("Centrality %d-%d%%", centrality[iCent][0], centrality[iCent][1]));
        }

        can->SaveAs(outPdf);

        // ---- Page 2: heatmaps (Mean + Sigma) ----
        can->cd();
        can->Clear();
        can->SetCanvasSize(2000, 1000);
        can->Divide(2,1);

        TLatex tex_title;
        tex_title.SetNDC();
        tex_title.SetTextSize(0.045);

        gStyle->SetOptStat(0);
        gStyle->SetPaintTextFormat("4.3f");
        gStyle->SetPalette(kBird);

        // Mean map
        can->cd(1);
        gPad->SetLeftMargin(0.18);
        gPad->SetRightMargin(0.15);
        gPad->SetBottomMargin(0.3);
        gPad->SetTopMargin(0.1);
        gPad->SetLogy(0);

        gStyle->SetTitleY(0.95);
        hMeanMap->SetTitle(Form("Mean resolution (var %d);; ", iVar));
        hMeanMap->GetXaxis()->LabelsOption("v"); // vertical x labels
        hMeanMap->Draw("COLZ TEXT");
        tex_title.DrawLatex(0.25,0.96,Form("Mean resolution (var %s)", varTitle[iVar].Data()));

        // Sigma map
        can->cd(2);
        gPad->SetLeftMargin(0.18);
        gPad->SetRightMargin(0.15);
        gPad->SetBottomMargin(0.3);
        gPad->SetTopMargin(0.1);
        gPad->SetLogy(0);

        gStyle->SetTitleY(0.95);
        hSigmaMap->SetTitle(Form("Sigma (RMS) resolution (var %d);; ", iVar));
        hSigmaMap->GetXaxis()->LabelsOption("v");
        hSigmaMap->Draw("COLZ TEXT");
        tex_title.DrawLatex(0.25,0.96,Form("Sigma (RMS) resolution (var %s)",  varTitle[iVar].Data()));


        can->SaveAs(outPdf);

        // ---- cleanup ----
        delete hMeanMap;
        delete hSigmaMap;
        delete tex;
        delete tex2;
    }

    // hotovo
}
inline bool Inside(double low, double up, double minV, double maxV, double eps=1e-12)
{
    return (low >= minV - eps) && (up <= maxV + eps);
}

inline double Center(double low, double up)
{
    return 0.5*(low + up);
}
void BuildResponse1DFromCache(int iCent, TH2D *hMatchPt, TH1D *hMissPt, TH1D *hFakePt)
{
    if (!hMatchPt || !hMissPt || !hFakePt) {
        cout << "Missing 1D cache histograms for cent " << iCent << endl;
        return;
    }

    const double recoMin = ptRecoBinsVec[iCent].front();
    const double recoMax = ptRecoBinsVec[iCent].back();
    const double trueMin = ptMcBinsVecCustom[iCent].front();
    const double trueMax = ptMcBinsVecCustom[iCent].back();

    // ---------- matched cache -> Fill / Miss / Fake ----------
    for (int ix = 1; ix <= hMatchPt->GetNbinsX(); ++ix) {
        const double recoLow = hMatchPt->GetXaxis()->GetBinLowEdge(ix);
        const double recoUp  = hMatchPt->GetXaxis()->GetBinUpEdge(ix);
        const bool recoIn = Inside(recoLow, recoUp, recoMin, recoMax);
        const double recoPt = Center(recoLow, recoUp);

        for (int iy = 1; iy <= hMatchPt->GetNbinsY(); ++iy) {
            const double w = hMatchPt->GetBinContent(ix, iy);
            if (w == 0.0) continue;

            const double trueLow = hMatchPt->GetYaxis()->GetBinLowEdge(iy);
            const double trueUp  = hMatchPt->GetYaxis()->GetBinUpEdge(iy);
            const bool trueIn = Inside(trueLow, trueUp, trueMin, trueMax);
            const double truePt = Center(trueLow, trueUp);

            if (recoIn && trueIn) {
                rurResponse[iCent].Fill(recoPt, truePt, w);
            } else if (!recoIn && trueIn && MissingJets) {
                rurResponse[iCent].Miss(truePt, w);
            } else if (recoIn && !trueIn && FakeJets) {
                rurResponse[iCent].Fake(recoPt, w);
            }
        }
    }

    // ---------- pure misses (true exists, reco not found / not in eta / ...) ----------
    for (int iy = 1; iy <= hMissPt->GetNbinsX(); ++iy) {
        const double w = hMissPt->GetBinContent(iy);
        if (w == 0.0) continue;

        const double trueLow = hMissPt->GetXaxis()->GetBinLowEdge(iy);
        const double trueUp  = hMissPt->GetXaxis()->GetBinUpEdge(iy);

        if (Inside(trueLow, trueUp, trueMin, trueMax) && MissingJets) {
            rurResponse[iCent].Miss(Center(trueLow, trueUp), w);
        }
    }

    // ---------- pure fakes (reco exists, true not found / not in eta / ...) ----------
    for (int ix = 1; ix <= hFakePt->GetNbinsX(); ++ix) {
        const double w = hFakePt->GetBinContent(ix);
        if (w == 0.0) continue;

        const double recoLow = hFakePt->GetXaxis()->GetBinLowEdge(ix);
        const double recoUp  = hFakePt->GetXaxis()->GetBinUpEdge(ix);

        if (Inside(recoLow, recoUp, recoMin, recoMax) && FakeJets) {
            rurResponse[iCent].Fake(Center(recoLow, recoUp), w);
        }
    }
}
void BuildResponse2DZFromCache(int iCent, THnSparseD *hMatchPtZ, TH2D *hMissPtZ, TH2D *hFakePtZ)
{
    if (!hMatchPtZ || !hMissPtZ || !hFakePtZ) {
        cout << "Missing 2D pT-z cache histograms for cent " << iCent << endl;
        return;
    }

    const double rPtMin = ptRecoBinsVec[iCent].front();
    const double rPtMax = ptRecoBinsVec[iCent].back();
    const double rZMin  = zRecoBinsVec[iCent].front();
    const double rZMax  = zRecoBinsVec[iCent].back();

    const double tPtMin = ptMcBinsVecCustom[iCent].front();
    const double tPtMax = ptMcBinsVecCustom[iCent].back();
    const double tZMin  = zMcBinsVecCustom[iCent].front();
    const double tZMax  = zMcBinsVecCustom[iCent].back();

    TAxis* axRPT = hMatchPtZ->GetAxis(0);
    TAxis* axRZ  = hMatchPtZ->GetAxis(1);
    TAxis* axTPT = hMatchPtZ->GetAxis(2);
    TAxis* axTZ  = hMatchPtZ->GetAxis(3);

    Int_t coord[4];
    const Long64_t n = hMatchPtZ->GetNbins();

    // ---------- matched cache -> Fill / Miss / Fake ----------
    for (Long64_t i = 0; i < n; ++i) {
        const double w = hMatchPtZ->GetBinContent(i, coord);
        if (w == 0.0) continue;

        const double rPtLow = axRPT->GetBinLowEdge(coord[0]);
        const double rPtUp  = axRPT->GetBinUpEdge (coord[0]);
        const double rZLow  = axRZ ->GetBinLowEdge(coord[1]);
        const double rZUp   = axRZ ->GetBinUpEdge (coord[1]);

        const double tPtLow = axTPT->GetBinLowEdge(coord[2]);
        const double tPtUp  = axTPT->GetBinUpEdge (coord[2]);
        const double tZLow  = axTZ ->GetBinLowEdge(coord[3]);
        const double tZUp   = axTZ ->GetBinUpEdge (coord[3]);

        const bool recoIn =
            Inside(rPtLow, rPtUp, rPtMin, rPtMax) &&
            Inside(rZLow,  rZUp,  rZMin,  rZMax);

        const bool trueIn =
            Inside(tPtLow, tPtUp, tPtMin, tPtMax) &&
            Inside(tZLow,  tZUp,  tZMin,  tZMax);

        const double rPt = Center(rPtLow, rPtUp);
        const double rZ  = Center(rZLow,  rZUp);
        const double tPt = Center(tPtLow, tPtUp);
        const double tZ  = Center(tZLow,  tZUp);

        if (recoIn && trueIn) {
            rurResponse2D[iCent][0].Fill(rPt, rZ, tPt, tZ, w);
        } else if (!recoIn && trueIn && MissingJets) {
            rurResponse2D[iCent][0].Miss(tPt, tZ, w);
        } else if (recoIn && !trueIn && FakeJets) {
            rurResponse2D[iCent][0].Fake(rPt, rZ, w);
        }
    }

    // ---------- pure misses ----------
    for (int ix = 1; ix <= hMissPtZ->GetNbinsX(); ++ix) {
        const double tPtLow = hMissPtZ->GetXaxis()->GetBinLowEdge(ix);
        const double tPtUp  = hMissPtZ->GetXaxis()->GetBinUpEdge (ix);
        const bool tPtIn = Inside(tPtLow, tPtUp, tPtMin, tPtMax);
        const double tPt = Center(tPtLow, tPtUp);

        for (int iy = 1; iy <= hMissPtZ->GetNbinsY(); ++iy) {
            const double w = hMissPtZ->GetBinContent(ix, iy);
            if (w == 0.0) continue;

            const double tZLow = hMissPtZ->GetYaxis()->GetBinLowEdge(iy);
            const double tZUp  = hMissPtZ->GetYaxis()->GetBinUpEdge (iy);

            if (tPtIn && Inside(tZLow, tZUp, tZMin, tZMax) && MissingJets) {
                rurResponse2D[iCent][0].Miss(tPt, Center(tZLow, tZUp), w);
            }
        }
    }

    // ---------- pure fakes ----------
    for (int ix = 1; ix <= hFakePtZ->GetNbinsX(); ++ix) {
        const double rPtLow = hFakePtZ->GetXaxis()->GetBinLowEdge(ix);
        const double rPtUp  = hFakePtZ->GetXaxis()->GetBinUpEdge (ix);
        const bool rPtIn = Inside(rPtLow, rPtUp, rPtMin, rPtMax);
        const double rPt = Center(rPtLow, rPtUp);

        for (int iy = 1; iy <= hFakePtZ->GetNbinsY(); ++iy) {
            const double w = hFakePtZ->GetBinContent(ix, iy);
            if (w == 0.0) continue;

            const double rZLow = hFakePtZ->GetYaxis()->GetBinLowEdge(iy);
            const double rZUp  = hFakePtZ->GetYaxis()->GetBinUpEdge (iy);

            if (rPtIn && Inside(rZLow, rZUp, rZMin, rZMax) && FakeJets) {
                rurResponse2D[iCent][0].Fake(rPt, Center(rZLow, rZUp), w);
            }
        }
    }
}
void BuildResponse2DAngFromCache(int iCent, int iAng,
                                 THnSparseD *hMatchPtAng,
                                 TH2D *hMissPtAng,
                                 TH2D *hFakePtAng)
{
    if (!hMatchPtAng || !hMissPtAng || !hFakePtAng) {
        cout << "Missing 2D pT-lambda cache histograms for cent "
             << iCent << " ang " << iAng << endl;
        return;
    }

    const double rPtMin  = ptRecoBinsVec[iCent].front();
    const double rPtMax  = ptRecoBinsVec[iCent].back();
    const double rAMin   = angRecoBinsVec[iCent][iAng].front();
    const double rAMax   = angRecoBinsVec[iCent][iAng].back();

    const double tPtMin  = ptMcBinsVecCustom[iCent].front();
    const double tPtMax  = ptMcBinsVecCustom[iCent].back();
    const double tAMin   = angMcBinsVecCustom[iCent][iAng].front();
    const double tAMax   = angMcBinsVecCustom[iCent][iAng].back();

    TAxis* axRPT = hMatchPtAng->GetAxis(0);
    TAxis* axRA  = hMatchPtAng->GetAxis(1);
    TAxis* axTPT = hMatchPtAng->GetAxis(2);
    TAxis* axTA  = hMatchPtAng->GetAxis(3);

    Int_t coord[4];
    const Long64_t n = hMatchPtAng->GetNbins();

    // ---------- matched cache -> Fill / Miss / Fake ----------
    for (Long64_t i = 0; i < n; ++i) {
        const double w = hMatchPtAng->GetBinContent(i, coord);
        if (w == 0.0) continue;

        const double rPtLow = axRPT->GetBinLowEdge(coord[0]);
        const double rPtUp  = axRPT->GetBinUpEdge (coord[0]);
        const double rALow  = axRA ->GetBinLowEdge(coord[1]);
        const double rAUp   = axRA ->GetBinUpEdge (coord[1]);

        const double tPtLow = axTPT->GetBinLowEdge(coord[2]);
        const double tPtUp  = axTPT->GetBinUpEdge (coord[2]);
        const double tALow  = axTA ->GetBinLowEdge(coord[3]);
        const double tAUp   = axTA ->GetBinUpEdge (coord[3]);

        const bool recoIn =
            Inside(rPtLow, rPtUp, rPtMin, rPtMax) &&
            Inside(rALow,  rAUp,  rAMin,  rAMax);

        const bool trueIn =
            Inside(tPtLow, tPtUp, tPtMin, tPtMax) &&
            Inside(tALow,  tAUp,  tAMin,  tAMax);

        const double rPt = Center(rPtLow, rPtUp);
        const double rA  = Center(rALow,  rAUp);
        const double tPt = Center(tPtLow, tPtUp);
        const double tA  = Center(tALow,  tAUp);

        if (recoIn && trueIn) {
            rurResponse2D[iCent][iAng + 1].Fill(rPt, rA, tPt, tA, w);
        } else if (!recoIn && trueIn && MissingJets) {
            rurResponse2D[iCent][iAng + 1].Miss(tPt, tA, w);
        } else if (recoIn && !trueIn && FakeJets) {
            rurResponse2D[iCent][iAng + 1].Fake(rPt, rA, w);
        }
    }

    // ---------- pure misses ----------
    for (int ix = 1; ix <= hMissPtAng->GetNbinsX(); ++ix) {
        const double tPtLow = hMissPtAng->GetXaxis()->GetBinLowEdge(ix);
        const double tPtUp  = hMissPtAng->GetXaxis()->GetBinUpEdge (ix);
        const bool tPtIn = Inside(tPtLow, tPtUp, tPtMin, tPtMax);
        const double tPt = Center(tPtLow, tPtUp);

        for (int iy = 1; iy <= hMissPtAng->GetNbinsY(); ++iy) {
            const double w = hMissPtAng->GetBinContent(ix, iy);
            if (w == 0.0) continue;

            const double tALow = hMissPtAng->GetYaxis()->GetBinLowEdge(iy);
            const double tAUp  = hMissPtAng->GetYaxis()->GetBinUpEdge (iy);

            if (tPtIn && Inside(tALow, tAUp, tAMin, tAMax) && MissingJets) {
                rurResponse2D[iCent][iAng + 1].Miss(tPt, Center(tALow, tAUp), w);
            }
        }
    }

    // ---------- pure fakes ----------
    for (int ix = 1; ix <= hFakePtAng->GetNbinsX(); ++ix) {
        const double rPtLow = hFakePtAng->GetXaxis()->GetBinLowEdge(ix);
        const double rPtUp  = hFakePtAng->GetXaxis()->GetBinUpEdge (ix);
        const bool rPtIn = Inside(rPtLow, rPtUp, rPtMin, rPtMax);
        const double rPt = Center(rPtLow, rPtUp);

        for (int iy = 1; iy <= hFakePtAng->GetNbinsY(); ++iy) {
            const double w = hFakePtAng->GetBinContent(ix, iy);
            if (w == 0.0) continue;

            const double rALow = hFakePtAng->GetYaxis()->GetBinLowEdge(iy);
            const double rAUp  = hFakePtAng->GetYaxis()->GetBinUpEdge (iy);

            if (rPtIn && Inside(rALow, rAUp, rAMin, rAMax) && FakeJets) {
                rurResponse2D[iCent][iAng + 1].Fake(rPt, Center(rALow, rAUp), w);
            }
        }
    }
}
void LoadDataCache()
{
    TFile *fCache = TFile::Open(CacheRMFileName, "READ");
    if (!fCache || fCache->IsZombie()) {
        cout << "Cannot open cache file: " << CacheRMFileName << endl;
        return;
    }

    CheckAllCacheCompatibility();

    for (int iCent = 0; iCent < nCentralityBins; iCent++) {

        // ---------- 1D ----------
        TH2D *hMatchPt = (TH2D*) fCache->Get(Form("hCacheMatchPt_cent%d", iCent));
        TH1D *hMissPt  = (TH1D*) fCache->Get(Form("hCacheMissPt_cent%d", iCent));
        TH1D *hFakePt  = (TH1D*) fCache->Get(Form("hCacheFakePt_cent%d", iCent));

        BuildResponse1DFromCache(iCent, hMatchPt, hMissPt, hFakePt);

        // ---------- 2D pT-z ----------
        THnSparseD *hMatchPtZ = (THnSparseD*) fCache->Get(Form("hCacheMatchPtZ_cent%d", iCent));
        TH2D *hMissPtZ        = (TH2D*) fCache->Get(Form("hCacheMissPtZ_cent%d", iCent));
        TH2D *hFakePtZ        = (TH2D*) fCache->Get(Form("hCacheFakePtZ_cent%d", iCent));

        BuildResponse2DZFromCache(iCent, hMatchPtZ, hMissPtZ, hFakePtZ);

        // ---------- 2D pT-lambda ----------
        for (int iAng = 0; iAng < nAngularities; iAng++) {
            THnSparseD *hMatchPtAng = (THnSparseD*) fCache->Get(Form("hCacheMatchPtAng_cent%d_ang%d", iCent, iAng));
            TH2D *hMissPtAng        = (TH2D*) fCache->Get(Form("hCacheMissPtAng_cent%d_ang%d", iCent, iAng));
            TH2D *hFakePtAng        = (TH2D*) fCache->Get(Form("hCacheFakePtAng_cent%d_ang%d", iCent, iAng));

            BuildResponse2DAngFromCache(iCent, iAng, hMatchPtAng, hMissPtAng, hFakePtAng);
        }
    }

    fCache->Close();
}
void LoadDataMC() {

    TRandom3 randSplit(12345);
    TFile *treeFile;
    treeFile = new TFile(McJetsFileData, "READ");
    if (!treeFile || treeFile->IsZombie()) return;

    TTree *jetTree = (TTree *) treeFile->Get("jets");

    TString preFix = (Method == "ICS" ? "ICS_" : "");

    // 3) TTreeCache (u vzdálených/velkých souborů často gamechanger)
    //jetTree->SetBranchStatus("*", 0);
    jetTree->SetCacheSize(128*1024*1024);          // klidně 256 MB
    jetTree->AddBranchToCache("*", kTRUE);
    jetTree->SetCacheLearnEntries(1000);


    // Event
    jetTree->SetBranchStatus("centrality", 1);
    jetTree->SetBranchStatus("centralityAlt", 1);
    jetTree->SetBranchStatus("weightCentrality", 1);
    jetTree->SetBranchStatus("eventMaxPtTrack", 1);


    // MC
    jetTree->SetBranchStatus("mcJetPt", 1);
    jetTree->SetBranchStatus("mcJetEta", 1);
    jetTree->SetBranchStatus("mcJetNConst", 1);
    jetTree->SetBranchStatus("mcJetD0Z", 1);
    jetTree->SetBranchStatus("mcD0Pt", 1);
    jetTree->SetBranchStatus("mcD0Eta", 1);
    jetTree->SetBranchStatus("gRefMult", 1);

    // MC lambdas
    jetTree->SetBranchStatus("mcJetLambda1_1", 1);
    jetTree->SetBranchStatus("mcJetLambda1_1_5", 1);
    jetTree->SetBranchStatus("mcJetLambda1_2", 1);
    jetTree->SetBranchStatus("mcJetLambda1_3", 1);
    jetTree->SetBranchStatus("mcJetLambda1_0_5", 1);
    jetTree->SetBranchStatus("mcJetMomDisp", 1);

    // RECO
    jetTree->SetBranchStatus(preFix + "recoJetPt" + (preFix == "ICS_" ? "" : "Corr"), 1); // nebo recoJetPt podle prefixu
    jetTree->SetBranchStatus(preFix + "recoJetEta", 1);
    jetTree->SetBranchStatus(preFix + "recoJetNConst", 1);
    jetTree->SetBranchStatus(preFix + "recoJetD0Z", 1);

    jetTree->SetBranchStatus("mcSmearedD0Pt", 1);
    jetTree->SetBranchStatus("mcSmearedD0Eta", 1);

    jetTree->SetBranchStatus(preFix + "recoJetArea", 1);
    jetTree->SetBranchStatus(preFix + "recoJetRho", 1);

    // RECO lambdas
    jetTree->SetBranchStatus(preFix + "recoJetLambda1_1", 1);
    jetTree->SetBranchStatus(preFix + "recoJetLambda1_1_5", 1);
    jetTree->SetBranchStatus(preFix + "recoJetLambda1_2", 1);
    jetTree->SetBranchStatus(preFix + "recoJetLambda1_3", 1);
    jetTree->SetBranchStatus(preFix + "recoJetLambda1_0_5", 1);
    jetTree->SetBranchStatus(preFix + "recoJetMomDisp", 1);

    StJetTreeStruct mcJet, recoJet;
    assignTree(jetTree, mcJet, recoJet);
    Double_t nEntries = jetTree->GetEntries();

 


    nEntries /= DividedMcDataBy;
    cout << "nEntries = " << (Float_t) nEntries / 1000. << "k" << endl
         << endl;


    TH1D *hCentr = new TH1D("hCentr", "hCentr", 9, -0.5, 8.5);

    TH2D* hMigRef[3];
    TH2D* hMigWgt[3];
    TH2D* hMigRefZ[3];
    TH2D* hMigWgtZ[3];
    for (int c = 0; c < 3; ++c) {
        hMigRef[c] = (TH2D*) hRespZ[0][c]->Clone(Form("hMigRef_c%i", c));
        hMigWgt[c] = (TH2D*) hRespZ[0][c]->Clone(Form("hMigWgt_c%i", c));
        hMigRef[c]->Reset("ICES");  // vynulovat obsah + sumw2
        hMigWgt[c]->Reset("ICES");
        hMigRefZ[c] = (TH2D*) hRespZ[1][c]->Clone(Form("hMigRefZ_c%i", c));
        hMigWgtZ[c] = (TH2D*) hRespZ[1][c]->Clone(Form("hMigWgtZ_c%i", c));
        hMigRefZ[c]->Reset("ICES");  // vynulovat obsah + sumw2
        hMigWgtZ[c]->Reset("ICES");
    }


    //new file for ekuivalent binning
    TFile *fEkvi = new TFile("Tests/EkviBinning.root", "RECREATE");
    //New tree for ekuivalent binning
    TTree *tEkvi = new TTree("tEkvi", "tEkvi");
    //větev ptReco 
    Float_t ptRecoEkvi;
    tEkvi->Branch("ptRecoEkvi", &ptRecoEkvi, "ptRecoEkvi/F");
    //větev zReco
    Float_t zRecoEkvi;
    tEkvi->Branch("zRecoEkvi", &zRecoEkvi, "zRecoEkvi/F");
    //větev váha
    Float_t weightEkvi;
    tEkvi->Branch("weightEkvi", &weightEkvi, "weightEkvi/F");
    //větev centBin
    Int_t centBinEkvi;
    tEkvi->Branch("centBinEkvi", &centBinEkvi, "centBinEkvi/I");
    Double_t lam[6];
    tEkvi->Branch("lam0", &lam[0], "lam0/D");
    tEkvi->Branch("lam1", &lam[1], "lam1/D");
    tEkvi->Branch("lam2", &lam[2], "lam2/D");
    tEkvi->Branch("lam3", &lam[3], "lam3/D");
    tEkvi->Branch("lam4", &lam[4], "lam4/D");
    tEkvi->Branch("lam5", &lam[5], "lam5/D");

    
    for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {

        Float_t progress = 0.;
        progress = (Float_t) iEntry / (1. * nEntries);

        if (iEntry % 100000 == 0) cout << "Training: \r (" << (progress * 100.0) << "%)" << std::flush;

        jetTree->GetEntry(iEntry);
        Int_t centBin = getCentralityBin99(recoJet.centrality);


        /**************************************************** */
        //Event
        double _centrality = recoJet.centrality;
        double _centralityAlt = mcJet.centralityAlt;
        double _weightCentrality = recoJet.weight;
        double _eventMaxPtTrack = recoJet.eventMaxPtTrack;

        if (_eventMaxPtTrack > 30) continue;

        //D0 meson
        double _mcD0Pt = mcJet.d0pt;
        double eta = mcJet.d0eta;

        double _mcSmearedD0Pt = recoJet.RecoD0Pt;
        double _mcSmearedD0Eta = recoJet.RecoD0Eta;

        double mD0 = 1.864; // GeV

        double pz = _mcD0Pt* sinh(eta);
        double p  = _mcD0Pt* cosh(eta);
        double E  = sqrt(p*p + mD0*mD0);

        double y = 0.5 * log((E + pz) / (E - pz));

        if (abs(y) > 0.6) continue; //D0 meson rapidity cut

        //double _mcSmearedD0Pt = recoJet.RecoD0Pt;
        //double _mcSmearedD0Eta = recoJet.RecoD0Eta;


        //D0-jet
        double _mcJetEta = mcJet.mcJetEta;
        double _recoJetEta = recoJet.recoJetEta;
        double _mcJetPt = mcJet.jetpt;
        double _recoJetPt = recoJet.jetpt;
        double _mcNConst = mcJet.numberofconstituents;

        double _mcD0Z = mcJet.d0z;
        double _recoD0Z = recoJet.d0z;
        double _mcLambda[6] =   {mcJet.lambda[0], mcJet.lambda[1], mcJet.lambda[2],
                                mcJet.lambda[3], mcJet.lambda[4], mcJet.lambda[5]};
        double _recoLambda[6] = {recoJet.lambda[0], recoJet.lambda[1], recoJet.lambda[2],
                                recoJet.lambda[3], recoJet.lambda[4], recoJet.lambda[5]};
        double _recoNConst = recoJet.numberofconstituents;                        
        /**************************************************** */


        if (centBin < 0) continue;
        if (_mcD0Pt < minPtD0Cut) continue;
        if (_mcD0Pt >= maxPtD0Cut) continue;
        if (_mcSmearedD0Pt < minPtD0Cut) continue;
        if (_mcSmearedD0Pt >= maxPtD0Cut) continue;
        //if((mcJet.gRefMultMc*0.1150315604 -3.0) >  recoJet.recoJetRho) continue;
        //if((mcJet.gRefMultMc*0.1742999236 -6) <  recoJet.recoJetRho) continue;

        if (abs(_mcJetEta) > 0.6 && abs(_recoJetEta) > 0.6) continue; //not found but not relevant

        if (deleteOneConstituentJets && (_mcNConst == 1 && _recoNConst == 1)) continue;

        if (_mcJetPt > 30) continue;

        bool goesToTrain = true;
        bool goesToTest  = false;

        if (ClosureTest) {
            if (UseTheSameSample) {
                goesToTrain = true;
                goesToTest  = true;
            } else {
                goesToTrain = (randSplit.Rndm() < TrainToTestRatio);
                goesToTest  = !goesToTrain;
            }
        }

        double recoJetPtMin = ptRecoBinsVec[centBin][0];
        double recoJetPtMax = ptRecoBinsVec[centBin][ptRecoBinsVec[centBin].size() - 1];
        double trueJetPtMin = ptMcBinsVecCustom[centBin][0];
        double trueJetPtMax = ptMcBinsVecCustom[centBin][ptMcBinsVecCustom[centBin].size() - 1];
        double recoJetZMin = zRecoBinsVec[centBin][0];
        double recoJetZMax = zRecoBinsVec[centBin][zRecoBinsVec[centBin].size() - 1];
        double trueJetZMin = zMcBinsVecCustom[centBin][0];
        double trueJetZMax = zMcBinsVecCustom[centBin][zMcBinsVecCustom[centBin].size() - 1];
        double recoJetLambdaMin[6];
        double recoJetLambdaMax[6];
        double trueJetLambdaMin[6];
        double trueJetLambdaMax[6];


        for (int i = 0; i < 6; i++) {
            recoJetLambdaMin[i] = angRecoBinsVec[centBin][i][0];
            recoJetLambdaMax[i] = angRecoBinsVec[centBin][i][angRecoBinsVec[centBin][i].size() - 1];
            trueJetLambdaMin[i] = angMcBinsVecCustom[centBin][i][0];
            trueJetLambdaMax[i] = angMcBinsVecCustom[centBin][i][angMcBinsVecCustom[centBin][i].size() - 1];
        }

        //Jet pT
        Bool_t isRecoJetFound = _recoNConst > 0;
        Bool_t isRecoJetMiss = !isRecoJetFound;
        Bool_t isFakeJetMiss = false; //It cannot happen for D0 jets
        Bool_t isTrueJetFound = true; //It is always true for D0 jets

        //Eta
        Bool_t isTrueJetInEtaFound = abs(_mcJetEta) < 0.6;
        Bool_t isRecoJetInEtaFound =  abs(_recoJetEta) < 0.6;
        Bool_t isMissJetInEta = isTrueJetInEtaFound && !isRecoJetInEtaFound;
        Bool_t isFakeJetInEta = isRecoJetInEtaFound && !isTrueJetInEtaFound;
        Bool_t isMatchedInEta = isTrueJetInEtaFound && isRecoJetInEtaFound;

        //NConst
        Bool_t isTrueNConst = true;
        Bool_t isRecoNConst = true;
        Bool_t isMissNConst = false;
        Bool_t isFakeNConst = false;
        Bool_t isMatchedNConst = true;

        if(deleteOneConstituentJets){
            isTrueNConst = _mcNConst > 1;
            isRecoNConst = _recoNConst > 1;
            isMissNConst = isTrueNConst && !isRecoNConst;
            isFakeNConst = isRecoNConst && !isTrueNConst;
            isMatchedNConst = isTrueNConst && isRecoNConst;
        }





        //JetPt
        Bool_t isRecoInPt = _recoJetPt >= recoJetPtMin && _recoJetPt < recoJetPtMax;
        Bool_t isTrueInPt = _mcJetPt >= trueJetPtMin && _mcJetPt < trueJetPtMax;
        Bool_t isMissInPt = isTrueInPt && !isRecoInPt;
        Bool_t isFakeInPt = isRecoInPt && !isTrueInPt;
        Bool_t isMatchedInPt = isRecoInPt && isTrueInPt;

        //1D unfolding (pT summary)
        Bool_t isMatchedPt = isRecoJetFound && isMatchedInEta && isMatchedInPt && isMatchedNConst;
        Bool_t isTruePt = isTrueJetInEtaFound && isTrueNConst && isTrueInPt;
        Bool_t isRecoPt = isRecoJetFound && isRecoJetInEtaFound && isRecoNConst && isRecoInPt;
        Bool_t isMissPt = isTruePt && (isRecoJetMiss || isMissJetInEta || isMissNConst || isMissInPt);
        Bool_t isFakePt = isRecoPt && (isFakeJetMiss || isFakeJetInEta || isFakeNConst || isFakeInPt);

        bool isMatchedPtCache = isRecoJetFound && isMatchedInEta;

        bool isMissPtCache =
            isTrueJetFound &&
            isTrueJetInEtaFound &&
            !(isRecoJetFound && isRecoJetInEtaFound);

        bool isFakePtCache =
            isRecoJetFound &&
            isRecoJetInEtaFound &&
            !(isTrueJetFound && isTrueJetInEtaFound);

        if (int(isMatchedPt) + int(isMissPt) + int(isFakePt) > 1) {
            cerr << "Mistake in logic (1D)" << endl;
            exit(1);
        }

        //if(!isMatchedPt && !isMissPt && !isFakePt) continue; //Do not care about them


        //Z value
        Bool_t isRecoInZ = _recoD0Z >= recoJetZMin && _recoD0Z < recoJetZMax;
        Bool_t isTrueInZ = _mcD0Z >= trueJetZMin && _mcD0Z < trueJetZMax;
        Bool_t isMissInZ =  isTrueInZ && !isRecoInZ;
        Bool_t isFakeInZ =  isRecoInZ && !isTrueInZ;
        Bool_t isMatchedInZ = isTrueInZ && isRecoInZ;

        //2D unfolding (pT, z) summary
        Bool_t isMatchedPtZ = isMatchedPt && isMatchedInZ;
        Bool_t isTruePtZ = isTruePt && isTrueInZ;
        Bool_t isRecoPtZ = isRecoPt && isRecoInZ;
        Bool_t isMissPtZ = isTruePtZ && (isMissPt || isMissInZ);
        Bool_t isFakePtZ = isRecoPtZ && (isFakePt || isFakeInZ);

        bool isMatchedPtZCache = isMatchedPtCache;

        bool isMissPtZCache = isMissPtCache;

        bool isFakePtZCache = isFakePtCache;

        if (int(isMatchedPtZ) + int(isMissPtZ) + int(isFakePtZ) > 1) {
            cerr << "Mistake in logic (2D)" << endl;
            exit(1);
        }

        //Angularity
        Bool_t isRecoInLambda[6];
        Bool_t isTrueInLambda[6];
        Bool_t isMissInLambda[6];
        Bool_t isFakeInLambda[6];
        Bool_t isMatchedInLambda[6];

        Bool_t isMatchedPtLambda[6];
        Bool_t isTruePtLambda[6];
        Bool_t isRecoPtLambda[6];
        Bool_t isMissPtLambda[6];
        Bool_t isFakePtLambda[6];

        bool isMatchedPtLambdaCache[6];

        bool isMissPtLambdaCache[6];
        bool isFakePtLambdaCache[6];

        for (int i = 0; i < 6; i++) {
            
            isMatchedPtLambdaCache[i] = isMatchedPtCache;

            isMissPtLambdaCache[i] = isMissPtCache;

            isFakePtLambdaCache[i] = isFakePtCache;

            //Lambda
            isRecoInLambda[i] = _recoLambda[i] >= recoJetLambdaMin[i] && _recoLambda[i] < recoJetLambdaMax[i];
            isTrueInLambda[i] = _mcLambda[i] >= trueJetLambdaMin[i] && _mcLambda[i] < trueJetLambdaMax[i];
            isMissInLambda[i] = isTrueInLambda[i] && !isRecoInLambda[i];
            isFakeInLambda[i] = isRecoInLambda[i] && !isTrueInLambda[i];
            isMatchedInLambda[i] = isRecoInLambda[i] && isTrueInLambda[i];

            //2D unfolding (pT, lam) summary
            isMatchedPtLambda[i] = isMatchedPt && isMatchedInLambda[i];
            isTruePtLambda[i] = isTruePt && isTrueInLambda[i];
            isRecoPtLambda[i] = isRecoPt && isRecoInLambda[i];
            isMissPtLambda[i] = isTruePtLambda[i] && (isMissPt || isMissInLambda[i]);
            isFakePtLambda[i] = isRecoPtLambda[i] && (isFakePt || isFakeInLambda[i]);

            if (int(isMatchedPtLambda[i]) + int(isMissPtLambda[i]) + int(isFakePtLambda[i]) > 1) {
                cerr << "Mistake in logic (pT, lambda), i = " << i << endl;
                exit(1);
            }

        }

        double cweight = 1;

        cweight = 1.0 * _weightCentrality //centr. weight
                  * nCentralityNumbers->GetBinContent(
                nCentralityNumbers->FindBin(_centrality)) //Number of events in real event
                  / nCentralityNumbersMC->GetBinContent(
                nCentralityNumbersMC->FindBin(_centrality)); //Number of events in MC event (Normalization)

        double weightPt = 1.;


        if (FONLLjet)
            weightPt = 1. / McPTRawD0Jet[centBin]->GetBinContent(McPTRawD0Jet[centBin]->FindBin(_mcJetPt)) *
                       pureFonll->Eval(_mcJetPt) *100000;        //pure FONLL //_mcJetPt
        else
            weightPt = 1. / McPTRawD0[centBin]->GetBinContent(McPTRawD0[centBin]->FindBin(_mcD0Pt)) *
            pureFonllD0meson->Eval(_mcD0Pt) * 100000 ;        //pure FONLL

        double weight = 1;
        weight = weightPt * cweight;// * vaha;

        //skip if inf
        if (!isfinite(weight)) {
            cerr << "error, wrong weight" << endl;
            exit(1);
        }
        
        if (weight < 0) cerr << weight << endl;

        hCentr->Fill(_centrality, weight);


        double flatweight = 1. / (McPTRawD0[centBin]->GetBinContent(McPTRawD0[centBin]->FindBin(recoJet.RecoD0Pt))) *
                            _weightCentrality;
        gRefMultMcCorr->Fill(mcJet.gRefMultMc, weight);

        //if (iEntry < (ClosureTest ? (TrainToTestRatio * nEntries) : 1.0 * nEntries)) { //Check later!!!


        double pSWeight[8];
        

        
        pSWeight[0] = getPriorShapeWeight(_usePriorShapeWeighting, _mcJetPt, 0); //jetPt
        pSWeight[1] = getPriorShapeWeight(_usePriorShapeWeighting, _mcD0Z, 1); //z
        for (int i = 0; i < 6; i++) {
            pSWeight[2 + i] = getPriorShapeWeight(_usePriorShapeWeighting, _mcLambda[i], 2 + i); //lambdas
        }
/*
        cout << "pSWeight: ";
        for (int i = 0; i < 8; i++) {
            cout << pSWeight[i] << " ";
        }
        cout << endl;
        */
        

        if (FillStandardRM) {
            if (!ClosureTest || goesToTrain) {

                D0MesonPtMcReco[centBin]->Fill(_mcSmearedD0Pt, weight);
                D0MesonPtMcTrue[centBin]->Fill(_mcD0Pt, weight);
                D0JetPtMcReco[centBin]->Fill(_recoJetPt, weight);
                D0JetPtMcTrue[centBin]->Fill(_mcJetPt, weight);

                if (isMatchedPt) {

                    int pTJetBinCut = _mcJetPt < 5 ? 0 : (_mcJetPt < 10 ? 1 : (_mcJetPt < 15 ? 2 : (_mcJetPt < 20 ? 3 : 4)));

                    hResVar[centBin][0][pTJetBinCut]->Fill((_mcJetPt - _recoJetPt)/_mcJetPt, weight);

                    double x[4] = { _recoJetPt, _recoD0Z, _mcJetPt, _mcD0Z};

                    //response check
                    hMigRef[centBin]->Fill(_recoJetPt, _mcJetPt, cweight*pSWeight[0]);
                    hMigWgt[centBin]->Fill(_recoJetPt, _mcJetPt, weight*pSWeight[0]);
                    hMigRefZ[centBin]->Fill(_recoD0Z, _mcD0Z, cweight*pSWeight[1]*pSWeight[0]);
                    hMigWgtZ[centBin]->Fill(_recoD0Z, _mcD0Z, weight*pSWeight[1]*pSWeight[0]);
                    rurResponse[centBin].Fill(_recoJetPt, _mcJetPt, weight*pSWeight[0]);

                    jetPtCheck[centBin]->Fill(_mcJetPt);
                    jetPtCheckScaled2[centBin]->Fill(_mcJetPt, weightPt);
                    jetPtRecoCheckScaled[centBin]->Fill(_recoJetPt, weightPt);
                    jetPtRecoCheckScaled2[centBin]->Fill(_recoJetPt, weightPt);
                    jetPtRecoCheck[centBin]->Fill(_recoJetPt);
                    hRespZ[0][centBin]->Fill(_recoJetPt, _mcJetPt, weight*pSWeight[1]*pSWeight[0]);

                    hRespZHighRes[0][centBin]->Fill(_recoJetPt, _mcJetPt);

                    ptRecoEkvi = _recoJetPt;
                    zRecoEkvi = _recoD0Z;
                    weightEkvi = weight;
                    centBinEkvi = centBin;
                    for (int i = 0; i < 6; i++) {
                        lam[i] = _recoLambda[i];
                    }
                    tEkvi->Fill();



                } else if (isMissPt && MissingJets) {

                    rurResponse[centBin].Miss(_mcJetPt, weight*pSWeight[0]);

                } else if (isFakePt && FakeJets) {

                    rurResponse[centBin].Fake(_recoJetPt, weight*pSWeight[0]);

                }

                if (Unfold2D) {

                    int pTJetBinCut = _mcJetPt < 5 ? 0 : (_mcJetPt < 10 ? 1 : (_mcJetPt < 15 ? 2 : (_mcJetPt < 20 ? 3 : 4)));





                    if (isMatchedPtZ) {
                        hRespZ[1][centBin]->Fill(_recoD0Z, _mcD0Z, weight*pSWeight[1]*pSWeight[0]);
                        
                        Double_t coord4D[4] = {_recoJetPt, _recoD0Z, _mcJetPt, _mcD0Z};
			hResponseFine4D[centBin][0]->Fill(coord4D, weight * pSWeight[1] * pSWeight[0]);
                        //hResponseFine4D[centBin][0]->Fill(_recoJetPt, _recoD0Z, _mcJetPt, _mcD0Z, weight*pSWeight[1]*pSWeight[0]);
                        if (_mcJetPt > 5 && _mcJetPt < 20) jetZ[centBin]->Fill(_mcD0Z);
                        hResVar[centBin][1][pTJetBinCut]->Fill((_mcD0Z - _recoD0Z)/_mcD0Z, weight*pSWeight[1]*pSWeight[0]);
                        hRespZHighRes[1][centBin]->Fill(_recoD0Z, _mcD0Z);

                        rurResponse2D[centBin][0].Fill(_recoJetPt, _recoD0Z, _mcJetPt, _mcD0Z, weight*pSWeight[1]*pSWeight[0]);

                        rurResponse2DTest[centBin].Fill(_recoJetPt, _recoD0Z, _mcJetPt, _mcD0Z, weight*pSWeight[1]*pSWeight[0]);
                        rurResponse2DTestW[centBin].Fill(_recoJetPt, _recoD0Z, _mcJetPt, _mcD0Z, cweight*pSWeight[1]*pSWeight[0]);

                    }

                    if (isMissPtZ && MissingJets) {

                        rurResponse2D[centBin][0].Miss(_mcJetPt, _mcD0Z, weight*pSWeight[1]*pSWeight[0]);

                    } else if (isFakePtZ && FakeJets) {

                        rurResponse2D[centBin][0].Fake(_recoJetPt, _recoD0Z, weight*pSWeight[1]*pSWeight[0]);

                    }

                    for (int iAng = 0; iAng < 6; iAng++){


                        if (isMatchedPtLambda[iAng]) {

			    Double_t coord4D[4] = {_recoJetPt, _recoLambda[iAng], _mcJetPt, _mcLambda[iAng]};
    			    hResponseFine4D[centBin][1 + iAng]->Fill(coord4D, weight * pSWeight[1] * pSWeight[0]);

                            //hResponseFine4D[centBin][1 + iAng]->Fill(_recoJetPt, _recoLambda[iAng], _mcJetPt,
                            //                               _mcLambda[iAng], weight*pSWeight[2 + iAng]*pSWeight[0]);

                            rurResponse2D[centBin][iAng + 1].Fill(_recoJetPt, _recoLambda[iAng], _mcJetPt,
                                                                _mcLambda[iAng], weight*pSWeight[2 + iAng]*pSWeight[0]);
                            hRespZHighRes[2 + iAng][centBin]->Fill(_recoLambda[iAng],_mcLambda[iAng]);

                            hResVar[centBin][iAng + 2][pTJetBinCut]->Fill((_mcLambda[iAng] - _recoLambda[iAng])/_mcLambda[iAng], weight*pSWeight[2 + iAng]*pSWeight[0]);
                            


                        } else if (isMissPtLambda[iAng] && MissingJets){

                            rurResponse2D[centBin][iAng + 1].Miss(_mcJetPt, _mcLambda[iAng], weight*pSWeight[2 + iAng]*pSWeight[0]);

                        } else if (isFakePtLambda[iAng] && FakeJets){

                            rurResponse2D[centBin][iAng + 1].Fake(_recoJetPt, _recoLambda[iAng], weight*pSWeight[2 + iAng]*pSWeight[0]);

                        }





                    }
                }

            }
    }

    // ---------- cache ----------
    if (FillCacheRM){

        // ---------- 1D pT ----------
        if (isMatchedPtCache) {
            hCacheMatchPt[centBin]->Fill(_recoJetPt, _mcJetPt, weight);
        } else if (isMissPtCache && MissingJets) {
            hCacheMissPt[centBin]->Fill(_mcJetPt, weight);
        } else if (isFakePtCache && FakeJets) {
            hCacheFakePt[centBin]->Fill(_recoJetPt, weight);
        }

        // ---------- 2D pT-z ----------
        if (isMatchedPtZCache) {
    	    Double_t coord4D[4] = {_recoJetPt, _recoD0Z, _mcJetPt, _mcD0Z};
 	    hCacheMatchPtZ[centBin]->Fill(coord4D, weight);
            //hCacheMatchPtZ[centBin]->Fill(_recoJetPt, _recoD0Z, _mcJetPt, _mcD0Z, weight);
        } else if (isMissPtZCache && MissingJets) {
            hCacheMissPtZ[centBin]->Fill(_mcJetPt, _mcD0Z, weight);
        } else if (isFakePtZCache && FakeJets) {
            hCacheFakePtZ[centBin]->Fill(_recoJetPt, _recoD0Z, weight);
        }

        // ---------- 2D pT-lambda ----------
        for (int iAng = 0; iAng < 6; iAng++) {

            if (isMatchedPtLambdaCache[iAng]) {
                Double_t coord4D[4] = {_recoJetPt, _recoLambda[iAng], _mcJetPt, _mcLambda[iAng]};
   		hCacheMatchPtAng[centBin][iAng]->Fill(coord4D, weight);
             //   hCacheMatchPtAng[centBin][iAng]->Fill(_recoJetPt, _recoLambda[iAng],
             //                                       _mcJetPt, _mcLambda[iAng], weight);
            } else if (isMissPtLambdaCache[iAng] && MissingJets) {
                hCacheMissPtAng[centBin][iAng]->Fill(_mcJetPt, _mcLambda[iAng], weight);
            } else if (isFakePtLambdaCache[iAng] && FakeJets) {
                hCacheFakePtAng[centBin][iAng]->Fill(_recoJetPt, _recoLambda[iAng], weight);
            }
        }
    }

        //iEntry < (ClosureTest ? (TrainToTestRatio * nEntries) : 1.0 * nEntries)
      /*  if (ClosureTest && (
            (UseTheSameSample && iEntry < (TrainToTestRatio * nEntries))
            || (!UseTheSameSample && iEntry > ((TrainToTestRatio * nEntries)) && iEntry < (nEntries))
            )   )
        {*/
        if (ClosureTest && goesToTest){
            if (isRecoPt) hRealData[centBin].Fill(_recoJetPt, weight);
            if (Unfold2D) {
                if (isRecoPtZ) hRealData2D[centBin][0].Fill(_recoJetPt, _recoD0Z, weight);
                for (int ang = 0; ang < nAngularities; ang++) {
                    if (isRecoPtLambda[ang]) hRealData2D[centBin][ang + 1].Fill(_recoJetPt, _recoLambda[ang], weight);
                }
            }
        }




    }
    tEkvi->Write();
    fEkvi->Close();

    if (FillCacheRM)
    {
    
        TFile *fCache = new TFile("Output/CacheRM.root", "RECREATE");
        fCache->cd();

        // --- 1D pT ---
        for (int ic = 0; ic < nCentralityBins; ic++) {

            if (hCacheMatchPt[ic]) hCacheMatchPt[ic]->Write();
            if (hCacheMissPt[ic])  hCacheMissPt[ic]->Write();
            if (hCacheFakePt[ic])  hCacheFakePt[ic]->Write();

            // --- pT-z ---
            if (hCacheMatchPtZ[ic]) hCacheMatchPtZ[ic]->Write();
            if (hCacheMissPtZ[ic])  hCacheMissPtZ[ic]->Write();
            if (hCacheFakePtZ[ic])  hCacheFakePtZ[ic]->Write();

            // --- pT-lambda ---
            for (int iAng = 0; iAng < 6; iAng++) {

                if (hCacheMatchPtAng[ic][iAng]) hCacheMatchPtAng[ic][iAng]->Write();
                if (hCacheMissPtAng[ic][iAng])  hCacheMissPtAng[ic][iAng]->Write();
                if (hCacheFakePtAng[ic][iAng])  hCacheFakePtAng[ic][iAng]->Write();
            }
        }

        fCache->Close();
    }

    //save rurResponse2D[0-3][0] into file
    TFile *outFile = new TFile("Tests/rurResponse2D.root", "RECREATE");
    for (int c = 0; c < 3; c++) {
        rurResponse2D[c][0].Write();
    }
    outFile->Close();


    //kontrola stability response matrix
    TH2D* hMigRefScaled[3];
    TH2D* hMigRefScaledZ[3];
    TH2D* hMigRefScaledZ4D[3];

    TCanvas *cMig = new TCanvas("cMig", "cMig", 1200, 800);
    cMig->Divide(3,2);


    for (int iCent = 0; iCent < 3; iCent++){
        cMig->cd(iCent+1);
        //right margin
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.15);
        hMigRefScaled[iCent] = (TH2D*) hMigRef[iCent]->Clone(Form("hMigRefScaled_c%i", iCent));
        hMigRefScaled[iCent]->Divide(hMigRef[iCent],hMigWgt[iCent]);
        gStyle->SetPaintTextFormat(".1f%%");

        //vydělík každý řádek histogramu jeho průměrem
        for (int j = 1; j <= hMigRefScaled[iCent]->GetNbinsY(); j++) {
            double rowSum = 0;
            int rowCount = 0;
            for (int i = 1; i <= hMigRefScaled[iCent]->GetNbinsX(); i++) {
                double binContent = hMigRefScaled[iCent]->GetBinContent(i, j);
                if (binContent != 0) {
                    rowSum += binContent;
                    rowCount++;
                }
            }
            double rowAvg = (rowCount > 0) ? (rowSum / rowCount) : 1.0; // zabránit dělení nulou
            for (int i = 1; i <= hMigRefScaled[iCent]->GetNbinsX(); i++) {
                double binContent = hMigRefScaled[iCent]->GetBinContent(i, j);
                if (binContent != 0) {
                    double deviation = (binContent - rowAvg) / rowAvg * 100.0;
                    hMigRefScaled[iCent]->SetBinContent(i, j, deviation);
                } else {
                    hMigRefScaled[iCent]->SetBinContent(i, j, NAN);  // <–– prázdný bin
                    hMigRefScaled[iCent]->SetBinError(i, j, -1);

                }
            }


        }




        hMigRefScaled[iCent]->GetXaxis()->SetTitle("Reco pT");
        hMigRefScaled[iCent]->GetYaxis()->SetTitle("True pT");
        Stejn2(*hMigRefScaled[iCent],Form("check_%.d",iCent));

        //hMigRefScaled[iCent]->Draw("COLZTEXT");

        cMig->cd(iCent+4);
        //right margin
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.15);
        hMigRefScaledZ[iCent] = (TH2D*) hMigRefZ[iCent]->Clone(Form("hMigRefScaledZ_c%i", iCent));
        hMigRefScaledZ[iCent]->Divide(hMigRefZ[iCent],hMigWgtZ[iCent]);
        gStyle->SetPaintTextFormat(".1f%%");
        //vydělík každý řádek histogramu jeho průměrem
        for (int j = 1; j <= hMigRefScaledZ[iCent]->GetNbinsY(); j++) {
            double rowSum = 0;
            int rowCount = 0;
            for (int i = 1; i <= hMigRefScaledZ[iCent]->GetNbinsX(); i++) {
                double binContent = hMigRefScaledZ[iCent]->GetBinContent(i, j);
                if (binContent != 0) {
                    rowSum += binContent;
                    rowCount++;
                }
            }
            double rowAvg = (rowCount > 0) ? (rowSum / rowCount) : 1.0; // zabránit dělení nulou
            for (int i = 1; i <= hMigRefScaledZ[iCent]->GetNbinsX(); i++) {
                double binContent = hMigRefScaledZ[iCent]->GetBinContent(i, j);
                if (binContent != 0) {
                    double deviation = (binContent - rowAvg) / rowAvg * 100.0;
                    hMigRefScaledZ[iCent]->SetBinContent(i, j, deviation);
                } else {
                    hMigRefScaledZ[iCent]->SetBinContent(i, j, NAN);  // <–– prázdný bin
                    hMigRefScaledZ[iCent]->SetBinError(i, j, -1);
                }
            }
        }

        hMigRefScaledZ[iCent]->GetXaxis()->SetTitle("Reco z");
        hMigRefScaledZ[iCent]->GetYaxis()->SetTitle("True z");
        Stejn2(*hMigRefScaledZ[iCent],Form("checkZ_%.d",iCent));
    }
    cMig->SaveAs("OutputPdf/MigrationMatrixStabilityCheck.pdf[");

    cMig->SaveAs("OutputPdf/MigrationMatrixStabilityCheck.pdf");

    cMig->Clear();

    cMig->SaveAs("OutputPdf/MigrationMatrixStabilityCheck.pdf");
    cMig->Clear();
        cMig->SetCanvasSize(1200, 400);

    cMig->Divide(3,1);
    //size
    TH2D*ProjectionTestUnw[3];
    TH2D*ProjectionTestW[3];
    for (int iCent = 0; iCent < 3; iCent++) {
        cMig->cd(iCent + 1);
        //right margin
        gPad->SetRightMargin(0.15);
        gPad->SetLeftMargin(0.15);
        hMigRefScaledZ4D[iCent] = (TH2D *) hMigRef[iCent]->Clone(Form("hMigRefScaled_c%i", iCent));
        hMigRefScaledZ4D[iCent]->Divide(hMigRef[iCent], hMigWgt[iCent]);
        gStyle->SetPaintTextFormat(".1f%%");
        ProjectionTestUnw[iCent]  = (TH2D *) rurResponse2DTest[iCent].Hresponse()->Clone(Form("hResponseTruthTest22D_%i_%i", iCent,0 ));
        ProjectionTestW[iCent]    = (TH2D *) rurResponse2DTestW[iCent].Hresponse()->Clone(Form("hResponseTruthTestW22D_%i_%i", iCent,0 ));

        //ratio
        ProjectionTestW[iCent]->Divide(ProjectionTestUnw[iCent]);
        ProjectionTestW[iCent]->Draw("colz");

        for (int j = 1; j <= ProjectionTestW[iCent]->GetNbinsY(); j++) {
            double rowSum = 0;
            int rowCount = 0;
            for (int i = 1; i <= ProjectionTestW[iCent]->GetNbinsX(); i++) {
                double binContent = ProjectionTestW[iCent]->GetBinContent(i, j);
                if (binContent != 0) {
                    rowSum += binContent;
                    rowCount++;
                }
            }
            double rowAvg = (rowCount > 0) ? (rowSum / rowCount) : 1.0; // zabránit dělení nulou
            for (int i = 1; i <= ProjectionTestW[iCent]->GetNbinsX(); i++) {
                double binContent = ProjectionTestW[iCent]->GetBinContent(i, j);
                if (binContent != 0) {
                    double deviation = (binContent - rowAvg) / rowAvg * 100.0;
                    ProjectionTestW[iCent]->SetBinContent(i, j, deviation);
                } else {
                    ProjectionTestW[iCent]->SetBinContent(i, j, NAN);  // <–– prázdný bin
                    ProjectionTestW[iCent]->SetBinError(i, j, -1);
                }
            }
        }

        ProjectionTestW[iCent]->GetZaxis()->SetRangeUser(-100, 100);

        int nxbinsZ = hMigRefScaledZ[iCent]->GetXaxis()->GetNbins();
        int nybinsZ = hMigRefScaledZ[iCent]->GetYaxis()->GetNbins();
        int nxbinsPt = hMigRefScaled[iCent]->GetXaxis()->GetNbins();
        int nybinsPt = hMigRefScaled[iCent]->GetYaxis()->GetNbins();

        //draw vertical lines
        for (int i = 0; i < nybinsZ + 1; i++) {
            TLine *line = new TLine(0, nybinsPt * i, nxbinsZ * nxbinsPt, nybinsPt * i);

            line->SetLineColor(kBlack);
            line->SetLineStyle(1);
            line->SetLineWidth(1);
            line->Draw(i == 0 ? "" : "same");


        }

        for (int i = 0; i < nxbinsZ + 1; i++) {

            TLine *line2 = new TLine(nxbinsPt * i, 0, nxbinsPt * i, nybinsZ * nybinsPt);

            line2->SetLineColor(kBlack);
            line2->SetLineStyle(1);
            line2->SetLineWidth(1);
            line2->Draw(i == 0 ? "" : "same");
        }
/*
        ProjectionTestW[iCent]->GetXaxis()->SetLabelSize(0);
        ProjectionTestW[iCent]->GetYaxis()->SetLabelSize(0);
        ProjectionTestW[iCent]->GetXaxis()->SetTickLength(0);
        ProjectionTestW[iCent]->GetYaxis()->SetTickLength(0);
        */
        TH2D hResponseClone = *(TH2D *) ProjectionTestW[iCent]->Clone(Form("hResponseClone_%d", iCent));
        //xaxis title
        ProjectionTestW[iCent]->GetXaxis()->SetTitle("Reco z (Reco p_{T})");
        ProjectionTestW[iCent]->GetYaxis()->SetTitle("True z (True p_{T})");
        //offset
        ProjectionTestW[iCent]->GetXaxis()->SetTitleOffset(1.2);


TString var = "z";
   

    }


    cMig->Update();

    cMig->SaveAs("OutputPdf/MigrationMatrixStabilityCheck.pdf");
    cMig->SaveAs("OutputPdf/MigrationMatrixStabilityCheck.pdf]");

}


void LoadDataRealParams() {



    TFile *filea = TFile::Open(D0SpectraBeforeShift);
    if (!filea || filea->IsZombie()) {
        cout << "File not found or is corrupted: " << PaperD0Spectrum << endl;
        return;
    }
    hPaperD0PtBeforeShift[0] = (TGraphErrors *) filea->Get("gD0_err_0_10")->Clone("gD0_err_0_10_copy");
    hPaperD0PtBeforeShift[1] = (TGraphErrors *) filea->Get("gD0_err_10_40")->Clone("gD0_err_10_40_copy");
    hPaperD0PtBeforeShift[2] = (TGraphErrors *) filea->Get("gD0_err_40_80")->Clone("gD0_err_40_80_copy");

    filea->Close();

    //load spectra from PaperD0Spectrum
    TFile *file = TFile::Open(PaperD0Spectrum);
    //otevřu složku D^0 spectra in AuAu collisions
    if (!file || file->IsZombie()) {
        cout << "File not found or is corrupted: " << PaperD0Spectrum << endl;
        return;
    }
    //file->cd("D^0 spectra in AuAu collisions");
    //load PaperD0Spectrum
    hPaperD0Pt[0] = (TGraphAsymmErrors *) file->Get("D^0 spectra in AuAu collisions/Graph1D_y1")->Clone("hPaperD0Pt_y1");
    hPaperD0Pt[1] = (TGraphAsymmErrors *) file->Get("D^0 spectra in AuAu collisions/Graph1D_y6")->Clone("hPaperD0Pt_y6");
    hPaperD0Pt[2] = (TGraphAsymmErrors *) file->Get("D^0 spectra in AuAu collisions/Graph1D_y7")->Clone("hPaperD0Pt_y7");

    //close
    file->Close();

    //new canvas
    TCanvas *can2d = new TCanvas("can2d", "can2d", 800, 800);
    can2d->cd();


    //********************************/
    /********D0 jet spectra **********/
    /********************************/
    TFile *realDataFileNew = new TFile("./Output/" + TString(outputFile) + TString(Method)+ TString(_sys) + "2.root", "READ");
    if (!realDataFileNew || realDataFileNew->IsZombie()) {
        cout << "????" << endl;
        return;
    }
    TTree *treeReal = (TTree *) realDataFileNew->Get("jets");
    StJetTreeStruct2 measured;
    assignTree2(treeReal, measured);

    Double_t nEntries2 = treeReal->GetEntries();
    cout << "nEntries = " << (Float_t) nEntries2 / 1000. << "k" << endl
    << endl;

    treeReal->SetBranchStatus("*", 0);  // Vypni všechny větve
    treeReal->SetBranchStatus("d0Pt", 1);
    treeReal->SetBranchStatus("centralityAlt", 1);
    treeReal->SetBranchStatus("centrality", 1);
    treeReal->SetBranchStatus("gRefMult", 1);
    treeReal->SetBranchStatus("sWeightSignal", 1);
    treeReal->SetBranchStatus("weightD0Efficiency", 1);
    treeReal->SetBranchStatus("weightCentrality", 1);
    treeReal->SetBranchStatus("weightDoubleCount", 1);

    //Filling 2D histograms
    for (Int_t iEntry = 0; iEntry < nEntries2; iEntry++) {

        treeReal->GetEntry(iEntry);

        double _d0Pt = measured.vd0pt;
        double _centralityAlt = measured.centralAlt;
        double _centrality = measured.central;
        double _gRefMult = measured.gRefMultVal;

        double _sWeightSignal = measured.s_weight_value;
        double _weightD0Efficiency = measured.eff_weight_value;
        double _weightCentrality = measured.centr_weight_value;
        double _weightDoubleCount = measured.doubleCount;
        double _weight = _sWeightSignal * _weightD0Efficiency * _weightCentrality * _weightDoubleCount;


        if (_d0Pt < minPtD0Cut) continue;
        if (_d0Pt >= maxPtD0Cut) continue;
        
        int iCent = _centralityAlt < 10 ? 0 : _centralityAlt < 40 ? 1 : _centralityAlt < 80 ? 2 : -9;
        
        if (iCent == -9) {
            cout << "Centrality out of range: " << _centralityAlt << endl;
            exit(0);
        }

            nCentralityNumbers->Fill(_centrality, _weight);
            gRefMult->Fill(_gRefMult, _weight);
            D0MesonPtReal[iCent]->Fill(_d0Pt, _weight);
    }
    
    realDataFileNew->Close();


    //********************************/
    /******D0 meson spectra *********/
    /********************************/

    TFile *realDataFileNewD0 = new TFile("./Output/" + TString(outputFile) + "_D0" + TString(_sys)  + "2.root", "READ");

    if (!realDataFileNewD0 || realDataFileNewD0->IsZombie()) {
        cout << "????" << endl;
        exit(1);
    }

    TTree *treeRealD0 = (TTree *) realDataFileNewD0->Get("jets");

    StJetTreeStruct2 measuredD0;
    assignTree2(treeRealD0, measuredD0);

    Double_t nEntries3 = treeRealD0->GetEntries();
    cout << "nEntries = " << (Float_t) nEntries2 / 1000. << "k" << endl
         << endl;

        treeRealD0->SetBranchStatus("*", 0);  // Vypni všechny větve
        treeRealD0->SetBranchStatus("d0Pt", 1);
        treeRealD0->SetBranchStatus("centralityAlt", 1);
        treeRealD0->SetBranchStatus("centrality", 1);
        treeRealD0->SetBranchStatus("gRefMult", 1);
        treeRealD0->SetBranchStatus("sWeightSignal", 1);
        treeRealD0->SetBranchStatus("weightD0Efficiency", 1);
        treeRealD0->SetBranchStatus("weightCentrality", 1);
        treeRealD0->SetBranchStatus("weightDoubleCount", 1);

    for (Int_t iEntry = 0; iEntry < nEntries3; iEntry++) {

        treeRealD0->GetEntry(iEntry);

        double _d0Pt = measuredD0.vd0pt;
        double _centralityAlt = measuredD0.centralAlt;
        double _centrality = measuredD0.central;

        double _sWeightSignal = measuredD0.s_weight_value;
        double _weightD0Efficiency = measuredD0.eff_weight_value;
        double _weightCentrality = measuredD0.centr_weight_value;
        double _weightDoubleCount = measuredD0.doubleCount;
        double _weight = _sWeightSignal * _weightD0Efficiency * _weightCentrality * _weightDoubleCount;

        int iCent = _centralityAlt < 10 ? 0 : _centralityAlt < 40 ? 1 : _centralityAlt < 80 ? 2 : -9;
        if (iCent == -9) {
            cout << "Centrality out of range: " << _centralityAlt << endl;
            exit(1);
        }

        hMeasuredD0MesonPt[iCent]->Fill(_d0Pt, _weight);

    }


    realDataFileNewD0->Close();



    //********************************/
    /******MC D0 jet spectra *********/
    /********************************/

    TFile *treeFile; // Open the file containing the tree.
    treeFile = new TFile(McJetsFileData, "READ");
    if (!treeFile || treeFile->IsZombie()) return;

    TTree *jetTree = (TTree *) treeFile->Get("jets");
/*
    // 1) vypnout všechno
    jetTree->SetBranchStatus("*", 0);

    // 2) povolit jen to, co fakt používáš v
    // (přesné názvy větví uprav podle assignTree / struktury)
    jetTree->SetBranchStatus("centrality", 1);
    jetTree->SetBranchStatus("centralityAlt", 1);
    jetTree->SetBranchStatus("weightCentrality", 1);
    jetTree->SetBranchStatus("gRefMult", 1);
    jetTree->SetBranchStatus("mcJetPt", 1);
    jetTree->SetBranchStatus("recoJetPt", 1);
    jetTree->SetBranchStatus("mcJetEta", 1);
    jetTree->SetBranchStatus("*recoJetEta", 1);
    jetTree->SetBranchStatus("mcD0Pt", 1);
    jetTree->SetBranchStatus("mcSmearedD0Pt", 1);
    jetTree->SetBranchStatus("mcJetD0Z", 1);
    jetTree->SetBranchStatus("*recoJetD0Z", 1);
    jetTree->SetBranchStatus("mcJetLambda*", 1);
    jetTree->SetBranchStatus("*recoJetLambda*", 1);
    jetTree->SetBranchStatus("mcJetNConst", 1);
    jetTree->SetBranchStatus("*recoJetNConst", 1);
*/
    // 3) TTreeCache (u vzdálených/velkých souborů často gamechanger)
    jetTree->SetCacheSize(128*1024*1024);          // klidně 256 MB
    jetTree->AddBranchToCache("*", kTRUE);
    jetTree->SetCacheLearnEntries(1000);

    StJetTreeStruct mcJet, recoJet;
    assignTree(jetTree, mcJet, recoJet);
    Double_t nEntries = jetTree->GetEntries();
    nEntries /= DividedMcDataBy;
    cout << "nEntries = " << (Float_t) nEntries / 1000. << "k" << endl
         << endl;

    TH1D *hCentr = new TH1D("hCentr", "hCentr", 9, -0.5, 8.5);

    for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {

        Float_t progress = 0.;
        progress = (Float_t) iEntry / nEntries;


        if (iEntry % 10000 == 0) cout << "Training: \r (" << (progress * 100.0) << "%)" << std::flush;


        jetTree->GetEntry(iEntry);

        //Event
        double _centrality = recoJet.centrality;
        double _centralityAlt = mcJet.centralityAlt;
        double _weightCentrality = recoJet.weight;
        double _gRefMult = mcJet.gRefMultMc;
        double _eventMaxPtTrack = recoJet.eventMaxPtTrack;

        if (_eventMaxPtTrack > 30) continue;

        //D0 meson
        double _mcD0Pt = mcJet.d0pt;
        double eta = mcJet.d0eta;

        double _mcSmearedD0Pt = recoJet.RecoD0Pt;
        double _mcSmearedD0Eta = recoJet.RecoD0Eta;

        double mD0 = 1.864; // GeV

        double pz = _mcD0Pt* sinh(eta);
        double p  = _mcD0Pt* cosh(eta);
        double E  = sqrt(p*p + mD0*mD0);

        double y = 0.5 * log((E + pz) / (E - pz));

        if (abs(y) > 0.6) continue; //D0 meson rapidity cut


        //D0-jet
        double _mcJetEta = mcJet.mcJetEta;
        double _recoJetEta = recoJet.recoJetEta;
        double _mcJetPt = mcJet.jetpt;
        double _recoJetPt = recoJet.jetpt;
        double _mcNConst = mcJet.numberofconstituents;

        ////double _recoJetPt = recoJet.jetpt_nocorr - (recoJet.recoJetRho+0.0) * recoJet.recoJetArea;
        double _mcD0Z = mcJet.d0z;
        double _recoD0Z = recoJet.d0z;
        double _mcLambda[6] =   {mcJet.lambda[0], mcJet.lambda[1], mcJet.lambda[2],
                                mcJet.lambda[3], mcJet.lambda[4], mcJet.lambda[5]};
        double _recoLambda[6] = {recoJet.lambda[0], recoJet.lambda[1], recoJet.lambda[2],
                                recoJet.lambda[3], recoJet.lambda[4], recoJet.lambda[5]};
        double _recoNConst = recoJet.numberofconstituents;                        


        Int_t centBin9 = getCentralityBin99(recoJet.centrality); //0, 1, 2
        Int_t centBin = recoJet.centrality;


        if (_centrality < 0) continue;
        if (_mcD0Pt < minPtD0Cut) continue;
        if (_mcD0Pt >= maxPtD0Cut) continue;
        if (_mcSmearedD0Pt < minPtD0Cut) continue;
        if (_mcSmearedD0Pt >= maxPtD0Cut) continue;

        if (abs(_mcJetEta) > 0.6 && abs(_recoJetEta) > 0.6) continue; //tyhle vůbec nechci
        if (_mcJetPt > 30) continue;

        if (deleteOneConstituentJets && _mcNConst == 1 && _recoNConst == 1) continue;

        JetPtZMc[centBin9]->Fill(_mcJetPt, _mcD0Z, _weightCentrality);
        nCentralityNumbersMC->Fill(centBin, _weightCentrality);
        gRefMultMc->Fill(_gRefMult, _weightCentrality);
        McPTRawD0[centBin9]->Fill(_mcD0Pt, _weightCentrality);
        McPTRawD0Jet[centBin9]->Fill(_mcJetPt, _weightCentrality);
        McPTRawD0JetD0Meson[centBin9]->Fill(_mcJetPt, _mcD0Pt, _weightCentrality);


        //kin eff:
        double recoJetPtMin = ptRecoBinsVec[centBin9][0];
        double recoJetPtMax = ptRecoBinsVec[centBin9][ptRecoBinsVec[centBin9].size() - 1];

        //eta KinEffEta
        if (_recoNConst > 0 && (abs(_recoJetEta) > 0.6 && (abs(_mcJetEta) < 0.6))) {
            KinEffEta[centBin9].Fill(true, _mcJetPt); //cout << "Filling KinEffEta false, _mcJetEta = " << _mcJetEta << endl;
    
        } else if (_recoNConst > 0 && (abs(_mcJetEta) < 0.6)&& (abs(_recoJetEta) < 0.6)) {
            KinEffEta[centBin9].Fill(false, _mcJetPt); //cout << "Filling KinEffEta true, _mcJetEta = " << _mcJetEta << endl;
        }

        //Real (1-fake)
        if (_recoNConst > 0 && (abs(_recoJetEta) < 0.6 && abs(_mcJetEta) > 0.6)){
            FakeEffEta[centBin9].Fill(true, _mcJetPt);
        } else if (_recoNConst > 0 && (abs(_mcJetEta) < 0.6 && (abs(_recoJetEta) < 0.6))){
            FakeEffEta[centBin9].Fill(false, _mcJetPt);
        }
        
        //1D pT //Měly by tu být váhy??? TO DO!!!
        if (_recoNConst > 0 && ((_recoJetPt) < recoJetPtMin || (_recoJetPt) >= recoJetPtMax)){
           KinEff1D[centBin9].Fill(false, _mcJetPt);
        } else if (_recoNConst > 0){
            KinEff1D[centBin9].Fill(true, _mcJetPt);
        }

        double recoJetZMin = zRecoBinsVec[centBin9][0];
        double recoJetZMax = zRecoBinsVec[centBin9][zRecoBinsVec[centBin9].size() - 1];
        double trueJetZMin = zMcBinsVecCustom[centBin9][0];
        double trueJetZMax = zMcBinsVecCustom[centBin9][zMcBinsVecCustom[centBin9].size() - 1];
        //2D pT z
        if (_recoNConst > 0 && ((_recoJetPt) < recoJetPtMin || (_recoJetPt) >= recoJetPtMax || (_recoD0Z) < recoJetZMin || (_recoD0Z) > recoJetZMax   )){

            KinEff2DpTZ[0][centBin9].Fill(false, _mcJetPt);
            KinEff2DpTZ[1][centBin9].Fill(false, _mcD0Z);
            KinEff2DpTZ[2][centBin9].Fill(false, _mcJetPt, _mcD0Z);

            if (_mcJetPt >= 5 && _mcJetPt < 20) KinEff2DpTZCut[centBin9].Fill(false, _mcD0Z);

        } else if (_recoNConst > 0){
            KinEff2DpTZ[0][centBin9].Fill(true, _mcJetPt);
            KinEff2DpTZ[1][centBin9].Fill(true, _mcD0Z);
            KinEff2DpTZ[2][centBin9].Fill(true, _mcJetPt, _mcD0Z);
            if (_mcJetPt >= 5 && _mcJetPt < 20)     KinEff2DpTZCut[centBin9].Fill(true, _mcD0Z);
        }
        if (_recoNConst > 0 &&
        (   _recoJetPt < recoJetPtMin    || _recoJetPt >= recoJetPtMax ||
            _recoD0Z < recoJetZMin       || _recoD0Z > recoJetZMax     ||
            _mcD0Z < trueJetZMin         || _mcD0Z > trueJetZMax)     ){

            KinEff2DpTZZCut[centBin9].Fill(false, _mcJetPt);

        } else if (_recoNConst > 0) KinEff2DpTZZCut[centBin9].Fill(true, _mcJetPt);


            for (int iAng = 0; iAng < nAngularities; iAng++) {


            double recoJetVarMin = angRecoBinsVec[centBin9][iAng][0];
            double recoJetVarMax = angRecoBinsVec[centBin9][iAng][angRecoBinsVec[centBin9][iAng].size() - 1];
            double trueJetVarMin = angMcBinsVecCustom[centBin9][iAng][0];
            double trueJetVarMax = angMcBinsVecCustom[centBin9][iAng][angMcBinsVecCustom[centBin9][iAng].size() - 1];
            //2D pT lambda
            if ((_recoNConst > 0 && ((_recoJetPt) < recoJetPtMin || (_recoJetPt) >= recoJetPtMax || (_recoLambda[iAng]) < recoJetVarMin || (_recoLambda[iAng]) >= recoJetVarMax))) {
                KinEff2DAng[iAng][0][centBin9].Fill(false, _mcJetPt);
                KinEff2DAng[iAng][1][centBin9].Fill(false, _mcLambda[iAng]);
                KinEff2DAng[iAng][2][centBin9].Fill(false, _mcJetPt, _mcLambda[iAng]);
                if (_mcJetPt >= 5 && _mcJetPt < 20)     KinEff2DAngCut[iAng][centBin9].Fill(false, _mcLambda[iAng]);
            } else if (_recoNConst > 0) {
                KinEff2DAng[iAng][0][centBin9].Fill(true, _mcJetPt);
                KinEff2DAng[iAng][1][centBin9].Fill(true, _mcLambda[iAng]);
                KinEff2DAng[iAng][2][centBin9].Fill(true, _mcJetPt, _mcLambda[iAng]);
                if (_mcJetPt >= 5 && _mcJetPt < 20)     KinEff2DAngCut[iAng][centBin9].Fill(true, _mcLambda[iAng]);
            }

                if (((_recoNConst > 0 && ((_recoJetPt) < recoJetPtMin ||
                (_recoJetPt) >= recoJetPtMax || (_recoLambda[iAng]) < recoJetVarMin ||
                (_recoLambda[iAng]) >= recoJetVarMax)) || _mcLambda[iAng] < trueJetVarMin ||
                _mcLambda[iAng] >= trueJetVarMax)) {

                    KinEff2DAngPtCut[iAng][centBin9].Fill(false, _mcJetPt);
                } else if (_recoNConst > 0) {

                    KinEff2DAngPtCut[iAng][centBin9].Fill(true, _mcJetPt);
                }


        }
    }

    //close file
    treeFile->Close();

    //logy
    gPad->SetLogy();
    nCentralityNumbers->Draw("histtext");
    nCentralityNumbers->SetLineColor(kRed);
    nCentralityNumbersMC->Draw("histtext same");

    //nastavím maximum na yové ose, aby byly vidět oba histogramy:
    nCentralityNumbers->GetYaxis()->SetRangeUser(min(nCentralityNumbers->GetMinimum(), nCentralityNumbersMC->GetMinimum())*0.8, max(nCentralityNumbers->GetMaximum(), nCentralityNumbersMC->GetMaximum())*1.2);

    realDataFileNew->Close();

    double bins[10] = {0,5,10,20,30,40,50,60,70,80};
    TH1D *nCentralityNumbersAlt = new TH1D("nCentralityNumbersAlt", "nCentralityNumbersAlt;Centrality [%];Counts",
             9, bins);
    for (int i = 1; i <= 9; ++i) {
        nCentralityNumbersAlt->SetBinContent(i, nCentralityNumbers->GetBinContent(9-i+1));
        nCentralityNumbersAlt->SetBinError  (i, nCentralityNumbers->GetBinError(9-i+1));
    }
    //new file
    TFile *treeFile_A; // Open the file containing the tree.
    treeFile_A= new TFile("./Output/CentralityDistr.root", "RECREATE");
    treeFile_A->cd();
    nCentralityNumbers->Write();
    nCentralityNumbersAlt->Write();
    treeFile_A->Close();

}

void LoadDataRealNeil() {
    NumberOfWEvents[0] =1.1248969e+08;
    NumberOfWEvents[1] =3.5062620e+08;
    NumberOfWEvents[2] =4.7563863e+08;
    LoadEfficiency1DPaper();

    // Načtení DoubleCounting histogramů
    TFile* fileDC = TFile::Open(DoubleCounting);
    if (!fileDC) return;

    TH1D* histDC[5];
    for (int i = 0; i < 5; ++i) {
        histDC[i] = (TH1D*)fileDC->Get(Form("DC_cent%d", i));
        if (!histDC[i]) {
            cout << "Error: DC histogram " << i << " not found!" << endl;
            return;
        }
    }


    //  TFile *treeFile = new TFile("./Data/Output_Step_0_IterParam_0_njptbin_12_nzbin_7.root", "READ");
    //  TFile *treeFile = new TFile("./Data/RM_dw.root", "READ");
    TFile *treeFile = new TFile("../Neil_Unfold/Data/Histograms_D01_10GeV_RecoJetPt_1_1000.root", "READ");

    if (!treeFile || treeFile->IsZombie()) {
        std::cout << "Error opening file: Response_Step_0_IterParam_0_njptbin_12_nzbin_7.root" << std::endl;
        return;
    }

    TString Names[3] = {"0_10", "10_40", "40_80"};
    for (int iCent = 0; iCent < 3; iCent++) {
        //  hRealDataCopy[iCent] = *(TH1D *)treeFile->Get(Form("Measured1D_%i", iCent));
        //  hRealData2DCopy[iCent][0] = *(TH2D *)treeFile->Get(Form("MeasuredWide_%i", iCent));
        hRealDataCopy[iCent] = *(TH1D *)treeFile->Get(Form("JetPt_Area_Wide_%s", Names[iCent].Data()));
        hRealData2DCopy[iCent][0] = *(TH2D *)treeFile->Get(Form("ZPt_Area_Wide_%s", Names[iCent].Data()));
    }

    treeFile->Close();


    //return;



    TFile *realDataFileNew = new TFile("../Neil_Unfold/Data/TestMerge2.root", "READ");

    if (!realDataFileNew || realDataFileNew->IsZombie()) {
        cout << "File " << outputFile << "2.root not found" << endl;
        exit(0);
    }

    TTree *treeReal = (TTree *) realDataFileNew->Get("Signal_sw");

    //Načtu branches
    Float_t mM = -69.;
    Float_t mR = -69;
    Float_t mZ = -69;
    Float_t mZArea = -69;
    Float_t mPt = -420.;
    Float_t mEta = -69;
    Float_t mJPt = -69;
    Float_t mJPtArea = -69;
    Int_t mCen = -69;
    Float_t mgRefMultCorr = -69;
    Float_t mJEta = -69;
    Float_t mWeight = -69;
    Float_t mJetArea = -69;
    Int_t mConstArea = -69;
    Int_t mConstCS = -69;
    Float_t mJetRhoVal = -69;
    Float_t sWeight = -69; // S-weight for the event

    //načtu je z treeReal
    treeReal->SetBranchAddress("mM", &mM);
    treeReal->SetBranchAddress("mR", &mR);
    treeReal->SetBranchAddress("mZ", &mZ);
    treeReal->SetBranchAddress("mZArea", &mZArea);
    treeReal->SetBranchAddress("sWeight", &sWeight);
    treeReal->SetBranchAddress("mPt", &mPt);
    treeReal->SetBranchAddress("mEta", &mEta);
    treeReal->SetBranchAddress("mJPt", &mJPt);
    treeReal->SetBranchAddress("mJPtArea", &mJPtArea);
    treeReal->SetBranchAddress("mCen", &mCen);
    treeReal->SetBranchAddress("mGRefMult", &mgRefMultCorr);
    treeReal->SetBranchAddress("mJEta", &mJEta);
    treeReal->SetBranchAddress("mWeight", &mWeight);
    treeReal->SetBranchAddress("mJetArea", &mJetArea);
    treeReal->SetBranchAddress("mConstArea", &mConstArea);
    treeReal->SetBranchAddress("mConstCS", &mConstCS);
    treeReal->SetBranchAddress("mJetRhoVal", &mJetRhoVal);


    Double_t nEntries = treeReal->GetEntries(); // Ensure nEntries is initialized

    //Filling 2D histograms
    for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {


        treeReal->GetEntry(iEntry);

//        treeReal->GetEntry(iEntry);



        int iCent = -9;
        int aCent = mCen;
        if (aCent < 10) {
            iCent = 0;
        } else if (aCent >= 10 && aCent < 40) {
            iCent = 1;
        } else if (aCent >= 40 && aCent < 80) {
            iCent = 2;
        } else {
            cout << "Centrality out of range: " << aCent << endl;
            continue;
        }

        //????
        if (mZ >= 1.) mZ = 0.999;
        if (mJPt < mPt) mJPt = mPt; //Nikdy se nestane
        if (mJPtArea <= nbinsjetpt_wide[0]) mJPtArea = nbinsjetpt_wide[0] + 0.001;
        if (mJPtArea >= nbinsjetpt_wide[njpt_bins_wide]) mJPtArea = nbinsjetpt_wide[njpt_bins_wide] - 0.001;
        if (mZArea <= nbinsz_wide[0]) mZArea = nbinsz_wide[0] + 0.001;
        if (mZArea >= nbinsz_wide[nz_bins_wide]) mZArea = nbinsz_wide[nz_bins_wide] - 0.001;

        //????????

        // double rec_eff = 1./D0_efficiency_1DPaper(mPt, mCen);/////
        double rec_eff = 1./getEff(mCen, mPt);

        double d0_double_count = (1-getDoubleCount(mCen, mPt));
        // double d0_double_count = 1 - D0_DoubleCounting(histDC, mPt, aCent);

        double weight =  mWeight * sWeight * d0_double_count * rec_eff;



        int d0PtBin = -1;
        if (mPt >= 1 && mPt < 2) d0PtBin = 0;
        else if (mPt >= 2 && mPt < 3) d0PtBin = 1;
        else if (mPt >= 3 && mPt < 4) d0PtBin = 2;
        else if (mPt >= 4 && mPt < 5) d0PtBin = 3;
        else if (mPt >= 5 && mPt < 10) d0PtBin = 4;


        hRealData[iCent].Fill(mJPtArea, weight);
        if (Unfold2D) {
            // hRealData2D[iCent][0].Fill(mJPtArea, mZArea, weight);
            hRealData2DD0Pt[iCent][d0PtBin]->Fill(mJPtArea, mZArea, weight);
            for (int ang = 0; ang < 6; ang++) {
                //hRealData2D[iCent][ang+1].Fill(mJPt, measured.lambda_value[ang], weight);

            }
        }


    }

    //Vykreslím na canvas s pady 3x2
/*
    //divide all hRealData[i] by bin width
    for (int i = 0; i < nCentralityBins; i++) {
        for (int j = 1; j <= hRealData[i].GetNbinsX(); j++) {
            double binWidth = hRealData[i].GetBinWidth(j);
            if (binWidth > 0) {
                hRealData[i].SetBinContent(j, hRealData[i].GetBinContent(j) / binWidth);
                hRealData[i].SetBinError(j, hRealData[i].GetBinError(j) / binWidth);
            }
        }
        for (int ang = 0; ang < (nAngularities + 1); ang++) {
            for (int j = 1; j <= hRealData2D[i][ang].GetNbinsX(); j++) {
                for (int k = 1; k <= hRealData2D[i][ang].GetNbinsY(); k++) {
                    double binWidthX = hRealData2D[i][ang].GetXaxis()->GetBinWidth(j);
                    double binWidthY = hRealData2D[i][ang].GetYaxis()->GetBinWidth(k);
                    if (binWidthX > 0 && binWidthY > 0) {
                        hRealData2D[i][ang].SetBinContent(j, k, hRealData2D[i][ang].GetBinContent(j, k) / (binWidthX * binWidthY));
                        hRealData2D[i][ang].SetBinError(j, k, hRealData2D[i][ang].GetBinError(j, k) / (binWidthX * binWidthY));
                    }
                }
            }
        }
    }
*/

    if (true) {
        for (int i = 0; i < nCentralityBins; i++) {
            for (int j = 1; j <= hRealData[i].GetNbinsX(); j++) {
                if (hRealData[i].GetBinContent(j) < 0) {
                    cout << "Histogram " << hRealData[i].GetName() << " has negative bin content: "
                         << hRealData[i].GetBinContent(j) << endl;
                    hRealData[i].SetBinContent(j, 0);
                }
            }
            for (int ang = 0; ang < 1; ang++) {
                for (int j = 1; j <= hRealData2D[i][ang].GetNbinsX(); j++) {
                    for (int k = 1; k <= hRealData2D[i][ang].GetNbinsY(); k++) {

                        /*  if (hRealData2D[i][ang].GetBinContent(j, k) < 0) {
                              cout << "Histogram " << hRealData2D[i][ang].GetName() << " has negative bin content: "
                                   << hRealData2D[i][ang].GetBinContent(j, k) << endl;
                              hRealData2D[i][ang].SetBinContent(j, k, 0);
                          }*/
                        for (int d0pt = 0; d0pt < 5; d0pt++) {
                            if (hRealData2DD0Pt[i][d0pt]->GetBinContent(j, k) < 0) {
                                cout << "Histogram " << hRealData2DD0Pt[i][d0pt]->GetName()
                                     << " has negative bin content: "
                                     << hRealData2DD0Pt[i][d0pt]->GetBinContent(j, k) << endl;
                                hRealData2DD0Pt[i][d0pt]->SetBinContent(j, k, 0);
                            }
                        }
                    }
                }
            }
        }
    }


    //write down binning of hRealData2DD0Pt[0][0] to cout


    //merge hRealData2DD0Pt[i][all] to hRealData2D[i][0]
    for (int i = 0; i < nCentralityBins; i++) {
        for (int d0pt = 0; d0pt < 5; d0pt++) {
            hRealData2D[i][0].Add(hRealData2DD0Pt[i][d0pt]);
            //number of new entries
            std::cout << "Merged hRealData2DD0Pt[" << i << "][" << d0pt << "] into hRealData2D[" << i
                      << "][0]. New entries: " << hRealData2DD0Pt[i][d0pt]->GetEntries() << std::endl;
        }
    }

    //clear hRealData
    for (int i = 0; i < nCentralityBins; i++) {
        //Make it projection of hRealData2D[i][0]
        hRealDataXXX[i] = *(TH1D *) hRealData2D[i][0].ProjectionX(
                Form("%s_ProjX", hRealData2D[i][0].GetName()))->Clone(
                Form("%s_ProjX_Clone", hRealData2D[i][0].GetName()));
        hRealData[i].Reset();
        hRealData[i] = hRealDataXXX[i];
    }
    TCanvas *can = new TCanvas("can", "can", 1200, 800);
    can->Divide(3, 2);
    can->cd(1);
    hRealDataXXX[0].SetTitle("Real Data Centrality 0");
    hRealDataXXX[0].Draw("hist");
    //set range
    hRealDataXXX[0].GetYaxis()->SetRangeUser(0.1, 1e7);
    gPad->SetLogy();
    hRealDataCopy[0].SetLineColor(kRed);
    hRealDataCopy[0].Draw("histsame");

    // hRealData2DD0Pt[0][0].SetLineColor(kGreen);
    // hRealData2DD0Pt[0][0].Draw("histsame");

    can->cd(2);
    hRealDataXXX[1].SetTitle("Real Data Centrality 1");
    hRealDataXXX[1].Draw("hist");
    hRealDataXXX[1].GetYaxis()->SetRangeUser(0.1, 1e7);
    gPad->SetLogy();
    hRealDataCopy[1].SetLineColor(kRed);
    hRealDataCopy[1].Draw("histsame");


    can->cd(3);
    hRealDataXXX[2].SetTitle("Real Data Centrality 2");
    hRealDataXXX[2].Draw("hist");
    hRealDataXXX[2].GetYaxis()->SetRangeUser(0.1, 1e7);
    gPad->SetLogy();
    hRealDataCopy[2].SetLineColor(kRed);
    hRealDataCopy[2].Draw("histsame");

    can->cd(4);
    hRealData2D[0][0].SetTitle("Real Data Centrality 0 2D");
    hRealData2D[0][0].Draw("colz");
    gPad->SetLogz();
    hRealData2D[0][0].GetYaxis()->SetRangeUser(-2, 2);

    can->cd(5);
    hRealData2D[1][0].SetTitle("Real Data Centrality 1 2D");
    hRealData2D[1][0].Draw("colz");
    gPad->SetLogz();
    hRealData2D[1][0].GetYaxis()->SetRangeUser(-2, 2);
    can->cd(6);
    hRealData2D[2][0].SetTitle("Real Data Centrality 2 2D");
    hRealData2D[2][0].Draw("colz");
    gPad->SetLogz();
    hRealData2D[2][0].GetYaxis()->SetRangeUser(-2, 2);

    //y range
    can->SaveAs("./OutputPdf/RealDataNeil.pdf");

    //pause

    //exit(0);
}

void LoadDataReal() {



    TFile *realDataFileNew = new TFile("./Output/" + TString(outputFile) + "_" +TString(Method) + TString(_sys) + "2.root", "READ");

    if (!realDataFileNew || realDataFileNew->IsZombie()) {
        cout << "File " << outputFile << "2.root not found" << endl;
        exit(0);
    }

    TTree *treeReal = (TTree *) realDataFileNew->Get("jets");


    StJetTreeStruct2 measured;
    assignTree2(treeReal, measured);

    Double_t nEntries = treeReal->GetEntries();
        cout << "nEntries = " << (Float_t) nEntries / 1000. << "k" << endl
             << endl;

        //Filling 2D histograms
        for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {
            treeReal->GetEntry(iEntry);


            double _d0Pt = measured.vd0pt;
            double _z = measured.z_value;
            double _jetPt = measured.pT_value;
            double _lambda[6] = {measured.lambda_value[0], measured.lambda_value[1], measured.lambda_value[2],
                                  measured.lambda_value[3], measured.lambda_value[4], measured.lambda_value[5]};

            double _nJetConst = measured.nJetConst;

            double _centralityAlt = measured.centralAlt;
            double _centrality = measured.central;
            double _gRefMult = measured.gRefMultVal;

            double _sWeightSignal = measured.s_weight_value;
            double _weightD0Efficiency = measured.eff_weight_value;
            double _weightCentrality = measured.centr_weight_value;
            double _weightDoubleCount = measured.doubleCount;
            double _weight = _sWeightSignal * _weightD0Efficiency * _weightCentrality * _weightDoubleCount;

            if (_d0Pt < minPtD0Cut) continue;
            if (_d0Pt >= maxPtD0Cut) continue;
            if (deleteOneConstituentJets && _nJetConst == 1) continue;

            int iCent = _centralityAlt < 10 ? 0 : _centralityAlt < 40 ? 1 : _centralityAlt < 80 ? 2 : -9;

            int d0PtBin = -1;
            if (_d0Pt >= 1 && _d0Pt < 2) d0PtBin = 0;
            else if (_d0Pt >= 2 && _d0Pt < 3) d0PtBin = 1;
            else if (_d0Pt >= 3 && _d0Pt < 4) d0PtBin = 2;
            else if (_d0Pt >= 4 && _d0Pt < 5) d0PtBin = 3;
            else if (_d0Pt >= 5 && _d0Pt < 10) d0PtBin = 4;

            hRealData[iCent].Fill(_jetPt, _weight);
            hRealFine[iCent][0]->Fill(_jetPt,_z, _weight);
            hRealData1D0Pt[iCent][d0PtBin].Fill(_jetPt, _weight);
            if (Unfold2D) {
                hRealData2D[iCent][0].Fill(_jetPt, _z, _weight);
                hRealData2DD0Pt[iCent][d0PtBin][0].Fill(_jetPt, _z, _weight);
                for (int ang = 0; ang < 6; ang++) {
                    hRealData2D[iCent][ang+1].Fill(_jetPt, _lambda[ang], _weight);
                    hRealFine[iCent][1 + ang]->Fill(_jetPt, _lambda[ang], _weight);
                    hRealData2DD0Pt[iCent][d0PtBin][1 + ang].Fill(_jetPt, _lambda[ang], _weight);
                }
            }

        }

    if (CutOfNegative) {
        int CutNegativeMode = 2; // 0 none, 1 post-sum, 2 pre-sum

        // --- (A) PRE-SUM CUT: vynuluj negativní biny na slice histogramech dřív, než je sečteš
        if (CutNegativeMode == 2) {
            for (int i = 0; i < nCentralityBins; i++) {
                for (int d0pt = 0; d0pt < 5; d0pt++) {
                    //1D
                    for (int j = 1; j <= hRealData1D0Pt[i][d0pt].GetNbinsX(); j++) {
                        if (hRealData1D0Pt[i][d0pt].GetBinContent(j) < 0)
                            hRealData1D0Pt[i][d0pt].SetBinContent(j, 0.0);
                    }

                    //2D
                    for (int ang = 0; ang < (nAngularities + 1); ang++) {
                        for (int j = 1; j <= hRealData2DD0Pt[i][d0pt][ang].GetNbinsX(); j++) {
                            for (int k = 1; k <= hRealData2DD0Pt[i][d0pt][ang].GetNbinsY(); k++) {
                                if (hRealData2DD0Pt[i][d0pt][ang].GetBinContent(j,k) < 0)
                                    hRealData2DD0Pt[i][d0pt][ang].SetBinContent(j,k, 0.0);
                            }
                        }
                    }
                }
            }
        }

        // --- (B) vždycky přepočítat hRealData2D jako sumu D0pT-sliců
        for (int i = 0; i < nCentralityBins; i++) {
            //1D
            hRealData[i].Reset("ICES");
            for (int d0pt = 0; d0pt < 5; d0pt++) {
                hRealData[i].Add(&hRealData1D0Pt[i][d0pt]);
            }

            //2D
            for (int ang = 0; ang < (nAngularities + 1); ang++) {
                hRealData2D[i][ang].Reset("ICES");
                for (int d0pt = 0; d0pt < 5; d0pt++) {
                    hRealData2D[i][ang].Add(&hRealData2DD0Pt[i][d0pt][ang]);
                }
            }
        }

        // --- (C) POST-SUM CUT: vynuluj negativní biny až na sečteném hRealData2D
        if (CutNegativeMode == 1) {
            for (int i = 0; i < nCentralityBins; i++) {
                //1D
                for (int j = 1; j <= hRealData[i].GetNbinsX(); j++) {
                    if (hRealData[i].GetBinContent(j) < 0)
                        hRealData[i].SetBinContent(j, 0.0);
                }

                //2D
                for (int ang = 0; ang < (nAngularities + 1); ang++) {
                    for (int j = 1; j <= hRealData2D[i][ang].GetNbinsX(); j++) {
                        for (int k = 1; k <= hRealData2D[i][ang].GetNbinsY(); k++) {
                            if (hRealData2D[i][ang].GetBinContent(j,k) < 0)
                                hRealData2D[i][ang].SetBinContent(j,k, 0.0);
                        }
                    }
                }
            }
        }
    }

}


void LoadDataFactors() {


    //otevřít soubor s FONLL
    TFile *fileFONLL2 = TFile::Open("./Output/FONLL_D0Pt_0_31.root");
    
    //TFile *fileFONLL2 = TFile::Open("./Output/FONLL_cPty0_6cut_0_31.root");
    //check if file is opened
    if (!fileFONLL2 || fileFONLL2->IsZombie()) {
        cout << "File FONLL_D0Pt_0_31.root not found" << endl;
        exit(1);
    }

    pureFonll = (TF1*) fileFONLL2->Get("fSpline");

    //check if pureFonll is loaded
    if (!pureFonll) {
        cout << "pureFonll not found in file" << endl;
        exit(1);
    }

    fileFONLL2->Close();

    //otevřít soubor s FONLL
    TFile *fileFONLL3 = TFile::Open("./Output/FONLL_D0Pt_0_10.root");
    //check if file is opened
    if (!fileFONLL3 || fileFONLL3->IsZombie()) {
        cout << "File FONLL_D0Pt_0_10.root not found" << endl;
        exit(1);
    }

    pureFonllD0meson = (TF1*) fileFONLL3->Get("fSpline");

    //check if pureFonll is loaded
    if (!pureFonllD0meson) {
        cout << "pureFonll not found in file" << endl;
        exit(1);
    }

    fileFONLL3->Close();

}

void LoadDataNeilRM() {
    TFile *treeFile = new TFile("../Neil_Unfold/Data/Response_Step_0_IterParam_0_njptbin_12_nzbin_7.root", "READ");

    // TFile *treeFile = new TFile("./Data/RM_dw.root", "READ");
    if (!treeFile || treeFile->IsZombie()) {
        std::cout << "Error opening file: Response_Step_0_IterParam_0_njptbin_12_nzbin_7.root" << std::endl;
        return;
    }

    for (int iCent = 0; iCent < 3; iCent++) {
        rurResponse[iCent] = *(RooUnfoldResponse *)treeFile->Get(Form("Resp1D_%i", iCent));
        rurResponse2D[iCent][0] = *(RooUnfoldResponse *)treeFile->Get(Form("RespWide_%i", iCent));
      //  PrintResponseMissFake(rurResponse[iCent], Form("1D_cent%d", iCent));
    }


    treeFile->Close();
    // exit(0);
}
void ApplyOverrideMacro(const char* overrideMacro)
{
    if (!overrideMacro || overrideMacro[0] == '\0') return;  // nic nezadané => defaulty z config.h

    // existuje soubor?
    if (gSystem->AccessPathName(overrideMacro)) {
        cout << "[override] file not found: " << overrideMacro << " -> using defaults from config.h" << endl;
        return;
    }

    cout << "[override] applying: " << overrideMacro << endl;

    // Spustí makro jako skript: v něm jen přepíšeš globální proměnné/vektory
    gROOT->Macro(overrideMacro);
}

void Machine(bool _fonllJet = true, bool _CutOfNegative = true, double _minJetPtRecoCut = -30,
             int _savedIter = 4,
             char *InputFileIn = 0,
             const char *OutputFile = "Output",
             double _minPtD0Cut = 1, double _maxPtD0Cut = 10, //min not applied
             const char *OverrideMacro = "",
             const char *ScanDir = "DefaultScanDir",
             int usePriorShapeWeighting = 0, //0 none, 1X JetpT, 2X second variable // X = 0 +20%; X = 1 -20%
             int systematicSPlot = 0)
{
    ApplyOverrideMacro(OverrideMacro);

    TString InputFile;     
    if (InputFileIn) {
        InputFile = InputFileIn;
    } else {
        InputFile = RealJetsFileData;
    }

    // _systematicSPlot = 0; // nominal: Gaussian + Exponential
    // _systematicSPlot = 1; // background variation: Gaussian + Chebychev 2nd order
    // _systematicSPlot = 2; // signal variation: Double Gaussian + Exponential
    // _systematicSPlot = 3  // Student-t signal + Exponential background
    // _systematicSPlot = 4; // narrower fit range
    // _systematicSPlot = 5; // wider fit range
    // _systematicSPlot = 6; // keep negative bins (no cut of negative values)

    TString sys = "";
    _systematicSPlot = systematicSPlot;
    if (_systematicSPlot == 1) sys = "_ChebysevBkg";
    if (_systematicSPlot == 2) sys = "_DoubleGaussSgn";
    if (_systematicSPlot == 3) sys = "_StudentTSignal";
    if (_systematicSPlot == 4) sys = "_NarrowFit";
    if (_systematicSPlot == 5) sys = "_WideFit";
    if (_systematicSPlot == 6) sys = "_KeepNegative"; 
    //Paper uncertainties
    if (_systematicSPlot == 7) sys = "_paperTPCTrackUp";
    if (_systematicSPlot == 8) sys = "_paperTPCTrackDown";
    if (_systematicSPlot == 9) sys = "_paperPIDUp";
    if (_systematicSPlot == 10) sys = "_paperPIDDown";
    if (_systematicSPlot == 11) sys = "_paperSingleTrackPtUp";
    if (_systematicSPlot == 12) sys = "_paperSingleTrackPtDown";
    if (_systematicSPlot == 13) sys = "_paperTopoEffUp";
    if (_systematicSPlot == 14) sys = "_paperTopoEffDown";
    if (_systematicSPlot == 15) sys = "_paperDoubleCountingUp";
    if (_systematicSPlot == 16) sys = "_paperDoubleCountingDown";
    if (_systematicSPlot == 17) sys = "_paperVertexCorrUp";
    if (_systematicSPlot == 18) sys = "_paperVertexCorrDown";
    if (_systematicSPlot == 19) sys = "_paperSecondaryTrackUp";
    if (_systematicSPlot == 20) sys = "_paperSecondaryTrackDown";
    //Jet reconstruction
    if (_systematicSPlot == 21){
        sys = "_jetRecEfficiency";
        RealJetsFileData = "./Data/Output_real_final_01022026.root";
        InputFile = RealJetsFileData;
        McJetsFileData = "./Data/Output_MC_MidLow_trackEff_05052026.root"; 
    }
    if (_systematicSPlot == 22){
        sys = "_jetnHitsFit13";
        RealJetsFileData = "./Data/Output_real_final_01022026.root";
        InputFile = RealJetsFileData;
        McJetsFileData = "./Data/Output_MC_MidLow_nHitsFit13_06052026.root"; 
    }
    if (_systematicSPlot == 23){
        sys = "_jetnHitsFit17";
        RealJetsFileData = "./Data/Output_real_final_01022026.root";
        InputFile = RealJetsFileData;
        McJetsFileData = "./Data/Output_MC_MidLow_nHitsFit17_06052026.root"; 
    }
    if (_systematicSPlot == 24){
        sys = "_jetKTDrop";
        RealJetsFileData = "./Data/Output_real_final_01022026.root";
        InputFile = RealJetsFileData;
        McJetsFileData = "./Data/Output_MC_MidLow_kTDrop_07052026.root"; 
    }
    if (_systematicSPlot == 25){
        sys = "_jetDCA2_8";
        RealJetsFileData = "./Data/Output_real_final_01022026.root";
        InputFile = RealJetsFileData;
        McJetsFileData = "./Data/Output_MC_MidLow_DCA2_8_09052026.root"; 
    }
    if (_systematicSPlot == 26){
        sys = "_jetDCA3_2";
        RealJetsFileData = "./Data/Output_real_final_01022026.root";
        InputFile = RealJetsFileData;
        McJetsFileData = "./Data/Output_MC_MidLow_DCA3_2_10052026.root"; 
    }
    if (_systematicSPlot == 27){
        sys = "_jetHadronicCorr";
        RealJetsFileData = "./Data/Output_real_final_01022026.root";
        InputFile = RealJetsFileData;
        McJetsFileData = "./Data/Output_MC_MidLow_hadrCorr_08052026.root"; 
    }

    _sys=sys;

    if (_systematicSPlot >= 21 && _systematicSPlot <= 27) {
    UseCachedRM = false;
    FillStandardRM = true;

    cout << "[Jet reco systematic]" << endl;
    cout << "  sys              = " << sys << endl;
    cout << "  RealJetsFileData = " << RealJetsFileData << endl;
    cout << "  McJetsFileData   = " << McJetsFileData << endl;
    cout << "  InputFile        = " << InputFile << endl;
    cout << "  UseCachedRM      = " << UseCachedRM << endl;
}

    FONLLjet = _fonllJet; //true = FONLL jet pT, false = FONLL D0 pT
    CutOfNegative = _CutOfNegative; //negative values in histograms will be set to 0
    minJetPtRecoCut = _minJetPtRecoCut;
    savedIter = _savedIter; //number of iterations to save
    GivenIter = _savedIter; //number of iterations to use for unfolding, can be different from savedIter, but usually it's the same
    _usePriorShapeWeighting = usePriorShapeWeighting;

 runId = gSystem->BaseName(OutputFile);   // "r000005"
 OpenStabilityFile(ScanDir);
outputFileMachine = OutputFile;
    outputFile = "Output";


    if( usePriorShapeWeighting != 0 && UseCachedRM) {
        cout << "Using prior shape weighting with cached RM. This doesnt work as expected!" << endl;
        exit(1);
    }

    std::string outputF;
    maxPtD0Cut = _maxPtD0Cut;
    SetMinD0Pt = minPtD0Cut;
    //Output PDF file

    TString _outPdf = "./OutputPdf/" + TString(OutputFile) + TString(runId) + TString(Method)+"Rcp.pdf";
    outPdf = _outPdf;

    gSystem->Load("libRooUnfold");

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);

    //Load prepared distributions
    LoadDataFactors();
    LoadPriorShapeWeights();

    //Histogram initialization
    HistogramInit();

    //Load data to rescale the MC
    if(!UseCachedRM) LoadDataRealParams();
    //Load Mc data to fill RM and input data for closure test
    if (UseCachedRM) {
        cout << "Loading RM from cache histograms..." << endl;
        LoadDataCache();
    } else {
        cout << "Loading RM from TTree..." << endl;
        LoadDataMC();
        //LoadDataNeilRM();
    }

  
                    
    //SaveFine();

    //Recalculate the centrality ranges
    std::vector <std::vector<int>> centralityRange(centrality.size());
    std::transform(centrality.begin(), centrality.end(), centralityRange.begin(), CentrRangeTransf);
    Centr(centralityRange, InputFile);

    //Load real data
    if (!ClosureTest) LoadDataReal();
   //if (!ClosureTest){LoadFiles(); LoadDataRealNeil();};


    //Main canvas for unfolding and RCP
    TCanvas *can = new TCanvas("Main Canvas", "Unfolding and RCP", 1200, 400);

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.055);

    can->Clear();

    can->SetCanvasSize(900, 600);
    can->SaveAs(outPdf + "[");

    can->Clear();
    DrawAllBinning();
    can->SaveAs(outPdf);
    can->Clear();

    can->Divide(3, 2);
    std::vector<TGraphErrors*> gD0ShiftedPt(nCentralityBins);

    if (!UseCachedRM){    
        for (int i = 0; i < nCentralityBins; i++) {

        can->cd(i + 1);
        //left margin
        gPad->SetLeftMargin(0.15);
        gPad->SetLogy();
        hMeasuredD0MesonPt[i]->Draw();

        double BranchingRatio = 0.0389;

        for (int j = 1; j <= hMeasuredD0MesonPt[i]->GetNbinsX(); j++) {
            hMeasuredD0MesonPt[i]->SetBinContent(j, 1.*hMeasuredD0MesonPt[i]->GetBinContent(j) / hMeasuredD0MesonPt[i]->GetBinCenter(j)); // /pT
        } // pT

        hMeasuredD0MesonPt[i]->Scale(1. / NumberOfWEvents[i]); // /(Nev)
        hMeasuredD0MesonPt[i]->Scale(1. / (2 * TMath::Pi())); // /Delta 2pi
        //hMeasuredD0MesonPt[i]->Scale(1. / 2.); // /Delta y
        hMeasuredD0MesonPt[i]->Scale(1. / 1.2); // /Delta y
        hMeasuredD0MesonPt[i]->Scale(1. / 2.); // (D0+antiD0/2)
        hMeasuredD0MesonPt[i]->Scale(1. / BranchingRatio); // /D0
        NormalizeByBinWidth(hMeasuredD0MesonPt[i],2000); // /delta pT

        hMeasuredD0MesonPt[i]->SetLineColor(kRed);
        hMeasuredD0MesonPt[i]->SetMarkerColor(kRed);
        hMeasuredD0MesonPt[i]->GetXaxis()->SetTitle("p_{T}^{D^{0}} [GeV/c]");
        //set range
        hMeasuredD0MesonPt[i]->GetYaxis()->SetRangeUser(1e-8,1);

        hMeasuredD0MesonPt[i]->GetYaxis()->SetTitle("1/(N_{ev} B.R.) d^{2}N/(2#pidp_{T} dy)");
        //offset
        hMeasuredD0MesonPt[i]->GetYaxis()->SetTitleOffset(1.3);

        //fPt[i]->Draw("same");
        hPaperD0PtBeforeShift[i]->Draw("same");
        can->cd(4 + i);

        gPad->SetLeftMargin(0.15);

        hMeasuredD0MesonPtRatio[i] = (TH1D *) hMeasuredD0MesonPt[i]->Clone(Form("hMeasuredD0MesonPtRatio_%d", i));
        hMeasuredD0MesonPtRatio[i]->SetTitle("Ratio to FONLL");
        //hMeasuredD0MesonPtRatio[i]->Divide(fPt[i]);
        TH1D graphToHist = TH1D("converted", "converted", BinyVl.size() - 1, &BinyVl[0]);
        for (int j = 0; j < hPaperD0PtBeforeShift[i]->GetN(); ++j) {
            double x, y;
            hPaperD0PtBeforeShift[i]->GetPoint(j, x, y);
            //get uncertainties
            double ey = hPaperD0PtBeforeShift[i]->GetErrorYhigh(j);
            int bin = graphToHist.FindBin(x);
            graphToHist.SetBinContent(bin, y);
            graphToHist.SetBinError(bin, ey);
        }
        //hMeasuredD0MesonPtRatio[i]->Divide(&graphToHist, 1.0, 1.0, "B");
        hMeasuredD0MesonPtRatio[i]->Divide(hMeasuredD0MesonPt[i], &graphToHist, 1.0, 1.0, "B");

        hMeasuredD0MesonPtRatio[i]->SetLineColor(kBlue);
        hMeasuredD0MesonPtRatio[i]->SetMarkerColor(kBlue);
        hMeasuredD0MesonPtRatio[i]->Draw("same");
        hMeasuredD0MesonPtRatio[i]->GetYaxis()->SetTitle("Paper / Measured D^{0} ratio");

        hMeasuredD0MesonPtRatio[i]->GetYaxis()->SetRangeUser(0.5, 1.5);
        hMeasuredD0MesonPtRatio[i]->GetYaxis()->SetTitleOffset(1.3);
        DrawLineOne2(hMeasuredD0MesonPtRatio[i]->GetXaxis()->GetXmin(),
                     hMeasuredD0MesonPtRatio[i]->GetXaxis()->GetXmax());

    }

    can->SaveAs(outPdf);

    can->Clear();
    can->SetCanvasSize(1200, 800);

    can->Divide(3, 2);

    for (int i = 0; i < nCentralityBins; i++) {

        can->cd(i + 1);
        gPad->SetLeftMargin(0.15);

        //logscale
        gPad->SetLogy();
        D0MesonPtMcTrue[i]->Draw();
        D0MesonPtMcTrue[i]->Scale(1. / D0MesonPtMcTrue[i]->Integral());
        D0MesonPtMcTrue[i]->GetXaxis()->SetTitle("p_{T}^{D^{0}} [GeV/c]");
        D0MesonPtMcTrue[i]->GetYaxis()->SetTitle("Counts");
        D0MesonPtMcTrue[i]->SetLineColor(kRed);
        D0MesonPtMcTrue[i]->SetMarkerColor(kRed);
        D0MesonPtMcReco[i]->Draw("same");

        D0MesonPtMcReco[i]->Scale(1. / D0MesonPtMcReco[i]->Integral());
        D0MesonPtMcReco[i]->SetLineColor(kBlue);
        D0MesonPtMcReco[i]->SetMarkerColor(kBlue);

        D0MesonPtReal[i]->Draw("same");
        D0MesonPtReal[i]->SetLineColor(kGreen);
        D0MesonPtReal[i]->SetMarkerColor(kGreen);
        D0MesonPtReal[i]->Scale(1. / D0MesonPtReal[i]->Integral());
        D0MesonPtMcTrue[i]->GetYaxis()->SetRangeUser(0.00001, 1);

        TLegend *leg = new TLegend(0.6, 0.6, 0.9, 0.9);
        leg->AddEntry(D0MesonPtMcTrue[i], "MC True", "l");

        pureFonllD0meson->SetRange(1, 10);
// Spočítáš integrál přes rozsah, který tě zajímá (velmi nižší tolerance pro spline)
        double integral = pureFonllD0meson->Integral(1, 10, 1e-2);

// Normalizuješ
        if (integral > 0)
            pureFonllD0meson->SetNormalized(true);  // ROOT 6+
        // nebo ručně přepočítat parametr scale:

        pureFonllD0meson->SetParameter(0, pureFonllD0meson->GetParameter(0) / integral);


// A pak kreslit
        pureFonllD0meson->SetLineColor(kBlack);
        pureFonllD0meson->Draw("same");
        //transparency
        leg->SetFillColorAlpha(0, 0);
        leg->SetBorderSize(0);
        leg->AddEntry(D0MesonPtMcReco[i], "MC Reco", "l");
        leg->AddEntry(D0MesonPtReal[i], "Real Data", "l");
        leg->AddEntry(pureFonllD0meson, "D^{0} FONLL", "l");

        leg->Draw();


        can->cd(4+i);
        gPad->SetLeftMargin(0.15);
        gPad->SetLogy();

        D0JetPtMcTrue[i]->Draw();
        D0JetPtMcTrue[i]->SetLineColor(kRed);
        D0JetPtMcTrue[i]->SetMarkerColor(kRed);
        D0JetPtMcTrue[i]->Scale(1. / D0JetPtMcTrue[i]->Integral());

        TLegend *leg2 = new TLegend(0.6, 0.6, 0.9, 0.9);
        leg2->AddEntry(D0JetPtMcTrue[i], "MC True", "l");
        leg2->SetFillColorAlpha(0, 0);
        leg2->SetBorderSize(0);
      //  leg2->AddEntry(D0JetPtMcTrue[i], "MC Reco", "l");
        //leg2->AddEntry(D0MesonPtReal[i], "Real Data", "l");


        double integral2 = pureFonll->Integral(1, 30, 1e-2);  // Velmi nižší tolerance pro spline

// Normalizuješ
       if (integral2 > 0)
        pureFonll->SetNormalized(true);  // ROOT 6+
        // nebo ručně přepočítat parametr scale:

        pureFonll->SetParameter(0, pureFonll->GetParameter(0) / integral2);

        pureFonll->Draw("same");
        leg2->AddEntry(pureFonll, "c-quark FONLL", "l");
        leg2->Draw();
    }
    can->SaveAs(outPdf);
    can->Clear();
    can->SetCanvasSize(1200, 400);

    can->Divide(3, 1);


    for (int iCent = 0; iCent < nCentralityBins; iCent++) {

        can->cd(iCent + 1);
        gPad->SetLeftMargin(0.15);

        KinEff1D[iCent].Draw("E1");
        gPad->Update();  // důležité

// správné xmin/xmax: první a poslední HRANA binů
        const double xmin = ptMcBinsVecCustom[iCent].front();
        const double xmax = ptMcBinsVecCustom[iCent].back();
        const double xminCut = ptRecoBinsVec[iCent].front();
        const double xmaxCut = ptRecoBinsVec[iCent].back();

        auto gr = KinEff1D[iCent].GetPaintedGraph();
        gr->GetYaxis()->SetRangeUser(0.0, 1.2);

// KLÍČOVÉ: pro TGraph použij SetLimits (ne SetRangeUser)
        gr->GetXaxis()->SetLimits(xmin, xmax);

        gr->GetYaxis()->SetTitleOffset(1.3);

        DrawLineOne2(xmin, xmax);

        gPad->Modified();
        gPad->Update();

        TLegend *leg = new TLegend(0.2, 0.1, 0.7, 0.4);
        leg->AddEntry(&KinEff1D[iCent], Form("#frac{N(|#eta_{Jet}^{reco}|<0.6 && %.f < p_{T,Jet}^{reco} < %.f GeV/c)}{N(|#eta_{Jet}^{reco}|<0.6)}",xminCut,xmaxCut), "lp");
        leg->SetFillColorAlpha(0, 0);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.035);
        leg->Draw();

    }

    double trueJetVarZMin[3] = {zMcBinsVecCustom[0].front(), zMcBinsVecCustom[1].front(), zMcBinsVecCustom[2].front()};
    double trueJetVarZMax[3] = {zMcBinsVecCustom[0].back(), zMcBinsVecCustom[1].back(), zMcBinsVecCustom[2].back()};
    double recoJetVarZMin[3] = {zRecoBinsVec[0].front(), zRecoBinsVec[1].front(), zRecoBinsVec[2].front()};
    double recoJetVarZMax[3] = {zRecoBinsVec[0].back(), zRecoBinsVec[1].back(), zRecoBinsVec[2].back()};

    can->SaveAs(outPdf);
    can->Clear();
    DrawKinEff2D(KinEff2DpTZ[2], KinEff2DpTZ[0], KinEff2DpTZ[1], KinEff2DpTZCut,KinEff2DpTZZCut, trueJetVarZMin, trueJetVarZMax, recoJetVarZMin,recoJetVarZMax, can);

    double trueJetVarMin[6][3];
    double trueJetVarMax[6][3];
    double recoJetVarMin[6][3];
    double recoJetVarMax[6][3];
    for (int iCent = 0; iCent < nCentralityBins; iCent++) {
        for (int iAng = 0; iAng < nAngularities; iAng++) {
            trueJetVarMin[iAng][iCent] = angMcBinsVecCustom[iCent][iAng][0];
            trueJetVarMax[iAng][iCent] = angMcBinsVecCustom[iCent][iAng][angMcBinsVecCustom[iCent][iAng].size() - 1];
            recoJetVarMin[iAng][iCent] = angRecoBinsVec[iCent][iAng][0];
            recoJetVarMax[iAng][iCent] = angRecoBinsVec[iCent][iAng][angRecoBinsVec[iCent][iAng].size() - 1];
        }
    }

    for (int iAng = 0; iAng < nAngularities; iAng++) {

         DrawKinEff2D(KinEff2DAng[iAng][2], KinEff2DAng[iAng][0], KinEff2DAng[iAng][1], KinEff2DAngCut[iAng],KinEff2DAngPtCut[iAng], trueJetVarMin[iAng] , trueJetVarMax[iAng],recoJetVarMin[iAng],recoJetVarMax[iAng], can);
    }



    can->Clear();

    can->SetCanvasSize(1200, 800);

    can ->Divide(3, 2);

    for (int iCent = 0; iCent < nCentralityBins; iCent++) {
        can->cd(iCent + 1);
        gPad->SetLeftMargin(0.15);

        KinEffEta[iCent].Draw();
        gPad->Update();  // DŮLEŽITÉ – vytvoří se interní graf

        auto graph = KinEffEta[iCent].GetPaintedGraph();
        graph->GetYaxis()->SetTitleOffset(1.4); 

        //yaxis range
        gPad->Update();  // důležité!

        KinEffEta[iCent].GetPaintedGraph()->SetMinimum(0.0);
        KinEffEta[iCent].GetPaintedGraph()->SetMaximum(0.2);

        can->cd(iCent + 4);
        gPad->SetLeftMargin(0.15);
        FakeEffEta[iCent].Draw();
        gPad->Update();
        auto graph2 = FakeEffEta[iCent].GetPaintedGraph();
        graph2->GetYaxis()->SetTitleOffset(1.4); 

        gPad->Update();  // důležité!

        FakeEffEta[iCent].GetPaintedGraph()->SetMinimum(0.0);
        FakeEffEta[iCent].GetPaintedGraph()->SetMaximum(0.2);

    }
    can->SaveAs(outPdf);
}
    can->Clear();

    can->SetCanvasSize(1200, 1000);

    //Unfolding
    for (int iCent = 0; iCent < nCentralityBins; iCent++) {

        plotComparison1D(can, iCent);

        if (Unfold2D) {
            plotComparison2D(can, iCent, "z");
                for (int iLambda = 0; iLambda < nAngularities; iLambda++) {
                    plotComparison2D(can, iCent, AngNames[iLambda]);
                }
            }
        }


    plotRcp1D(can);
    plotRcp2D(can);
    plotFinalComp(can, ScanDir);
    if(!UseCachedRM) plotResolution(can);

    can->SaveAs(outPdf+"]");

    //Změním rozměr canvasu
    can->SetCanvasSize(1200, 800);
    can->cd();
    can->Clear();


   // TH1D *hRcpPt_[NSuperIter][nCentralityBins][nCentralityBins][nIter];
    //TH2D *hRcpPTZ_[nCentralityBins][nCentralityBins][nIter];
   // TH2D *hRcpPTang[nCentralityBins][nCentralityBins][nIter][nAngularities];

    can->Clear();

    cout << NumberOfWEvents[0] << " " << NumberOfWEvents[1] << " " << NumberOfWEvents[2] << endl;


}


