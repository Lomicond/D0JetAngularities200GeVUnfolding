#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
#include "TLine.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TNtuple.h"
//TGraphErrors
#include "TGraphErrors.h"
//TPaletteAxis
#include "TPaletteAxis.h"
//TLatex
#include "TLatex.h"
//TGaxis
#include "TGaxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "RooUnfoldResponse.h"
//TObjString
//TTreeReader
#include "TTreeReader.h"
//TKey
#include "TKey.h"
//TProfile
#include "TProfile.h"
#include "TObjString.h"
#include "TEfficiency.h"
#include <fstream>

                    //Centrality ranges
std::vector<std::vector<int>> centrality = {{0, 10}, {10, 40}, {40, 80}}; 		                        //Centrality ranges
std::vector<std::vector<double>> momenta = {{0,1},{1,1.5},{1.5,2},{2,3},{3,5},{5,10}};


	                        //jetPT, z, lam_1_1, lam_1_1.5, lam_1_2, lam_1_3
////Int_t NSuperIter = 1;
Int_t GivenIter = 4; //(3 => 3-1 = 2 v poli, převedeno v kodu)
//const Int_t nRecoBins = 8;
//const Int_t nMcBins = 6;
//int NumOfBin[6] = {nRecoBins,nRecoBins,nRecoBins,nRecoBins,nRecoBins,nRecoBins};
const Int_t DividedMcDataBy = 1;
//bool RewriteWeights = true;

const char *outputFile;
const char *outputFileMachine;

double minPtD0Cut = 1;
double maxPtD0Cut = 10;

TString Method = "ICS"; //ICS vs AREA

bool FillStandardRM = true;
bool FillCacheRM = false; //Be carefoul!!!
//bool LoadResponseMatrix = false;
bool UseCachedRM = false;
TString CacheRMFileName = "./Output/CacheRM.root";

//Histograms
TH1D Unfolded2D_X[7][3];
TH1D Unfolded2D_Y[3];
TH1D Unfolded1D[3];
TH1D Unfolded2DRatio[7][3];
TH1D Unfolded1DRatio[3];
TH1D Unfolded2DRatioVar[7][3];

TH2D* hCacheMatchPt[3];
TH1D* hCacheMissPt[3];
TH1D* hCacheFakePt[3];

THnSparseD* hCacheMatchPtZ[3];
TH2D* hCacheMissPtZ[3];
TH2D* hCacheFakePtZ[3];

THnSparseD* hCacheMatchPtAng[3][6];
TH2D* hCacheMissPtAng[3][6];
TH2D* hCacheFakePtAng[3][6];

//TString McJetsFileData = "./Data/Output_sim_final_01022026_slimfix.root"; //inclusive
//TString McJetsFileData = "./Data/Output_sim_cf_1902206_slim.root"; //inclusive
////TString McJetsFileData = "./Data/Output_sim_konecny_3032026.root"; //inclusive

TString McJetsFileData = "/home/ondrej/Desktop/Pracovna/Koš/testHadrCorrNew28082026.root"; //inclusive

TString NeilFONLL = "./Data/FONLL_Pt_1_30.root";
TString D0SpectraBeforeShift = "./Data/new/D0_Spectra_Run14_HFT_beforePtShift.root"; //inclusive


TString RealJetsFileData = "./Data/Output_real_final_01022026.root"; //inclusive
//TString RealJetsFileData = "./Data/Output_real_cfactor.root";

//TString RealJetsFileData = "./Data/Output_pionmassreal_21092025.root"; //inclusive

TString DoubleCounting = "./Data/DoubleCounting2014.root";

//TString NeilFile = "./Data/SmallSampleForSasha.root";
TString PaperD0Spectrum = "./Data/new/HEPData-ins1711377-v3-D^0_spectra_in_AuAu_collisions.root";


const TString AngNames[6] = {"#lambda^{1}_{1}", "#lambda^{1}_{1.5}", "#lambda^{1}_{2}", "#lambda^{1}_{3}", "#lambda^{1}_{0.5}","p_{T}^{D}"};

const Double_t TruthJetPtMin = 1;
const Double_t TruthJetPtMax = 20;
const Double_t RecoJetPtMin[3] = {1, 1, 1};

const Double_t RecoJetPtMax[3] = {20, 20, 20};

const Double_t TrainToTestRatio = 0.5;
//graphs in createResponseMatrix
bool PearsonCoeff = true;
bool ResponseMatrix = true;
bool BinMigration = true;
bool MissingJets = true;
bool FakeJets = true;

bool deleteOneConstituentJets = false;



//2D unfolding
bool Unfold2D = true;
//SVD
bool SVD = true;
//
bool FlatPrior = false;
//Use superiteration
bool SuperIteration = false;


//Relative to previous iteration
bool UseRelativeP = true;
//ClosureTest
bool ClosureTest = false;
bool UseTheSameSample = false; ///!!!

const Int_t nKterm = 3;
int kterm[nKterm] = {1,2, 3};
int ChosenKterm = 0;
int SVDindex = 0;

//Weihgted prior
bool WeightedPrior = true;
double McPtShift[3] = {0, 0, 0};


const vector<Int_t> PlotIterations = {1,2,3, 4, 5};

//const vector<Int_t> PlotIterations = {1, 2, 3, 4, 5, 10, 15, 20};
const Int_t nIter = PlotIterations.size();

const vector<Int_t> sPlotIterations = {0,1,2};
//const vector<Int_t> sPlotIterations = {0,1,2,3,4,5,10,15,20};//,25,30,40,50,60,70,80,90,100,130,190,200};

Int_t nsIter = sPlotIterations.size();

//RM
bool UseOverflow = false;
double underflowplot = -35.;
double overflowplot = 60.;
const Int_t nCentralityBins = 3;
const Int_t nAngularities = 6;
double NumberOfWEvents[nCentralityBins] = {1.0, 1.0, 1.0};
vector<Double_t> ptMcBinsVec[nCentralityBins];
vector<Double_t> pTMcBinsVecTest[nCentralityBins];
vector<Double_t> zMcBinsVec[nCentralityBins];
vector<Double_t> angularityMcBinsVec[nCentralityBins][nAngularities];

bool useCustomPtMcBins = true;
bool CompOfDifferentPt = true;
TF1 *levy[7];
TF1 *levy2[3];

// https://www.star.bnl.gov/protected/lfsupc/tdrk/Centrality/Run19AuAu200/top20_tables/table_Ncoll_vs_centrality_systematicerror.txt
TString centralityTitles[nCentralityBins] = {"0-10%", "10-40%", "40-80%"};
TString centralityNames[nCentralityBins] = {"0_10", "10_40", "40_80"};
TString RcpTitles[3] = {"0-10/40-80", "0-10/10-40", "10-40/40-80"};
TString RcpTitles2[3] = {"0-10", "10-40", "40-80"};

TH1D *hMcWeight[3];
TH1D *hMcWeightJetPtReco[3];
TH1D *jetPtCheckScaled[3];
TH1D *jetPtCheckScaled2[3];
TH1D *jetPtRecoCheckScaled[3];
TH1D *jetPtRecoCheckScaled2[3];
TH1D *jetPtRecoCheck[3];

TH1D *D0MesonPtReal[3];
TH1D *D0MesonPtMcReco[3];
TH1D *D0MesonPtMcTrue[3];
TH1D *D0JetPtMcReco[3];
TH1D *D0JetPtMcTrue[3];
TH2D *JetPtZMc[3];

TH2F *Efficiency2D[5];
TEfficiency *Efficiency2D_TEff[5];
TH1D *Efficiency1D[5];
TH1F *Efficiency1DPaper[5];

TString outPdf;
TString outRoot = "Output";

TH1D hVarReduced_Y_ang[7][3][6];

TString _sys;
//static std::ofstream fout("stability.txt", std::ios::app);
TString runId;
//std::ofstream fout("scan_pTD/stability.tsv", std::ios::app);
std::ofstream fout;
////
//1D
TH1D hRealData[3];
RooUnfoldResponse rurResponse[3]= {
RooUnfoldResponse("dummy1D0", "dummy1D0"),
RooUnfoldResponse("dummy1D1", "dummy1D1"),
RooUnfoldResponse("dummy1D2", "dummy1D2")
};
TH1D hUnfoldedPt[3][15]; //NIter
TH1D hBackfoldedPt[3][15]; //NIter
////
//Cent + Variables
TH2D hRealData2D[3][7];
TH2D hMcReco2D[3][7];
TH2D hMcTrue2D[3][7];
TH2D hUnfolded2D[3][7][15]; //NIter
//TH1D hRcpPtVarReduced_Y[numOfComb][nVar][nIter];
TH1D hRcpPtVarReduced_Y[3][7][15]; //NIter
TH1D hVarReduced_Y[3][7][15]; //NIter

THnSparseF* hResp4D[3];
THnSparseF* hResp4DWeighted[3];

//Simple_splot
Int_t _systematicSPlot = 0;

RooUnfoldResponse rurResponse2D[3][8] = {
    {RooUnfoldResponse("dummy0", "dummy0"), RooUnfoldResponse("dummy1", "dummy1"), RooUnfoldResponse("dummy2", "dummy2"),
     RooUnfoldResponse("dummy3", "dummy3"), RooUnfoldResponse("dummy4", "dummy4"), RooUnfoldResponse("dummy5", "dummy5"),
     RooUnfoldResponse("dummy6", "dummy6"), RooUnfoldResponse("dummy7", "dummy7")},
    {RooUnfoldResponse("dummy0_2", "dummy0_2"), RooUnfoldResponse("dummy1_2", "dummy1_2"), RooUnfoldResponse("dummy2_2", "dummy2_2"),
     RooUnfoldResponse("dummy3_2", "dummy3_2"), RooUnfoldResponse("dummy4_2", "dummy4_2"), RooUnfoldResponse("dummy5_2", "dummy5_2"),
     RooUnfoldResponse("dummy6_2", "dummy6_2"), RooUnfoldResponse("dummy7_2", "dummy7_2")},
    {RooUnfoldResponse("dummy0_3", "dummy0_3"), RooUnfoldResponse("dummy1_3", "dummy1_3"), RooUnfoldResponse("dummy2_3", "dummy2_3"),
     RooUnfoldResponse("dummy3_3", "dummy3_3"), RooUnfoldResponse("dummy4_3", "dummy4_3"), RooUnfoldResponse("dummy5_3", "dummy5_3"),
     RooUnfoldResponse("dummy6_3", "dummy6_3"), RooUnfoldResponse("dummy7_3", "dummy7_3")}
};
RooUnfoldResponse rurResponse2DTest[3] = { RooUnfoldResponse("dummy2D0", "dummy2D0"), RooUnfoldResponse("dummy2D1", "dummy2D1"), RooUnfoldResponse("dummy2D2", "dummy2D2")};
RooUnfoldResponse rurResponse2DTestW[3] = { RooUnfoldResponse("dummy2D0W", "dummy2D0W"), RooUnfoldResponse("dummy2D1W", "dummy2D1W"), RooUnfoldResponse("dummy2D2W", "dummy2D2W")};

//kin efficiency
TEfficiency KinEff1D[3];
TEfficiency KinEff2DpTZ[3][3];
TEfficiency KinEff2DpTZCut[3];
TEfficiency KinEff2DpTZZCut[3];
TEfficiency KinEff2DAng[6][3][3];
TEfficiency KinEff2DAngCut[6][3];
TEfficiency KinEff2DAngPtCut[6][3];

TEfficiency KinEffEta[3];
TEfficiency FakeEffEta[3];

//resolutions
//cent + variable + cuts
// second dimension needs to hold indices up to 7 (2 + iLambda, iLambda in [0..5])
TH1D *hResVar[3][8][5];

TH1D *HistReal_pTTemp[3];
TH1D *hTruthPtTemp[3];
TH2D *HistReal_pTZTemp[3];
TH2D *hTruthPTZTemp[3];
TH2D *HistReal_ZLam[3];
TH2D *hTruthZLam[3];
TH2D *HistReal_pTAngTemp[3][6];
TH2D *hTruthAngTemp[3][6];

TH2D hRealData2DD0Pt[3][5][7];
TH1D hRealData1D0Pt[3][5];
TGraphAsymmErrors *hPaperD0Pt[3];
TGraphErrors *hPaperD0PtBeforeShift[3];


TH1D *hMeasuredPtRealTest[nCentralityBins];
TGraphAsymmErrors *grPaperYield[7];


TH1D * NeilsFONLL;

//TString angularityTitle[nAngularities] = {"#lambda_{1}^{1}", "#lambda_{1}{3/2}", "#lambda_{1}^{2}", "#lambda_{1}^{3}", "#lambda_{1}^{0.5}"};
//TString angularityNames[nAngularities] = {"lambda_1_1", "lambda_1_1half", "lambda_1_2", "lambda_1_3"};

//vector <Double_t> BinyVl = {0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5 ,6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10};
//vector <Double_t> BinyVl = {0, 0.484, 0.972, 1.468, 1.952, 2.468, 2.952, 3.768, 4.952, 5.788, 7.392, 9.888};
vector <Double_t> BinyVl = {0,0.5,1,1.5,2,2.5,3,4,5,6,8,10};



vector<Double_t> ptMcBinsVecCustom[nCentralityBins] = {
        //Neil's
/*
        {1, 2, 3, 4, 5, 7, 9,11,13, 15, 20,25},
        {1, 2, 3, 4, 5, 7, 9,11,13, 15, 20,25},
        {1, 2, 3, 4, 5, 7, 9,11,13, 15, 20,25}*/
       /* {1, 2, 3, 4, 5, 7, 9,11,13, 15, 20,25},
        {1, 2, 3, 4, 5, 7, 9,11,13, 15, 20,25},
        {1, 2, 3, 4, 5, 7, 9,11,13, 15, 20,25}*/

        {1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 20},
        {1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 20},
        {1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 20}

    //1,1.5,2,2.5,3,3.5,4
 };

vector<Double_t> zMcBinsVecCustom[nCentralityBins] = {
        {0,0.2,0.4,0.6,0.7,0.8,0.9,1.01},
        {0,0.2,0.4,0.6,0.7,0.8,0.9,1.01},
        {0,0.2,0.4,0.6,0.7,0.8,0.9,1.01}
};

bool useCustomPtRecoBins = true;


vector<Double_t> angMcBinsVecCustom[nCentralityBins][nAngularities] = {
/*
        {//Centrality 0-10%
                {0,0.1,0.2,0.3,0.4,0.5, 0.7},
                {0,0.1,0.2,0.3,0.4,0.5, 0.7},
                {0,0.1,0.2,0.3,0.4,0.5, 0.7},
                {0,0.1,0.2,0.3,0.4,0.5, 0.75,1},
                {0,0.1,0.2,0.3,0.4,0.5, 0.7},
                {0,0.1,0.2,0.3,0.4,0.5, 0.7,0.9},

        },
        {//Centrality 10-40%
                {0,0.1,0.2,0.3,0.4,0.5, 0.7},
                {0,0.1,0.2,0.3,0.4,0.5, 0.7},
                {0,0.1,0.2,0.3,0.4,0.5, 0.7},
                {0,0.1,0.2,0.3,0.4,0.5, 0.75,1},
                {0,0.1,0.2,0.3,0.4,0.5, 0.7},
                {0,0.1,0.2,0.3,0.4,0.5, 0.7,0.9},

        },
        {//Centrality 40-80%
                {0,0.1,0.2,0.3,0.4,0.5, 0.7},
                {0,0.1,0.2,0.3,0.4,0.5, 0.7},
                {0,0.1,0.2,0.3,0.4,0.5, 0.7},
                {0,0.1,0.2,0.3,0.4,0.5, 0.75,1},
                {0,0.1,0.2,0.3,0.4,0.5, 0.7},
                {0,0.1,0.2,0.3,0.4,0.5, 0.7,0.9},
        }*/

        {//Centrality 0-10%
                {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 1}, //l11
                {0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.3, 1}, //l11.5
                {0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 1}, //l12
                {0, 0.0125, 0.025, 0.0375, 0.05, 0.075, 0.1, 0.15, 1},
                {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1}, //l10.5
                {0, 0.3, 0.5, 0.65, 0.75, 0.85, 1.01} //pTD

        },
        {//Centrality 10-40%
                {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 1}, //l11
                {0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.3, 1}, //l11.5
                {0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 1}, //l12
                {0, 0.0125, 0.025, 0.0375, 0.05, 0.075, 0.1, 0.15, 1},
                {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1}, //l10.5
                {0, 0.3, 0.5, 0.65, 0.75, 0.85, 1.01} //pTD
        },
        {//Centrality 40-80%
                {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 1}, //l11
                {0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.3, 1}, //l11.5
                {0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 1}, //l12
                {0, 0.0125, 0.025, 0.0375, 0.05, 0.075, 0.1, 0.15, 1},
                {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1}, //l10.5
                {0, 0.3, 0.5, 0.65, 0.75, 0.85, 1.01} //pTD
        }
        
/*
        {//Centrality 0-10%
                {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.05},
                {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.05},
                {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.05},
                {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.05},
                {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.05},
                {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.05},

        },
        {//Centrality 0-10%
                {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.05},
                {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.05},
                {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.05},
                {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.05},
                {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.05},
                {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.05},

        },
        {//Centrality 0-10%
                {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.05},
                {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.05},
                {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.05},
                {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.05},
                {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.05},
                {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.05},

        },*/
        
};





vector<Double_t> ptRecoBinsVec[nCentralityBins] = {

  //best four
/*
        {1, 1.4, 1.8, 2.2, 2.6, 4, 6.4, 8.8, 11.2, 13.6, 16, 20,30},
        {1, 1.4, 1.8, 2.2, 2.6, 4, 5.4, 6.8, 8.2, 9.6, 14, 20,30},
        {1, 1.4, 1.8, 2.2, 2.6, 3, 3.4, 3.8, 4.2, 8.6, 14, 20,30}
*/

//ICS
/*
{ 1, 2, 3,4,5, 6, 9, 12, 15, 17, 21, 30,50},
{ 1, 2, 3,4,5, 6, 9, 12, 15, 17, 21, 30,50},
{ 1,1.5, 2, 2.5, 3,3.5,4,4.5,5, 6, 9,12,15,20,30,50},
*/
/*
{ 1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18, 21, 25},
{ 1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18, 21, 25},
{ 1,1.5,2,2.5, 3,3.5,4,4.5,5, 6, 9,15,20},*/
/*{ 1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 26},
{ 1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18, 21, 25},
{ 1,1.5,2,2.5, 3,3.5,4,4.5,5, 6, 9,15,20}*/

/*
{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 16, 25},
{1, 1.5, 2, 2.5, 3, 3.5, 4.5, 7, 10, 13, 16, 19, 25},
{1, 1.5, 2, 2.5, 3, 4, 5, 7, 9, 11, 13, 18, 25}
*/
{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 16, 25},
{1, 1.5, 2, 2.5, 3, 3.5, 4.5, 7, 10, 13, 16, 19, 25},
{1, 1.5, 2, 2.5, 3, 4, 5, 7, 9, 11, 13, 18, 25}



/*
{0,1,1.5,2,2.5,3,6,9,12,15,18,21,24},
{0,1,1.5,2,2.5,3,6,9,12,15,18,21,24},
{0,1,1.5,2,2.5,3,6,9,12,15,18}
*//*
      {-10,-5,0, 1, 2, 3, 6, 9, 12, 15, 18, 21, 25,30},
        {-5,0, 1, 2, 3, 6, 9, 12, 15, 18, 21, 25,30},
        {-1,0,0.5, 1,1.5, 2, 3, 6, 10, 20,25},
  */
//AREA
/*
        {-15,-10,-5,0, 1, 2, 3, 6, 9, 12, 15, 18, 21, 25,30},
        {-15,-10,-5,0, 1, 2, 3, 6, 9, 12, 15, 18, 21, 25,30},
        {-15,-10,-5,0, 1, 2, 3, 6, 9, 12, 15, 18, 21, 25,30},
  
*/
     /*
            {3, 6, 9, 12, 15, 18, 21, 25,30},
        {3, 6, 9, 12, 15, 18, 21, 25,30},
        {3, 6, 9, 12, 15, 18, 21, 25,30},*/
/*
            {-9.608856201, -3.215456247, -1.229485631, 0.2109038979, 1.479412079, 2.835394859, 4.332457066, 6.057146549, 8.205479622, 11.30413914,16, 21.65945053,50},
        {-6.220077991, -1.574356556, -0.2933921516, 0.6501793265, 1.555423975, 2.5460639, 3.646336079, 4.945100307, 6.599853992, 9.087589264,14, 18.13246155,25,50},
        {-1.638304114, 0.1897594035, 0.7939822674, 1.278977633, 1.7433635, 2.229038239, 2.77445507, 3.428291559, 4.282679558, 5.633323669, 11.00477123,20,50},
*/
    };
/*
vector<Double_t> DG[3] = {
    {0.5, 0.461538, 0.522695, 0.53076, 0.520245, 0.483003, 0.519231, 0.514634, 0.413534, 0.556886, 0.536304, 0.537415},
    {0.0769231, 0.428025, 0.45825, 0.459542, 0.450811, 0.445956, 0.479452, 0.442177, 0.492647, 0.45768, 0.557971},
    {0, 0.48337, 0.537399, 0.461538, 0.340909, 0.857143}
};*/


vector<Double_t> zRecoBinsVec[nCentralityBins] = {
/*
        {0,0.3,0.35,0.4,0.5,0.6,0.8,0.99,1.01},
        {0,0.3,0.35,0.4,0.5,0.6,0.8,0.99,1.01},
        {0,0.3,0.35,0.4,0.5,0.6,0.8,0.99,1.01}
*/

/*
        {-100,0,0.3,0.35,0.4,0.5,0.6,0.8,1,100,1000},
        {-100,0,0.3,0.35,0.4,0.5,0.6,0.8,1,100,1000},
        {-100,0,0.3,0.35,0.4,0.5,0.6,0.8,1,100,1000}
*/
/*
        {-1000, -500, -100, -50, -20, -10, -5, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 5, 10, 20, 50, 100, 500, 1000},
        {-1000, -500, -100, -50, -20, -10, -5, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 5, 10, 20, 50, 100, 500, 1000},
        {-1000, -500, -100, -50, -20, -10, -5, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 5, 10, 20, 50, 100, 500, 1000}
*/
//ICS

{0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01},
{0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01},
{0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01}

//AREA
/*
      {-1000, -500, -100, -50, -20, -10, -5, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 5, 10, 20, 50, 100, 500, 1000},
        {-1000, -500, -100, -50, -20, -10, -5, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 5, 10, 20, 50, 100, 500, 1000},
        {-1000, -500, -100, -50, -20, -10, -5, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 5, 10, 20, 50, 100, 500, 1000}
*/
/*
        {-100, 0, 0.3, 0.35, 0.4, 0.5, 0.6, 0.8, 1, 100, 1000},
        {-100, 0, 0.3, 0.35, 0.4, 0.5, 0.6, 0.8, 1, 100, 1000},
        {-100, 0, 0.3, 0.35, 0.4, 0.5, 0.6, 0.8, 1, 100, 1000}
*/
/*
      {-1000, -100, -20, -5, -2.0, -0.5, 0, 0.5, 1.0, 2.0,  5, 10,  50, 500, 1000},
      {-1000, -100, -20, -5, -2.0, -0.5, 0, 0.5, 1.0, 2.0,  5, 10,  50, 500, 1000},
      {-1000, -100, -20, -5, -2.0, -0.5, 0, 0.5, 1.0, 2.0,  5, 10,  50, 500, 1000},
      */
        /*
        {0,0.2,0.4,0.6,0.7,0.8,0.9,0.99,1.01},
        {0,0.2,0.4,0.6,0.7,0.8,0.9,0.99,1.01},
        {0,0.2,0.4,0.6,0.7,0.8,0.9,0.99,1.01}
*/

};

vector<Double_t> angRecoBinsVec[nCentralityBins][nAngularities] = {
/*
        {//Centrality 0-10%
                {-3,-2,-1,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.01}, //l11
                {-100,0.001, 0.2,0.25,0.3,0.35,0.4, 0.45,0.8,1,100}, //l11.5
                {-100,0.001,0.1,0.15, 0.2,0.25,0.3,0.4,0.8,1,100}, //l12
                {-5,-1,0.001,0.05,0.10, 0.15,0.20,0.25,0.3,1,2}, //l13
                {-100,0.25,0.4,0.45,0.5,0.6,0.7,0.8,1,100},
                {-1000,0.4,0.5,0.6, 0.7,0.75,0.8,0.85,0.95,0.99,1,1000},

        },
        {//Centrality 0-10%
                {-100,0.001, 0.2,0.25,0.3,0.35,0.4, 0.45,0.8,1,100},
                {-100,0.001, 0.2,0.25,0.3,0.35,0.4, 0.45,0.8,1,100},
                {-100,0.001,0.1,0.15, 0.2,0.25,0.3,0.4,0.8,1,100},
                {-5,-1,0.001,0.05,0.10, 0.15,0.20,0.25,0.3,1,2}, //l13
                {-100,0.001,0.25,0.4,0.45,0.5,0.6,0.7,0.8,1,100},
                {-1000,0.2,0.4,0.5,0.6, 0.7,0.75,0.8,0.85,1,1000}
//                {-100,0.2,0.4,0.5,0.6, 0.7,0.75,0.8,0.85,0.95,0.99,1.01,1,100},
        },
        {//Centrality 40-80%
                {-100,0.001, 0.2,0.25,0.3,0.35,0.4, 0.45,0.8,1,100},
                {-100,0.001, 0.2,0.25,0.3,0.35,0.4, 0.45,0.8,1,100},
                {-100,0.001,0.1,0.15, 0.2,0.25,0.3,0.4,0.8,1,100},
                {-15,-10,-5,0.001,0.05,0.10, 0.15,0.20,0.25,0.3,1,2}, //l13
                {-100,0.001,0.25,0.4,0.45,0.5,0.6,0.7,0.8,1,100},
                {-1000,0.2,0.4,0.5,0.6, 0.7,0.75,0.8,0.85,0.95,1,1000},

        }
*/
/*
        {//Centrality 0-10%
                {-3,-0.001,0.001, 0.2,0.25,0.3,0.35,0.4, 0.45,0.8},
                {-3,-0.001,0.001, 0.2,0.25,0.3,0.35,0.4, 0.45,0.8},
                {-3,-0.001,0.001,0.1,0.15, 0.2,0.25,0.3,0.4,0.8,1},
                {-3,-0.001,0.001,0.05,0.10, 0.15,0.20,0.25,0.3},
                {-3,-0.001,0.25,0.4,0.45,0.5,0.6,0.7,0.8},
                {0.2,0.4,0.5,0.6, 0.7,0.75,0.8,0.85,0.95,0.99,1.01},

        },
        {//Centrality 0-10%
                {-3,-0.001,0.001, 0.2,0.25,0.3,0.35,0.4, 0.45,0.8},
                {-3,-0.001,0.001, 0.2,0.25,0.3,0.35,0.4, 0.45,0.8},
                {-3,-0.001,0.001,0.1,0.15, 0.2,0.25,0.3,0.4,0.8},
                {-3,-0.001,0.001,0.05,0.10, 0.15,0.20,0.25,0.3},
                {-3,-0.001,0.001,0.25,0.4,0.45,0.5,0.6,0.7,0.8},
                {0.2,0.4,0.5,0.6, 0.7,0.75,0.8,0.85,0.95,0.99,1.01},

        },
        {//Centrality 40-80%
                {-3,-0.001,0.001,0.001, 0.2,0.25,0.3,0.35,0.4, 0.45,0.8},
                {-3,-0.001,0.001,0.001, 0.2,0.25,0.3,0.35,0.4, 0.45,0.8},
                {-3,-0.001,0.001,0.001,0.1,0.15, 0.2,0.25,0.3,0.4},
                {-3,-0.001,0.001,0.001,0.05,0.10, 0.15,0.20,0.25,0.3},
                {-3,-0.001,0.001,0.25,0.4,0.45,0.5,0.6,0.7,0.8},
                {0.2,0.4,0.5,0.6, 0.7,0.75,0.8,0.85,0.95,0.99,1.01},

        },
*/

//ICS
/*
        {//Centrality 0-10%
                {0,0.05,0.1,0.15,0.2,0.25,0.3,0.4}, //l11
                {0,0.05,0.1,0.15,0.2,0.25,0.3,0.4}, //l11
                {0,0.025,0.05,0.075,0.1,0.15,0.2,0.3,0.4}, //l12
                {0,0.025,0.05,0.075,0.1,0.15,0.2,0.3}, //l12
                {0,0.1,0.2,0.3,0.4,0.5,0.6, 0.7,0.8}, 
                //{0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1.01}, 
{0,0.30, 0.50, 0.65, 0.75, 0.85, 1.01}

        },
        {//Centrality 0-10%
                {0,0.05,0.1,0.15,0.2,0.25,0.3,0.4}, //l11
                {0,0.05,0.1,0.15,0.2,0.25,0.3,0.4}, //l11
                {0,0.025,0.05,0.075,0.1,0.15,0.2,0.3,0.4}, //l12
                {0,0.025,0.05,0.075,0.1,0.15,0.2,0.3}, //l12
                {0,0.1,0.2,0.3,0.4,0.5,0.6, 0.7,0.8}, 
               // {0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1.01}, 
               {0,0.30, 0.50, 0.65, 0.75, 0.85, 1.01}

        },
        {//Centrality 40-80%
                {0,0.05,0.1,0.15,0.2,0.25,0.3,0.4}, //l11
                {0,0.05,0.1,0.15,0.2,0.25,0.3,0.4}, //l11
                {0,0.025,0.05,0.075,0.1,0.15,0.2,0.3,0.4}, //l12
                {0,0.025,0.05,0.075,0.1,0.15,0.2,0.3}, //l12
                {0,0.1,0.2,0.3,0.4,0.5,0.6, 0.7,0.8}, 
               // {0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1.01}, 
               {0,0.30, 0.50, 0.65, 0.75, 0.85, 1.01}


        },
     
        */

        {//Centrality 0-10%
                {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1}, //l11
                {0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.225, 0.325, 0.425, 0.625, 1}, //l11.5
                {0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.25, 0.4}, //l12
                {0, 0.0125, 0.025, 0.0375, 0.05, 0.0625, 0.075, 0.175, 0.475, 1}, //l13
                {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1}, //l10.5
                {0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01} //pTD

        },
        {//Centrality 0-10%
                {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1}, //l11
                {0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.25, 0.35, 0.45, 1}, //l11.5
                {0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.225, 0.325, 0.425, 0.525, 0.625, 0.8}, //l12
                {0, 0.0125, 0.025, 0.0375, 0.05, 0.0625, 0.2625, 0.4625, 0.8}, //l13
                {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1}, 
                {0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01}

        },
        {//Centrality 40-80%
                {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1}, //l11
                {0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.225, 0.325, 0.425, 0.6}, //l11.5
                {0, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.5, 0.7, 1}, //l12
                {0, 0.0125, 0.025, 0.0375, 0.05, 0.0625, 0.2625, 0.4625, 0.8}, //l13
                {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1}, 
                {0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.01}


        },

//area
/*
        {//Centrality 0-10%
                {-21.18506622, -0.6395395994, -0.5046246648, -0.418797195, -0.3511375487, -0.2919230461, -0.2349690795, -0.1770468652, -0.1116337478, -0.02537871152, 21.35414314},
                {-13.62477875, -0.5922381878, -0.4722691476, -0.3958807588, -0.3350563645, -0.2822803557, -0.231336087, -0.1789815277, -0.1207230911, -0.0422225818, 57.80092621},
                {-12.06307983, -0.5578217506, -0.4477470815, -0.3773155212, -0.3217560649, -0.2726988494, -0.2258463949, -0.1778436452, -0.123432219, -0.05130774528, 1.143356919},
                {-13.47061443, -0.5116429329, -0.4126211405, -0.3490117788, -0.299590081, -0.2556153834, -0.2134719491, -0.1700161546, -0.1205456853, -0.05266328529, 1.836455226},
                {-710.6211548, -0.7026255131, -0.5457248688, -0.4457566738, -0.3666979074, -0.2976997197, -0.2319573462, -0.1640931368, -0.08789813519, 0.01284991857, 218.945755},
                {0.1147195771, 0.220450148, 0.2433865815, 0.2626351118, 0.2814859748, 0.3017706573, 0.3254807293, 0.3564285338, 0.404173702, 0.5114161372, 6.806921959},

        },
        {//Centrality 0-10%
        {-35.37883759, -0.6014803648, -0.4147776663, -0.3034613132, -0.2202934623, -0.1495978534, -0.0846452862, -0.02037858963, 0.04972814769, 0.1398187727, 29.20513535},
        {-18.4844799, -0.5919837952, -0.4225452244, -0.3213024139, -0.2458240986, -0.1818170398, -0.1229432672, -0.06470581889, -0.001473854529, 0.08076714724, 73.70556641},
        {-21.36901855, -0.5826363564, -0.4241725504, -0.3307829499, -0.2604392469, -0.2012131214, -0.1467386782, -0.09315121174, -0.03503075242, 0.03999764472, 1.593991756},
        {-42.45981598, -0.5618024468, -0.4176850021, -0.3330792189, -0.2702974081, -0.2171404958, -0.1682988405, -0.1202346608, -0.06787889451, -0.0008985317545, 1.7289536},
        {-1278.444824, -0.5906765461, -0.3831316531, -0.2595567107, -0.1655131727, -0.08582697809, -0.01186532062, 0.06170696393, 0.1422400028, 0.2461509109, 461.817627},
        {0.1004553065, 0.2809029222, 0.3194768131, 0.3535116911, 0.3878872693, 0.4260223806, 0.4716630578, 0.531858027, 0.6245649457, 0.8164517283, 7.23605442},

        },
        {//Centrality 40-80%
        {-23.94165611, -0.536305666, -0.292716682, -0.1624126434, -0.06986171007, 0.003951221704, 0.07076773047, 0.1346552372, 0.2036074251, 0.2906728089, 60.21204376},
        {-8.322198868, -0.5552617311, -0.342225641, -0.2242382169, -0.1399384737, -0.07176597416, -0.01053231116, 0.0475221239, 0.1115479618, 0.1946247071, 27.42455292},
        {-15.2677002, -0.5932496786, -0.3833314478, -0.2695368528, -0.1884868592, -0.1230270267, -0.06508649141, -0.009672820568, 0.04760535806, 0.1226373911, 1.477837682},
        {-30.20765305, -0.594188869, -0.4000194967, -0.2959260345, -0.2233205587, -0.1649475992, -0.1135926694, -0.06481643766, -0.013707228, 0.04964916781, 1.757585645},
        {-653.9255371, -0.3933562934, -0.1376113743, -0.004758452065, 0.0878502205, 0.1651063412, 0.2340113074, 0.301397264, 0.3733809292, 0.4649389982, 412.9133911},
        {0.1103229746, 0.4280516505, 0.4995870292, 0.5627983809, 0.6263791919, 0.6952294111, 0.7760039568, 0.8766238093, 1.020272136, 1.282564282, 5.115240574},


        },*/
        
/*

        {//Centrality 0-10%
                {-21.18506622, -0.6395395994, -0.5046246648, -0.418797195, -0.3511375487, -0.2919230461, -0.2349690795, -0.1770468652, -0.1116337478, -0.02537871152, 21.35414314},
                {-13.62477875, -0.5922381878, -0.4722691476, -0.3958807588, -0.3350563645, -0.2822803557, -0.231336087, -0.1789815277, -0.1207230911, -0.0422225818, 57.80092621},
                {-12.06307983, -0.5578217506, -0.4477470815, -0.3773155212, -0.3217560649, -0.2726988494, -0.2258463949, -0.1778436452, -0.123432219, -0.05130774528, 1.143356919},
                {-13.47061443, -0.5116429329, -0.4126211405, -0.3490117788, -0.299590081, -0.2556153834, -0.2134719491, -0.1700161546, -0.1205456853, -0.05266328529, 1.836455226},
                {-710.6211548, -0.7026255131, -0.5457248688, -0.4457566738, -0.3666979074, -0.2976997197, -0.2319573462, -0.1640931368, -0.08789813519, 0.01284991857, 218.945755},
                {0.1147195771, 0.220450148, 0.2433865815, 0.2626351118, 0.2814859748, 0.3017706573, 0.3254807293, 0.3564285338, 0.404173702, 0.5114161372, 6.806921959},

        },
        {//Centrality 0-10%
        {-35.37883759, -0.6014803648, -0.4147776663, -0.3034613132, -0.2202934623, -0.1495978534, -0.0846452862, -0.02037858963, 0.04972814769, 0.1398187727, 29.20513535},
        {-18.4844799, -0.5919837952, -0.4225452244, -0.3213024139, -0.2458240986, -0.1818170398, -0.1229432672, -0.06470581889, -0.001473854529, 0.08076714724, 73.70556641},
        {-21.36901855, -0.5826363564, -0.4241725504, -0.3307829499, -0.2604392469, -0.2012131214, -0.1467386782, -0.09315121174, -0.03503075242, 0.03999764472, 1.593991756},
        {-42.45981598, -0.5618024468, -0.4176850021, -0.3330792189, -0.2702974081, -0.2171404958, -0.1682988405, -0.1202346608, -0.06787889451, -0.0008985317545, 1.7289536},
        {-1278.444824, -0.5906765461, -0.3831316531, -0.2595567107, -0.1655131727, -0.08582697809, -0.01186532062, 0.06170696393, 0.1422400028, 0.2461509109, 461.817627},
        {0.1004553065, 0.2809029222, 0.3194768131, 0.3535116911, 0.3878872693, 0.4260223806, 0.4716630578, 0.531858027, 0.6245649457, 0.8164517283, 7.23605442},

        },
        {//Centrality 40-80%
        {-23.94165611, -0.536305666, -0.292716682, -0.1624126434, -0.06986171007, 0.003951221704, 0.07076773047, 0.1346552372, 0.2036074251, 0.2906728089, 60.21204376},
        {-8.322198868, -0.5552617311, -0.342225641, -0.2242382169, -0.1399384737, -0.07176597416, -0.01053231116, 0.0475221239, 0.1115479618, 0.1946247071, 27.42455292},
        {-15.2677002, -0.5932496786, -0.3833314478, -0.2695368528, -0.1884868592, -0.1230270267, -0.06508649141, -0.009672820568, 0.04760535806, 0.1226373911, 1.477837682},
        {-30.20765305, -0.594188869, -0.4000194967, -0.2959260345, -0.2233205587, -0.1649475992, -0.1135926694, -0.06481643766, -0.013707228, 0.04964916781, 1.757585645},
        {-653.9255371, -0.3933562934, -0.1376113743, -0.004758452065, 0.0878502205, 0.1651063412, 0.2340113074, 0.301397264, 0.3733809292, 0.4649389982, 412.9133911},
        {0.1103229746, 0.4280516505, 0.4995870292, 0.5627983809, 0.6263791919, 0.6952294111, 0.7760039568, 0.8766238093, 1.020272136, 1.282564282, 5.115240574},


        },
     */
};


TH1D *Test[3];
TH1D* KinEfficiency[3][11];

TF1*fPt[3];
TH1D *gPt[3]
;/*
vector<Double_t> ptMcBinsVecCustom[nCentralityBins] = {
        {1, 2, 3, 4, 5, 6, 7, 8},
        {1, 2, 3, 4, 5, 6, 7, 8},
        {1, 2, 3, 4, 5, 6, 7, 8}
};*/

///////////////////////////////////////
TH1D *nCentralityNumbers = new TH1D("nCentralityNumbers", "nCentralityNumbers", 9, -0.5, 8.5);
TH1D *nCentralityNumbersMC = new TH1D("nCentralityNumbersMC", "nCentralityNumbersMC", 9, -0.5, 8.5);
TH1D *gRefMult = new TH1D("gRefMult", "gRefMult", 200, 0, 1000);
TH1D *gRefMultMc = new TH1D("gRefMultMc", "gRefMultMc", 200, 0, 1000);
TH1D *gRefMultMcCorr = new TH1D("gRefMultMcCorr", "gRefMultMcCorr", 200, 0, 1000);
TH1D *McPTRawD0[3];
TH1D *McPTRawD0Jet[3];
TH2D *McPTRawD0JetD0Meson[3];
TH1D *hMeasuredD0MesonPt[3];
TH1D *hMeasuredD0MesonPtRatio[3];
TH1D* HistD0PT_reweighted2[3];


TH1D * jetPtCheck[3]; //biny podle ptMcBinsVecCustom
TH1D * jetZ[3];
//Jemné binování
TH2D *hRealFine[3][7];

Int_t nbins_[4]    = {250, 150, 250, 150};
Double_t xmin_[4]  = {1, 0, 1, 0};
Double_t xmax_[4]  = {25, 1.5, 25, 1.5};


THnSparseD *hResponseFine4D[3][7];

TF1 * pureFonll;
TF1 * pureFonllD0meson;

//TH1D hRealData[3];
TH1D hRealDataCopy[3];
TH1D hRealDataXXX[3];
const int njpt_bins_wide = 16;
double nbinsjetpt_wide[njpt_bins_wide + 1] = {-10, -5, 0, 1, 2, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 40, 50};

const int nz_bins_wide = 24;
double nbinsz_wide[nz_bins_wide + 1] =  {-1000, -500, -100, -50, -20, -10, -5, -2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5, 5, 10, 20, 50, 100, 500, 1000};
TH2D hRealData2DCopy[3][7];

////////////////////////////////////////
TGraph *gDoub0, *gDoub1, *gDoub2, *gDoub3;
const int NumberOfR1 = 5;
const int NumberOfR2 = 5;
const int NumberOfC = 3;
TH1F *hEff0_Sys, *hEff1_Sys, *hEff2_Sys, *hEff3_Sys, *hEff4_Sys;
TH1F *hEff0, *hEff1, *hEff2, *hEff3, *hEff4;



TH1D *PriorShapeWeights[8][2]; //[jetPt, z, ang0, ang1, ang2, ang3, ang4, ang5] [+20%/-20%]
void LoadPriorShapeWeights(){

    TFile *fileJetPt = TFile::Open("./Data/PriorShape/jetPtTiltWeights.root");
    if (!fileJetPt) {
        cout << "File jetPtTiltWeights.root not found" << endl;
        exit(1);
    }

    PriorShapeWeights[0][0] = (TH1D*)fileJetPt->Get("hWeightFONLLPlus20");
    PriorShapeWeights[0][0]->SetDirectory(0);
    PriorShapeWeights[0][1] = (TH1D*)fileJetPt->Get("hWeightFONLLMinus20");
    PriorShapeWeights[0][1]->SetDirectory(0);
    fileJetPt->Close();

    TFile *fileZ = TFile::Open("./Data/PriorShape/zTiltWeights.root");
    if (!fileZ) {
        cout << "File zTiltWeights.root not found" << endl;
        exit(1);
    }

    PriorShapeWeights[1][0] = (TH1D*)fileZ->Get("hWeightZPlus20");
    PriorShapeWeights[1][0]->SetDirectory(0);
    PriorShapeWeights[1][1] = (TH1D*)fileZ->Get("hWeightZMinus20");
    PriorShapeWeights[1][1]->SetDirectory(0);
    fileZ->Close();

    TFile *fileAng = TFile::Open("./Data/PriorShape/angTiltWeights.root");
    if (!fileAng) {
        cout << "File angTiltWeights.root not found" << endl;
        exit(1);
    }

    for (int i = 0; i < 6; i++) {
        PriorShapeWeights[i+2][0] = (TH1D*)fileAng->Get(Form("hWeightAngPlus20_ang%d", i));
        PriorShapeWeights[i+2][0]->SetDirectory(0);
        PriorShapeWeights[i+2][1] = (TH1D*)fileAng->Get(Form("hWeightAngMinus20_ang%d", i));
        PriorShapeWeights[i+2][1]->SetDirectory(0);
    }
     fileAng->Close();

}
int _usePriorShapeWeighting = 0;

//int usePriorShapeWeighting = 0 //0 none, 1X JetpT, 2X second observable // X = 0 +20%; X = 1 -20%
//iSecondVar = 0 for jetPt, 1 for z, 2-7 for ang0-ang5
double getPriorShapeWeight(int usePriorShapeWeighting, double value, int iVariable){

    int variableIndex = usePriorShapeWeighting / 10; //0-2
    int direction = usePriorShapeWeighting % 10; //0-1
    int variableId = variableIndex - 1;
    int variableType = (iVariable == 0) ? 0 /*JetPt*/ : (iVariable >= 1 && iVariable <= 7) ? 1 /*z or lambda*/ : -1;

    bool isWeightNotUsed = (usePriorShapeWeighting==0);
    bool isWeightNotUsedForVar = (variableId != variableType);

    if (isWeightNotUsedForVar || isWeightNotUsed ) return 1.0; //no weighting

    if (variableIndex < 1 || variableIndex > 2 || direction < 0 || direction > 1) {
        cerr << "Invalid usePriorShapeWeighting = "
             << usePriorShapeWeighting << endl;
        exit(1);
    }

    if (iVariable < 0 || iVariable > 7) {
        cerr << "Invalid iVariable = " << iVariable << endl;
        exit(1);
    }


    double weight = PriorShapeWeights[iVariable][direction]->GetBinContent(PriorShapeWeights[iVariable][direction]->FindBin(value));

    if (weight != weight || weight < 0) { //check for NaN or negative values
        cerr << "Invalid weight for usePriorShapeWeighting = " << usePriorShapeWeighting
             << ", iVariable = " << iVariable
             << ", value = " << value << endl;
        exit(1);
    }
    return weight;

}

TH1D * Efficiency1DCut[5];

void LoadEfficiency1DCut() {
    //TFile *file2 = TFile::Open("../Unfolding/Output/Efficiency2D_weightedRecoTrue.root");
   TFile *file2 = TFile::Open("./Data/zKodu2014_y06.root");
//TFile *file2 = TFile::Open("Files/zKodu2014.root");
    //check if file is open
    if (!file2) {
        cout << "File Efficiency1D.root not found" << endl;
        exit(1);
    }
    for (int i = 0; i < 5; i++) {
        Efficiency1DCut[i] = (TH1D*)file2->Get(Form("heffBinD0_%d", i));
        //Efficiency1DCut[i] = (TH1D*)file2->Get(Form("Efficiency1D_original_%d", i));
        //check
        if (!Efficiency1DCut[i]) {
            cout << "Histogram hEffPt_y0.6_cent" << i << " not found" << endl;
            exit(1);
        }
        //set directory 0
        Efficiency1DCut[i]->SetDirectory(0);
    }
    file2->Close();
}

TFile *fPaperD0Sys = NULL;

// indices: [centrality group][error source]
// centrality group: 0 = 0-10, 1 = 10-40, 2 = 40-80
// error source according to paper plot:
// 0 TPC track
// 1 PID
// 2 raw yield extract
// 3 single track pT
// 4 topo. eff.
// 5 double count
// 6 vertex correction
// 7 secondary track
// 8 TPC luminosity
// 9 general sys.
TGraph *gPaperD0Sys[9][10];

bool LoadPaperD0Systematics(const char *fileName = "Data/Paper/SysErrorsPaper.root")
{
    if (fPaperD0Sys) return true;

    fPaperD0Sys = new TFile(fileName, "READ");

    if (!fPaperD0Sys || fPaperD0Sys->IsZombie()) {
        cout << "[LoadPaperD0Systematics] ERROR: Cannot open " << fileName << endl;
        return false;
    }

    for (int ic = 0; ic < 9; ic++) {
        for (int ie = 0; ie < 10; ie++) {
            gPaperD0Sys[ic][ie] = (TGraph*)fPaperD0Sys->Get(
                Form("gSys_cent%i_err%i", ic, ie)
            );

            if (!gPaperD0Sys[ic][ie]) {
                cout << "[LoadPaperD0Systematics] WARNING: Missing "
                     << Form("gSys_cent%i_err%i", ic, ie)
                     << endl;
            }
        }
    }

    cout << "[LoadPaperD0Systematics] Loaded paper D0 systematics from "
         << fileName << endl;

    return true;
}

int GetPaperD0SysCentIndex(double cent)
{
    if (cent >= 0  && cent < 10) return 0; // 0-10
    if (cent >= 10 && cent < 20) return 1; // 10-20
    if (cent >= 20 && cent < 40) return 2; // 20-40
    if (cent >= 40 && cent < 60) return 3; // 40-60
    if (cent >= 60 && cent < 80) return 4; // 60-80

    return -1;
}

double GetPaperD0SysWeight(int systematicSPlot, double centAlt, double d0pt)
{
    // nominal or non-paper variations
    if (systematicSPlot < 7 || systematicSPlot > 20) return 1.0;

    int centIndex = GetPaperD0SysCentIndex(centAlt);
    if (centIndex < 0) return 1.0;

    int errIndex = -1;
    int sign = 0;

    if (systematicSPlot == 7)  { errIndex = 0; sign = +1; } // TPC track up
    if (systematicSPlot == 8)  { errIndex = 0; sign = -1; } // TPC track down

    if (systematicSPlot == 9)  { errIndex = 1; sign = +1; } // PID up
    if (systematicSPlot == 10) { errIndex = 1; sign = -1; } // PID down

    if (systematicSPlot == 11) { errIndex = 3; sign = +1; } // single track pT up
    if (systematicSPlot == 12) { errIndex = 3; sign = -1; } // single track pT down

    if (systematicSPlot == 13) { errIndex = 4; sign = +1; } // topo eff up
    if (systematicSPlot == 14) { errIndex = 4; sign = -1; } // topo eff down

    if (systematicSPlot == 15) { errIndex = 5; sign = +1; } // double count up
    if (systematicSPlot == 16) { errIndex = 5; sign = -1; } // double count down

    if (systematicSPlot == 17) { errIndex = 6; sign = +1; } // vertex corr up
    if (systematicSPlot == 18) { errIndex = 6; sign = -1; } // vertex corr down

    if (systematicSPlot == 19) { errIndex = 7; sign = +1; } // secondary track up
    if (systematicSPlot == 20) { errIndex = 7; sign = -1; } // secondary track down

    if (errIndex < 0) return 1.0;

    TGraph *g = gPaperD0Sys[centIndex][errIndex];
    if (!g) return 1.0;

    double relSys = fabs(g->Eval(d0pt));

    double weight = 1.0 + sign * relSys;

    // safety protection
    if (weight < 0.0) weight = 0.0;

    return weight;
}

double D0_efficiency_1DCut(int centrality, double D0_pT){

    int centr_range_ =   (centrality >= 0 &&  centrality < 10 ) ? 0 :
                         (centrality >= 10 &&  centrality < 20) ? 1 :
                         (centrality >= 20 &&  centrality < 40) ? 2 :
                         (centrality >= 40 &&  centrality < 60) ? 3 :
                         4;

    double efficiency = Efficiency1DCut[centr_range_]->GetBinContent(Efficiency1DCut[centr_range_]->FindBin(D0_pT));
    return efficiency;
}


void LoadFiles() {
    TFile *doub = new TFile("../Neil_Unfold/Data/EssentialFiles/MisPID_SB_Final.root", "READ");


    //check
    if (!doub || doub->IsZombie()) {
        cout << "Error opening file Data/MisPID_SB_Final.root" << endl;
        exit(1);
    }
    gDoub0 = (TGraph *) doub->Get("DoubleCounting_Cen_0_SB");
    gDoub1 = (TGraph *) doub->Get("DoubleCounting_Cen_1_SB");
    gDoub2 = (TGraph *) doub->Get("DoubleCounting_Cen_2_SB");
    gDoub3 = (TGraph *) doub->Get("DoubleCounting_Cen_3_SB");

    TFile *eff = new TFile("../Neil_Unfold/Data/EssentialFiles/D0Eff.root", "READ");
    if (!eff || eff->IsZombie()) {
        cout << "Error opening file Data/D0Eff.root" << endl;
        exit(1);
    }

    hEff0 = (TH1F *)eff->Get("hEff0");
    hEff1 = (TH1F *)eff->Get("hEff1");
    hEff2 = (TH1F *)eff->Get("hEff2");
    hEff3 = (TH1F *)eff->Get("hEff3");
    hEff4 = (TH1F *)eff->Get("hEff4");

    hEff0_Sys = (TH1F *)hEff0->Clone("hEff0_Sys");
    hEff1_Sys = (TH1F *)hEff1->Clone("hEff1_Sys");
    hEff2_Sys = (TH1F *)hEff2->Clone("hEff2_Sys");
    hEff3_Sys = (TH1F *)hEff3->Clone("hEff3_Sys");
    hEff4_Sys = (TH1F *)hEff4->Clone("hEff4_Sys");
}

//MAHCINE

bool FONLLjet = true;
bool CutOfNegative = true;
double minJetPtRecoCut = -999;
int savedIter = -999;

//////////////////////////////////////////
// Centralities = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80}
// Centralities = { 8, 7,  6,  5,  4,  3,  2,  1,  0}
double centBins[nCentralityBins + 1] = {0, 10, 40, 80}; // in icreasing order
double getEff(int mCen, float mPt)
{
    double eff = 0;

    // Note I am using the efficiency histograms with "_Sys" here. If running applyWeights() only, this doesn't matter
    // If running doSys(), this will shuffle the efficinecy with a Gaussian each iteration in the _Sys histograms.
    // To not have to change the code here, _Sys is used all the time but is equivelent to original histograms in nominal running
    if (mCen >= 0 && mCen < 10)
        eff = hEff0_Sys->GetBinContent(hEff0_Sys->FindBin(mPt));
    else if (mCen >= 10 && mCen < 20)
        eff = hEff1_Sys->GetBinContent(hEff1_Sys->FindBin(mPt));
    else if (mCen >= 20 && mCen < 40)
        eff = hEff2_Sys->GetBinContent(hEff2_Sys->FindBin(mPt));
    else if (mCen >= 40 && mCen < 60)
        eff = hEff3_Sys->GetBinContent(hEff3_Sys->FindBin(mPt));
    else if (mCen >= 60 && mCen < 80)
        eff = (2. / 3.) * hEff4_Sys->GetBinContent(hEff4_Sys->FindBin(mPt)); // efficiency is from paper which has scale factor
    return eff;
}
void SaveFine(){

    //uložim soubor s jemným binováním do Output/FineBinning.root
    TFile *file = TFile::Open("Output/FineBinning.root", "RECREATE");
    if (!file || file->IsZombie()) {
        cout << "Error opening file Output/FineBinning.root" << endl;
        return;
    }
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 7; j++) {
            hResponseFine4D[i][j]->Write();
            hRealFine[i][j]->Write();
        }
    }
    file->Close();

}
double ComputeTVShapeDistance(const TH1D& hUnfolded,
                              const TH1D& hPrior,
                              double xMin = -1e99,
                              double xMax =  1e99)
{
    double sumU = 0.0;
    double sumP = 0.0;

    // integraly v danem rozsahu
    for (int iBin = 1; iBin <= hUnfolded.GetNbinsX(); iBin++) {
        double x = hUnfolded.GetBinCenter(iBin);
        if (x <= xMin || x >= xMax) continue;

        double u = hUnfolded.GetBinContent(iBin);
        double p = hPrior.GetBinContent(iBin);
        double bw = hUnfolded.GetBinWidth(iBin);

        if (!std::isfinite(u) || !std::isfinite(p)) continue;
        if (u < 0.0) u = 0.0;
        if (p < 0.0) p = 0.0;

        sumU += u * bw;
        sumP += p * bw;
    }

    if (sumU <= 0.0 || sumP <= 0.0) return -1.0;

    double tv = 0.0;

    for (int iBin = 1; iBin <= hUnfolded.GetNbinsX(); iBin++) {
        double x = hUnfolded.GetBinCenter(iBin);
        if (x <= xMin || x >= xMax) continue;

        double u = hUnfolded.GetBinContent(iBin);
        double p = hPrior.GetBinContent(iBin);
        double bw = hUnfolded.GetBinWidth(iBin);

        if (!std::isfinite(u) || !std::isfinite(p)) continue;
        if (u < 0.0) u = 0.0;
        if (p < 0.0) p = 0.0;

        double uProb = (u * bw) / sumU;
        double pProb = (p * bw) / sumP;

        tv += std::abs(uProb - pProb);
    }

    return 0.5 * tv;
}
double getDoubleCount(int mCen, float mPt)
{
    double dcount = 1;
    if (mCen == 0 && mCen == 80)
        dcount = gDoub0->Eval(mPt);
    else if (mCen >= 0 && mCen < 10)
        dcount = gDoub1->Eval(mPt);
    else if (mCen >= 10 && mCen < 40)
        dcount = gDoub2->Eval(mPt);
    else if (mCen >= 40 && mCen < 80)
        dcount = gDoub3->Eval(mPt);
    return dcount;
}

Int_t getCentralityBin(const Float_t &centrality)
{
    for (Int_t i = 0; i < nCentralityBins; i++)
    {
        if (centrality >= centBins[i] && centrality < centBins[i + 1])
            return i;
    }
    //cout << "Error: centrality not in range" << endl;
    return -1;
}

Int_t getCentralityBin9(const Float_t &centrality)
{
    return (centrality == 75 || centrality == 70) ? 0 :
           (centrality == 65 || centrality == 60) ? 1 :
           (centrality == 55 || centrality == 50) ? 2 :
           (centrality == 45 || centrality == 40) ? 3 :
           (centrality == 35 || centrality == 30) ? 4 :
           (centrality == 25 || centrality == 20) ? 5 :
           (centrality == 15 || centrality == 10) ? 6 :
           (centrality == 5)                     ? 7 :
           (centrality == 0)                     ? 8 : -1;

}
long factorial(const int n)
{
    long f = 1;
    for (int i=1; i<=n; ++i)
        f *= i;
    return f;
}
int GetCombIndex(int i, int j, int n) {
    return i * (n - 1) - (i * (i - 1)) / 2 + (j - i - 1);
}
Int_t getCentralityBin99(const Float_t &centrality)
{
    return (centrality >= 7) ? 0 :
           (centrality >= 4) ? 1 :
           2;
}

Int_t getVariable(TString var) {
    int variable = 0;
    if (var == "z") variable = 1;
    if (var == "#lambda^{1}_{1}") variable = 2;
    if (var == "#lambda^{1}_{1.5}") variable = 3;
    if (var == "#lambda^{1}_{2}") variable = 4;
    if (var == "#lambda^{1}_{3}") variable = 5;
    if (var == "#lambda^{1}_{0.5}") variable = 6;
    if (var == "p_{T}^{D}") variable = 7;

    return  variable;

}

TH1D *NeilEff[5];

//načtu soubour Output/Efficiency2D.root
void LoadEfficiency2D() {
    TFile *file = TFile::Open("Output/Efficiency2D.root");

    //check if file is open
    if (!file) {
        cout << "File Efficiency2D.root not found" << endl;
        exit(1);
    }
    for (int i = 0; i < 5; i++) {
        Efficiency2D[i] = (TH2F*)file->Get(Form("Efficiency2D_%d", i));
        Efficiency2D_TEff[i] = (TEfficiency*)file->Get(Form("Efficiency2D_TEfficiency_%d", i));

        //check
        if (!Efficiency2D[i]) {
            cout << "Histogram Efficiency2D_" << i << " not found" << endl;
            exit(1);
        }
        //set directory 0
        Efficiency2D[i]->SetDirectory(0);
    }
   // file->Close();

    TFile *NeilEffF = TFile::Open("./Data/D0Eff.root");
    if (!NeilEffF || NeilEffF->IsZombie()) {
        cout << "File not found or is corrupted: NeilEff.root" << endl;
        return;
    }
    //Do NeilEff načtu histogramy s názvy "hEff0", "hEff1", "hEff2"
    for (int i = 0; i < 5; i++) {
        NeilEff[i] = (TH1D *) NeilEffF->Get(Form("hEff%i", i))->Clone(Form("hEff%i", i));
        if (!NeilEff[i]) {
            cout << "Histogram hEff" << i << " not found in NeilEff.root" << endl;
            return;
        }
        NeilEff[i]->SetDirectory(0);
    }

   NeilEffF->Close();

}


double BrRatioD0 = 0.0395;
void NormalizeFinalSpectraPt(TH1D *hist, const Int_t color) {

for (int i = 1; i <= hist->GetNbinsX(); i++) {
        double binContent = hist->GetBinContent(i);
        double binError = hist->GetBinError(i);
        double binCenter = hist->GetBinCenter(i);
        double binWidth = hist->GetBinWidth(i);

        double normalizedContent = binContent;
        normalizedContent /= (2 * TMath::Pi()); // 2 * Pi for the azimuthal angle
        normalizedContent /= 1.2; // deta_jet from -0.6 to 0.6
        normalizedContent /= binWidth; // dpT_jet
        normalizedContent /= binCenter; // pT_jet
        normalizedContent /= 2.; // 2 for the two D0 mesons (D0 and D0bar)
        normalizedContent /= BrRatioD0; // Branching ratio for D0 -> K pi

        double normalizedError = binError;
        normalizedError /= (2 * TMath::Pi()); // 2 * Pi for the azimuthal angle
        normalizedError /= 1.2; // deta_jet from -0.6 to 0.6
        normalizedError /= binWidth; // dpT_jet
        normalizedError /= binCenter; // pT_jet
        normalizedError /= 2.; // 2 for the two D0 mesons (D0 and D0bar)
        normalizedError /= BrRatioD0; // Branching ratio for D0 -> K pi

        // Set the normalized content and error back to the histogram
        hist->SetBinContent(i, normalizedContent);
        hist->SetBinError(i, normalizedError);

        hist->SetLineColor(color);
        hist->SetMarkerColor(color);
        hist->SetMarkerStyle(20);
    }


}
void NormalizeFinalSpectra(TH1D *hist, const Int_t color) {

for (int i = 1; i <= hist->GetNbinsX(); i++) {
        double binContent = hist->GetBinContent(i);
        double binError = hist->GetBinError(i);
        double binCenter = hist->GetBinCenter(i);
        double binWidth = hist->GetBinWidth(i);

        double normalizedContent = binContent;
        //normalizedContent /= (2 * TMath::Pi()); // 2 * Pi for the azimuthal angle
        normalizedContent /= 1.2; // deta_jet from -0.6 to 0.6
        normalizedContent /= binWidth; // dpT_jet
        //normalizedContent /= binCenter; // pT_jet
        normalizedContent /= 2.; // 2 for the two D0 mesons (D0 and D0bar)
        normalizedContent /= BrRatioD0; // Branching ratio for D0 -> K pi

        double normalizedError = binError;
        //normalizedError /= (2 * TMath::Pi()); // 2 * Pi for the azimuthal angle
        normalizedError /= 1.2; // deta_jet from -0.6 to 0.6
        normalizedError /= binWidth; // dpT_jet
        //normalizedError /= binCenter; // pT_jet
        normalizedError /= 2.; // 2 for the two D0 mesons (D0 and D0bar)
        normalizedError /= BrRatioD0; // Branching ratio for D0 -> K pi

        // Set the normalized content and error back to the histogram
        hist->SetBinContent(i, normalizedContent);
        hist->SetBinError(i, normalizedError);

        hist->SetLineColor(color);
        hist->SetMarkerColor(color);
        hist->SetMarkerStyle(20);
    }


}
void LoadEfficiency1DPaper() {

    TFile *file = TFile::Open("./Data/new/effFromPaper.root");

    if (!file || file->IsZombie()) {
        cout << "Error opening file ./Data/new/effFromPaper.root" << endl;
        return;
    }

    for (int i = 0; i < 5; i++) {
        Efficiency1DPaper[i] = (TH1F*)file->Get(Form("heffD0_%d", i));
        //check
        if (!Efficiency1DPaper[i]) {
            cout << "Histogram Efficiency1D_" << i << " not found" << endl;
            exit(1);
        }
        //set directory 0
        Efficiency1DPaper[i]->SetDirectory(0);
    }





}
//načtu soubour Output/Efficiency2D.root
void LoadEfficiency1D() {
    TFile *file = TFile::Open("Output/Efficiency2D_weightedRecoTrue.root");

    //check if file is open
    if (!file) {
        cout << "File Efficiency1D.root not found" << endl;
        exit(1);
    }
    for (int i = 0; i < 5; i++) {
       // Efficiency1D[i] = (TH1D*)file->Get(Form("hEffPt_y0.6_cent%d", i));
        //D0_efficiency_1D
        Efficiency1D[i] = (TH1D*)file->Get(Form("Efficiency1D_original_%d", i));
        //check
        if (!Efficiency1D[i]) {
            cout << "Histogram hEffPt_y0.6_cent" << i << " not found" << endl;
            exit(1);
        }
        //set directory 0
        Efficiency1D[i]->SetDirectory(0);
    }
    file->Close();
}

double D0_efficiency_2D(double D0_pT, double D0_y, int centrality){
    int centr_range_ =   (centrality >= 7) ? 0 ://0-10%
                        (centrality == 6) ? 1 ://10-20%
                        (centrality >= 4) ? 2 ://20-40%
                        (centrality >= 2) ? 3 ://40-60%
                        4; //60-80%

    double efficiency = Efficiency2D[centr_range_]->GetBinContent(Efficiency2D[centr_range_]->FindBin( D0_y,D0_pT));
    TH2F* hist = (TH2F*) const_cast<TH1*>(Efficiency2D_TEff[centr_range_]->GetTotalHistogram());
    int bin = hist->FindBin(D0_y, D0_pT);
    double efficiency_TEff = Efficiency2D_TEff[centr_range_]->GetEfficiency(bin);
    return efficiency_TEff;
}

//D0Eff.root
double D0_efficiency_Neil(double D0_pT, int centrality){
    int centr_range_ =   (centrality >= 7) ? 0 :
                         (centrality == 6) ? 1 :
                         (centrality >= 4) ? 2 :
                         (centrality >= 2) ? 3 :
                         4;

    double efficiency = NeilEff[centr_range_]->GetBinContent(NeilEff[centr_range_]->FindBin(D0_pT));
    //!!!!//
    if (centr_range_ == 4)   efficiency = efficiency * (2./3.);
    return efficiency;

}
double D0_efficiencyError_Neil(double D0_pT, int centrality) {
    int centr_range_ = (centrality >= 7) ? 0 :
                       (centrality == 6) ? 1 :
                       (centrality >= 4) ? 2 :
                       (centrality >= 2) ? 3 :
                       4;

    double efficiency = NeilEff[centr_range_]->GetBinError(NeilEff[centr_range_]->FindBin(D0_pT));
    return efficiency;

}
double D0_efficiencyError_2D(double D0_pT, double D0_y, int centrality){
    int centr_range_ =   (centrality >= 7) ? 0 :
                         (centrality == 6) ? 1 :
                         (centrality >= 4) ? 2 :
                         (centrality >= 2) ? 3 :
                         4;

    double efficiency = Efficiency2D[centr_range_]->GetBinError(Efficiency2D[centr_range_]->FindBin( D0_y,D0_pT));
    //TH2F* hist = (TH2F*) const_cast<TH1*>(Efficiency2D_TEff[centr_range_]->GetTotalHistogram());
    //int bin = hist->FindBin(D0_y, D0_pT);
    //double efficiency_TEff = Efficiency2D_TEff[centr_range_]->GetEfficiency(bin);
    return efficiency;
}

TGraphAsymmErrors* HistToGraphAsymmErrors(TH1D* hist) {
    int n = hist->GetNbinsX();
    TGraphAsymmErrors* graph = new TGraphAsymmErrors(n);

    for (int i = 1; i <= n; ++i) {
        double x = hist->GetBinCenter(i);
        double y = hist->GetBinContent(i);
        double ex = hist->GetBinWidth(i) / 2.0;
        double ey = hist->GetBinError(i);

        graph->SetPoint(i - 1, x, y);
        graph->SetPointError(i - 1, ex, ex, ey, ey);  // symetrické chyby
    }

    return graph;
}
double D0_efficiency_1D(double D0_pT, int centrality){
    int centr_range_ =   (centrality >= 7) ? 0 :
                         (centrality == 6) ? 1 :
                         (centrality >= 4) ? 2 :
                         (centrality >= 2) ? 3 :
                         4;

    double efficiency = Efficiency1D[centr_range_]->GetBinContent(Efficiency1D[centr_range_]->FindBin(D0_pT));
    return efficiency;
}


double D0_efficiency_1DPaper(double D0_pT, int centrality){

    //otevřu soubor v ,/Data/new/effFromPaper.root
   // TFile *file = TFile::Open("./Data/new/effFromPaper.root");


    int centr_range_ =   (centrality >= 7) ? 0 :
                         (centrality == 6) ? 1 :
                         (centrality >= 4) ? 2 :
                         (centrality >= 2) ? 3 :
                         4;

    double efficiency = Efficiency1DPaper[centr_range_]->GetBinContent(Efficiency1DPaper[centr_range_]->FindBin(D0_pT));
    if (centr_range_ == 4)   efficiency = efficiency * (3./2.); //efficiency is from paper which has scale factor

    // file->Close();
    return efficiency;

}

std::vector<int> CentrRangeTransf(const std::vector<int>& centrality) {
    std::vector<int> range(2);
    range[1] = (centrality[0] == 70) ? 0 :
               (centrality[0] == 60) ? 1 :
               (centrality[0] == 50) ? 2 :
               (centrality[0] == 40) ? 3 :
               (centrality[0] == 30) ? 4 :
               (centrality[0] == 20) ? 5 :
               (centrality[0] == 10) ? 6 :
               (centrality[0] == 5)  ? 7 :
               (centrality[0] == 0)  ? 8 :
               ([]() -> Int_t {
                   std::cout << "Nepodařilo se nastavit dolní hranici centrality" << std::endl;
                   exit(1);
               })();

    range[0] = (centrality[1] == 80) ? 0 :
               (centrality[1] == 70) ? 1 :
               (centrality[1] == 60) ? 2 :
               (centrality[1] == 50) ? 3 :
               (centrality[1] == 40) ? 4 :
               (centrality[1] == 30) ? 5 :
               (centrality[1] == 20) ? 6 :
               (centrality[1] == 10) ? 7 :
               (centrality[1] == 5)  ? 8 :
               ([]() -> Int_t {
                   std::cout << "Nepodařilo se nastavit horní hranici centrality" << std::endl;
                   exit(1);
               })();


    return range;
}
void DrawLineOne()
{
    TLine *lineOne = new TLine(ptMcBinsVecCustom[0][0], 1, ptMcBinsVecCustom[0][ptMcBinsVecCustom[0].size()-1], 1);
    lineOne->SetLineColor(kBlack);
    lineOne->SetLineStyle(2);
    lineOne->SetLineWidth(1);
    lineOne->Draw();
}
void DrawLineZero()
{
    TLine *lineOne = new TLine(ptMcBinsVecCustom[0][0], 0, ptMcBinsVecCustom[0][ptMcBinsVecCustom[0].size()-1], 0);
    lineOne->SetLineColor(kBlack);
    lineOne->SetLineStyle(2);
    lineOne->SetLineWidth(1);
    lineOne->Draw();
}
void DrawLineVertical(double x)
{
    TLine *lineOne = new TLine(x, -1, x, 1);
    lineOne->SetLineColor(kBlack);
    lineOne->SetLineStyle(2);
    lineOne->SetLineWidth(1);
    lineOne->Draw();
}
void DrawLineOne2(double min, double max)
{
    TLine *lineOne = new TLine(min, 1, max, 1);
    lineOne->SetLineColor(kBlack);
    lineOne->SetLineStyle(2);
    lineOne->SetLineWidth(1);
    lineOne->Draw();
}
void DrawLineZero2(double min, double max)
{
    TLine *lineOne = new TLine(min, 0, max, 0);
    lineOne->SetLineColor(kBlack);
    lineOne->SetLineStyle(2);
    lineOne->SetLineWidth(1);
    lineOne->Draw();
}
void DrawKinEff2D(TEfficiency pTVar[], TEfficiency pT[], TEfficiency var[], TEfficiency varCut[], TEfficiency varPtCut[], double min[], double max[],double minR[], double maxR[], TCanvas *can){
    can->Clear();
    can->SetCanvasSize(1200, 1200);
    can->Divide(3, 3);

    const TH1* htot = var[0].GetTotalHistogram(); // nebo GetCopyTotalHisto() a pak delete
    const char* xTitle = htot ? htot->GetXaxis()->GetTitle() : "";
    std::string xTitleReco(xTitle);
    size_t pos = xTitleReco.find("true");
    if (pos != std::string::npos)
        xTitleReco.replace(pos, 4, "reco");


    for (int iCent = 0; iCent < nCentralityBins; iCent++) {

        const double xmin = ptMcBinsVecCustom[iCent].front();
        const double xmax = ptMcBinsVecCustom[iCent].back();
        const double xminCut = ptRecoBinsVec[iCent].front();
        const double xmaxCut = ptRecoBinsVec[iCent].back();


        can->cd(iCent + 1);
            gPad->SetLeftMargin(0.10);
            gPad->SetRightMargin(0.22);

            pTVar[iCent].Draw("COLZ");
        gPad->Update();  // jinak palette ještě neexistuje

// 1) Osa Z přes painted histogram (TH2*)
        TH2 *h2 = dynamic_cast<TH2*>(pTVar[iCent].GetPaintedHistogram());
        if (h2) {
            // nový titulek Z
            h2->GetZaxis()->SetTitle(Form("#frac{N(|#eta_{Jet}^{reco}|<0.6 && %.1f < p_{T,Jet}^{reco} < %.1f GeV/c && %.2f < %s < %.2f)}{N(|#eta_{Jet}^{reco}|<0.6)}",
                                          xminCut, xmaxCut, minR[iCent],xTitleReco.c_str(), maxR[iCent]));
            h2->GetZaxis()->SetTitleOffset(2.8);
            h2->GetZaxis()->SetTitleSize(0.03);

            // sundej starou paletu
            if (auto pal = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject("palette")) {
                h2->GetListOfFunctions()->Remove(pal);
                delete pal; // bezpečné – paleta je jen kreslící objekt
            }

            // znovu ji vytvoř (vezme nový Z-title)
            h2->Draw("COLZ");
        }

        gPad->Modified();
        gPad->Update();
        can->cd(iCent + 4);

        gPad->SetLeftMargin(0.15);
            pT[iCent].Draw("E1");
            gPad->Update();              // <-- důležité, jinak GetPaintedGraph() vrací nullptr
            //set range
            pT[iCent].GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.2);
            pT[iCent].GetPaintedGraph()->GetXaxis()->SetRangeUser(ptMcBinsVecCustom[iCent][0], ptMcBinsVecCustom[iCent][ptMcBinsVecCustom[iCent].size() - 1]);
            //offset
            pT[iCent].GetPaintedGraph()->GetYaxis()->SetTitleOffset(1.3);

            varPtCut[iCent].Draw("E1same");
            gPad->Update();              // <-- důležité, jinak GetPaintedGraph
            //color
            varPtCut[iCent].SetLineColor(kRed);
            varPtCut[iCent].SetMarkerColor(kRed);
            DrawLineOne2(ptMcBinsVecCustom[iCent][0], ptMcBinsVecCustom[iCent][ptMcBinsVecCustom[iCent].size() - 1]);

            //legend
            TLegend *legPt = new TLegend(0.2, 0.2, 0.3, 0.3);
            legPt->AddEntry(pT[iCent].GetPaintedGraph(), Form("(-#infty < %s < +#infty)", xTitle), "lp");
            legPt->AddEntry(varPtCut[iCent].GetPaintedGraph(), Form("(%.2f < %s < %.2f)", min[iCent], xTitle,max[iCent]), "lp");
           /* legPt->AddEntry((TH1D*)0, Form("#frac{N(|#eta_{Jet}^{reco}|<0.6 && %.1f < p_{T,Jet}^{reco} < %.1f GeV/c && %.2f < %s < %.2f)}{N(|#eta_{Jet}^{reco}|<0.6)}",
                                           xminCut, xmaxCut, minR[iCent],xTitleReco.c_str(), maxR[iCent]),"");*/
            legPt->SetFillColorAlpha(0, 0);
            //text size
            legPt->SetTextSize(0.04);
            legPt->SetBorderSize(0);
            legPt->Draw();


        TLatex tex;
        tex.SetNDC();
        tex.SetTextFont(42);
        tex.SetTextSize(0.028);
        tex.DrawLatex(0.13, 0.94, Form("#frac{N(|#eta_{Jet}^{reco}|<0.6 && %.1f < p_{T,Jet}^{reco} < %.1f GeV/c && %.2f < %s < %.2f)}{N(|#eta_{Jet}^{reco}|<0.6)}",
                                       xminCut, xmaxCut, minR[iCent],xTitleReco.c_str(), maxR[iCent]));

        can->cd(iCent + 7);
        gPad->SetLeftMargin(0.15);
        var[iCent].Draw("E1");
        gPad->Update();              // <-- důležité, jinak GetPaintedGraph
        //set range
        var[iCent].GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.2);
        var[iCent].GetPaintedGraph()->GetXaxis()->SetRangeUser(min[iCent], max[iCent]);
        //offset
        var[iCent].GetPaintedGraph()->GetYaxis()->SetTitleOffset(1.3);
        DrawLineOne2(min[iCent], max[iCent]);

        varCut[iCent].Draw("E1same");
        gPad->Update();              // <-- důležité, jinak GetPaintedGraph() vrací nullptr
        //set range
        varCut[iCent].GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.2);
        varCut[iCent].GetPaintedGraph()->GetXaxis()->SetRangeUser(min[iCent], max[iCent]);
        //offset
        varCut[iCent].GetPaintedGraph()->GetYaxis()->SetTitleOffset(1.3);
        DrawLineOne2(min[iCent], max[iCent]);
        //color
        varCut[iCent].SetLineColor(kRed);
        varCut[iCent].SetMarkerColor(kRed);
        varCut[iCent].SetMarkerStyle(20);
        varCut[iCent].SetMarkerSize(1.5);
        varCut[iCent].SetLineWidth(2);

        //tlegend
        TLegend *leg = new TLegend(0.2, 0.2, 0.3, 0.3);
        leg->AddEntry(var[iCent].GetPaintedGraph(), "(-#infty < p_{T}^{true} < +#infty GeV/c)", "lp");
        leg->AddEntry(varCut[iCent].GetPaintedGraph(), "(5 < p_{T}^{true} < 20 GeV/c)", "lp");
        leg->SetFillColorAlpha(0, 0);
        //text size
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->Draw();

        tex.DrawLatex(0.13, 0.94, Form("#frac{N(|#eta_{Jet}^{reco}|<0.6 && %.1f < p_{T,Jet}^{reco} < %.1f GeV/c && %.2f < %s < %.2f)}{N(|#eta_{Jet}^{reco}|<0.6)}",
                                       xminCut, xmaxCut, minR[iCent],xTitleReco.c_str(), maxR[iCent]));


    }

    can->SaveAs(outPdf);
    can->Clear();
}

const Double_t Ncoll[nCentralityBins] = {952., 397., 58.};
//double taa[3] = {959.42547, 401.44672, 59.01580}; Neils


//centralities, variables, iterations
double PrintCheckNumbers[3][8][11];

TH2D *hRespZ[8][3];
TH2D *hRespZHighRes[8][3];
void Stejn2(TH2D &h2,const char* name = "Transf") {

    vector <Double_t> x_edges;
    vector <Double_t> y_edges;


    for (int i = 1; i <= h2.GetNbinsX() + 1; i++) {
        x_edges.push_back(h2.GetXaxis()->GetBinLowEdge(i));
        // cout << x_edges[i-1] << endl;
    }
    //cout << "--" << endl;

    for (int i = 1; i <= h2.GetNbinsY() + 1; i++) {
        y_edges.push_back(h2.GetYaxis()->GetBinLowEdge(i));
        // cout << y_edges[i-1] << endl;
    }



    TH2D * Transf = new TH2D(name, "Transformed 2D Histogram", x_edges.size() - 1, 0,10*(x_edges.size()-1), y_edges.size() - 1,0,10*(y_edges.size()-1));

    //přepíšu názvy os
    Transf->GetXaxis()->SetTitle(h2.GetXaxis()->GetTitle());
    Transf->GetYaxis()->SetTitle(h2.GetYaxis()->GetTitle());

    // Vyplnění nového histogramu s původními hodnotami
    for (int i = 1; i <= h2.GetNbinsX(); i++) {
        for (int j = 1; j <= h2.GetNbinsY(); j++) {
            // Přenos obsahu původního histogramu
            Transf->SetBinContent(i, j, h2.GetBinContent(i, j));
        }
    }

    Transf->GetXaxis()->SetLabelSize(0);
    Transf->GetYaxis()->SetLabelSize(0);
    Transf->GetXaxis()->SetTickSize(0);
    Transf->GetYaxis()->SetTickSize(0);

    Transf->Draw("colztext");
    //keep ratio
    TGaxis *custom_axis = new TGaxis(0, 0, 10*(x_edges.size()-1), 0, 0, 10*(x_edges.size()-1), 510+(h2.GetNbinsX() > 10?10:0), "");
    TGaxis *custom_axis2 = new TGaxis(0, 0,0, 10*(y_edges.size()-1), 0, 10*(y_edges.size()-1), 510+(h2.GetNbinsY() > 10?10:0), "");

    for (int i = 1; i <= h2.GetNbinsX()+1; i++) {
        custom_axis->ChangeLabel(i, -1, -1, -1, -1, -1, Form("%.3g", x_edges[i-1]));
        //cout << x_edges[i-1] << endl;
    }




    //custom_axis->ChangeLabel(1, -1, -1, -1, -1, -1, "-#infty");
    //custom_axis->ChangeLabel(h2.GetNbinsX()+1, -1, -1, -1, -1, -1, "+#infty");
    for (int i = 1; i <= h2.GetNbinsY()+1; i++) {
        custom_axis2->ChangeLabel(i, -1, -1, -1, -1, -1, Form("%.3g", y_edges[i-1]));
        // cout << y_edges[i-1] << endl;
    }
    //custom_axis2->ChangeLabel(1, -1, -1, -1, -1, -1, "-#infty");
    //custom_axis2->ChangeLabel(h2.GetNbinsY()+1, -1, -1, -1, -1, -1, "+#infty");
    custom_axis->Draw();
    custom_axis2->Draw();
    //set label titles from h2
    Transf->GetXaxis()->SetTitle(h2.GetXaxis()->GetTitle());
    Transf->GetYaxis()->SetTitle(h2.GetYaxis()->GetTitle());
    //y label offset
    Transf->GetYaxis()->SetTitleOffset(1.3);
    //center y label
    Transf->GetYaxis()->CenterTitle();
    //Transf->GetZaxis()->SetRangeUser(0.0001, 10);
    Transf->GetZaxis()->SetRangeUser(-100, 100);

    //set logz
    //gPad->SetLogz();

    //nahradim h2 Transf
    h2 = *Transf;
}
void Stejn(TH2D *h2,const char* name = "Transf") {

    vector <Double_t> x_edges;
    vector <Double_t> y_edges;


    for (int i = 1; i <= h2->GetNbinsX() + 1; i++) {
        x_edges.push_back(h2->GetXaxis()->GetBinLowEdge(i));
       // cout << x_edges[i-1] << endl;
    }
    //cout << "--" << endl;

    for (int i = 1; i <= h2->GetNbinsY() + 1; i++) {
        y_edges.push_back(h2->GetYaxis()->GetBinLowEdge(i));
        // cout << y_edges[i-1] << endl;
    }



    TH2D * Transf = new TH2D(name, "Transformed 2D Histogram", x_edges.size() - 1, 0,10*(x_edges.size()-1), y_edges.size() - 1,0,10*(y_edges.size()-1));

    //přepíšu názvy os
    Transf->GetXaxis()->SetTitle(h2->GetXaxis()->GetTitle());
    Transf->GetYaxis()->SetTitle(h2->GetYaxis()->GetTitle());

    // Vyplnění nového histogramu s původními hodnotami
    for (int i = 1; i <= h2->GetNbinsX(); i++) {
        for (int j = 1; j <= h2->GetNbinsY(); j++) {
            // Přenos obsahu původního histogramu
            Transf->SetBinContent(i, j, h2->GetBinContent(i, j));
        }
    }

    Transf->GetXaxis()->SetLabelSize(0);
    Transf->GetYaxis()->SetLabelSize(0);
    Transf->GetXaxis()->SetTickSize(0);
    Transf->GetYaxis()->SetTickSize(0);

    Transf->Draw("colz");
    //keep ratio
    TGaxis *custom_axis = new TGaxis(0, 0, 10*(x_edges.size()-1), 0, 0, 10*(x_edges.size()-1), 510+(h2->GetNbinsX() > 10?10:0), "");
    TGaxis *custom_axis2 = new TGaxis(0, 0,0, 10*(y_edges.size()-1), 0, 10*(y_edges.size()-1), 510+(h2->GetNbinsY() > 10?10:0), "");

    for (int i = 1; i <= h2->GetNbinsX()+1; i++) {
        custom_axis->ChangeLabel(i, -1, -1, -1, -1, -1, Form("%.3g", x_edges[i-1]));
        //cout << x_edges[i-1] << endl;
    }




    //custom_axis->ChangeLabel(1, -1, -1, -1, -1, -1, "-#infty");
    //custom_axis->ChangeLabel(h2->GetNbinsX()+1, -1, -1, -1, -1, -1, "+#infty");
    for (int i = 1; i <= h2->GetNbinsY()+1; i++) {
        custom_axis2->ChangeLabel(i, -1, -1, -1, -1, -1, Form("%.3g", y_edges[i-1]));
        // cout << y_edges[i-1] << endl;
    }
    //custom_axis2->ChangeLabel(1, -1, -1, -1, -1, -1, "-#infty");
    //custom_axis2->ChangeLabel(h2->GetNbinsY()+1, -1, -1, -1, -1, -1, "+#infty");
    custom_axis->Draw();
    custom_axis2->Draw();
    //set label titles from h2
    Transf->GetXaxis()->SetTitle(h2->GetXaxis()->GetTitle());
    Transf->GetYaxis()->SetTitle(h2->GetYaxis()->GetTitle());
    //y label offset
    Transf->GetYaxis()->SetTitleOffset(1.3);
    //center y label
    Transf->GetYaxis()->CenterTitle();
    //Transf->GetZaxis()->SetRangeUser(0.0001, 10);
    Transf->GetZaxis()->SetRangeUser(0.0001, 10000000);

    //set logz
    //gPad->SetLogz();
}

#include "RooUnfoldResponse.h"
#include "TH2D.h"

TH2D* ProjectRecoTrue2D(const RooUnfoldResponse* response,
                        const TString& axisGroupX, const TString& axisNameX,
                        const TString& axisGroupY, const TString& axisNameY)
{
    // Získání os
    const TAxis* recoPt = response->Hmeasured()->GetXaxis();
    const TAxis* recoZ  = response->Hmeasured()->GetYaxis();
    const TAxis* truePt = response->Htruth()->GetXaxis();
    const TAxis* trueZ  = response->Htruth()->GetYaxis();

    const int nRecoPt = recoPt->GetNbins();
    const int nRecoZ  = recoZ->GetNbins();
    const int nTruePt = truePt->GetNbins();
    const int nTrueZ  = trueZ->GetNbins();

    // Výběr os podle vstupu
    const TAxis *axisX = nullptr, *axisY = nullptr;
    int nBinsX = 0, nBinsY = 0;

    auto setAxis = [&](const TString& group, const TString& axis, const TAxis*& out, int& n) {
        if (group == "reco") {
            if (axis == "x") { out = recoPt; n = nRecoPt; }
            else             { out = recoZ;  n = nRecoZ;  }
        } else {
            if (axis == "x") { out = truePt; n = nTruePt; }
            else             { out = trueZ;  n = nTrueZ;  }
        }
    };

    setAxis(axisGroupX, axisNameX, axisX, nBinsX);
    setAxis(axisGroupY, axisNameY, axisY, nBinsY);

    TH2D* hOut = new TH2D(Form("h%s%s_vs_%s%s", axisGroupY.Data(), axisNameY.Data(),
                               axisGroupX.Data(), axisNameX.Data()),
                          Form("%s%s vs %s%s", axisGroupY.Data(), axisNameY.Data(),
                               axisGroupX.Data(), axisNameX.Data()),
                          nBinsX, axisX->GetXbins()->GetArray(),
                          nBinsY, axisY->GetXbins()->GetArray());

    // Zdroje binů
    TH2D* hResp = (TH2D*) response->Hresponse();

    for (int iReco = 1; iReco <= hResp->GetNbinsX(); ++iReco) {
        for (int iTrue = 1; iTrue <= hResp->GetNbinsY(); ++iTrue) {
            double content = hResp->GetBinContent(iReco, iTrue);

            // Rozbalení z flat indexu – pT je vnitřní index
            int iRecoPt = (iReco - 1) % nRecoPt + 1;
            int iRecoZ  = (iReco - 1) / nRecoPt + 1;

            int iTruePt = (iTrue - 1) % nTruePt + 1;
            int iTrueZ  = (iTrue - 1) / nTruePt + 1;

            //cout << "x: " << iReco << " y: " << iTrue << " iRecoZ: " << iRecoZ << " iTrueZ: " << iTrueZ << endl;
           // cout << "iTrueZ: " << iTrueZ << " iRecoZ: " << iRecoZ << endl;
            //cout << "(" << iRecoPt << "," << iRecoZ << ") vs (" << iTruePt << "," << iTrueZ << ") = " << content << endl;

            // Výpočet hodnot
            double valX = 0, valY = 0;

            if (axisGroupX == "reco") {
                valX = (axisNameX == "x") ? recoPt->GetBinCenter(iRecoPt)
                                          : recoZ ->GetBinCenter(iRecoZ);
            } else {
                valX = (axisNameX == "x") ? truePt->GetBinCenter(iTruePt)
                                          : trueZ ->GetBinCenter(iTrueZ);
            }

            if (axisGroupY == "reco") {
                valY = (axisNameY == "x") ? recoPt->GetBinCenter(iRecoPt)
                                          : recoZ ->GetBinCenter(iRecoZ);
            } else {
                valY = (axisNameY == "x") ? truePt->GetBinCenter(iTruePt)
                                          : trueZ ->GetBinCenter(iTrueZ);
            }

            //cout << "iTrueZ: " << iTrueZ << " iRecoZ: " << iRecoZ << endl;

            hOut->Fill(valX, valY, content);
        }
    }

    return hOut;
}












#endif // CONFIG_H
