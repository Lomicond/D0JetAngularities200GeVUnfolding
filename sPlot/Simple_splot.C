#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "RooGenericPdf.h"
#include <iomanip>
#include "../config.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "RooUnfoldResponse.h"

// use this order for safety on library loading
using namespace RooFit;
using namespace RooStats;

Bool_t reject = true;

Double_t funkce(Double_t *x, Double_t *par)
{
    if (reject && x[0] > (1.8) && x[0] < (1.91)) {
        TF1::RejectPoint();
        return 0;
    }
    return x[0]*x[0]*par[0]+x[0]*par[1]+par[2];
}

Double_t funkce2(Double_t *x, Double_t *par)
{
    if (reject && x[0] > (1.81) && x[0] < (1.91)) {
        TF1::RejectPoint();
        return 0;
    }
    return x[0]*par[0]+par[1];
}

TString CentrRange(int C[2]){
    //Centrality cut
//	 	8 = 0-5%; 7 = 5-10%; 6 = 10-20%; 5 = 20-30%; 4 = 30-40%; 3 = 40-50%; 2 = 50-60%; 1 = 60-70%; 0 = 70-80%
//!bin:     9          8          7           6          5             4         3             2          1
    if (C[0]>9 || C[0]<0 || C[1]>9 || C[1]<0) return ("Wrong centrality range");
    TString CentrRange1[9] = {"70","60","50","40","30","20","10","5","0"};
    TString CentrRange2[9] = {"80","70","60","50","40","30","20","10","5"};
    return (CentrRange1[C[0]]+"-"+CentrRange2[C[1]]+"%");
}

TString CentrRangeshort(int C[2]){
    //Centrality cut
//	 	8 = 0-5%; 7 = 5-10%; 6 = 10-20%; 5 = 20-30%; 4 = 30-40%; 3 = 40-50%; 2 = 50-60%; 1 = 60-70%; 0 = 70-80%
//!bin:     9          8          7           6          5             4         3             2          1
    if (C[0]>9 || C[0]<0 || C[1]>9 || C[1]<0) return ("Wrong centrality range");
    TString CentrRange1[9] = {"70","60","50","40","30","20","10","5","0"};
    TString CentrRange2[9] = {"80","70","60","50","40","30","20","10","5"};
    return (CentrRange1[C[0]]+"-"+CentrRange2[C[1]]);
}
/*
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
               })();  // Volání lambda funkce v případě, že žádná z podmínek neplatí

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
*/
//Vytvořím funkci, která otevře daný root soubor načte daný tntuple a vrátí ho (vstupy budou název souboru a název tntuple)
TNtuple *openFileAndNtuple(const char *fileName, const char *ntupleName)
{
    // Otevření souboru fileName
    TFile *file = TFile::Open(fileName);

    // Kontrola, zda se soubor podařilo otevřít
    if (!file || file->IsZombie())
    {
        std::cout << "Nepodařilo se otevřít soubor " << fileName << std::endl;
        return nullptr;
    }

    // Načtení stromu TNtuple s názvem ntupleName
    TNtuple *ntuple = (TNtuple*)file->Get(ntupleName);

    // Kontrola, zda se podařilo načíst TNtuple
    if (!ntuple)
    {
        std::cout << "Nepodařilo se načíst TNtuple s názvem " << ntupleName << std::endl;
        file->Close();
        return nullptr;
    }

    return ntuple;
}

TTree *openFileAndNtuple2(const char *fileName, const char *ntupleName)
{
    // Otevření souboru fileName
    TFile *file = TFile::Open(fileName);

    // Kontrola, zda se soubor podařilo otevřít
    if (!file || file->IsZombie())
    {
        std::cout << "Nepodařilo se otevřít soubor " << fileName << std::endl;
        return nullptr;
    }

    // Načtení stromu TNtuple s názvem ntupleName
    TTree *ntuple = (TTree*)file->Get(ntupleName);

    // Kontrola, zda se podařilo načíst TNtuple
    if (!ntuple)
    {
        std::cout << "Nepodařilo se načíst TNtuple s názvem " << ntupleName << std::endl;
        file->Close();
        return nullptr;
    }

    return ntuple;
}

struct StJetTreeStruct2
{
    Double_t z_value;
    Double_t pT_value;
    Double_t lambda_value[4];
    Double_t s_weight_value;
    Double_t centr_weight_value;
    Double_t eff_weight_value;
};

void assignTree2(TTree *jetTree, StJetTreeStruct2 &measured_, TString name)
{
    jetTree->SetBranchAddress("z", &measured_.z_value);
    jetTree->SetBranchAddress("pT", &measured_.pT_value);
    jetTree->SetBranchAddress("lambda_1_1", &measured_.lambda_value[0]);
    jetTree->SetBranchAddress("lambda_1_1half", &measured_.lambda_value[1]);
    jetTree->SetBranchAddress("lambda_1_2", &measured_.lambda_value[2]);
    jetTree->SetBranchAddress("lambda_1_3", &measured_.lambda_value[3]);
    jetTree->SetBranchAddress("n_signal"+name+"_sw", &measured_.s_weight_value);
    jetTree->SetBranchAddress("centr_weight", &measured_.centr_weight_value);
    jetTree->SetBranchAddress("rev_weight_ef", &measured_.eff_weight_value);
}

double D0_DoubleCounting(TH1D *histDC[], double D0_pT, int centrality){

    int centr_range =   (centrality >= 7) ? 0 :
                        (centrality == 6) ? 1 :
                        (centrality >= 4) ? 2 :
                        (centrality >= 2) ? 3 :
                        4;

    return histDC[centr_range]->GetBinContent(histDC[centr_range]->FindBin(D0_pT));
}

double D0_efficiency(double D0_pT, int centrality){
    //https://inspirehep.net/literature/2051708
    // https://journals.aps.org/prc/abstract/10.1103/PhysRevC.99.034908

//D0 2014 efficiency                                  0-10%       10-20%         20-40%      40-60%     60-80%
 /*   const double efficiencies[11][5]={     {0.000622, 0.000740, 0.00088, 0.00108, 0.00140},         //pT = 0-0.5
                                           {0.000717,0.000717, 0.00110, 0.00143, 0.00160},          //pT = 0.5-1.0
                                           {0.001025, 0.001300, 0.00195, 0.00304, 0.00289},         //pT = 1.0-1.5
                                           {0.002248, 0.003014, 0.00411, 0.00610, 0.00528},         //pT = 1.5-2.0
                                           {0.003794,0.006104, 0.00738, 0.01072, 0.00997},          //pT = 2.0-2.5
                                           {0.006107, 0.009671, 0.01006, 0.01606, 0.01447},         //pT = 2.5-3.0
                                           {0.011249, 0.016457, 0.02087, 0.02565, 0.02482},         //pT = 3.0-4.0
                                           {0.015830, 0.023905, 0.02984, 0.03553, 0.03663},         //pT = 4.0-5.0
                                           {0.026520, 0.034723,	0.04134, 0.04656,0.04032},          //pT = 5.0-6.0
                                           {0.029904, 0.040095,0.04888,	0.05334,0.05080},           //pT = 6.0-8.0
                                           {0.032680,0.045229,	0.05514,0.06017,0.05685} };         //pT = 8.0-10.0
*/
 //D0 2014 efficiency    new
/*
    const double efficiencies[20][5]={
            {0.000285808, 0.000295828, 0.000279927, 0.000596018, 0.000790277},
            {0.000421541, 0.00047328, 0.000699687, 0.000715744, 0.000884337},
            {0.00092226, 0.00109118, 0.00168991, 0.00183595, 0.00180872},
            {0.0021721, 0.00288962, 0.00386534, 0.00415397, 0.0037352},
            {0.00438539, 0.00584226, 0.00684054, 0.00824008, 0.00896017},
            {0.00763724, 0.0093116, 0.011732, 0.0143702, 0.0154568},
            {0.0126991, 0.0164649, 0.0189316, 0.022335, 0.0226727},
            {0.0162577, 0.022248, 0.0255858, 0.0295719, 0.0300738},
            {0.0198151, 0.0266207, 0.0299858, 0.0346729, 0.0363673},
            {0.0232576, 0.0294244, 0.034567, 0.0390242, 0.0422618},
            {0.0302067, 0.0399938, 0.042766, 0.0465076, 0.0471732},
            {0.0320288, 0.0417054, 0.0464353, 0.0503547, 0.0484834},
            {0.0346935, 0.0450093, 0.0493515, 0.0534783, 0.0536006},
            {0.0369348, 0.0471141, 0.0520191, 0.0584171, 0.0565219},
            {0.0375094, 0.0486814, 0.0546515, 0.0588733, 0.0568737},
            {0.0379939, 0.0486559, 0.0564554, 0.0613643, 0.0580419},
            {0.0392913, 0.0529509, 0.0577245, 0.0630573, 0.059229},
            {0.0398821, 0.0509655, 0.0589518, 0.0630945, 0.0604374},
            {0.0408044, 0.0534084, 0.0606451, 0.0655055, 0.061268},
            {0.0410468, 0.0552335, 0.0605044, 0.0678948, 0.0616327},
    };*/
/*
    const double efficiencies[20][5]={
            {0.000244545, 0.000367257, 0.000427767, 0.000454979, 0.000549783},
            {0.000380307, 0.00049307, 0.000597279, 0.000707326, 0.000791917},
            {0.000837471, 0.00110467, 0.00136865, 0.00178568, 0.00199086},
            {0.00214439, 0.00272176, 0.00337829, 0.00423206, 0.0046667},
            {0.00422396, 0.00559107, 0.00680265, 0.00854194, 0.00912033},
            {0.00702908, 0.00890747, 0.010636, 0.0142595, 0.0148756},
            {0.0117397, 0.0152797, 0.0186788, 0.0222518, 0.0231999},
            {0.015694, 0.020488, 0.0242975, 0.0285326, 0.0298792},
            {0.0189075, 0.0249128, 0.029717, 0.0348436, 0.0367783},
            {0.0216719, 0.0283274, 0.0338012, 0.039493, 0.0414567},
            {0.0278901, 0.0358069, 0.0416066, 0.0460792, 0.0470259},
            {0.0301887, 0.0383621, 0.0451416, 0.0503925, 0.0509475},
            {0.0325233, 0.0413149, 0.0483517, 0.0538907, 0.0539527},
            {0.0341851, 0.044189, 0.0510151, 0.0559377, 0.0559185},
            {0.0351183, 0.0454528, 0.0532184, 0.0582579, 0.0585084},
            {0.0360937, 0.0470084, 0.0549074, 0.0602768, 0.0604317},
            {0.0369486, 0.0481914, 0.0560226, 0.0619409, 0.061743},
            {0.0372465, 0.04901, 0.0577869, 0.0633726, 0.0625597},
            {0.037978, 0.0502736, 0.0586718, 0.0646702, 0.0642836},
            {0.0386916, 0.0511218, 0.0600533, 0.0658832, 0.0653667},
    };*/
    double efficiencies[5][5] = {
            {0.00110748, 0.00167546, 0.00245666, 0.0033165, 0.00366217},
            {0.00257974, 0.00351865, 0.00485276, 0.00659388, 0.00710769},
            {0.0056708, 0.00762648, 0.00958756, 0.0130879, 0.0142395},
            {0.0146233, 0.0193535, 0.0238065, 0.0291359, 0.0308294},
            {0.0300003, 0.0385753, 0.0449201, 0.0500194, 0.0506072}
    };


/*
    //https://drupal.star.bnl.gov/STAR/system/files/2018_1109_D0spectra_Note.pdf
    //D0 2016 efficiency                                0-10%       10-20%         20-40%        40-60%          60-80%
    const double efficiencies2016[11][5] = {  {0.000669, 0.000741, 0.000845, 0.000981, 0.001094},    //pT = 0-0.5
                                              {0.000741, 0.000707, 0.001018, 0.001216, 0.001305},    //pT = 0.5-1.0
                                              {0.001037, 0.001285, 0.001957, 0.002952, 0.002222},    //pT = 1.0-1.5
                                              {0.002358, 0.002978, 0.004169, 0.005517, 0.004192},    //pT = 1.5-2.0
                                              {0.004128, 0.006111, 0.007096, 0.009658, 0.008442},    //pT = 2.0-2.5
                                              {0.006282, 0.009654, 0.010210, 0.014426, 0.012034},    //pT = 2.5-3.0
                                              {0.012071, 0.016277, 0.020178, 0.024096, 0.020480},    //pT = 3.0-4.0
                                              {0.018021, 0.023850, 0.028217, 0.033696, 0.029730},    //pT = 4.0-5.0
                                              {0.026160, 0.034300, 0.042521, 0.047122, 0.033852},    //pT = 5.0-6.0
                                              {0.029773, 0.039774, 0.047944, 0.051664, 0.038528},    //pT = 6.0-8.0
                                              {0.032934, 0.044408, 0.054539, 0.058225, 0.043018} }; //pT = 8.0-10.0
*/
/*
    int pT_range =  (D0_pT < 0.5) ? 0 :
                    (D0_pT < 1.0) ? 1 :
                    (D0_pT < 1.5) ? 2 :
                    (D0_pT < 2.0) ? 3 :
                    (D0_pT < 2.5) ? 4 :
                    (D0_pT < 3.0) ? 5 :
                    (D0_pT < 3.5) ? 6 :
                    (D0_pT < 4.0) ? 7 :
                    (D0_pT < 4.5) ? 8 :
                    (D0_pT < 5.0) ? 9 :
                    (D0_pT < 5.5) ? 10 :
                    (D0_pT < 6.0) ? 11 :
                    (D0_pT < 6.5) ? 12 :
                    (D0_pT < 7.0) ? 13 :
                    (D0_pT < 7.5) ? 14 :
                    (D0_pT < 8.0) ? 15 :
                    (D0_pT < 8.5) ? 16 :
                    (D0_pT < 9.0) ? 17 :
                    (D0_pT < 9.5) ? 18 :
                    19;
*/
    int pT_range =  (D0_pT < 1.0) ? -1 :
                    (D0_pT < 1.5) ? 0 :
                    (D0_pT < 2.0) ? 1 :
                    (D0_pT < 3.0) ? 2 :
                    (D0_pT < 5.0) ? 3 :
                    4;
    int centr_range =   (centrality >= 7) ? 0 :
                        (centrality == 6) ? 1 :
                        (centrality >= 4) ? 2 :
                        (centrality >= 2) ? 3 :
                        4;

    return efficiencies[pT_range][centr_range];
}



void PrintText(int padN,double D0_pt_min,TCanvas* cdata,int centrA, int centrB, double momentA, double momentB ){
    cdata->cd(padN);
    //latex text
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.055);
    tex->DrawLatex(0.2, 0.85, "Run 14, Au+Au, #sqrt{s_{NN}} = 200 GeV");
    tex->DrawLatex(0.2, 0.75, Form("%d - %d %%, " , centrA, centrB));
    //tex->DrawLatex(0.2, 0.65, Form("%.1f GeV/c < p_{T}^{jet} < %.1f GeV/c",momentA, momentB));
    tex->DrawLatex(0.2, 0.65, Form("p_{T,D^{0}} > %.1f GeV/c", D0_pt_min));
    tex->DrawLatex(0.2, 0.55, "anti-kt algorithm, ICS");
    tex->DrawLatex(0.2, 0.45, "R = 0.4, |#eta| < 1 - R");

}
void PrintHistogram(int padN, double D0_pt_min, string observable, RooRealVar z, TCanvas* cdata, RooDataSet* data, RooStats::SPlot sData, RooDataSet* dataw_signal, RooDataSet* dataw_background,int centrA, int centrB, double momentA, double momentB){

    cdata->cd(padN);
    gPad->SetLeftMargin(0.15);
    RooPlot *frame3 = z.frame();

    data->plotOn(frame3, RooFit::LineColor(kBlack), Name(observable.c_str()));
    dataw_signal->plotOn(frame3, LineStyle(kDashed),DataError(RooAbsData::SumW2), LineColor(kRed), MarkerColor(kRed), Name("signal_sw"));
    dataw_background->plotOn(frame3, LineStyle(kDashed),DataError(RooAbsData::SumW2), LineColor(kGray), MarkerColor(kGray), Name("background_sw"));

    //lambda_1_1    lambda_1_1half  lambda_1_2  lambda_1_3 jet_pt_corr z
    string x_axis;
    string equation;
    if (observable == "z") {
        x_axis = observable;
        equation = "z = #frac{#vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} }{|#vec{p}_{T, jet}|^{2}}";
        frame3->GetXaxis()->SetRangeUser(-6,18);
        frame3->GetYaxis()->SetRangeUser(0.01,100000);
        gPad->SetLogy();

    }
    else if (observable == "jet_pt_corr"){
        x_axis = "p_{T,jet}^{sub} [GeV/c]";
        equation = "p_{T, jet}^{sub} = p_{T, jet}^{raw}  - #rho A_{jet}";
        frame3->GetXaxis()->SetRangeUser(-15,40);
        frame3->GetYaxis()->SetRangeUser(0.01,10000);
        gPad->SetLogy();
    }
    else if (observable == "d0Pt") {
        x_axis = "p_{T}^{D^{0}} [GeV/c]";
        equation = "";
        frame3->GetXaxis()->SetRangeUser(0,10);
        frame3->GetYaxis()->SetRangeUser(0.01,1000000);
        gPad->SetLogy();
    }
    else if (observable == "lambda_1_1") {
        x_axis = "#lambda_{1}^{1}";
        equation = "#lambda_{1}^{1} = #sum_{i #in jet}#left(#frac{p_{T,i}}{p_{T,jet}}#right)#left(#frac{#Delta R_{Jet,i}}{R}#right)";
        frame3->GetXaxis()->SetRangeUser(-5,13);
        frame3->GetYaxis()->SetRangeUser(0.01,10000);
        gPad->SetLogy();
    }
    else if (observable == "lambda_1_1half") {
        x_axis = "#lambda_{1.5}^{1}";
        equation = "#lambda_{1.5}^{1} = #sum_{i #in jet}#left(#frac{p_{T,i}}{p_{T,jet}}#right)#left(#frac{#Delta R_{Jet,i}}{R}#right)^{1.5}";
        frame3->GetXaxis()->SetRangeUser(-5,13);
        frame3->GetYaxis()->SetRangeUser(0.01,10000);
        gPad->SetLogy();
    }
    else if (observable == "lambda_1_2") {
        x_axis = "#lambda_{2}^{1}";
        equation = "#lambda_{2}^{1} = #sum_{i #in jet}#left(#frac{p_{T,i}}{p_{T,jet}}#right)#left(#frac{#Delta R_{Jet,i}}{R}#right)^{2}";
        frame3->GetXaxis()->SetRangeUser(-5,13);
        frame3->GetYaxis()->SetRangeUser(0.01,10000);
        gPad->SetLogy();

    }
    else if (observable == "lambda_1_3") {
        x_axis = "#lambda_{3}^{1}";
        equation = "#lambda_{3}^{1} = #sum_{i #in jet}#left(#frac{p_{T,i}}{p_{T,jet}}#right)#left(#frac{#Delta R_{Jet,i}}{R}#right)^{3}";
        frame3->GetXaxis()->SetRangeUser(-5,13);
        frame3->GetYaxis()->SetRangeUser(0.01,10000);
        gPad->SetLogy();
    }
    else if (observable == "lambda_1_0half") {
        x_axis = "#lambda_{0.5}^{1}";
        equation = "#lambda_{0.5}^{1} = #sum_{i #in jet}#left(#frac{p_{T,i}}{p_{T,jet}}#right)#left(#frac{#Delta R_{Jet,i}}{R}#right)^{0.5}";
        frame3->GetXaxis()->SetRangeUser(-5,13);
        frame3->GetYaxis()->SetRangeUser(0.01,10000);
        gPad->SetLogy();

    }
    else if (observable == "lambda_2_0") {
        x_axis = "#lambda_{2}^{0}";
        equation = "#lambda_{2}^{0} = #sum_{i #in jet}#left(#frac{p_{T,i}}{p_{T,Jet}}#right)^{2}";
        frame3->GetXaxis()->SetRangeUser(-5,13);
        frame3->GetYaxis()->SetRangeUser(0.01,10000);

        gPad->SetLogy();
    }
    else if (observable == "D0_pt") {
        x_axis = "p_{T}^{D^{0}} [GeV/c]";
        equation = "";
        frame3->GetXaxis()->SetRangeUser(-5,13);
    }
    else {
        x_axis = observable;
        equation = "";
    }

    frame3->GetXaxis()->SetTitle(x_axis.c_str());
    frame3->SetTitle("");
    frame3->GetYaxis()->SetTitle("Counts");
    //Offset
    frame3->GetYaxis()->SetTitleOffset(1.5);

    TLatex latex;
    latex.SetTextSize(0.02);
    latex.SetTextAlign(12); // Left aligned

    //Větší rozestup vlevo
    gPad->SetLeftMargin(0.15);
    //TLegend leg2(0.4008759,0.501533,0.798658,0.877719);
    TLegend leg2(0.448941,0.707053,0.84765,0.850068);
    leg2.SetTextSize(0.03);
    //leg2.AddEntry((TObject *) nullptr, "run14, Au-Au, #sqrt{s_{NN}} = 200 GeV", "");
    //leg2.AddEntry((TObject *) nullptr, Form("%d - %d %%, " , centrA, centrB), "");
    //leg2.AddEntry((TObject *) nullptr, Form("%.1f GeV/c < p_{T}^{jet} < %.1f GeV/c",momentA, momentB), "");
    //leg2.AddEntry((TObject *) nullptr, "R = 0.4, |#eta| < 1 - R", "");
    //leg2.AddEntry((TObject *) nullptr, Form("p_{T,D^{0}} > %.1f GeV/c", D0_pt_min), "");
    //leg2.AddEntry((TObject *) nullptr, "", "");
    //leg2.AddEntry(frame3->findObject(observable.c_str()), equation.c_str(), "PL"); //#lambda^{1}_{2}
    leg2.AddEntry(frame3->findObject(observable.c_str()), "raw signal", "PL"); //#lambda^{1}_{2}

    //leg2.AddEntry((TObject *) nullptr, "", "");
    leg2.AddEntry(frame3->findObject("signal_sw"), "signal (s-weighted)", "PL");
    leg2.AddEntry(frame3->findObject("background_sw"), "background (s-weighted)", "PL");
    leg2.SetBorderSize(0);
    leg2.SetFillStyle(0);
    frame3->Draw();
    //set x range
    frame3->GetXaxis()->SetRangeUser(0,10);
    //y range
    frame3->GetYaxis()->SetRangeUser(1,10000000);
    leg2.DrawClone();


//next legend
TLegend leg3(0.148941,0.807053,0.21765,0.930068);
leg3.SetTextSize(0.05);
leg3.AddEntry((TObject *) nullptr, "Raw Data", "");
//no border
leg3.SetBorderSize(0);
leg3.SetFillStyle(0);
leg3.DrawClone();

}






void splot(const char* FileInput, const char* FileOutput, double D0_pt_min, double maxPtD0Cut,  string observable, std::vector<std::vector<int>> centralityRange,std::vector<std::vector<int>> centrality,  std::vector<std::vector<double>> momenta, bool setJets, TString methodUsed){


    // Načtení stromu TNtuple s názvem "Jets_2"
    TTree* tree = openFileAndNtuple2(FileInput, "jets");
    if (!tree) return;

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

    // Definice rozsahů proměnných
    double z_range[]             = {-10000, 10000};
    double pT_range[]            = {-10000, 10000};
    double lam_1_1_range[]       = {-10000, 10000};
    double lam_1_1half_range[]   = {-10000, 10000};
    double lam_1_2_range[]       = {-10000, 10000};
    double lam_1_3_range[]       = {-10000, 10000};
    double lam_1_0half_range[]   = {-10000, 10000};
    double lam_2_0_range[]      = {-10000, 10000};
    double D0_pt_range[]         = {0, 10};

    // Načtení 2D efektivity
    LoadEfficiency2D();
    LoadEfficiency1D();
    LoadEfficiency1DPaper();
    LoadPaperD0Systematics();

    // Vytvoření readeru pro strom
    TTreeReader reader(tree);

    // Nastavení canvasu a zahájení PDF
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    TString outPdf = "./OutputPdf/" + TString(FileOutput) + ".pdf";
    c1->SaveAs(outPdf + "[");
    c1->SaveAs(outPdf);

    // Inicializace počitadla
    int c = 0;
    int d = 0;

    double mass_min = 1.75;
    double mass_max = 2.02;

    if (_systematicSPlot == 4) {
        mass_min = 1.80;
        mass_max = 1.95;
        cout << "[sPlot systematic] Using narrower fit range: "
            << mass_min << " - " << mass_max << " GeV/c^2" << endl;
    }

    if (_systematicSPlot == 5) {
        mass_min = 1.70;
        mass_max = 2.10;
        cout << "[sPlot systematic] Using wider fit range: "
            << mass_min << " - " << mass_max << " GeV/c^2" << endl;
    }

    // Smyčka přes centrality
    for (auto it_centrality = centralityRange.begin(); it_centrality != centralityRange.end(); ++it_centrality) {
        cout << "--------------------------- Centralita: " << c << " ---------------------------" << endl;
        d = 0;

        // Smyčka přes momenta
        for (auto it_momenta = momenta.begin(); it_momenta != momenta.end(); ++it_momenta) {
            if ((*it_momenta)[0] < 1 && setJets) continue;
            // Vytvoření názvu histogramu
            TString name;
            name += observable;
            name += "_hist2D_";
            name += Form("%.1f-%.1fGeV_", (*it_momenta)[0], (*it_momenta)[1]);
            name += Form("%d.0-%d.0%%", centrality[c][0], centrality[c][1]);

            // Název datového objektu
            //TString name2 = Form("Data_%.0f%d%d", (*it_momenta)[1], centrality[c][0], centrality[c][1]);
            TString name2 = Form("jets%02d_%02d_%.0f_%.0f", centrality[c][0], centrality[c][1],(*it_momenta)[0],(*it_momenta)[1]);


            // ---------------------------
            // Deklarace RooRealVar proměnných
            // ---------------------------

            // Hlavní proměnné
            RooRealVar mass("d0Mass", "Mass of D0", mass_min, mass_max);                  mass.setBins(50); // 0.000666 GeV/bin při 300 bin/GeV
            RooRealVar z("z", "z", z_range[0], z_range[1]);                     z.setBins(100);
            RooRealVar jetEta("jetEta", "jetEta", -10, 10);                jetEta.setBins(100);
            RooRealVar jetPhi("jetPhi", "jetPhi", -10, 10);                jetPhi.setBins(100);
            RooRealVar jetRapidity("jetRapidity", "jetRapidity", -10, 10);                jetRapidity.setBins(100);
            RooRealVar jetArea("jetArea", "jetArea", -10, 10);                jetArea.setBins(100);
            RooRealVar pT("jetPt", "jetPt", pT_range[0], pT_range[1]);                pT.setBins(100);
            RooRealVar jetD0DeltaR("jetD0DeltaR", "jetD0DeltaR", 0, 200);                jetD0DeltaR.setBins(100);
            RooRealVar jetNeutralPtFrac("jetNeutralPtFraction", "jetNeutralPtFraction", 0, 1);                jetNeutralPtFrac.setBins(100);
            RooRealVar nJetConst("nJetConst", "nJetConst", 0, 200);                nJetConst.setBins(100);
            RooRealVar nJetsInEvent("nJetsInEvent", "nJetsInEvent", 0, 200);                nJetsInEvent.setBins(100);


            RooRealVar gRefMult("gRefMult", "gRefMult", 0, 10000);
            RooRealVar runId("runId", "runId", 15107008, 18167014);  runId.setBins(10000); // Rozsah runId, může být upraven podle potřeby
            RooRealVar rho("rho", "rho", 0, 1000);                          rho.setBins(100000); // Rozsah rho, může být upraven podle potřeby
            RooRealVar eventId("eventId", "eventId", 0, 10000000); eventId.setBins(10000); // Rozsah eventId, může být upraven podle potřeby

            gRefMult.setBins(10000); // Globální referenční multiplicita, může být použita pro další analýzy
            // Globální referenční multiplicita, může být použita pro další analýzy

            // Lambda observables
            RooRealVar lam_1_1("lambda1_1", "lambda_{1}^{1}", lam_1_1_range[0], lam_1_1_range[1]);                     lam_1_1.setBins(100);
            RooRealVar lam_1_1half("lambda1_1_5", "lambda_{1.5}^{1}", lam_1_1half_range[0], lam_1_1half_range[1]); lam_1_1half.setBins(100);
            RooRealVar lam_1_2("lambda1_2", "lambda_{2}^{1}", lam_1_2_range[0], lam_1_2_range[1]);                     lam_1_2.setBins(100);
            RooRealVar lam_1_3("lambda1_3", "lambda_{3}^{1}", lam_1_3_range[0], lam_1_3_range[1]);                     lam_1_3.setBins(100);
            RooRealVar lam_1_0half("lambda1_0_5", "lambda_{0.5}^{1}", lam_1_0half_range[0], lam_1_0half_range[1]); lam_1_0half.setBins(100);
            RooRealVar lam_2_0("momDisp", "Momentum dispersion", lam_2_0_range[0], lam_2_0_range[1]);                     lam_2_0.setBins(100);


            // Kinematické a pomocné proměnné
            RooRealVar vD0_pT("d0Pt", "d0Pt", 0, 10);                           vD0_pT.setBins(10);
            RooRealVar D0_rap("d0Rapidity", "d0Rapidity", -2, 2);                                                  D0_rap.setBins(10000);
            RooRealVar D0_eta("d0Eta", "d0Eta", -10, 10);                                                D0_eta.setBins(10000);
            RooRealVar jet_area("jetArea", "jetArea", 0, 1);

            // Váhové proměnné
            RooRealVar w_c("weightCentrality", "weightCentrality", 0., 5e8);                                      w_c.setBins(10000);
            RooRealVar w_cErr("weightCentralityErr", "weightCentralityErr", 0., 5e8);                                      w_cErr.setBins(10000);
            RooRealVar w_ef("weightD0Efficiency", "weightD0Efficiency (1/recEff)", -1, 5e10);                                  w_ef.setBins(100000);
            RooRealVar w_dc("weightDoubleCount", "weightDoubleCount", 0, 5e8);                                        w_dc.setBins(10000);
            RooRealVar central("centrality", "centrality", 0, 10);
            RooRealVar centralityAlt("centralityAlt", "Aletrantive centrality", 0, 100);


            // ---------------------------
            // Deklarace RooArgSet proměnných
            // ---------------------------
            RooArgSet vars(runId,
                           rho,
                           eventId,
                           central,
                           centralityAlt,
                           gRefMult,
                           mass,
                           vD0_pT,
                           D0_rap,
                           D0_eta,
                           pT,
                           jetEta,
                           jetPhi,
                           jetRapidity,
                           jetArea,
                           jetD0DeltaR,
                           jetNeutralPtFrac,
                           nJetConst,
                           nJetsInEvent,
                           lam_1_0half,
                           lam_1_1,
                           lam_1_1half,
                           lam_1_2,
                           lam_1_3,
                           lam_2_0,
                           z,
                           w_c,
                           w_ef,
                           w_cErr,
                           w_dc);

            // ---------------------------
            // Vytvoření RooDataSetu (zatím bez vážení)
            // ---------------------------

            TString dsName = "data" + name2;
            RooDataSet data(dsName, "Data from TNtuple " + name2, vars);
            // Alternativa s vážením: RooDataSet data(dsName, "Data", vars, WeightVar(w));

            // ---------------------------
            // Reader pro načítání stromu
            // ---------------------------

            TTreeReader reader(tree);
            // TTreeReader reader2(ntuple2);  // zakomentováno, nepoužívá se

            // ---------------------------
            // Načítání hodnot z větví stromu
            // ---------------------------
/*
            // Hlavní fyzikální proměnné
            TTreeReaderValue<Float_t> massValue(reader, "D0mass");
            TTreeReaderValue<Float_t> zValue(reader, "z");
            TTreeReaderValue<Float_t> D0_ptValue(reader, "D0_pT");
            TTreeReaderValue<Float_t> centralityValue(reader, "centrality");
            TTreeReaderValue<Float_t> centr_weight(reader, "centr_weight");
            TTreeReaderValue<Float_t> JetPtValue(reader, "jet_pt_corr");
            TTreeReaderValue<Float_t> NpTfraction(reader, "NpTfraction");

            // Lambda observables
            TTreeReaderValue<Float_t> lam1Value(reader,      "lambda_1_1");
            TTreeReaderValue<Float_t> lam1halfValue(reader,  "lambda_1_1half");
            TTreeReaderValue<Float_t> lam2Value(reader,      "lambda_1_2");
            TTreeReaderValue<Float_t> lam3Value(reader,      "lambda_1_3");
            TTreeReaderValue<Float_t> lam0halfValue(reader,  "lambda_1_0half");
            TTreeReaderValue<Float_t> lam20Value(reader,     "lambda_2_0");

            // Další pomocné a geometrické proměnné
            TTreeReaderValue<Float_t> NJet(reader,       "NJet");
            TTreeReaderValue<Float_t> JetPtRaw(reader,   "jet_pt");
            TTreeReaderValue<Float_t> Area(reader,       "jet_area");
            TTreeReaderValue<Float_t> pseudorapidity(reader, "jet_eta");
            TTreeReaderValue<Float_t> D0_eta2(reader,    "D0_eta");
            TTreeReaderValue<Float_t> D0_rap2(reader,    "D0_rap");


*/

            // Hlavní fyzikální proměnné
                       // TString branchPrefix = "ICS_";
// Prefix pro větve (např. "ICS_" nebo "")
TString branchPrefix;
if (methodUsed == "ICS") branchPrefix = "ICS_";
else branchPrefix = "";

LoadFiles() ;
LoadEfficiency1DCut();
// --- persistentní názvy větví s prefixem ---
const TString brZ               = branchPrefix + "z";
const TString brJetPtCorr       = branchPrefix + "jetPt" + (branchPrefix == "ICS_" ? "" : "Corr");
const TString brJetNeutralPtFrac= branchPrefix + "jetNeutralPtFrac";

const TString brLam11           = branchPrefix + "lambda1_1";
const TString brLam115          = branchPrefix + "lambda1_1_5";
const TString brLam12           = branchPrefix + "lambda1_2";
const TString brLam13           = branchPrefix + "lambda1_3";
const TString brLam105          = branchPrefix + "lambda1_0_5";
const TString brMomDisp         = branchPrefix + "momDisp";

const TString brNJetConst       = branchPrefix + "nJetConst";
const TString brJetArea         = branchPrefix + "jetArea";
const TString brJetEta          = branchPrefix + "jetEta";
const TString brJetPhi          = branchPrefix + "jetPhi";
const TString brJetRapidity     = branchPrefix + "jetRapidity";
const TString brJetD0DeltaR     = branchPrefix + "jetD0DeltaR";

// --- reader values (bez dočasných TString) ---
TTreeReaderValue<Float_t> massValue(reader, "d0Mass");
TTreeReaderValue<Float_t> d0EtaValue(reader, "d0Eta");

TTreeReaderValue<Float_t> zValue(reader, brZ.Data());
TTreeReaderValue<Float_t> D0_ptValue(reader, "d0Pt");
TTreeReaderValue<Int_t>   centralityValue(reader, "centrality");
TTreeReaderValue<Int_t>   gRefMultValue(reader, "gRefMult");
TTreeReaderValue<Int_t>   runIdValue(reader, "runId");
TTreeReaderValue<Float_t> rhoValue(reader, "backgroundDensity");

TTreeReaderValue<Int_t>   eventIdValue(reader, "eventId");
TTreeReaderValue<Float_t> centr_weight(reader, "weightCentrality");

TTreeReaderValue<Float_t> JetPtValue(reader, brJetPtCorr.Data());
TTreeReaderValue<Float_t> NpTfraction(reader, brJetNeutralPtFrac.Data());

// Lambda observables
TTreeReaderValue<Float_t> lam1Value(reader,     brLam11.Data());
TTreeReaderValue<Float_t> lam1halfValue(reader, brLam115.Data());
TTreeReaderValue<Float_t> lam2Value(reader,     brLam12.Data());
TTreeReaderValue<Float_t> lam3Value(reader,     brLam13.Data());
TTreeReaderValue<Float_t> lam0halfValue(reader, brLam105.Data());
TTreeReaderValue<Float_t> lam20Value(reader,    brMomDisp.Data());

// Další pomocné a geometrické proměnné
TTreeReaderValue<Int_t>   NJet(reader, brNJetConst.Data());
TTreeReaderValue<Float_t> Area(reader, brJetArea.Data());
TTreeReaderValue<Float_t> D0_eta2(reader, "d0Eta");
TTreeReaderValue<Float_t> D0_rap2(reader, "d0Rapidity");

// Extra větve (pokud je opravdu chceš duplicitně číst zvlášť)
TTreeReaderValue<Float_t> jetEtaValue(reader, brJetEta.Data());
TTreeReaderValue<Float_t> jetPhiValue(reader, brJetPhi.Data());
TTreeReaderValue<Float_t> jetRapidityValue(reader, brJetRapidity.Data());
TTreeReaderValue<Float_t> jetAreaValue(reader, brJetArea.Data());
TTreeReaderValue<Float_t> jetD0DeltaRValue(reader, brJetD0DeltaR.Data());
TTreeReaderValue<Float_t> jetNeutralPtFracValue(reader, brJetNeutralPtFrac.Data());
TTreeReaderValue<Int_t>   nJetConstValue(reader, brNJetConst.Data());
TTreeReaderValue<Int_t>   nJetsInEventValue(reader, "nJetsInEvent");



            // ---------------------------
            // Plnění RooDataSet z TNtuple s aplikací řezů
            // ---------------------------
            while (reader.Next()) {

                // -------------------
                // Výběrové podmínky (cuts)
                // -------------------
                if (*massValue < mass_min || *massValue > mass_max)                      continue;
                if (setJets && (*D0_ptValue < D0_pt_min || *D0_ptValue > maxPtD0Cut))        continue;
                if (*centralityValue < (*it_centrality)[0] || *centralityValue > (*it_centrality)[1])                     continue;
                if (*D0_ptValue < (*it_momenta)[0] || *D0_ptValue > (*it_momenta)[1]) continue;  // TODO: vrátit zpět
                if (setJets && fabs(*jetEtaValue) > 0.6)                                 continue;          // TODO: potvrdit, že má zůstat
                if (setJets && *NpTfraction>0.95) continue;
                if (abs(*D0_rap2) > 0.6) continue; 
               // if (setJets && abs(*d0Eta)>1) continue; //DELETE

                int iCent = -1;
                if (*centralityValue < 4) iCent = 2;
                else if (*centralityValue < 7) iCent = 1;
                else iCent = 0;
                //if (setJets && (*JetPtValue < RecoJetPtMin[iCent] || *JetPtValue > RecoJetPtMax[iCent])) continue;  // TODO: vrátit zpět

                // -------------------
                // Plnění proměnných
                // -------------------
                //runid
                runId.setVal(*runIdValue);
                rho.setVal(*rhoValue);
                eventId.setVal(*eventIdValue);
                gRefMult.setVal(*gRefMultValue);
                mass.setVal(*massValue);
                z.setVal(*zValue);
                pT.setVal(*JetPtValue);

                lam_1_1.setVal(*lam1Value);
                lam_1_1half.setVal(*lam1halfValue);
                lam_1_2.setVal(*lam2Value);
                lam_1_3.setVal(*lam3Value);
                lam_1_0half.setVal(*lam0halfValue);
                lam_2_0.setVal(*lam20Value);

                jetEta.setVal(*jetEtaValue);
                jetPhi.setVal(*jetPhiValue);
                jetRapidity.setVal(*jetRapidityValue);
                jetArea.setVal(*jetAreaValue);
                jetD0DeltaR.setVal(*jetD0DeltaRValue);
                jetNeutralPtFrac.setVal(*jetNeutralPtFracValue);
                nJetConst.setVal(*nJetConstValue);
                nJetsInEvent.setVal(*nJetsInEventValue);

                vD0_pT.setVal(*D0_ptValue);
                w_c.setVal(*centr_weight);
                //w_ef.setVal(1./D0_efficiency_Neil(*D0_ptValue, *centralityValue));


              // w_ef.setVal(1./D0_efficiency_1D(*D0_ptValue,*centralityValue));
                //w_cErr.setVal(1./D0_efficiencyError_Neil(*D0_ptValue, *centralityValue));

               // w_ef.setVal(1./D0_efficiency_2D(*D0_ptValue, *D0_rap2, *centralityValue));
               // w_cErr.setVal(1./D0_efficiencyError_2D(*D0_ptValue, *D0_rap2, *centralityValue));



                //w_ef.setVal(1./D0_efficiency_1DPaper(*D0_ptValue, *centralityValue));
                /////w_ef.setVal(1./D0_efficiency(*D0_ptValue, *centralityValue));

                //if (1./D0_efficiency(*D0_ptValue, *centralityValue) < 0) exit(0);
               // cout << "D0 pT: " << *D0_ptValue << ", centrality: " << *centralityValue << ", efficiency: " << 1./D0_efficiency(*D0_ptValue, *centralityValue)<< endl;

                w_dc.setVal(1.-D0_DoubleCounting(histDC, *D0_ptValue, *centralityValue));

                central.setVal(*centralityValue);
                //přepočítám C_ID na alternativní centralitu
                //
                double centralityAltValue = -999;
                if (*centralityValue == 8) {
                    centralityAltValue = 2.5; // 0-5%
                } else if (*centralityValue == 7) {
                    centralityAltValue = 7.5; // 5-10%
                } else if (*centralityValue == 6) {
                    centralityAltValue = 15; // 10-20%
                } else if (*centralityValue == 5) {
                    centralityAltValue = 25; // 20-30%
                } else if (*centralityValue == 4) {
                    centralityAltValue = 35; // 30-40%
                } else if (*centralityValue == 3) {
                    centralityAltValue = 45; // 40-50%
                } else if (*centralityValue == 2) {
                    centralityAltValue = 55; // 50-60%
                } else if (*centralityValue == 1) {
                    centralityAltValue = 65; // 60-70%
                } else if (*centralityValue == 0) {
                    centralityAltValue = 75; // 70-80%
                }

                double rec_eff = 1./D0_efficiency_1DCut(centralityAltValue,*D0_ptValue);

                //rec_eff includes also systematics
                if (_systematicSPlot > 6 && _systematicSPlot < 21) {
                    rec_eff *= GetPaperD0SysWeight(_systematicSPlot, centralityAltValue, *D0_ptValue);
                   // cout << "[sPlot systematic] Applying additional systematic variation to efficiency: " << GetPaperD0SysWeight(_systematicSPlot, centralityAltValue, *D0_ptValue) << endl;
                }

                w_ef.setVal(rec_eff);


                //pokud je inf nebo NaN, zastavím kód:
                if (std::isinf(w_ef.getVal()) || std::isnan(w_ef.getVal())) {
                  //  cout << "Error: Efficiency is infinite or NaN for D0 pT = " << *D0_ptValue << ", D0 rap = " << *D0_rap2 << ", centrality = " << *centralityValue << endl;
                  //  exit(0); // přeskočím tento záznam
                  w_ef.setVal(0);
                }

                centralityAlt.setVal(centralityAltValue);
                D0_rap.setVal(*D0_rap2);
                D0_eta.setVal(*D0_eta2);


                // Přidání do RooDataSet (neváženě)
                data.add(vars);

                // Vážená varianta (pokud použiješ add(vars, váha)):
                // data.add(vars, w.getVal());
            }


/*

            // Vytvořím model pro signál a pozadí
            RooRealVar mean("mean", "Mean of Gaussian", 1.864, mass_min, mass_max);
            RooRealVar sigma("sigma", "Sigma of Gaussian", 0.01, 0.001, 0.1);
            RooGaussian signal("signal", "Signal PDF", mass, mean, sigma);
            */
/*
            // Pozadí (model pozadí obsahuje polynom druhého řádu)
            RooRealVar coef2("par0", "Coefficient 0", 1, -10000.0, 10000.0);
            RooRealVar coef1("par1", "Coefficient 1", 0, -10000.0, 10000.0);
            RooRealVar coef0("par2", "Coefficient 2", 0, -10000.0, 10000.0);

            RooPolynomial background("background", "Background", mass, RooArgList(coef0, coef1, coef2));
            */
/*
            RooRealVar c0("c0", "c0", 0.3, -0.5, 0.6+l);
            RooRealVar c1("c1", "c1", 0.3, -0.5, 0.5+l);
            RooRealVar c2("c2", "c2", 0.1, -0.3, 0.3+l);  // můžeš i fixnout c2 na 0
            RooChebychev background("background", "Stabilní Chebyshev background", mass, RooArgList(c0, c1));
*/
/*
            RooRealVar lambda("lambda", "slope", -2.0, -10.0, 0.0);
            RooExponential background("background", "Combinatorial Background", mass, lambda);
            */


            // Nominal: single Gaussian
            RooRealVar mean("mean", "Mean of Gaussian", 1.864, mass_min, mass_max);
            RooRealVar sigma("sigma", "Sigma of Gaussian", 0.01, 0.001, 0.1);

            RooGaussian signal("signal", "Single Gaussian signal", mass, mean, sigma);

            // Systematic variation [2]: double Gaussian with common mean
            RooRealVar sigmaRatio("sigmaRatio", "sigma_{wide}/sigma_{core}", 2.0, 1.01, 5.0);
            RooFormulaVar sigmaWide("sigmaWide", "@0*@1", RooArgList(sigma, sigmaRatio));

            RooGaussian signalCore("signalCore", "Core Gaussian signal", mass, mean, sigma);

            RooGaussian signalWide("signalWide", "Wide Gaussian signal", mass, mean, sigmaWide);

            RooRealVar fracCore("fracCore", "Core fraction", 0.8, 0.5, 0.95);

            RooAddPdf signalDouble("signalDouble", "Double Gaussian signal", RooArgList(signalCore, signalWide), RooArgList(fracCore));


            // Systematic variation: Student-t signal
            RooRealVar nuStudent("nuStudent", "Student-t #nu", 5.0, 2.1, 50.0);

            // For the first test I would fix nu.
            // Otherwise the fit can simply send nu -> very large and recover a Gaussian.
            nuStudent.setConstant(kTRUE);

            RooGenericPdf signalStudent(
                "signalStudent",
                "Student-t signal",
                "pow(1.0 + ((@0-@1)*(@0-@1))/(@2*@3*@3), -0.5*(@2+1.0))",
                RooArgList(mass, mean, nuStudent, sigma)
            );

            // Select signal model
            RooAbsPdf* signalPdf = &signal;
            TString signalLabel = "Gaussian";

            if (_systematicSPlot == 2) {
                signalPdf = &signalDouble;
                signalLabel = "DoubleGaussian";

                cout << "[sPlot systematic] Using Double Gaussian + nominal background" << endl;
            }

            if (_systematicSPlot == 3) {
                signalPdf = &signalStudent;
                signalLabel = "Student-t";

                cout << "[sPlot systematic] Using Student-t + nominal background" << endl;
            }

            // Nominal: exponential background
            RooRealVar lambda("lambda" + name2, "slope", -2.0, -10.0, 0.0);
            RooExponential backgroundExp("background_exp" + name2, "Exponential background", mass, lambda);

            // Systematic variation [1]: Chebychev 2nd order background
            RooRealVar cheb_c1("cheb_c1", "Chebychev c1", 0.0, -0.95, 0.95);
            RooRealVar cheb_c2("cheb_c2", "Chebychev c2", 0.0, -0.95, 0.95);
            RooChebychev backgroundCheb2("backgroundCheb2", "Chebychev 2nd order background", mass, RooArgList(cheb_c1, cheb_c2));


            // Definuji váhovací koeficienty pro signál a pozadí
            RooRealVar n_signal("n_signal" + name2, "Number of signal events", 22100, 0, 1000000);



            // Select background model
            RooAbsPdf* backgroundPdf = &backgroundExp;
            TString backgroundName = backgroundExp.GetName();
            TString backgroundLabel = "Exponential";

            if (_systematicSPlot == 1) {
                backgroundPdf = &backgroundCheb2;
                backgroundName = backgroundCheb2.GetName();
                backgroundLabel = "Chebychev2";

                cout << "[sPlot systematic] Using Gaussian + Chebychev2 background" << endl;
            } else {
                cout << "[sPlot nominal] Using Gaussian + Exponential background" << endl;
            }


            RooRealVar n_background("n_background" + name2, "Number of background events", 205200, 0, 1000000);

            // Vytvořím seznam modelů a seznam váhovacích koeficientů
            RooArgList pdfs(*signalPdf, *backgroundPdf);
            RooArgList coeffs(n_signal, n_background);

            // Vytvořím vážený model
            RooAddPdf total_pdf_weighted("total_pdf_weighted" + name2, "Total PDF with weights", pdfs, coeffs);

            //Provedu fit
            total_pdf_weighted.fitTo(data, RooFit::Extended(kTRUE), RooFit::Save());
            //total_pdf_weighted.fitTo(data, RooFit::Save());

            mean.setConstant(kTRUE);
            sigma.setConstant(kTRUE);

            if (_systematicSPlot == 3) {
                nuStudent.setConstant(kTRUE);
            }

            if (_systematicSPlot == 2) {
                sigmaRatio.setConstant(kTRUE);
                fracCore.setConstant(kTRUE);
            }

            if (_systematicSPlot == 1) {
                cheb_c1.setConstant(kTRUE);
                cheb_c2.setConstant(kTRUE);
            } else {
                lambda.setConstant(kTRUE);
            }


            double mean_value = mean.getVal();
            double mean_error = mean.getError();
            double sigma_value = sigma.getVal();
            double sigma_error = sigma.getError();

            double lowerMass = mean_value - 3.0 * sigma_value;
            double upperMass = mean_value + 3.0 * sigma_value;

            // Vytvořím workspace pro práci se splotem
            RooWorkspace wspace{"myWS" + name2};

            //Vytvořím kopii data (proč?)
            RooDataSet data_copy = data;
/*
            RooDataSet dataNarrow("dataNarrow" + name2, "cut around peak", vars);

// Získám binning z RooRealVar
            const RooAbsBinning& binning = mass.getBinning("cache");

            int binLow = binning.binNumber(mean_value - 6.0 * sigma_value);
            int binHigh = binning.binNumber(mean_value + 6.0 * sigma_value);

            double lowerCut = binning.binLow(binLow);
            double upperCut = binning.binHigh(binHigh);  // už není třeba +1
            mass.setRange(lowerCut, upperCut);

            for (int i = 0; i < data.numEntries(); ++i) {
                const RooArgSet* row = data.get(i);
                double m = ((RooRealVar*)row->find("d0Mass"))->getVal();

                if (m >= lowerCut && m < upperCut)
                    dataNarrow.add(*row);
            }

            // Zafixuj shape parametry
            mass.setRange("fitrange_narrow", mean_value - 3*sigma_value, mean_value + 3*sigma_value);
*/


            /*
            // Refit jen výnosy
            total_pdf_weighted.fitTo(dataNarrow,
                                     RooFit::Extended(kTRUE),
                                     RooFit::Range("reduced"),
                                     RooFit::SumCoefRange("reduced"),
                                     RooFit::Save());


            RooStats::SPlot sData("sData" + name2, "SPlot result", dataNarrow, &total_pdf_weighted,
                                  RooArgList(n_signal, n_background));
            */
            // Provádím vážený splot
            RooStats::SPlot sData("sData" + name2, "SPlot result", data_copy, &total_pdf_weighted,
                                  RooArgList(n_signal, n_background));
            /*RooStats::SPlot sData("sData_cut", "SPlot on cut", *data_cut, &total_pdf_weighted,
                                  RooArgList(n_signal, n_background));*/

            std::cout << ">>> Check SWeights:\n\n";
            std::cout << "Yield of D0 is "    << n_signal.getVal()
                      << ". From sWeights it is " << sData.GetYieldFromSWeight("n_signal" + name2 + "_sw") << "\n";
            std::cout << "Yield of BKG is "   << n_background.getVal()
                      << ". From sWeights it is " << sData.GetYieldFromSWeight("n_background" + name2 + "_sw") << "\n\n";

            cout << "--------------------------------------------------"<< endl;
            sData.Print();
            // create weighted data sets

           RooDataSet dataw_signal{data_copy.GetName(), data_copy.GetTitle(), &data_copy, *data_copy.get(), nullptr,
                                    "n_signal" + name2 + "_sw"};
           /* RooDataSet dataw_signal{dataNarrow.GetName(), dataNarrow.GetTitle(), &dataNarrow, *dataNarrow.get(), nullptr,
                                    "n_signal" + name2 + "_sw"};*/
            //Rebinning(hist_signal_sw_th1d, observable);
            //--------------------
            RooDataSet* clonedData = dynamic_cast<RooDataSet*>(dataw_signal.Clone("clonedData"));





            // Convert cloned dataset to tree store
            dataw_signal.convertToTreeStore();

            // Create a new ROOT file for the tree
            TFile* treeFile = new TFile("./Output/"+TString(FileOutput)+".root", (c==0&&d==0)?"RECREATE":"UPDATE");

            //Prvně si ten strom načtu
            TTree* tree_unbin = const_cast<TTree*>(dataw_signal.tree());


            //Přejmenuji ten strom
            tree_unbin->SetName(name2);

            // Write the tree to the file
            // Write the tree to the file
            tree_unbin->Write();

/*
// Získám původní strom z RooDataSet
TTree* oldTree = const_cast<TTree*>(clonedData->tree());

// Připravím nový strom s novými větvemi
float observableVar, weight;
TTree* newTree = new TTree(name2, "Tree with renamed branch");
newTree->Branch("observable", &observableVar, "observable/F");
newTree->Branch("sWeightSignal", &weight, "sWeightSignal/F");

// Přiřazení adres původních větví
oldTree->SetBranchAddress("observable", &observableVar);
oldTree->SetBranchAddress(("n_signal" + name2 + "_sw").Data(), &weight);

// Naplním nový strom
Long64_t nentries = oldTree->GetEntries();
for (Long64_t i = 0; i < nentries; ++i) {
    oldTree->GetEntry(i);
    newTree->Fill();
}

// Uložení
newTree->Write();

 * */

            std::cout << ">>> Check 2 SWeights:\n\n";
            std::cout << "Yield of D0 is "    << n_signal.getVal()
                      << ". From sWeights it is " << sData.GetYieldFromSWeight("n_signal" + name2 + "_sw") << "\n";
            std::cout << "Yield of BKG is "   << n_background.getVal()
                      << ". From sWeights it is " << sData.GetYieldFromSWeight("n_background" + name2 + "_sw") << "\n\n";
            // Close the file
            treeFile->Close();
            //delete treeFile;
            delete clonedData;



            cout << "--------------------------------------------------"<< endl;

            RooDataSet dataw_background{data_copy.GetName(), data_copy.GetTitle(), &data_copy, *data_copy.get(),
                                        nullptr,
                                        "n_background" + name2 + "_sw"};


 //dataNarrow
/*
            RooDataSet dataw_background{dataNarrow.GetName(), dataNarrow.GetTitle(), &dataNarrow, *dataNarrow.get(),
                                        nullptr,
                                        "n_background" + name2 + "_sw"};*/
            cout << "Mean: " << mean_value << " +- " << mean_error << endl;
            cout << "Sigma: " << sigma_value << " +- " << sigma_error << endl;
            //Vypočítám signifikanci v oblasti píku mu +/- 5 sigma
            double n_signal_peak = dataw_signal.sumEntries(
                    Form("d0Mass > %f - 5*%f && d0Mass < %f + 5*%f", mean_value, sigma_value, mean_value, sigma_value));
            double n_background_peak = dataw_background.sumEntries(
                    Form("d0Mass > %f - 5*%f && d0Mass < %f + 5*%f", mean_value, sigma_value, mean_value, sigma_value));
            double significance = n_signal_peak / sqrt(n_signal_peak + n_background_peak);
            cout << "Significance: " << significance << endl;

            //------------------------------------------------------------------------
            TCanvas *cdata;
            bool dvoji = true;
            bool vse = false;
            if (vse) {
                cdata = new TCanvas("sPlot" + name, "sPlot_" + name, 2000, 2000);
                cdata->Divide(3, 3);
                cdata->cd(1);
            } else if (dvoji){
                cdata = new TCanvas("sPlot" + name, "sPlot_" + name, 2400, 1200);
                cdata->Divide(2, 1);
                cdata->cd(1);
            }else{
                cdata = new TCanvas("cdata", "cdata", 1600, 900);
            }

            RooDataSet* sWeightedData = sData.GetSDataSet();
            gPad->SetLeftMargin(0.17);
            RooPlot *frame_mass = mass.frame();

// Použij správný dataset pro plotování
            //plotData->plotOn(frame_mass, DataError(RooAbsData::SumW2), LineColor(kBlack));
            sWeightedData->plotOn(frame_mass, DataError(RooAbsData::SumW2), LineColor(kBlack));

            total_pdf_weighted.plotOn(frame_mass, Name("mass"));
            //frame_mass->GetYaxis()->SetRangeUser(0.01, 100); // Nastavení rozsahu y-ové osy
            frame_mass->GetYaxis()->SetRangeUser(0, 2 * frame_mass->GetMaximum());
            frame_mass->GetXaxis()->SetRangeUser(1.78, 2.00); // Nastavení rozsahu x-ové osy
            //logscale
            total_pdf_weighted.plotOn(frame_mass, Components(*signalPdf), LineStyle(kDashed), LineColor(kRed), Name("signal"));
            total_pdf_weighted.plotOn(frame_mass, Components(*backgroundPdf), LineStyle(kDashed), LineColor(kGreen), Name("background"));


            //Nastavím titulek

            // Nastavení titulku pro graf
            frame_mass->SetTitle("Mass distribution" + name);
            frame_mass->SetTitle("");
// Nastavení popisků os
            frame_mass->GetXaxis()->SetTitle("m_{D^{0}} [GeV/c^{2}]");
            frame_mass->GetYaxis()->SetTitle("Counts");
            //set offset

            //set offset
            frame_mass->GetYaxis()->SetTitleOffset(1.55);

            //Set x range, skipp first and last bin
            ////frame_mass->GetXaxis()->SetRangeUser(1.8, 1.96);

            //frame_mass->GetYaxis()->SetTitle("Count");
            TLegend leg(0.565181,0.543379,0.790338,0.886986);
            //Nastavím velikost písma
            leg.SetTextSize(0.025);
            //setmargin
            leg.SetMargin(0.25);
            //Přidám pouze text do první řádku, budu odkaovat na nulový pointer
            leg.AddEntry((TObject *) nullptr, "STAR Run 14, Au+Au", "");
            leg.AddEntry((TObject *) nullptr, Form("#sqrt{s_{NN}} = 200 GeV, %d - %d%%",centrality[c][0], centrality[c][1]), "");
            ////leg.AddEntry((TObject *) nullptr, Form("%d - %d %%, %.0f < p_{T, Jet}^{sub} [GeV/c] < %.0f" , centrality[c][0], centrality[c][1],momenta[d][0], momenta[d][1]), "");
           //// leg.AddEntry((TObject *) nullptr, Form("inclusive D^{0}-jets, %.0f < p_{T, D^{0}} [GeV/c] < 10", D0_pt_min ), "");
           // leg.AddEntry((TObject *) nullptr, Form("%d - %d %%, %.0f < p_{T, D^{0}} [GeV/c] < 10",centrality[c][0], centrality[c][1], D0_pt_min ), "");
          //  leg.AddEntry((TObject *) nullptr, Form("%.0f < p_{T}^{D^{0}} [GeV/c] < 10", 0 ), "");
            leg.AddEntry((TObject *) nullptr, Form("%.0f < p_{T}^{D^{0}} [GeV/c] < %.0f", (*it_momenta)[0], (*it_momenta)[1] ), "");

            //leg.AddEntry((TObject *) nullptr, "anti-k_{T}, R = 0.4, |#eta_{jet}| < 1 - R", "");
            //leg.AddEntry((TObject *) nullptr, Form("p_{T,D^{0}} > %.1f GeV/c", D0_pt_min), "");
            TString sig_name = "Gaussian";
            TString bg_name  = "Exponential";

            if (_systematicSPlot == 1) {
                bg_name = "Chebychev2";
            }

            if (_systematicSPlot == 2) {
                sig_name = "Double Gaussian";
            }
            if (_systematicSPlot == 3) sig_name = "Student-t";

            TString fullFitLabel = sig_name + " + " + bg_name;

            leg.AddEntry(frame_mass->findObject("mass"),fullFitLabel.Data(),"L");

            leg.AddEntry(frame_mass->findObject("signal"),sig_name.Data(),"L");

            leg.AddEntry(frame_mass->findObject("background"),bg_name.Data(),"L");

            //Dopíšu signifikanci
            leg.AddEntry((TObject *) nullptr, Form("#mu = (%.1f #pm %.1f) MeV/c^{2}", mean_value*1000, mean_error*1000), "");
            /*
            leg.AddEntry((TObject *) nullptr,
                         Form("#sigma = (%.1f #pm %.1f) MeV/c^{2} ", sigma_value * 1000, sigma_error * 1000), "");
                         */
            if (_systematicSPlot == 2) {
            leg.AddEntry((TObject *) nullptr,
                        Form("#sigma_{core} = (%.1f #pm %.1f) MeV/c^{2}",
                            sigma_value * 1000, sigma_error * 1000), "");

            leg.AddEntry((TObject *) nullptr,
                        Form("#sigma_{wide} = %.1f MeV/c^{2}",
                            sigmaWide.getVal() * 1000), "");

            leg.AddEntry((TObject *) nullptr,
                        Form("f_{core} = %.2f", fracCore.getVal()), "");
            } else if (_systematicSPlot == 3) {
                leg.AddEntry((TObject *) nullptr,
                            Form("#sigma_{scale} = (%.1f #pm %.1f) MeV/c^{2}",
                                sigma_value * 1000, sigma_error * 1000), "");

                leg.AddEntry((TObject *) nullptr,
                            Form("#nu = %.1f", nuStudent.getVal()), "");
            } else {
            leg.AddEntry((TObject *) nullptr,
                        Form("#sigma = (%.1f #pm %.1f) MeV/c^{2}",
                            sigma_value * 1000, sigma_error * 1000), "");
            }
                    //leg.AddEntry((TObject *) nullptr, Form("S = %.0f, B = %.0f", n_signal_peak, n_background_peak), "");
            leg.AddEntry((TObject *) nullptr, Form("Significance = %.2f", significance), "");
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
/*
            RooHist *hist = (RooHist*) frame_mass->getObject(0);  // 0 je index objektu
            if (hist) {
                for (int i = 0; i < hist->GetN(); i++) {
                    double x, y;
                    hist->GetPoint(i, x, y);
                    cout << "Bin center: " << x << ", obsahuje: " << y << endl;

            } else {
                cout << "RooHist nebyl nalezen!" << endl;
            }
*/
            frame_mass->Draw();
            leg.DrawClone();
            TLegend leg3(0.178941,0.807053,0.23765,0.930068);
            leg3.SetTextSize(0.05);
            leg3.AddEntry((TObject *) nullptr, "Raw Data", "");
//no border
            leg3.SetBorderSize(0);
            leg3.SetFillStyle(0);
            leg3.DrawClone();
            if (vse) {

               // PrintHistogram(2, D0_pt_min, "D0_pT", vD0_pT,cdata, &data, sData, &dataw_signal, &dataw_background,centrality[c][0],centrality[c][1],momenta[d][0], momenta[d][1]);
                PrintHistogram(2, D0_pt_min, "z", z,cdata, &data, sData, &dataw_signal, &dataw_background,centrality[c][0],centrality[c][1],momenta[d][0], momenta[d][1]);
            PrintHistogram(3, D0_pt_min, "jet_pt_corr", pT,cdata, &data, sData, &dataw_signal, &dataw_background,centrality[c][0],centrality[c][1],momenta[d][0], momenta[d][1]);
            //PrintText(4,D0_pt_min,cdata,centrality[c][0],centrality[c][1],momenta[d][0], momenta[d][1]);
            PrintHistogram(4, D0_pt_min, "lambda_1_1", lam_1_1,cdata, &data, sData, &dataw_signal, &dataw_background,centrality[c][0],centrality[c][1],momenta[d][0], momenta[d][1]);
            PrintHistogram(5, D0_pt_min, "lambda_1_1half", lam_1_1half,cdata, &data, sData, &dataw_signal, &dataw_background,centrality[c][0],centrality[c][1],momenta[d][0], momenta[d][1]);
            PrintHistogram(6, D0_pt_min, "lambda_1_2", lam_1_2,cdata, &data, sData, &dataw_signal, &dataw_background,centrality[c][0],centrality[c][1],momenta[d][0], momenta[d][1]);
            PrintHistogram(7, D0_pt_min, "lambda_1_3", lam_1_3,cdata, &data, sData, &dataw_signal, &dataw_background,centrality[c][0],centrality[c][1],momenta[d][0], momenta[d][1]);
            PrintHistogram(8, D0_pt_min, "lambda_1_0half", lam_1_0half,cdata, &data, sData, &dataw_signal, &dataw_background,centrality[c][0],centrality[c][1],momenta[d][0], momenta[d][1]);
            PrintHistogram(9, D0_pt_min, "lambda_2_0", lam_2_0,cdata, &data, sData, &dataw_signal, &dataw_background,centrality[c][0],centrality[c][1],momenta[d][0], momenta[d][1]);
            } if (dvoji){
                cdata->cd(2);
              //  PrintHistogram(1, D0_pt_min, "D0_pt", vD0_pT,cdata, data, sData, dataw_signal, dataw_background,centrality[c][0],centrality[c][1],momenta[d][0], momenta[d][1]);
              //  PrintHistogram(9, D0_pt_min, "lambda_2_0", lam_2_0,cdata, data, sData, dataw_signal, dataw_background,centrality[c][0],centrality[c][1],momenta[d][0], momenta[d][1]);
                PrintHistogram(9, D0_pt_min, "d0Pt", vD0_pT, cdata, &data, sData, &dataw_signal, &dataw_background, centrality[c][0],centrality[c][1],momenta[d][0], momenta[d][1]);

            }
            // Assuming you have already calculated sWeights using SPlot instance 'sPlot'
            ////RooDataSet* sWeightedData = sData.GetSDataSet();

// Create a new dataset for unbinarized data and their weights
            //RooDataSet unbinnedResult("unbinnedResult" + name, "Unbinned result from SPlot", vars);
            cdata->SaveAs("./OutputPdf/"+TString(FileOutput)+".pdf");
            //clean frame_mass
            frame_mass->Clear();
            d++;
            //delete cdata;
            //delete frame_mass;
        }
        c++;
    }
    c1->SaveAs("./OutputPdf/"+TString(FileOutput)+".pdf]");


}
// Pomocná rekurzivní funkce pro kopírování objektů (stromy, histogramy, složky atd.)
void CopyOnlyFoldersAndHistograms(TDirectory* source, TDirectory* dest) {
    TIter nextkey(source->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)nextkey())) {
        TObject* obj = key->ReadObj();

        if (obj->InheritsFrom("TDirectory")) {
            // Rekurzivně kopírujeme podadresář
            TDirectory* newDir = dest->mkdir(obj->GetName());
            source->cd(obj->GetName());
            CopyOnlyFoldersAndHistograms((TDirectory*)obj, newDir);
            dest->cd();
        }
        else if (obj->InheritsFrom("TH1")) {
            // Klonujeme histogram, abychom neměli dva vlastníky téhož objektu
            dest->cd();
            TH1* hClone = dynamic_cast<TH1*>(obj->Clone(obj->GetName()));
            if (hClone) {
                hClone->SetDirectory(nullptr);  // Odpojit od gDirectory
                hClone->Write();
            }
        }
        else if (obj->InheritsFrom("TProfile")) {
            // Stejně jako histogramy – clone + zápis
            dest->cd();
            TProfile* pClone = dynamic_cast<TProfile*>(obj->Clone(obj->GetName()));
            if (pClone) {
                pClone->SetDirectory(nullptr);
                pClone->Write();
            }
        }

        delete obj;  // Uvolníme objekt načtený ze vstupu (ne clone)
    }
}



void RenameTreeAndMerge(const char* FileInput, const char* FileOutput,
                const std::vector<std::vector<int>>& centralityRange,
                const std::vector<std::vector<int>>& centrality, TString methodUsed, TString sys)
{
    TFile* inputFile = TFile::Open("./Output/"+outRoot+methodUsed+sys+".root", "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error opening input file: " << FileInput << std::endl;
        return;
    }

    //inputFile->ls();

    TFile* outputFile = TFile::Open("./Output/"+outRoot+methodUsed+sys+"2.root", "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error creating output file: " << FileOutput << std::endl;
        inputFile->Close();
        return;
    }
    outputFile->cd();  // ⬅️ důležité!
    TTree* mergedTree = new TTree("jets", "All centralities merged");
    TFile* realInput = TFile::Open(RealJetsFileData, "READ");
    if (!realInput || realInput->IsZombie()) {
        std::cerr << "Error creating realInput file: " << FileOutput << std::endl;
        realInput->Close();
        return;
    }

    // Persistentní proměnné napříč všemi stromy
    std::map<TString, Double_t> doubleVars;
    std::map<TString, Float_t> floatVars;
    std::map<TString, Int_t> intVars;

    Double_t sWeightSigVal = 0;
    Double_t sWeightSigLVal = 0;
    Double_t sWeightBgVal = 0;
    Double_t sWeightBgLVal = 0;

    Double_t correctedWeight = 0;



    bool branchesDefined = false;

    for (size_t c = 0; c < centralityRange.size(); ++c) {
        int d = 0;
        for (auto it_momenta = momenta.begin(); it_momenta != momenta.end(); ++it_momenta) {

            //TString name2 = Form("jets%02d_%02d", centrality[c][0], centrality[c][1]);
            TString name2 = Form("jets%02d_%02d_%.0f_%.0f", centrality[c][0], centrality[c][1],(*it_momenta)[0],(*it_momenta)[1]);


            TTree *tree = (TTree *) inputFile->Get(name2);
            if (!tree) {
                std::cerr << "Error: Tree " << name2 << " not found in input file." << std::endl;
                continue;
            }

            // Připojit běžné větve (kromě těch, které přejmenujeme)
            TObjArray *branches = tree->GetListOfBranches();
            for (int i = 0; i < branches->GetEntries(); ++i) {
                TBranch *br = (TBranch *) branches->At(i);
                TString brName = br->GetName();

                // Přeskočíme ty, které budeme přejmenovávat
                if (brName.BeginsWith("n_sig") || brName.BeginsWith("L_n_sig") ||
                    brName.BeginsWith("n_bac") || brName.BeginsWith("L_n_bac")) {
                    continue;
                }

                TLeaf *leaf = br->GetLeaf(brName);
                if (!leaf) continue;

                TString type = leaf->GetTypeName();

                if (type == "Double_t") {
                    doubleVars[brName] = 0;
                    tree->SetBranchAddress(brName, &doubleVars[brName]);
                    if (!branchesDefined) mergedTree->Branch(brName, &doubleVars[brName], brName + "/D");
                } else if (type == "Float_t") {
                    floatVars[brName] = 0;
                    tree->SetBranchAddress(brName, &floatVars[brName]);
                    if (!branchesDefined) mergedTree->Branch(brName, &floatVars[brName], brName + "/F");
                } else if (type == "Int_t") {
                    intVars[brName] = 0;
                    tree->SetBranchAddress(brName, &intVars[brName]);
                    if (!branchesDefined) mergedTree->Branch(brName, &intVars[brName], brName + "/I");
                }
            }

            // Připojit sWeight větve
            TBranch *br_nsig = nullptr;
            TBranch *br_lnsig = nullptr;
            TBranch *br_nbac = nullptr;
            TBranch *br_lnbac = nullptr;

            TObjArray *allBranches = tree->GetListOfBranches();
            for (int i = 0; i < allBranches->GetEntries(); ++i) {
                TString bname = allBranches->At(i)->GetName();
                if (bname.BeginsWith("n_sig")) br_nsig = (TBranch *) allBranches->At(i);
                if (bname.BeginsWith("L_n_sig")) br_lnsig = (TBranch *) allBranches->At(i);
                if (bname.BeginsWith("n_bac")) br_nbac = (TBranch *) allBranches->At(i);
                if (bname.BeginsWith("L_n_bac")) br_lnbac = (TBranch *) allBranches->At(i);
            }

            if (br_nsig) tree->SetBranchAddress(br_nsig->GetName(), &sWeightSigVal);
            if (br_lnsig) tree->SetBranchAddress(br_lnsig->GetName(), &sWeightSigLVal);
            if (br_nbac) tree->SetBranchAddress(br_nbac->GetName(), &sWeightBgVal);
            if (br_lnbac) tree->SetBranchAddress(br_lnbac->GetName(), &sWeightBgLVal);


            if (!branchesDefined) {
                mergedTree->Branch("sWeightSignal", &sWeightSigVal, "sWeightSignal/D");
                mergedTree->Branch("sWeightSigLikelihood", &sWeightSigLVal, "sWeightSigLikelihood/D");
                mergedTree->Branch("sWeightBackground", &sWeightBgVal, "sWeightBackground/D");
                mergedTree->Branch("sWeightBgLikelihood", &sWeightBgLVal, "sWeightBgLikelihood/D");
               // mergedTree->Branch("correctedWeight", &correctedWeight, "correctedWeight/D");
                branchesDefined = true;
            }

            // Plnění stromu
            Long64_t nEntries = tree->GetEntries();
            for (Long64_t i = 0; i < nEntries; ++i) {
                tree->GetEntry(i);
                mergedTree->Fill();
            }
        }

        d++;
    }


    outputFile->cd();
    //mergedTree->Write();
    mergedTree->Write("", TObject::kOverwrite);

    CopyOnlyFoldersAndHistograms(realInput, outputFile);
    cout << "closing 1" << endl;
    outputFile->Close();
    //delete outputFile;
    cout << "closing 2" << endl;
    inputFile->Close();
    //delete inputFile;
    cout << "both closed" << endl;

    realInput->Close();

    //delete inputFile;
    //delete outputFile;
    //delete realInput;


}




    void Simple_splot(const char* InputFileIn = 0, const char* OutputFile = "Output", const char* OutputFile2 = "Output2", double minPtD0Cut =1, double maxPtD0Cut = 10, int systematicSPlot = 0)
    {

    TString InputFile;     
    if (InputFileIn) {
        InputFile = InputFileIn;
    } else {
        InputFile = RealJetsFileData;
    }

    _systematicSPlot = systematicSPlot;
    // _systematicSPlot = 0; // nominal: Gaussian + Exponential
    // _systematicSPlot = 1; // background variation: Gaussian + Chebychev 2nd order
    // _systematicSPlot = 2; // signal variation: Double Gaussian + Exponential
    // _systematicSPlot = 3  // Student-t signal + Exponential background
    // _systematicSPlot = 4; // narrower fit range
    // _systematicSPlot = 5; // wider fit range
    // _systematicSPlot = 6; // keep negative bins 

    //Transfering vector of vectors to vector of ranges
        std::vector<std::vector<int>> centrality = {{0, 10}, {10, 20}, {20, 30}, {30, 40}, {40, 80}}; 		                        //Centrality ranges
      //  std::vector<std::vector<int>> centrality = { {10, 20}}; 		                        //Centrality ranges
     //   std::vector<std::vector<int>> centrality = {{0, 10}, {10, 40}, {40,80}}; 		                        //Centrality ranges

        std::vector<std::vector<int>> centralityRange(centrality.size());
   std::transform(centrality.begin(), centrality.end(), centralityRange.begin(), CentrRangeTransf);

    TString sys = "";
    //sPlot variations
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

    splot(InputFile, outRoot+"_ICS"+sys, minPtD0Cut, maxPtD0Cut,"z",centralityRange, centrality, momenta, true, "ICS");
    RenameTreeAndMerge(OutputFile, OutputFile2, centralityRange, centrality,"_ICS",sys);

    splot(InputFile, outRoot+"_AREA"+sys, minPtD0Cut, maxPtD0Cut,"z",centralityRange, centrality, momenta, true, "AREA");
    RenameTreeAndMerge(OutputFile, OutputFile2, centralityRange, centrality,"_AREA",sys);



    //přidám _D0 k outRoot

    outRoot += "_D0";
    splot(InputFile, outRoot+sys, minPtD0Cut, maxPtD0Cut,"z",centralityRange, centrality, momenta, false, "");
    RenameTreeAndMerge(OutputFile, OutputFile2, centralityRange, centrality,"",sys);

        std::cout << "Hard exit" << std::endl;
        ::_exit(0);

        /*
    gROOT->Reset(); // brutální, ale efektivní reset paměťových struktur ROOTu
*/
   }
