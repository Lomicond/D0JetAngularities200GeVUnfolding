#ifndef CONFIG_HIST_H
#define CONFIG_HIST_H

#include <vector>
#include "TLine.h"
//include congif.h
#include "config.h"


#endif // CONFIG_H

struct StJetTreeStruct2
{
    Double_t z_value;
    Double_t pT_value;
    Double_t lambda_value[6];
    Double_t s_weight_value;
    Double_t s_weight_error_value;
    Double_t centr_weight_value;
    Double_t eff_weight_value_error;
    Double_t eff_weight_value;
    Double_t jetD0DeltaR;
    Double_t doubleCount;
    Double_t central;
    Double_t centralAlt;
    Double_t gRefMultVal;
    Double_t nJetConst;
    Double_t jet_area;
    Double_t bg_dens;
    Double_t pT_raw;
    Double_t mass;
    Double_t vd0pt;
    Double_t D0_eta;
    Double_t jetEta;
};

struct StJetTreeStruct
{
    float refmult;
    float grefmult;
    float centrality;
    float centralityAlt;
    float eventMaxPtTrack;
    float refcorr2;
    float mcrefmult;
    float weight;
    float jetpt;
    float jetpt_nocorr;
    float recoJetArea;
    float recoJetRho;
    float jetcorrectedpt;
    float mcJetEta;
    float jetphi;
    float jetarea;
    float jetradius;
    float jetenergy;
    float fRhoValforjet;
    int numberofconstituents;
    //int numberofconstituents;

    float recoJetEta;
    float RecoD0Pt;
    float RecoD0Eta;
    float McRecoJetEta;
    float RecoPionEta;
    float RecoKaonEta;
    float d0z;
    float d0mass;
    float d0pt;
    float d0eta;
    float gRefMultMc;
    float d0phi;
    float pionpt;
    float pioneta;
    float pionphi;
    float kaonpt;
    float kaoneta;
    float kaonphi;
    float lambda[6];

};
/*
void assignTree(TTree *jetTree, StJetTreeStruct &mcJet, StJetTreeStruct &recoJet)
{
    cout << "Reading tree" << endl;

    jetTree->SetBranchAddress("Centrality", &recoJet.centrality);
    jetTree->SetBranchAddress("Weight", &recoJet.weight);
    jetTree->SetBranchAddress("McJetPt", &mcJet.jetpt);
    jetTree->SetBranchAddress("McJetLambda_1_1", &mcJet.lambda[0]);
    jetTree->SetBranchAddress("McJetLambda_1_1half", &mcJet.lambda[1]);
    jetTree->SetBranchAddress("McJetLambda_1_2", &mcJet.lambda[2]);
    jetTree->SetBranchAddress("McJetLambda_1_3", &mcJet.lambda[3]);
    jetTree->SetBranchAddress("McJetLambda_1_half", &mcJet.lambda[4]);
    jetTree->SetBranchAddress("McJetDispersion", &mcJet.lambda[5]);
    jetTree->SetBranchAddress("McJetD0Z", &mcJet.d0z);
    jetTree->SetBranchAddress("RecoD0Pt", &recoJet.RecoD0Pt);
    jetTree->SetBranchAddress("RecoJetPt", &recoJet.jetpt);
    jetTree->SetBranchAddress("RecoJetNConst", &recoJet.numberofconstituents);
    jetTree->SetBranchAddress("RecoJetLambda_1_1", &recoJet.lambda[0]);
    jetTree->SetBranchAddress("RecoJetLambda_1_1half", &recoJet.lambda[1]);
    jetTree->SetBranchAddress("RecoJetLambda_1_2", &recoJet.lambda[2]);
    jetTree->SetBranchAddress("RecoJetLambda_1_3", &recoJet.lambda[3]);
    jetTree->SetBranchAddress("RecoJetLambda_1_half", &recoJet.lambda[4]);
    jetTree->SetBranchAddress("RecoJetDispersion", &recoJet.lambda[5]);
    jetTree->SetBranchAddress("RecoJetEta", &recoJet.recoJetEta);
    jetTree->SetBranchAddress("RecoJetD0Z", &recoJet.d0z);
    jetTree->SetBranchAddress("RecoD0Eta", &recoJet.RecoD0Eta);
    // RecoPionEta
    jetTree->SetBranchAddress("RecoPionEta", &recoJet.RecoPionEta);
    //RecoKaonEta
    jetTree->SetBranchAddress("RecoKaonEta", &recoJet.RecoKaonEta);

}
*/
//TESTING

void assignTree(TTree *jetTree, StJetTreeStruct &mcJet, StJetTreeStruct &recoJet)
{
    cout << "Reading tree" << endl;
//Int_t nnum = 2;
    jetTree->SetBranchAddress("centrality", &recoJet.centrality);
    jetTree->SetBranchAddress("centralityAlt", &recoJet.centralityAlt);
    jetTree->SetBranchAddress("eventMaxPtTrack", &recoJet.eventMaxPtTrack);

    jetTree->SetBranchAddress("weightCentrality", &recoJet.weight);
    jetTree->SetBranchAddress("mcJetPt", &mcJet.jetpt);
    jetTree->SetBranchAddress("mcJetLambda1_1", &mcJet.lambda[0]);
    jetTree->SetBranchAddress("mcJetLambda1_1_5", &mcJet.lambda[1]);
    jetTree->SetBranchAddress("mcJetLambda1_2", &mcJet.lambda[2]);
    jetTree->SetBranchAddress("mcJetLambda1_3", &mcJet.lambda[3]);
    jetTree->SetBranchAddress("mcJetLambda1_0_5", &mcJet.lambda[4]);
    jetTree->SetBranchAddress("mcJetMomDisp", &mcJet.lambda[5]);
    jetTree->SetBranchAddress("mcJetD0Z", &mcJet.d0z);
    jetTree->SetBranchAddress("mcD0Pt", &mcJet.d0pt);
    jetTree->SetBranchAddress("mcD0Eta", &mcJet.d0eta);
    jetTree->SetBranchAddress("gRefMult", &mcJet.gRefMultMc);
    jetTree->SetBranchAddress("mcJetEta", &mcJet.mcJetEta);

    jetTree->SetBranchAddress("mcSmearedD0Pt", &recoJet.RecoD0Pt);
    jetTree->SetBranchAddress("mcJetNConst", &mcJet.numberofconstituents);

    



/*
    jetTree->SetBranchAddress(Form("RecoJetPt_t_%d_0", nnum), &recoJet.jetpt);
    jetTree->SetBranchAddress(Form("RecoNOfAllConst_t_%d_0", nnum), &recoJet.numberofconstituents);
    jetTree->SetBranchAddress(Form("RecoAlpha11_t_%d_0", nnum), &recoJet.lambda[0]);
    jetTree->SetBranchAddress(Form("RecoAlpha11half_t_%d_0", nnum), &recoJet.lambda[1]);
    jetTree->SetBranchAddress(Form("RecoAlpha12_t_%d_0", nnum), &recoJet.lambda[2]);
    jetTree->SetBranchAddress(Form("RecoAlpha13_t_%d_0", nnum), &recoJet.lambda[3]);
    jetTree->SetBranchAddress(Form("RecoAlpha10half_t_%d_0", nnum), &recoJet.lambda[4]);
    jetTree->SetBranchAddress(Form("RecoMomDisp_t_%d_0", nnum), &recoJet.lambda[5]);
        jetTree->SetBranchAddress(Form("D0z_t%d_0",nnum), &recoJet.d0z);

*/
/*    float jetpt_nocorr;
    float recoJetArea;
    float recoJetRho;*/
/*
    jetTree->SetBranchAddress("recoJetArea", &recoJet.recoJetArea);
    jetTree->SetBranchAddress("recoJetRho", &recoJet.recoJetRho);

    jetTree->SetBranchAddress("recoJetPt", &recoJet.jetpt_nocorr); //ZDE +Corr
    jetTree->SetBranchAddress("recoJetPtCorr", &recoJet.jetpt); //ZDE +Corr
    jetTree->SetBranchAddress("recoJetNConst", &recoJet.numberofconstituents);
    jetTree->SetBranchAddress("recoJetLambda1_1", &recoJet.lambda[0]);
    jetTree->SetBranchAddress("recoJetLambda1_1_5", &recoJet.lambda[1]);
    jetTree->SetBranchAddress("recoJetLambda1_2", &recoJet.lambda[2]);
    jetTree->SetBranchAddress("recoJetLambda1_3", &recoJet.lambda[3]);
    jetTree->SetBranchAddress("recoJetLambda1_0_5", &recoJet.lambda[4]);
    jetTree->SetBranchAddress("recoJetMomDisp", &recoJet.lambda[5]);
    jetTree->SetBranchAddress("recoJetD0Z", &recoJet.d0z);
    jetTree->SetBranchAddress("recoJetEta", &recoJet.recoJetEta);
    jetTree->SetBranchAddress("mcSmearedD0Eta", &recoJet.RecoD0Eta);
    jetTree->SetBranchAddress("mcSmearedJetEta", &recoJet.McRecoJetEta);
*/
 //TString branchPrefix = "ICS_"; // "ICS_" or ""
  TString branchPrefix = (Method == "ICS" ? "ICS_" : "");
    jetTree->SetBranchAddress(branchPrefix+"recoJetArea", &recoJet.recoJetArea);
    jetTree->SetBranchAddress(branchPrefix + "recoJetRho", &recoJet.recoJetRho);

   // jetTree->SetBranchAddress(branchPrefix + "recoJetPt", &recoJet.jetpt_nocorr); //ZDE +Corr
    jetTree->SetBranchAddress(branchPrefix + "recoJetPt" + (branchPrefix == "ICS_" ? "" : "Corr"), &recoJet.jetpt); //ZDE +Corr
    jetTree->SetBranchAddress(branchPrefix + "recoJetNConst", &recoJet.numberofconstituents);
    jetTree->SetBranchAddress(branchPrefix + "recoJetLambda1_1", &recoJet.lambda[0]);
    jetTree->SetBranchAddress(branchPrefix + "recoJetLambda1_1_5", &recoJet.lambda[1]);
    jetTree->SetBranchAddress(branchPrefix + "recoJetLambda1_2", &recoJet.lambda[2]);
    jetTree->SetBranchAddress(branchPrefix + "recoJetLambda1_3", &recoJet.lambda[3]);
    jetTree->SetBranchAddress(branchPrefix + "recoJetLambda1_0_5", &recoJet.lambda[4]);
    jetTree->SetBranchAddress(branchPrefix + "recoJetMomDisp", &recoJet.lambda[5]);
    jetTree->SetBranchAddress(branchPrefix + "recoJetD0Z", &recoJet.d0z);
    jetTree->SetBranchAddress(branchPrefix + "recoJetEta", &recoJet.recoJetEta);
    jetTree->SetBranchAddress("mcSmearedD0Eta", &recoJet.RecoD0Eta);
    jetTree->SetBranchAddress("mcSmearedJetEta", &recoJet.McRecoJetEta);



    // RecoPionEta
    ////jetTree->SetBranchAddress("mcSmearedPionEta", &recoJet.RecoPionEta);
    //RecoKaonEta
    ////jetTree->SetBranchAddress("mcSmearedKaonEta", &recoJet.RecoKaonEta);

}
/*
void assignTree2(TTree *jetTree, StJetTreeStruct2 &measured_, TString name)
{
    jetTree->SetBranchAddress("z", &measured_.z_value);
    jetTree->SetBranchAddress("pT", &measured_.pT_value);
    jetTree->SetBranchAddress("lambda_1_1", &measured_.lambda_value[0]);
    jetTree->SetBranchAddress("lambda_1_1half", &measured_.lambda_value[1]);
    jetTree->SetBranchAddress("lambda_1_2", &measured_.lambda_value[2]);
    jetTree->SetBranchAddress("lambda_1_3", &measured_.lambda_value[3]);
    jetTree->SetBranchAddress("lambda_1_0half", &measured_.lambda_value[4]);
    jetTree->SetBranchAddress("lambda_2_0", &measured_.lambda_value[5]);
    jetTree->SetBranchAddress("n_signal"+name+"_sw", &measured_.s_weight_value);
    jetTree->SetBranchAddress("centr_weight", &measured_.centr_weight_value);
    jetTree->SetBranchAddress("rev_weight_ef", &measured_.eff_weight_value);
    jetTree->SetBranchAddress("doubleCount", &measured_.doubleCount);
    jetTree->SetBranchAddress("central", &measured_.central);
    jetTree->SetBranchAddress("jet_area", &measured_.jet_area);
    jetTree->SetBranchAddress("bg_dens", &measured_.bg_dens);
    jetTree->SetBranchAddress("pT_raw", &measured_.pT_raw);
    jetTree->SetBranchAddress("mass", &measured_.mass);
    jetTree->SetBranchAddress("vD0pT", &measured_.vd0pt);
    jetTree->SetBranchAddress("D0_eta", &measured_.D0_eta);



}*/
void assignTree2(TTree *jetTree, StJetTreeStruct2 &measured_)
{
    jetTree->SetBranchAddress("z", &measured_.z_value);
    jetTree->SetBranchAddress("jetPt", &measured_.pT_value);
    jetTree->SetBranchAddress("lambda1_1", &measured_.lambda_value[0]);
    jetTree->SetBranchAddress("lambda1_1_5", &measured_.lambda_value[1]);
    jetTree->SetBranchAddress("lambda1_2", &measured_.lambda_value[2]);
    jetTree->SetBranchAddress("lambda1_3", &measured_.lambda_value[3]);
    jetTree->SetBranchAddress("lambda1_0_5", &measured_.lambda_value[4]);
    jetTree->SetBranchAddress("momDisp", &measured_.lambda_value[5]);
    jetTree->SetBranchAddress("sWeightSignal", &measured_.s_weight_value);
    //jetTree->SetBranchAddress("correctedWeight", &measured_.s_weight_value);
    jetTree->SetBranchAddress("sWeightSigLikelihood", &measured_.s_weight_error_value);
    jetTree->SetBranchAddress("weightCentrality", &measured_.centr_weight_value);
    jetTree->SetBranchAddress("weightD0Efficiency", &measured_.eff_weight_value);
    jetTree->SetBranchAddress("weightCentralityErr", &measured_.eff_weight_value_error);
    jetTree->SetBranchAddress("weightDoubleCount", &measured_.doubleCount);
    jetTree->SetBranchAddress("centrality", &measured_.central);
    jetTree->SetBranchAddress("centralityAlt", &measured_.centralAlt);
    jetTree->SetBranchAddress("gRefMult", &measured_.gRefMultVal);
    jetTree->SetBranchAddress("nJetConst", &measured_.nJetConst);

    //jetTree->SetBranchAddress("jet_area", &measured_.jet_area);
    //jetTree->SetBranchAddress("bg_dens", &measured_.bg_dens);
    //jetTree->SetBranchAddress("pT_raw", &measured_.pT_raw);
    jetTree->SetBranchAddress("d0Mass", &measured_.mass);
    jetTree->SetBranchAddress("d0Pt", &measured_.vd0pt);
    jetTree->SetBranchAddress("d0Eta", &measured_.D0_eta);
    //jetTree->SetBranchAddress("bg_dens", &measured_.bg_dens);
    jetTree->SetBranchAddress("jetEta", &measured_.jetEta);
    jetTree->SetBranchAddress("jetD0DeltaR", &measured_.jetD0DeltaR);



}
void assignTree2(TTree *jetTree, StJetTreeStruct2 &measured_, TString name)
{
    jetTree->SetBranchAddress("z", &measured_.z_value);
    jetTree->SetBranchAddress("pT", &measured_.pT_value);
    jetTree->SetBranchAddress("lambda_1_1", &measured_.lambda_value[0]);
    jetTree->SetBranchAddress("lambda_1_1_5", &measured_.lambda_value[1]);
    jetTree->SetBranchAddress("lambda_1_2", &measured_.lambda_value[2]);
    jetTree->SetBranchAddress("lambda_1_3", &measured_.lambda_value[3]);
    jetTree->SetBranchAddress("lambda_1_0_5", &measured_.lambda_value[4]);
    jetTree->SetBranchAddress("lambda_2_0", &measured_.lambda_value[5]);
    jetTree->SetBranchAddress("n_signal"+name+"_sw", &measured_.s_weight_value);
    jetTree->SetBranchAddress("centr_weight", &measured_.centr_weight_value);
    jetTree->SetBranchAddress("rev_weight_ef", &measured_.eff_weight_value);
    jetTree->SetBranchAddress("doubleCount", &measured_.doubleCount);
    jetTree->SetBranchAddress("central", &measured_.central);
    jetTree->SetBranchAddress("jet_area", &measured_.jet_area);
    jetTree->SetBranchAddress("bg_dens", &measured_.bg_dens);
    jetTree->SetBranchAddress("pT_raw", &measured_.pT_raw);
    jetTree->SetBranchAddress("mass", &measured_.mass);
    jetTree->SetBranchAddress("vD0pT", &measured_.vd0pt);
    jetTree->SetBranchAddress("jetEta", &measured_.jetEta);
    jetTree->SetBranchAddress("jetD0DeltaR", &measured_.jetD0DeltaR);



}

// H: THnSparseF s osami [0]=reco pT, [1]=reco z, [2]=true pT, [3]=true z
TH2D* MakeFlat2D(THnSparseF* H, const char* name="RMflat") {
    auto axRpt = H->GetAxis(0);
    auto axRz  = H->GetAxis(1);
    auto axTpt = H->GetAxis(2);
    auto axTz  = H->GetAxis(3);

    const int nRpt = axRpt->GetNbins();
    const int nRz  = axRz ->GetNbins();
    const int nTpt = axTpt->GetNbins();
    const int nTz  = axTz ->GetNbins();

    // plocha: X = (true pT × true z), Y = (reco pT × reco z)
    TH2D* h = new TH2D(name, ";(p_{T}^{true}, z^{true});(p_{T}^{reco}, z^{reco})",
                       nTpt*nTz, 0, nTpt*nTz,
                       nRpt*nRz, 0, nRpt*nRz);
    h->Sumw2();

    // projdeme pouze obsazené sparse bity
    Long64_t nBins = H->GetNbins();
    std::vector<Int_t> idx(4,0);
    for (Long64_t ib = 0; ib < nBins; ++ib) {
        double val = H->GetBinContent(ib, idx.data());
        if (val == 0.0) continue;

        // bin indexy jsou 1..N → převedeme na 0-based
        int iRpt = idx[0]-1;
        int iRz  = idx[1]-1;
        int iTpt = idx[2]-1;
        int iTz  = idx[3]-1;

        // lineární indexy (stejné „ploché“ pořadí, jaké používá RooUnfold)
        int iY = iRpt + iRz * nRpt;     // reco kombinovaný (0..nRpt*nRz-1)
        int iX = iTpt + iTz * nTpt;     // true kombinovaný

        h->AddBinContent(h->GetBin(iX+1, iY+1), val);
    }
    return h;
}

    // ---------- cache binning ----------
    const Int_t nPtRecoCache  = 98;
    const Double_t ptRecoMin  = 1.0;
    const Double_t ptRecoMax  =  50.0;

    const Int_t nPtTrueCache  = 58;
    const Double_t ptTrueMin  =   1.0;
    const Double_t ptTrueMax  =  30.0;

    const Int_t nZRecoCache   = 101;
    const Double_t zRecoMin   = 0;
    const Double_t zRecoMax   =  1.01;

    const Int_t nZTrueCache   = 101;
    const Double_t zTrueMin   = 0;
    const Double_t zTrueMax   =  1.01;

    const Int_t nAngRecoCache = 1200;
    const Double_t angRecoMin = 0;
    const Double_t angRecoMax =  3;

    const Int_t nAngTrueCache = 1200;
    const Double_t angTrueMin = 0;
    const Double_t angTrueMax =  3;

void HistogramInit(){



    for (int ic = 0; ic < nCentralityBins; ic++) {

        // ---------- 1D pT ----------
        hCacheMatchPt[ic] = new TH2D(
            Form("hCacheMatchPt_cent%d", ic),
            Form("Matched pT cent %d; p_{T,jet}^{reco}; p_{T,jet}^{true}", ic),
            nPtRecoCache, ptRecoMin, ptRecoMax,
            nPtTrueCache, ptTrueMin, ptTrueMax
        );
        hCacheMatchPt[ic]->Sumw2();
        hCacheMatchPt[ic]->SetDirectory(0);

        hCacheMissPt[ic] = new TH1D(
            Form("hCacheMissPt_cent%d", ic),
            Form("Miss pT cent %d; p_{T,jet}^{true}; counts", ic),
            nPtTrueCache, ptTrueMin, ptTrueMax
        );
        hCacheMissPt[ic]->Sumw2();
        hCacheMissPt[ic]->SetDirectory(0);

        hCacheFakePt[ic] = new TH1D(
            Form("hCacheFakePt_cent%d", ic),
            Form("Fake pT cent %d; p_{T,jet}^{reco}; counts", ic),
            nPtRecoCache, ptRecoMin, ptRecoMax
        );
        hCacheFakePt[ic]->Sumw2();
        hCacheFakePt[ic]->SetDirectory(0);

        // ---------- 2D pT-z ----------
        {
            Int_t nbins[4] = {
                nPtRecoCache,
                nZRecoCache,
                nPtTrueCache,
                nZTrueCache
            };

            Double_t xmin[4] = {
                ptRecoMin,
                zRecoMin,
                ptTrueMin,
                zTrueMin
            };

            Double_t xmax[4] = {
                ptRecoMax,
                zRecoMax,
                ptTrueMax,
                zTrueMax
            };

            hCacheMatchPtZ[ic] = new THnSparseD(
                Form("hCacheMatchPtZ_cent%d", ic),
                Form("Matched pT-z cent %d; p_{T}^{reco}; z^{reco}; p_{T}^{true}; z^{true}", ic),
                4, nbins, xmin, xmax
            );
            hCacheMatchPtZ[ic]->Sumw2();
        }

        hCacheMissPtZ[ic] = new TH2D(
            Form("hCacheMissPtZ_cent%d", ic),
            Form("Miss pT-z cent %d; p_{T}^{true}; z^{true}", ic),
            nPtTrueCache, ptTrueMin, ptTrueMax,
            nZTrueCache, zTrueMin, zTrueMax
        );
        hCacheMissPtZ[ic]->Sumw2();
        hCacheMissPtZ[ic]->SetDirectory(0);

        hCacheFakePtZ[ic] = new TH2D(
            Form("hCacheFakePtZ_cent%d", ic),
            Form("Fake pT-z cent %d; p_{T}^{reco}; z^{reco}", ic),
            nPtRecoCache, ptRecoMin, ptRecoMax,
            nZRecoCache, zRecoMin, zRecoMax
        );
        hCacheFakePtZ[ic]->Sumw2();
        hCacheFakePtZ[ic]->SetDirectory(0);

        // ---------- 2D pT-lambda ----------
        for (int iAng = 0; iAng < nAngularities; iAng++) {

            Int_t nbinsAng[4] = {
                nPtRecoCache,
                nAngRecoCache,
                nPtTrueCache,
                nAngTrueCache
            };

            Double_t xminAng[4] = {
                ptRecoMin,
                angRecoMin,
                ptTrueMin,
                angTrueMin
            };

            Double_t xmaxAng[4] = {
                ptRecoMax,
                angRecoMax,
                ptTrueMax,
                angTrueMax
            };

            hCacheMatchPtAng[ic][iAng] = new THnSparseD(
                Form("hCacheMatchPtAng_cent%d_ang%d", ic, iAng),
                Form("Matched pT-lambda cent %d ang %d; p_{T}^{reco}; #lambda^{reco}; p_{T}^{true}; #lambda^{true}", ic, iAng),
                4, nbinsAng, xminAng, xmaxAng
            );
            hCacheMatchPtAng[ic][iAng]->Sumw2();

            hCacheMissPtAng[ic][iAng] = new TH2D(
                Form("hCacheMissPtAng_cent%d_ang%d", ic, iAng),
                Form("Miss pT-lambda cent %d ang %d; p_{T}^{true}; #lambda^{true}", ic, iAng),
                nPtTrueCache, ptTrueMin, ptTrueMax,
                nAngTrueCache, angTrueMin, angTrueMax
            );
            hCacheMissPtAng[ic][iAng]->Sumw2();
            hCacheMissPtAng[ic][iAng]->SetDirectory(0);

            hCacheFakePtAng[ic][iAng] = new TH2D(
                Form("hCacheFakePtAng_cent%d_ang%d", ic, iAng),
                Form("Fake pT-lambda cent %d ang %d; p_{T}^{reco}; #lambda^{reco}", ic, iAng),
                nPtRecoCache, ptRecoMin, ptRecoMax,
                nAngRecoCache, angRecoMin, angRecoMax
            );
            hCacheFakePtAng[ic][iAng]->Sumw2();
            hCacheFakePtAng[ic][iAng]->SetDirectory(0);
        }
    }

    for (int i = 0; i < 3; i++) {
        KinEff1D[i] = TEfficiency(Form("_ahRealData_%d", i),Form("_hRealData_%d;p_{T,Jet}^{true} [GeV/c];Kin. efficiency", i), ptMcBinsVecCustom[i].size()-1, &ptMcBinsVecCustom[i][0]);
        KinEff2DpTZ[0][i] = TEfficiency(Form("_haRealDatapTz1_%d", i),Form("_hRealData_%d;p_{T,Jet}^{true} [GeV/c];Kin. efficiency (p_{T,Jet},z)", i), ptMcBinsVecCustom[i].size()-1, &ptMcBinsVecCustom[i][0]);
        KinEff2DpTZ[1][i] = TEfficiency(Form("_hdRealDatapTz2_%d", i),Form("_hRealData_%d;z^{true};Kin. efficiency (p_{T,Jet},z)", i), zMcBinsVecCustom[i].size()-1, &zMcBinsVecCustom[i][0]);
        KinEff2DpTZ[2][i] = TEfficiency(Form("_hfRealDatapTz3_%d", i),Form("_hRealData_%d;p_{T,Jet}^{true} [GeV/c];z^{true};Kin. efficiency (p_{T,Jet},z)", i), ptMcBinsVecCustom[i].size()-1, &ptMcBinsVecCustom[i][0], zMcBinsVecCustom[i].size()-1, &zMcBinsVecCustom[i][0]);
        KinEff2DpTZCut[i] = TEfficiency(Form("_hdRsealDatapTz2_%d", i),Form("_hRealData_%d;z^{true};Kin. efficiency (p_{T,Jet},z)", i), zMcBinsVecCustom[i].size()-1, &zMcBinsVecCustom[i][0]);
        KinEff2DpTZZCut[i] = TEfficiency(Form("_hdRsealZDatapTz2_%d", i),Form("_hRealData_%d;p_{T,Jet}^{true} [GeV/c];Kin. efficiency (p_{T,Jet},z)", i), ptMcBinsVecCustom[i].size()-1, &ptMcBinsVecCustom[i][0]);
        KinEffEta[i] = TEfficiency(Form("_hadRealDataEta_%d", i), Form("_hRealData_%d;p_{T,Jet}^{true} [GeV/c];Missing jets in #eta", i), ptMcBinsVecCustom[i].size()-1, &ptMcBinsVecCustom[i][0]);
        FakeEffEta[i] = TEfficiency(Form("_hadFakeDataEta_%d", i), Form("_hRealData_%d;p_{T,Jet}^{true} [GeV/c];Fake jets in #eta", i), ptMcBinsVecCustom[i].size()-1, &ptMcBinsVecCustom[i][0]);
        for (int iAng = 0; iAng < nAngularities; iAng++) {
            KinEff2DAng[iAng][0][i] = TEfficiency(Form("_hadRealDataAng%d_%d_%d_",0, i, iAng), Form("_hRealData_%d;p_{T,Jet}^{true} [GeV/c];Kin. efficiency (p_{T,Jet}, #lambda)", i),  ptMcBinsVecCustom[i].size()-1, &ptMcBinsVecCustom[i][0]);
            KinEff2DAng[iAng][1][i] = TEfficiency(Form("_hddRealDataAng%d_%d_%d_",1, i, iAng), Form("_hRealData_%d;%s ^{true};Kin. efficiency (p_{T,Jet}, %s)", i, AngNames[iAng].Data(),AngNames[iAng].Data()),angMcBinsVecCustom[i][iAng].size()-1, &angMcBinsVecCustom[i][iAng][0]);
            KinEff2DAng[iAng][2][i] = TEfficiency(Form("_hfdRealDataAng%d_%d_%d_",2, i, iAng), Form("_hRealData_%d;p_{T,Jet}^{true} [GeV/c];%s ^{true};Kin. efficiency (p_{T,Jet}, #lambda)", i, AngNames[iAng].Data()),  ptMcBinsVecCustom[i].size()-1, &ptMcBinsVecCustom[i][0], angMcBinsVecCustom[i][iAng].size()-1, &angMcBinsVecCustom[i][iAng][0]);
            KinEff2DAngCut[iAng][i] = TEfficiency(Form("_hadRealDataAngCut%d_%d_", iAng, i), Form("_hRealData_%d;#lambda_{%d}^{true};Kin. efficiency (p_{T,Jet}, #lambda)", i, iAng+1), angMcBinsVecCustom[i][iAng].size()-1, &angMcBinsVecCustom[i][iAng][0]);
            KinEff2DAngPtCut[iAng][i] = TEfficiency(Form("_hadRealDataAngCutpT%d_%d_", iAng, i), Form("_hRealData_%d;p_{T,Jet}^{true} [GeV/c];Kin. efficiency (p_{T,Jet}, %s)", i, AngNames[iAng].Data()), ptMcBinsVecCustom[i].size()-1, &ptMcBinsVecCustom[i][0]);

        }

        for (int j = 0; j < 7; j++) {
            hRealFine[i][j] = new TH2D(Form("hRealFine_%d_%d", i, j), Form("hRealFine_%d_%d", i, j), 100, 0, 20, 150, 0, 1.5);
            hResponseFine4D[i][j] = new THnSparseD(Form("hResponseFine4D_%d_%d", i, j), Form("hResponseFine4D_%d_%d", i, j), 4, nbins_, xmin_, xmax_);
        }
    }

    for (int iCent = 0; iCent < 3; iCent++) {
        hRealData[iCent] = TH1D(Form("hRealData_%d", iCent), Form("hRealData_%d", iCent), ptRecoBinsVec[iCent].size()-1, &ptRecoBinsVec[iCent][0]);
        hRealData[iCent].GetXaxis()->SetTitle("p_{T,Jet} [GeV/c]");
        hRealData[iCent].GetYaxis()->SetTitle("Counts");
    }

    //hRealData2D[centfrom0to2][7 variables]
    for (int iCent = 0; iCent < 3; iCent++) {
            hRealData2D[iCent][0] = TH2D(Form("hRealData2D_%d_%d", iCent, 0), Form("hRealData2D_%d_%d", iCent, 0), ptRecoBinsVec[iCent].size()-1, &ptRecoBinsVec[iCent][0], zRecoBinsVec[iCent].size()-1, &zRecoBinsVec[iCent][0]);
            hRealData2D[iCent][0].GetXaxis()->SetTitle("p_{T,Jet} [GeV/c]");
            hRealData2D[iCent][0].GetYaxis()->SetTitle("z");
        for (int ptD0 = 0; ptD0 < 5; ptD0++) {
            hRealData2DD0Pt[iCent][ptD0][0] = TH2D(Form("0hRealData2D_PtD0_%d_%d", iCent, ptD0), Form("hRealData2D_PtD0_%d_%d", iCent, ptD0), ptRecoBinsVec[iCent].size()-1, &ptRecoBinsVec[iCent][0], zRecoBinsVec[iCent].size()-1, &zRecoBinsVec[iCent][0]);
            hRealData2DD0Pt[iCent][ptD0][0].GetXaxis()->SetTitle("p_{T,Jet} (GeV/c)");
            hRealData2DD0Pt[iCent][ptD0][0].GetYaxis()->SetTitle("z");
            hRealData1D0Pt[iCent][ptD0] = TH1D(Form("hRealData1D0Pt_%d_%d", iCent, ptD0), Form("hRealData1D0Pt_%d_%d", iCent, ptD0), ptRecoBinsVec[iCent].size()-1, &ptRecoBinsVec[iCent][0]);
            hRealData1D0Pt[iCent][ptD0].GetXaxis()->SetTitle("p_{T,Jet} (GeV/c)");
        }
        for (Int_t iLambda = 0; iLambda < nAngularities; iLambda++) {
            hRealData2D[iCent][iLambda+1] = TH2D(Form("hRealData2D_%d_%d", iCent, iLambda+1), Form("hRealData2D_%d_%d", iCent, iLambda+1), ptRecoBinsVec[iCent].size()-1, &ptRecoBinsVec[iCent][0], angRecoBinsVec[iCent][iLambda].size()-1, &angRecoBinsVec[iCent][iLambda][0]);
            for (int ptD0 = 0; ptD0 < 5; ptD0++) {
                hRealData2DD0Pt[iCent][ptD0][iLambda+1] = TH2D(Form("hRealData2D_PtD0_%d_%d_%d", iCent, ptD0,  iLambda+1), Form("hRealData2D_PtD0_%d_%d_%d", iCent, ptD0,  iLambda+1), ptRecoBinsVec[iCent].size()-1, &ptRecoBinsVec[iCent][0], angRecoBinsVec[iCent][iLambda].size()-1, &angRecoBinsVec[iCent][iLambda][0]);
                hRealData2DD0Pt[iCent][ptD0][iLambda+1].GetXaxis()->SetTitle("p_{T,Jet} (GeV/c)");
                hRealData2DD0Pt[iCent][ptD0][iLambda+1].GetYaxis()->SetTitle(Form("#lambda_%d", iLambda));
            }


            }

    }

    

    for (int iCent = 0; iCent < 3; ++iCent) {
          JetPtZMc[iCent] = new TH2D("PTZTrueTemp_" + centralityNames[iCent],
                                   ";p_{T, Jet} [GeV/c]; z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}",
                                   ptMcBinsVecCustom[iCent].size() - 1, &ptMcBinsVecCustom[iCent][0],
                                   zMcBinsVecCustom[iCent].size() - 1, &zMcBinsVecCustom[iCent][0]);

        HistReal_pTZTemp[iCent] = new TH2D("PTZMeasTemp" + centralityNames[iCent],
                                           ";p_{T, Jet} [GeV/c]; z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}",
                                           ptRecoBinsVec[iCent].size() - 1, &ptRecoBinsVec[iCent][0],
                                           zRecoBinsVec[iCent].size() - 1, &zRecoBinsVec[iCent][0]);
        hTruthPTZTemp[iCent] = new TH2D("PTZTrueTemp" + centralityNames[iCent],
                                        ";p_{T, Jet} [GeV/c]; z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}",
                                        ptMcBinsVecCustom[iCent].size() - 1, &ptMcBinsVecCustom[iCent][0],
                                        zMcBinsVecCustom[iCent].size() - 1, &zMcBinsVecCustom[iCent][0]);
        
        HistReal_ZLam[iCent] = new TH2D("LamZMeasTemp" + centralityNames[iCent],
                                           ";z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}; #lambda_{1}^{1, reco}",
                                           zRecoBinsVec[iCent].size() - 1, &zRecoBinsVec[iCent][0],
                                           angMcBinsVecCustom[iCent][0].size() - 1, &angMcBinsVecCustom[iCent][0][0]);
        hTruthZLam[iCent] = new TH2D("LamZTrueTemp" + centralityNames[iCent],
                                        ";z = #vec{p}_{T, jet}#dot#vec{p}_{T,D^{0}} /|#vec{p}_{T, jet}|^{2}; #lambda_{1}^{1, true}",
                                        zMcBinsVecCustom[iCent].size() - 1, &zMcBinsVecCustom[iCent][0],
                                        angMcBinsVecCustom[iCent][0].size() - 1, &angMcBinsVecCustom[iCent][0][0]);

        hRespZ[1][iCent] = new TH2D("respZ" + centralityNames[iCent],
                                    ";z^{reco}; z^{truth}",
                                    zRecoBinsVec[iCent].size() - 1, &zRecoBinsVec[iCent][0],
                                    zMcBinsVecCustom[iCent].size() - 1, &zMcBinsVecCustom[iCent][0]);
        hRespZ[0][iCent] = new TH2D("respPT" + centralityNames[iCent],
                                    ";p_{T}^{reco}; p_{T}^{truth}",
                                    ptRecoBinsVec[iCent].size() - 1, &ptRecoBinsVec[iCent][0],
                                    ptMcBinsVecCustom[iCent].size() - 1, &ptMcBinsVecCustom[iCent][0]);
        hRespZHighRes[1][iCent] = new TH2D("respZ2" + centralityNames[iCent],
                                           ";z^{reco}; z^{truth}",
                                           100, zRecoBinsVec[iCent][0],
                                           zRecoBinsVec[iCent][zRecoBinsVec[iCent].size() - 1],
                                           100, zMcBinsVecCustom[iCent][0],
                                           zMcBinsVecCustom[iCent][zMcBinsVecCustom[iCent].size() - 1]);
        for (int i = 0; i < 5; i++){                                                       
        hResVar[iCent][0][i] = new TH1D("respVarPT" + centralityNames[iCent] + Form("_PtJet_%d", i),
                                    ";(p_{T}^{truth} - p_{T}^{reco})/p_{T}^{truth} [GeV/c]",
                                    100, -5, 5);
        hResVar[iCent][0][i]->Sumw2();
        hResVar[iCent][1][i] = new TH1D("respVarD0Z" + centralityNames[iCent] + Form("_PtJet_%d", i),
                                    ";(z^{truth} - z^{reco})/z^{truth}",
                                    100, -5, 5);
        hResVar[iCent][1][i]->Sumw2();
        }
        hRespZHighRes[0][iCent] = new TH2D("respPT2" + centralityNames[iCent],
                                           ";p_{T,Jet}^{reco} [GeV/c]; p_{T,Jet}^{truth} [GeV/c]",
                                           100, ptRecoBinsVec[iCent][0],
                                           ptRecoBinsVec[iCent][ptRecoBinsVec[iCent].size() - 1],
                                           100, ptMcBinsVecCustom[iCent][0],
                                           ptMcBinsVecCustom[iCent][ptMcBinsVecCustom[iCent].size() - 1]);
        hMeasuredPtRealTest[iCent] = new TH1D(Form("MeasPtRealTest%s", centralityNames[iCent].Data()), ";p_{T, Jet} [GeV/c]",
                                              ptRecoBinsVec[iCent].size() - 1, &ptRecoBinsVec[iCent][0]);


        rurResponse2D[iCent][0].Setup(HistReal_pTZTemp[iCent], hTruthPTZTemp[iCent]);
        rurResponse2DTest[iCent].Setup(HistReal_pTZTemp[iCent], hTruthPTZTemp[iCent]);
        rurResponse2DTestW[iCent].Setup(HistReal_pTZTemp[iCent], hTruthPTZTemp[iCent]);

        //angularities
        for (Int_t iLambda = 0; iLambda < 6; iLambda++) {

            for (int ptJet = 0; ptJet < 5; ptJet++) {
                hResVar[iCent][2+iLambda][ptJet] = new TH1D(Form("resdsfpVar%d_%d_PtJet_%d", iLambda, iCent, ptJet),
                                            ";(" + AngNames[iLambda] + "^{truth} - " + AngNames[iLambda] + "^{reco}) / " + AngNames[iLambda] + "^{truth}",
                                            100, -5, 5);
                hResVar[iCent][2+iLambda][ptJet]->Sumw2();
            }
            hRespZ[iLambda + 2][iCent] = new TH2D(Form("respZ%s_%i", centralityNames[iCent].Data(), iLambda),
                                                  ";" + AngNames[iLambda] + "(reco); " + AngNames[iLambda] + "(truth)",
                                                  angRecoBinsVec[iCent][iLambda].size() - 1,
                                                  &angRecoBinsVec[iCent][iLambda][0],
                                                  angMcBinsVecCustom[iCent][iLambda].size() - 1,
                                                  &angMcBinsVecCustom[iCent][iLambda][0]);
            hRespZHighRes[iLambda + 2][iCent] = new TH2D(
                    Form("respZHighRes%s_%i", centralityNames[iCent].Data(), iLambda),
                    ";" + AngNames[iLambda] + " ^{,reco}; " + AngNames[iLambda] + " ^{,true}",
                    100,
                    angRecoBinsVec[iCent][iLambda][0],
                    angRecoBinsVec[iCent][iLambda][angRecoBinsVec[iCent][iLambda].size() - 1],
                    100,
                    angMcBinsVecCustom[iCent][iLambda][0],
                    angMcBinsVecCustom[iCent][iLambda][angMcBinsVecCustom[iCent][iLambda].size() - 1]);


            HistReal_pTAngTemp[iCent][iLambda] = new TH2D(Form("MeasuredAngTemp%s_%i", centralityNames[iCent].Data(), iLambda),
                                                          ";p_{T, Jet} [GeV/c];#lambda_{%i}", ptRecoBinsVec[iCent].size() - 1,
                                                          &ptRecoBinsVec[iCent][0],
                                                          angRecoBinsVec[iCent][iLambda].size() - 1,
                                                          &angRecoBinsVec[iCent][iLambda][0]);
            hTruthAngTemp[iCent][iLambda] = new TH2D(Form("TrueAngTemp%s_%i", centralityNames[iCent].Data(), iLambda),
                                                     ";p_{T, Jet} [GeV/c];#lambda_{%i}", ptMcBinsVecCustom[iCent].size() - 1,
                                                     &ptMcBinsVecCustom[iCent][0],
                                                     angMcBinsVecCustom[iCent][iLambda].size() - 1,
                                                     &angMcBinsVecCustom[iCent][iLambda][0]);

            rurResponse2D[iCent][1+iLambda].Setup(HistReal_pTAngTemp[iCent][iLambda], hTruthAngTemp[iCent][iLambda]);

        }

        //rurResponse2D[iCent][7].Setup(HistReal_pTAngTemp[iCent][0], hTruthAngTemp[iCent][0]);


        HistReal_pTTemp[iCent] = new TH1D(Form("MeasPtTemp%s", centralityNames[iCent].Data()),
                                          ";p_{T, Jet} [GeV/c]",
                                          ptRecoBinsVec[iCent].size() - 1, &ptRecoBinsVec[iCent][0]);
        hTruthPtTemp[iCent] = new TH1D(Form("TruePtTemp%s", centralityNames[iCent].Data()),
                                       ";p_{T, Jet} [GeV/c]",
                                       ptMcBinsVecCustom[iCent].size() - 1, &ptMcBinsVecCustom[iCent][0]);
        rurResponse[iCent].Setup(HistReal_pTTemp[iCent], hTruthPtTemp[iCent]);
        if (UseOverflow) rurResponse[iCent].UseOverflow();
    }
       

    jetPtCheck[0] = new TH1D("jetPtCheck0", "jetPtCheck0", ptMcBinsVecCustom[0].size()-1, &ptMcBinsVecCustom[0][0]);
    jetPtCheck[1] = new TH1D("jetPtCheck1", "jetPtCheck1", ptMcBinsVecCustom[1].size()-1, &ptMcBinsVecCustom[1][0]);
    jetPtCheck[2] = new TH1D("jetPtCheck2", "jetPtCheck2", ptMcBinsVecCustom[2].size()-1, &ptMcBinsVecCustom[2][0]);

    jetZ[0] = new TH1D("jetZ0", "jetZ0", zMcBinsVecCustom[0].size()-1, &zMcBinsVecCustom[0][0]);
    jetZ[1] = new TH1D("jetZ1", "jetZ1", zMcBinsVecCustom[1].size()-1, &zMcBinsVecCustom[1][0]);
    jetZ[2] = new TH1D("jetZ2", "jetZ2", zMcBinsVecCustom[2].size()-1, &zMcBinsVecCustom[2][0]);

    McPTRawD0[0] = new TH1D("McPTRawD0_0", "McPTRawD0_0", 10, 0, 10);
    McPTRawD0[1] = new TH1D("McPTRawD0_1", "McPTRawD0_1", 10, 0, 10);
    McPTRawD0[2] = new TH1D("McPTRawD0_2", "McPTRawD0_2", 10, 0, 10);

    McPTRawD0Jet[0] = new TH1D("McPTRawD0Jet_0", "McPTRawD0Jet_0", 60, 1, 31);
    McPTRawD0Jet[1] = new TH1D("McPTRawD0Jet_1", "McPTRawD0Jet_1", 60, 1, 31);
    McPTRawD0Jet[2] = new TH1D("McPTRawD0Jet_2", "McPTRawD0Jet_2", 60, 1, 31);

    //Proklety radek
/*
    McPTRawD0Jet[0] = new TH1D("McPTRawD0Jet_0", "McPTRawD0Jet_0", ptMcBinsVecCustom[0].size()-1, &ptMcBinsVecCustom[0][0]);
    McPTRawD0Jet[1] = new TH1D("McPTRawD0Jet_1", "McPTRawD0Jet_1", ptMcBinsVecCustom[0].size()-1, &ptMcBinsVecCustom[0][0]);
    McPTRawD0Jet[2] = new TH1D("McPTRawD0Jet_2", "McPTRawD0Jet_2", ptMcBinsVecCustom[0].size()-1, &ptMcBinsVecCustom[0][0]);
*/
    McPTRawD0JetD0Meson[0] = new TH2D("McPTRawD0JetD0Meson_1)","McPTRawD0JetD0Meson_1);p_{T, Jet} [GeV/c];p_{T}^{D}",60, 1, 30, 40, 0, 10);
    McPTRawD0JetD0Meson[1] = new TH2D("McPTRawD0JetD0Meson_2)","McPTRawD0JetD0Meson_1);p_{T, Jet} [GeV/c];p_{T}^{D}",60, 1, 30, 40, 0, 10);
    McPTRawD0JetD0Meson[2] = new TH2D("McPTRawD0JetD0Meson_3)","McPTRawD0JetD0Meson_1);p_{T, Jet} [GeV/c];p_{T}^{D}",60, 1, 30, 40, 0, 10);


    hMeasuredD0MesonPt[0] = new TH1D("hMeasuredD0Meson_0", "hMeasuredD0Meson_0",  BinyVl.size() - 1, &BinyVl[0]);
    hMeasuredD0MesonPt[1] = new TH1D("hMeasuredD0Meson_1", "hMeasuredD0Meson_1",  BinyVl.size() - 1, &BinyVl[0]);
    hMeasuredD0MesonPt[2] = new TH1D("hMeasuredD0Meson_2", "hMeasuredD0Meson_2",  BinyVl.size() - 1, &BinyVl[0]);
    hMeasuredD0MesonPtRatio[0] = new TH1D("hMeasuredD0MesonRatio_0", "hMeasuredD0MesonRatio_0", BinyVl.size() - 1, &BinyVl[0]);
    hMeasuredD0MesonPtRatio[1] = new TH1D("hMeasuredD0MesonRatio_1", "hMeasuredD0MesonRatio_1", BinyVl.size() - 1, &BinyVl[0]);
    hMeasuredD0MesonPtRatio[2] = new TH1D("hMeasuredD0MesonRatio_2", "hMeasuredD0MesonRatio_2", BinyVl.size() - 1, &BinyVl[0]);

    jetPtCheckScaled[0] = new TH1D("jetPtCheckScaled0", "jetPtCheckScaled0", ptMcBinsVecCustom[0].size()-1, &ptMcBinsVecCustom[0][0]);
    jetPtCheckScaled[1] = new TH1D("jetPtCheckScaled1", "jetPtCheckScaled1", ptMcBinsVecCustom[1].size()-1, &ptMcBinsVecCustom[1][0]);
    jetPtCheckScaled[2] = new TH1D("jetPtCheckScaled2", "jetPtCheckScaled2", ptMcBinsVecCustom[2].size()-1, &ptMcBinsVecCustom[2][0]);

    jetPtCheckScaled2[0] = new TH1D("jetPtCheckScaled02", "jetPtCheckScaled0", ptMcBinsVecCustom[0].size()-1, &ptMcBinsVecCustom[0][0]);
    jetPtCheckScaled2[1] = new TH1D("jetPtCheckScaled12", "jetPtCheckScaled1", ptMcBinsVecCustom[1].size()-1, &ptMcBinsVecCustom[1][0]);
    jetPtCheckScaled2[2] = new TH1D("jetPtCheckScaled22", "jetPtCheckScaled2", ptMcBinsVecCustom[2].size()-1, &ptMcBinsVecCustom[2][0]);

    jetPtRecoCheckScaled[0] = new TH1D("jetPtRecoCheckScaled0", "jetPtRecoCheckScaled0", ptRecoBinsVec[0].size()-1, &ptRecoBinsVec[0][0]);
    jetPtRecoCheckScaled[1] = new TH1D("jetPtRecoCheckScaled1", "jetPtRecoCheckScaled1", ptRecoBinsVec[1].size()-1, &ptRecoBinsVec[1][0]);
    jetPtRecoCheckScaled[2] = new TH1D("jetPtRecoCheckScaled2", "jetPtRecoCheckScaled2", ptRecoBinsVec[2].size()-1, &ptRecoBinsVec[2][0]);

    jetPtRecoCheckScaled2[0] = new TH1D("jetPtRecoCheckScaled02", "jetPtRecoCheckScaled0", ptRecoBinsVec[0].size()-1, &ptRecoBinsVec[0][0]);
    jetPtRecoCheckScaled2[1] = new TH1D("jetPtRecoCheckScaled12", "jetPtRecoCheckScaled1", ptRecoBinsVec[1].size()-1, &ptRecoBinsVec[1][0]);
    jetPtRecoCheckScaled2[2] = new TH1D("jetPtRecoCheckScaled22", "jetPtRecoCheckScaled2", ptRecoBinsVec[2].size()-1, &ptRecoBinsVec[2][0]);

    jetPtRecoCheck[0] = new TH1D("jetPtRecoCheck0", "jetPtRecoCheck0", ptRecoBinsVec[0].size()-1, &ptRecoBinsVec[0][0]);
    jetPtRecoCheck[1] = new TH1D("jetPtRecoCheck1", "jetPtRecoCheck1", ptRecoBinsVec[1].size()-1, &ptRecoBinsVec[1][0]);
    jetPtRecoCheck[2] = new TH1D("jetPtRecoCheck2", "jetPtRecoCheck2", ptRecoBinsVec[2].size()-1, &ptRecoBinsVec[2][0]);




    D0MesonPtReal[0] = new TH1D("D0MesonPtReal0", "D0MesonPtReal0", 10, 0, 10);
    D0MesonPtReal[1] = new TH1D("D0MesonPtReal1", "D0MesonPtReal1", 10, 0, 10);
    D0MesonPtReal[2] = new TH1D("D0MesonPtReal2", "D0MesonPtReal2", 10, 0, 10);

    D0MesonPtMcReco[0] = new TH1D("D0MesonPtMcReco0", "D0MesonPtMcReco0", 10, 0, 10);
    D0MesonPtMcReco[1] = new TH1D("D0MesonPtMcReco1", "D0MesonPtMcReco1", 10, 0, 10);
    D0MesonPtMcReco[2] = new TH1D("D0MesonPtMcReco2", "D0MesonPtMcReco2", 10, 0, 10);

    D0MesonPtMcTrue[0] = new TH1D("D0MesonPtMcTrue0", "D0MesonPtMcTrue0", 10, 0, 10);
    D0MesonPtMcTrue[1] = new TH1D("D0MesonPtMcTrue1", "D0MesonPtMcTrue1", 10, 0, 10);
    D0MesonPtMcTrue[2] = new TH1D("D0MesonPtMcTrue2", "D0MesonPtMcTrue2", 10, 0, 10);



    D0JetPtMcReco[0] = new TH1D("D0JetPtMcReco0", "D0JetPtMcReco0;p_{T,Jet}^{reco} [GeV/c]", 58, 1, 30);
    D0JetPtMcReco[1] = new TH1D("D0JetPtMcReco1", "D0JetPtMcReco1;p_{T,Jet}^{reco} [GeV/c]", 58, 1, 30);
    D0JetPtMcReco[2] = new TH1D("D0JetPtMcReco2", "D0JetPtMcReco2;p_{T,Jet}^{reco} [GeV/c]", 58, 1, 30);

    D0JetPtMcTrue[0] = new TH1D("D0JetPtMcTrue0", "D0JetPtMcTrue0;p_{T,Jet}^{true} [GeV/c]", 58, 1, 30);
    D0JetPtMcTrue[1] = new TH1D("D0JetPtMcTrue1", "D0JetPtMcTrue1;p_{T,Jet}^{true} [GeV/c]", 58, 1, 30);
    D0JetPtMcTrue[2] = new TH1D("D0JetPtMcTrue2", "D0JetPtMcTrue2;p_{T,Jet}^{true} [GeV/c]", 58, 1, 30);

}
