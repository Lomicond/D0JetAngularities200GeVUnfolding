// Place in Unfolding/, beside Machine.C. Run from the repository root:
// root -l -b -q 'Unfolding/PlotPriorShapes.C(5,10)'
// Uses ../config.h for MC path, Method, truth bins and constituent selection.
// No unfolding, no response/cache changes, no data-derived prior tilts.
// Four weights on IDENTICAL selected truth events:
//   raw=1; centrality=c; before=c*F(pt)*1e10/H_selected(pt);
//   after=c*F(pt)*1e10/H_broad(pt).
// c=weightCentrality*Ndata(centrality)/Nmc_selected(centrality).
// Denominators are 60 bins [1,31], filled with weightCentrality (Machine2).
// All MC entries / DividedMcDataBy are used; no closure train/test split.
// Output ROOT histograms are unnormalised SUMS OF WEIGHTS, with Sumw2.
// PDF: each curve separately normalised in its displayed range, bin widths
// accounted for; ratios have no error bars (same correlated MC sample).
// Input defaults: config McJetsFileData, Output/Output_<Method>2.root,
// Output/FONLL_D0Pt_0_31.root (charm curve, evaluated at truth JET pt).
// Optional arguments: (d0min,d0max,outputStem,mcFile,dataFile,fonllFile,
//                      include5to7,broadMin,broadMax).
// Changes made ONLY inside Machine2() (e.g. a bin override) are not inherited:
// put the same truth bins in config.h if you need an identical binning.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>
#include <string>
#include <limits>
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TF1.h"
#include "TNamed.h"
using namespace std;
#include "../config.h"

namespace PriorShapesQA {
const int NV=9, NW=4, NC=3;
const char *keys[NV]={"jetPt","z","lambda1","lambda1p5","lambda2","lambda3","lambda0p5","pTD","D0Pt"};
const char *titles[NV]={"p_{T,jet}^{true} [GeV/c]","z^{true}","#lambda_{1}^{1}","#lambda_{1.5}^{1}","#lambda_{2}^{1}","#lambda_{3}^{1}","#lambda_{0.5}^{1}","p_{T}^{D}","p_{T,D^{0}}^{true} [GeV/c]"};
const char *labels[NW]={"Default MC (w = 1)","Default MC #times cweight","Cut before reshape","Cut after reshape"};
const char *cent[NC]={"0-10%","10-40%","40-80%"};
const int colors[NW]={kBlack,kAzure+2,kOrange+7,kMagenta+1};
const int styles[NW]={2,1,1,7};

void require(bool ok,const string &msg) { if(!ok) throw runtime_error(msg); }
unique_ptr<TFile> open(const TString &path) {
    unique_ptr<TFile> f(TFile::Open(path,"READ"));
    require(f && !f->IsZombie(),string("Cannot open ")+path.Data()); return f;
}
TTree *tree(TFile &f) {
    auto t=dynamic_cast<TTree*>(f.Get("jets"));
    require(t,string("Missing jets tree in ")+f.GetName()); return t;
}
TLeaf *leaf(TTree *t,const TString &name) {
    // Branch names and leaf names need not agree (e.g. ICS_recoJetEta
    // can contain a leaf named recoJetEta). SetBranchAddress in Machine2
    // resolves the branch, so resolve that exact branch first here too.
    TBranch *branch=t->GetBranch(name);
    TLeaf *l=nullptr;
    if(branch) {
        l=branch->GetLeaf(name);
        auto leaves=branch->GetListOfLeaves();
        if(!l && leaves && leaves->GetEntriesFast()==1)
            l=dynamic_cast<TLeaf*>(leaves->At(0));
    } else {
        l=t->GetLeaf(name);
        if(l)branch=l->GetBranch();
    }
    if(!l || !branch) {
        cerr<<"[PriorShapes] Cannot resolve "<<name<<"; available branches:"<<endl;
        auto branches=t->GetListOfBranches();
        for(int i=0;branches && i<branches->GetEntriesFast();++i)
            cerr<<"  "<<branches->At(i)->GetName()<<endl;
        throw runtime_error(string("Missing or ambiguous branch/leaf ")+name.Data()+" in "+t->GetCurrentFile()->GetName());
    }
    require(l->GetLeafCount()==nullptr && l->GetLenStatic()==1,
            string("Expected scalar leaf in branch ")+name.Data());
    t->SetBranchStatus(branch->GetName(),1);
    if(name!=l->GetName())cout<<"[PriorShapes] branch "<<branch->GetName()<<" -> leaf "<<l->GetName()<<endl;
    return l;
}
int central(double c) { return c>=0 && c<=8 && c==floor(c) ? (c>=7?0:c>=4?1:2) : -1; }
bool in(double x,double lo,double hi) { return x>=lo && x<hi; }
void checkEdges(const vector<double> &e) {
    require(e.size()>=2,"Empty binning in config.h");
    for(size_t i=0;i<e.size();++i) require(isfinite(e[i]) && (!i || e[i]>e[i-1]),"Invalid bin edges");
}
TH1D *make(const TString &name,const vector<double> &e) {
    checkEdges(e); auto h=new TH1D(name,"",e.size()-1,e.data());
    h->SetDirectory(nullptr); h->Sumw2(); return h;
}
void style(TH1 *h,int w) {
    h->SetStats(false); h->SetLineColor(colors[w]); h->SetLineStyle(styles[w]); h->SetLineWidth(3);
    h->GetXaxis()->SetTitleSize(.052); h->GetYaxis()->SetTitleSize(.05);
    h->GetXaxis()->SetLabelSize(.044); h->GetYaxis()->SetLabelSize(.044);
    h->GetYaxis()->SetTitleOffset(1.4);
}
void text(double x,double y,const TString &s,double size=.025) {
    TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(size); t.DrawLatex(x,y,s);
}

void page(TCanvas &can,const TString &pdf,TH1D *h[NC][NW],int v,
          const TString &selection) {
    can.Clear(); can.cd();
    text(.06,.97,TString("Prior shapes: ")+titles[v]+"   |   "+selection,.027);
    TLegend leg(.08,.885,.96,.94); leg.SetNColumns(4); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(.021);
    for(int w=0;w<NW;++w) { style(h[0][w],w); leg.AddEntry(h[0][w],labels[w],"l"); }
    leg.Draw();
    vector<unique_ptr<TH1D>> norm;
    for(int c=0;c<NC;++c) {
        TH1D *p[NW]; double ymax=0, ymin=numeric_limits<double>::max(); bool valid=true;
        for(int w=0;w<NW;++w) {
            norm.emplace_back(static_cast<TH1D*>(h[c][w]->Clone(Form("PS_norm_%d_%d",c,w))));
            p[w]=norm.back().get(); p[w]->SetDirectory(nullptr);
            double n=p[w]->Integral(); if(n<=0) {valid=false; continue;}
            p[w]->Scale(1./n,"width"); style(p[w],w);
            for(int b=1;b<=p[w]->GetNbinsX();++b) {double y=p[w]->GetBinContent(b); ymax=max(ymax,y); if(y>0)ymin=min(ymin,y);}
        }
        can.cd(); double left=.01+c/3., right=(c+1)/3.-.007;
        auto up=new TPad(Form("PS_up%d",c),"",left,.40,right,.875);
        up->SetLeftMargin(.17); up->SetRightMargin(.035); up->SetTopMargin(.08); up->SetBottomMargin(.04); up->SetTicks(1,1); up->Draw(); up->cd();
        text(.20,.92,cent[c],.055);
        if(!valid) {text(.2,.5,"Empty selection",.05); continue;}
        up->SetLogy();
        p[0]->SetTitle(""); p[0]->GetYaxis()->SetTitle("Unit-area density"); p[0]->GetXaxis()->SetLabelSize(0);
        p[0]->SetMinimum(max(1e-10,ymin*.45)); p[0]->SetMaximum(ymax*3.);
        p[0]->DrawCopy("HIST"); for(int w=1;w<NW;++w)p[w]->DrawCopy("HIST SAME");
        text(.20,.92,cent[c],.055);
        can.cd(); auto down=new TPad(Form("PS_down%d",c),"",left,.065,right,.385);
        down->SetLeftMargin(.17); down->SetRightMargin(.035); down->SetTopMargin(.03); down->SetBottomMargin(.24); down->SetTicks(1,1); down->Draw(); down->cd();
        TH1D frame(*p[1]); frame.SetDirectory(nullptr); frame.Reset(); frame.SetTitle("");
        double rmin=1,rmax=1;
        TGraph gr[NW];
        for(int w=0;w<NW;++w) {
            gr[w].SetLineColor(colors[w]); gr[w].SetMarkerColor(colors[w]); gr[w].SetMarkerStyle(20+w);
            gr[w].SetMarkerSize(.65); gr[w].SetLineWidth(2); gr[w].SetLineStyle(styles[w]);
            for(int b=1;b<=p[w]->GetNbinsX();++b) {
                double den=p[1]->GetBinContent(b); if(den<=0)continue;
                double r=p[w]->GetBinContent(b)/den;
                gr[w].SetPoint(gr[w].GetN(),p[w]->GetBinCenter(b),r);
                rmin=min(rmin,r); rmax=max(rmax,r);
            }
        }
        double pad=max(.08,(rmax-rmin)*.12);
        frame.SetMinimum(max(0.,rmin-pad)); frame.SetMaximum(rmax+pad);
        frame.GetYaxis()->SetTitle("Shape / (default #times c)"); frame.GetXaxis()->SetTitle(titles[v]);
        frame.GetXaxis()->SetLabelSize(.064); frame.GetXaxis()->SetTitleSize(.075);
        frame.GetYaxis()->SetLabelSize(.058); frame.GetYaxis()->SetTitleSize(.056); frame.GetYaxis()->SetTitleOffset(1.3); frame.GetYaxis()->SetNdivisions(505);
        frame.DrawCopy("AXIS");
        TLine line(frame.GetXaxis()->GetXmin(),1,frame.GetXaxis()->GetXmax(),1); line.SetLineStyle(2); line.SetLineColor(kGray+1); line.DrawClone();
        for(int w=0;w<NW;++w) if(w!=1)gr[w].DrawClone("P SAME");
    }
    can.cd(); text(.06,.023,"Same MC events for all curves; shapes normalised separately. Ratio errors omitted (correlated sample).",.020);
    can.Print(pdf); can.Clear();
}
} // namespace

void PlotPriorShapes(double d0min=5,double d0max=10,
    const char *outputStem="OutputPdf/PriorShapes_5_10",
    const char *mcFile="",const char *dataFile="",const char *fonllFile="Output/FONLL_D0Pt_0_31.root",
    bool include5to7=true,double broadMin=1,double broadMax=10) {
    using namespace PriorShapesQA;
    try {
        require(d0min>=broadMin && d0max<=broadMax && d0max>d0min,"Selected D0 range must lie inside broad range");
        require(useCustomPtMcBins,"This diagnostic expects useCustomPtMcBins=true, as in the supplied configuration");
        require(DividedMcDataBy>0,"DividedMcDataBy must be positive");
        TString mcPath=(*mcFile?TString(mcFile):McJetsFileData);
        TString dataPath=(*dataFile?TString(dataFile):"Output/Output_"+Method+"2.root");
        auto fm=open(mcPath); auto fd=open(dataPath); auto ff=open(fonllFile);
        auto f=dynamic_cast<TF1*>(ff->Get("fSpline")); require(f,"Missing TF1 fSpline in FONLL input");
        auto mt=tree(*fm); auto dt=tree(*fd);
        mt->SetBranchStatus("*",0); dt->SetBranchStatus("*",0);
        TLeaf *dc=leaf(dt,"centrality"), *dp=leaf(dt,"d0Pt");
        TLeaf *dw[4]={leaf(dt,"sWeightSignal"),leaf(dt,"weightD0Efficiency"),leaf(dt,"weightCentrality"),leaf(dt,"weightDoubleCount")};
        TH1D ndata("PS_Ndata","Data centrality sums",9,-.5,8.5), nmc("PS_Nmc","Selected MC centrality sums",9,-.5,8.5);
        ndata.SetDirectory(nullptr); nmc.SetDirectory(nullptr); ndata.Sumw2(); nmc.Sumw2();
        for(Long64_t i=0;i<dt->GetEntries();++i) {
            require(dt->GetEntry(i)>0,"Cannot read data entry");
            if(!in(dp->GetValue(),d0min,d0max))continue;
            double cc=dc->GetValue(); require(central(cc)>=0,"Invalid data centrality index (expected 0..8)");
            double w=1; for(auto l:dw)w*=l->GetValue(); require(isfinite(w),"Nonfinite data weight");
            ndata.Fill(cc,w); // signed sWeights are allowed here, as in Machine2
        }
        TString prefix=Method=="ICS"?"ICS_":"";
        TLeaf *cc=leaf(mt,"centrality"), *wc=leaf(mt,"weightCentrality"), *mx=leaf(mt,"eventMaxPtTrack");
        TLeaf *pt=leaf(mt,"mcJetPt"), *eta=leaf(mt,"mcJetEta"), *nc=leaf(mt,"mcJetNConst");
        TLeaf *re=leaf(mt,prefix+"recoJetEta"), *rn=leaf(mt,prefix+"recoJetNConst");
        TLeaf *dpt=leaf(mt,"mcD0Pt"), *deta=leaf(mt,"mcD0Eta"), *spt=leaf(mt,"mcSmearedD0Pt");
        TLeaf *obs[NV]={pt,leaf(mt,"mcJetD0Z"),leaf(mt,"mcJetLambda1_1"),leaf(mt,"mcJetLambda1_1_5"),leaf(mt,"mcJetLambda1_2"),leaf(mt,"mcJetLambda1_3"),leaf(mt,"mcJetLambda1_0_5"),leaf(mt,"mcJetMomDisp"),dpt};
        mt->SetCacheSize(64*1024*1024); mt->AddBranchToCache("*",true);
        Long64_t entries=mt->GetEntries()/DividedMcDataBy;
        TH1D *before[NC],*after[NC];
        for(int c=0;c<NC;++c) {
            before[c]=new TH1D(Form("PS_denomBefore_c%d",c),"",60,1,31);
            after[c]=new TH1D(Form("PS_denomAfter_c%d",c),"",60,1,31);
            for(auto h:{before[c],after[c]}) {h->SetDirectory(nullptr);h->Sumw2();}
        }
        auto base=[&]() {
            if(central(cc->GetValue())<0 || mx->GetValue()>30)return false;
            double p=dpt->GetValue(), et=deta->GetValue();
            if(!isfinite(p)||!isfinite(et))return false;
            double y=asinh(p*sinh(et)/sqrt(p*p+1.864*1.864));
            if(!isfinite(y)||abs(y)>.6)return false;
            if(abs(eta->GetValue())>.6 && abs(re->GetValue())>.6)return false;
            if(pt->GetValue()>30 || !isfinite(pt->GetValue()))return false;
            if(deleteOneConstituentJets && nc->GetValue()==1 && rn->GetValue()==1)return false;
            return true;
        };
        auto selected=[&](double lo,double hi) {return in(dpt->GetValue(),lo,hi)&&in(spt->GetValue(),lo,hi);};
        cout<<"[PriorShapes] Building both denominators from "<<entries<<" MC entries."<<endl;
        for(Long64_t i=0;i<entries;++i) {
            require(mt->GetEntry(i)>0,"Cannot read MC entry"); if(!base())continue;
            double w=wc->GetValue(); require(isfinite(w)&&w>=0,"Invalid MC weightCentrality");
            int c=central(cc->GetValue());
            if(selected(broadMin,broadMax))after[c]->Fill(pt->GetValue(),w);
            if(selected(d0min,d0max)) {before[c]->Fill(pt->GetValue(),w); nmc.Fill(cc->GetValue(),w);}
        }
        double factor[9];
        for(int j=0;j<9;++j) {
            double m=nmc.GetBinContent(j+1), d=ndata.GetBinContent(j+1);
            require(isfinite(d)&&d>=0,"Negative/nonfinite data centrality sum: not a positive MC prior");
            require(m>0 || d==0,"Data centrality bin has no MC support");
            factor[j]=m>0?d/m:0;
            cout<<"[PriorShapes] cent9="<<j<<" data/MC="<<factor[j]<<endl;
        }
        TH1D *h[2][NV][NC][NW]={};
        TH2D *joint[7][NC][NW]={};
        for(int c=0;c<NC;++c) {
            checkEdges(ptMcBinsVecCustom[c]);
            for(int v=0;v<NV;++v) {
                vector<double> edges;
                if(v==0)edges=ptMcBinsVecCustom[c];
                else if(v==1)edges=zMcBinsVecCustom[c];
                else if(v<8)edges=angMcBinsVecCustom[c][v-2];
                else {int n=int(ceil((d0max-d0min)/.25)); for(int b=0;b<=n;++b)edges.push_back(d0min+(d0max-d0min)*b/n);}
                for(int s=0;s<2;++s)for(int w=0;w<NW;++w)
                    h[s][v][c][w]=make(Form("PS_%s_s%d_c%d_w%d",keys[v],s,c,w),edges);
                if(v>=1&&v<=7)for(int w=0;w<NW;++w) {
                    auto &e=ptMcBinsVecCustom[c];
                    joint[v-1][c][w]=new TH2D(Form("PS_joint_%s_c%d_w%d",keys[v],c,w),"",e.size()-1,e.data(),edges.size()-1,edges.data());
                    joint[v-1][c][w]->SetDirectory(nullptr);joint[v-1][c][w]->Sumw2();
                }
            }
        }
        Long64_t count[NC]={};
        for(Long64_t i=0;i<entries;++i) {
            require(mt->GetEntry(i)>0,"Cannot read MC entry"); if(!base()||!selected(d0min,d0max))continue;
            int c=central(cc->GetValue()); double p=pt->GetValue();
            if(!(abs(eta->GetValue())<.6)|| (deleteOneConstituentJets && nc->GetValue()<=1))continue;
            if(!in(p,ptMcBinsVecCustom[c].front(),ptMcBinsVecCustom[c].back()))continue;
            double db=before[c]->GetBinContent(before[c]->FindFixBin(p));
            double da=after[c]->GetBinContent(after[c]->FindFixBin(p));
            require(wc->GetValue()==0 || (db>0&&da>0),"Empty reshape denominator for selected truth event");
            double F=f->Eval(p)*1e10;
            require(isfinite(F)&&F>=0,"Negative/nonfinite FONLL evaluation");
            double cw=wc->GetValue()*factor[int(cc->GetValue())];
            double w[NW]={1,cw,cw>0?cw*F/db:0.,cw>0?cw*F/da:0.};
            for(double x:w)require(isfinite(x)&&x>=0,"Invalid final weight");
            ++count[c];
            for(int v=0;v<NV;++v) {
                double x=obs[v]->GetValue(); require(isfinite(x),"Nonfinite truth observable");
                auto ax=h[0][v][c][0]->GetXaxis(); if(!in(x,ax->GetXmin(),ax->GetXmax()))continue;
                for(int k=0;k<NW;++k) {
                    if(v>=1&&v<=7)joint[v-1][c][k]->Fill(p,x,w[k]);
                    // Jet pt: full truth acceptance. Other observables: physical 5-20 range.
                    if(v==0||in(p,5,20))h[0][v][c][k]->Fill(x,w[k]);
                    if(in(p,5,7))h[1][v][c][k]->Fill(x,w[k]);
                }
            }
        }
        TString stem=outputStem; require(!stem.IsNull(),"Empty output stem");
        gSystem->mkdir(gSystem->DirName(stem),true);
        TString pdf=stem+".pdf", root=stem+".root";
        TFile out(root,"RECREATE"); require(!out.IsZombie(),"Cannot create output ROOT file");
        TNamed meta("PriorShapes_definition",Form("MC=%s; data=%s; FONLL=%s; Method=%s; D0=[%g,%g), broad=[%g,%g); all MC / DividedMcDataBy=%d; no split; w0=1,w1=c,w2=c*F*1e10/Hselected,w3=c*F*1e10/Hbroad; c=weightCentrality*Ndata/NmcSelected; denominator:60 bins [1,31]; nominal observables: truth jetPt [5,20); jetPt plot: full truth; s1: [5,7); raw sums stored; truth includes misses; no prior tilt; no ratio errors",mcPath.Data(),dataPath.Data(),fonllFile,Method.Data(),d0min,d0max,broadMin,broadMax,DividedMcDataBy));
        meta.Write(); ndata.Write();nmc.Write();
        for(int c=0;c<NC;++c) {
            before[c]->Write();after[c]->Write();
            for(int v=0;v<NV;++v)for(int w=0;w<NW;++w)for(int s=0;s<2;++s)h[s][v][c][w]->Write();
            for(int v=0;v<7;++v)for(int w=0;w<NW;++w)joint[v][c][w]->Write();
        }
        out.Close();
        int oldStat=gStyle->GetOptStat(); gStyle->SetOptStat(0);
        TCanvas can("PriorShapesCanvas","Prior shapes",1500,900); can.Print(pdf+"[");
        can.cd();text(.07,.92,"Prior-shape comparison",.045);
        text(.07,.85,Form("D^{0}: %.1f-%.1f GeV/c (truth AND smeared); %s; three centralities",d0min,d0max,Method.Data()));
        text(.07,.79,"Four variants on exactly the same selected truth events; no train/test split.");
        text(.07,.72,"Default = unit weights. Default #times cweight = centrality correction only.");
        text(.07,.66,"cweight = weightCentrality #times N_{data}(centrality) / N_{MC,selected}(centrality)");
        text(.07,.60,"N_{data} uses sWeightSignal #times weightD0Efficiency #times weightCentrality #times weightDoubleCount.",.022);
        text(.07,.53,Form("Before: cweight #times F(p_{T,jet}) / H_{%.1f-%.1f}(p_{T,jet})",d0min,d0max));
        text(.07,.47,Form("After:   cweight #times F(p_{T,jet}) / H_{%.1f-%.1f}(p_{T,jet})",broadMin,broadMax));
        text(.07,.40,"Both reshape weights include the original factor 10^{10}; denominators use weightCentrality.",.023);
        text(.07,.34,"F is the charm FONLL curve evaluated at jet p_{T}: a chosen shape, not a tagged-jet prediction.",.023);
        text(.07,.27,"PDF curves: unit area, bin-width corrected. Bottom panels: ratio to default #times cweight.",.023);
        text(.07,.21,"Truth includes missed reco jets; Machine2 D^{0}, rapidity, eta and constituent selections retained.",.023);
        text(.07,.15,"Raw weighted histograms and joint (jet p_{T}, observable) priors are saved in the companion ROOT file.",.023);
        text(.07,.09,Form("Selected truth entries: %lld / %lld / %lld. Binning from config.h; no Machine2-local overrides.",count[0],count[1],count[2]),.022);
        can.Print(pdf);can.Clear();
        for(int v=0;v<NV;++v)page(can,pdf,h[0][v],v,Form("D^{0} %.1f-%.1f; %s",d0min,d0max,v==0?"full truth jet-p_{T} range":"5 #leq p_{T,jet}^{true} < 20 GeV/c"));
        if(include5to7)for(int v=1;v<NV;++v)page(can,pdf,h[1][v],v,Form("D^{0} %.1f-%.1f; 5 #leq p_{T,jet}^{true} < 7 GeV/c",d0min,d0max));
        can.Print(pdf+"]");gStyle->SetOptStat(oldStat);
        cout<<"[PriorShapes] Saved "<<pdf<<" and "<<root<<endl;
        for(int c=0;c<NC;++c) {
            delete before[c];delete after[c];
            for(int v=0;v<NV;++v)for(int w=0;w<NW;++w)for(int s=0;s<2;++s)delete h[s][v][c][w];
            for(int v=0;v<7;++v)for(int w=0;w<NW;++w)delete joint[v][c][w];
        }
    } catch(const exception &e) {
        cerr<<"[PriorShapes ERROR] "<<e.what()<<endl;
        throw;
    }
}
