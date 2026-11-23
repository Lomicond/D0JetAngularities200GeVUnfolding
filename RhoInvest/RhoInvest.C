// RhoInvest: data / embedding background-density diagnostic.
// Place in <project>/RhoInvest/ and run FROM <project> (next to config.h):
//   root -l -b -q 'RhoInvest/RhoInvest.C(1,10)'
// Optional reconstructed jet-pT selection:
//   root -l -b -q 'RhoInvest/RhoInvest.C(1,10,5,20)'
// No arguments: D0 interval from config.h, no extra jet-pT cut.
// Data path defaults to ./Output/<outRoot>_<Method><_sys>2.root;
// MC path is McJetsFileData. The raw RealJetsFileData has no sWeights and is NOT used.
// Fifth argument: optional prepared sPlot file; sixth: optional output stem.
// Arguments 7-10: plotting bin counts/ranges. Argument 11: fitDegree (default 1).
// fitDegree=0 disables the new correction; 2 allows quadratic log-weight shapes.
// Example quadratic check (same selection as the default run):
// root -l -b -q 'RhoInvest/RhoInvest.C(1,5,-1e9,1e9,"","",120,600,200,100,2)'
// Fit weights are learned from INPUT TREES, not old RhoInvest ROOT outputs.
// Output in RhoInvest/output/: 13-page PDF, ROOT graphs/raw rho histograms,
// fine-centrality yields/factors in TSV, and a settings/provenance text file.
// PDF: mean/ratio/difference vs gRefMult, vs 9 fine centralities, vs 3 broad
// centralities, then normalized rho distributions in the 3 broad centralities.
// Nine extra pages show the gRefMult fit, ratio and extra weight per fine class.
// New default output stems contain _gRefFit1 (or requested order); old files survive.
// _grefFit.tsv records coefficients, fit decisions, effective counts and yields.
// ROOT TF1 gRefExtraWeight_codeN multiplies EXISTING centrality-reweighted MC.
// ROOT TF1 gRefTotalWeight_codeN instead multiplies BASE weightCentrality.
// Do not multiply both functions together. Validate these diagnostic fits before
// using them for a response matrix. No changes to Machine.C are made.
// The positive extra weight is exp(a*x+b*x*x)/Z_c; x is frozen outside MC 1-99%
// quantiles and Z_c is recomputed exactly in a second MC pass. The pre-cut
// normalization yields in all nine classes are preserved. Diagnostic yields may
// change if selection efficiencies depend on gRefMult. Negative sWeights remain.
// Fits retain a normalization nuisance (discarded before reweighting);
// chi2_fit refers to that fit. MC cell integration uses width 1 gRefMult;
// event-level application is continuous. Neff<40 retains unity, quadratic needs
// Neff>=150, unstable/boundary fits lower the degree. All decisions are logged.
// Fit covariance and data/MC histograms are saved. Displayed rho errors do not
// include fitted-weight uncertainty or training-sample correlations.
// Weighted histogram errors: https://root.cern.ch/doc/master/classTH1.html
// Defaults use gRefMult bins of width 5; set gRefBins=600 for the old width 1.
// MC "base" ALREADY includes weightCentrality. "Weighted" additionally uses
// Ndata(code)/Nmc(code), derived before the diagnostic reco-jet selections,
// exactly as in Machine::LoadDataRealParams. No FONLL/prior-shape weights.
// MC truth preselection is retained; data inherits all upstream sPlot selections.
// This is a candidate/jet-weighted comparison, not an unweighted MB-event study.
// Requires the same ROOT / RooUnfold setup as Unfolding/Machine.C.
// Reference: Lomicond/D0JetAngularities200GeVUnfolding,
// commit 9541f387d53d150b063875c9ee3a01fa646632da (2026-11-16).
// This reads STORED rho. It neither changes occupancy correction nor recomputes rho.
// Signed sWeights are retained. Errors are fixed-weight, independent-entry errors:
// no fitted-weight uncertainty, centrality-map uncertainty or repeated-event covariance.

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <map>
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TMatrixDSym.h"
#include "TF1.h"
#include <stdexcept>
#include <string>
#include <vector>
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TBranch.h"
#include "TObjArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TNamed.h"
using namespace std;
#include "../config.h"

namespace RhoInvestDetail {
constexpr int NC = 9;
const double centEdges[] = {0,5,10,20,30,40,50,60,70,80};
const char* broadLabels[] = {"0-10%", "10-40%", "40-80%"};
int broad(int c) { return c >= 7 ? 0 : c >= 4 ? 1 : 2; }
void require(bool ok, const string& msg) { if (!ok) throw runtime_error(msg); }

// S = sum w, R = sum w*rho, Qk = sum w^2*rho^k.
// These sufficient statistics allow exact reweighting without another tree scan.
struct Moment {
    double s=0, r=0, q0=0, q1=0, q2=0;
    void fill(double rho, double w) { s+=w; r+=w*rho; q0+=w*w; q1+=w*w*rho; q2+=w*w*rho*rho; }
    void add(const Moment& a, double f=1) { s+=f*a.s; r+=f*a.r; q0+=f*f*a.q0; q1+=f*f*a.q1; q2+=f*f*a.q2; }
    bool mean(double& value, double& error) const {
        if (!(s>0) || !isfinite(s)) return false;
        value=r/s;
        const double var=(q2-2*value*q1+value*value*q0)/(s*s);
        error=sqrt(max(0.,var));
        return isfinite(value) && isfinite(error);
    }
};
// A positive, slowly varying correction to the MC gRefMult template.
// x is frozen outside the central 98% of the MC distribution (no tail extrapolation).
struct GrefFit {
    int degree=0, ndf=0, status=-1, covStatus=0;
    double center=0, scale=1, a=0, b=0, normalization=1;
    double neff=0, chi2=0, chi2Before=0, chi2Normalized=0, amplitude=1, minWeight=1, maxWeight=1;
    double cov[9]={0}; // parameters (log amplitude, a, b), row-major
    string reason="disabled";
    vector<double> edges;
    double raw(double g) const {
        double x=max(-1.,min(1.,(g-center)/scale));
        return exp(a*x+b*x*x);
    }
    double eval(double g) const { return raw(g)/normalization; }
};
using GrefFits=array<GrefFit,NC>;
struct Reader {
    TTree* tree;
    vector<TLeaf*> leaves;
    explicit Reader(TTree* t): tree(t) { tree->SetBranchStatus("*",0); tree->SetCacheSize(64*1024*1024); }
    int bind(const TString& name) {
        // Branch and leaf names need not match. Resolve the scalar leaf ONLY
        // within the requested branch; an unqualified GetLeaf can miss ICS_
        // branches or select an identically named leaf from an AREA branch.
        auto* branch=tree->GetBranch(name);
        require(branch!=nullptr, "Missing required branch '"+string(name.Data())+"' in "+tree->GetCurrentFile()->GetName());
        auto* branchLeaves=branch->GetListOfLeaves();
        require(branchLeaves && branchLeaves->GetEntriesFast()==1,
                "Expected exactly one scalar leaf in branch '"+string(name.Data())+"'");
        auto* leaf=dynamic_cast<TLeaf*>(branchLeaves->At(0));
        require(leaf!=nullptr, "Cannot resolve scalar leaf in branch '"+string(name.Data())+"'");
        require(leaf->GetLeafCount()==nullptr && leaf->GetLenStatic()==1, "Expected scalar leaf: "+string(name.Data()));
        tree->SetBranchStatus(leaf->GetBranch()->GetName(),1);
        tree->AddBranchToCache(leaf->GetBranch()->GetName(),true);
        leaves.push_back(leaf); return int(leaves.size())-1;
    }
    double v(int i) const { return leaves.at(i)->GetValue(); }
    void get(Long64_t i) {
        require(tree->GetEntry(i)>0,"Failed reading tree entry "+to_string(i));
        for(auto* l:leaves) require(isfinite(l->GetValue()), "Non-finite value in "+string(l->GetName())+" at entry "+to_string(i));
    }
};
struct Sample {
    array<Moment,NC> norm, selected;
    array<map<int,Moment>,NC> normG; // unit-width gRefMult cells, before diagnostic cuts
    vector<vector<Moment>> ref, rho;
    int ng, nr; double gmax, rmax;
    Long64_t accepted=0, negative=0, badCent=0;
    Sample(int g, double gm, int r, double rm):ref(NC,vector<Moment>(g+2)),rho(NC,vector<Moment>(r+2)),ng(g),nr(r),gmax(gm),rmax(rm) {}
    static int bin(double x,int n,double xmax) { return x<0 ? 0 : x>=xmax ? n+1 : 1+int(x/xmax*n); }
    void fillNorm(int c,double g,double w) {
        require(g>=0 && g<1000000,"Invalid gRefMult in normalization sample");
        norm[c].fill(0,w); normG[c][int(floor(g))].fill(g,w);
    }
    void scaleClass(int c,double f) {
        auto scale=[f](Moment& m){ Moment r; r.add(m,f); m=r; };
        scale(norm[c]); scale(selected[c]);
        for(auto& kv:normG[c]) scale(kv.second);
        for(auto& m:ref[c]) scale(m);
        for(auto& m:rho[c]) scale(m);
    }
    void fill(int c,double g,double y,double w) {
        selected[c].fill(y,w); ref[c][bin(g,ng,gmax)].fill(y,w); rho[c][bin(y,nr,rmax)].fill(y,w);
        ++accepted; if(w<0) ++negative;
    }
};
int code(double x) {
    const int c=int(lround(x));
    return (c>=0 && c<NC && abs(x-c)<1e-5) ? c : -1;
}
bool in(double x,double lo,double hi) { return x>=lo && x<hi; }
unique_ptr<TFile> open(const TString& name) {
    unique_ptr<TFile> f(TFile::Open(name,"READ"));
    require(f && !f->IsZombie(),"Cannot open "+string(name.Data())); return f;
}
TTree* tree(TFile* f) {
    auto* t=dynamic_cast<TTree*>(f->Get("jets"));
    require(t!=nullptr,"Missing jets tree in "+string(f->GetName())); return t;
}
void progress(Long64_t i,Long64_t n,const char* label) {
    if(i%1000000==0) cout << label << ": " << i << "/" << n << endl;
}
void readData(Sample& a,const TString& file,double dmin,double dmax,double jmin,double jmax) {
    auto f=open(file); Reader rd(tree(f.get()));
    const int dc=rd.bind("d0Pt"), c=rd.bind("centrality"), ca=rd.bind("centralityAlt");
    const int sw=rd.bind("sWeightSignal"), eff=rd.bind("weightD0Efficiency"), cw=rd.bind("weightCentrality"), dw=rd.bind("weightDoubleCount");
    const int jpt=rd.bind("jetPt"), eta=rd.bind("jetEta"), nc=rd.bind("nJetConst"), r=rd.bind("rho"), g=rd.bind("gRefMult");
    const Long64_t n=rd.tree->GetEntries();
    for(Long64_t i=0;i<n;++i) {
        progress(i,n,"Data"); rd.get(i);
        if(!in(rd.v(dc),dmin,dmax)) continue;
        const int cc=code(rd.v(c));
        require(cc>=0,"Invalid data centrality code at entry "+to_string(i));
        require(in(rd.v(ca),0,80),"Invalid data centralityAlt at entry "+to_string(i));
        require(broad(cc)==(rd.v(ca)<10?0:rd.v(ca)<40?1:2),"Inconsistent data centrality / centralityAlt");
        const double w=rd.v(sw)*rd.v(eff)*rd.v(cw)*rd.v(dw);
        require(isfinite(w),"Non-finite product of data weights");
        a.fillNorm(cc,rd.v(g),w); // EXACT placement of nCentralityNumbers in Machine.
        if(abs(rd.v(eta))>0.6 || rd.v(nc)<1 || !in(rd.v(jpt),jmin,jmax)) continue;
        if(deleteOneConstituentJets && rd.v(nc)==1) continue;
        a.fill(cc,rd.v(g),rd.v(r),w);
    }
}
void readMC(Sample& a,const TString& file,double dmin,double dmax,double jmin,double jmax,const GrefFits* fit=nullptr) {
    auto f=open(file); Reader rd(tree(f.get()));
    const TString prefix=Method=="ICS"?"ICS_":"";
    const int c=rd.bind("centrality"), w=rd.bind("weightCentrality"), maxtrack=rd.bind("eventMaxPtTrack");
    const int dpt=rd.bind("mcD0Pt"), deta=rd.bind("mcD0Eta"), rpt=rd.bind("mcSmearedD0Pt");
    const int jt=rd.bind("mcJetPt"), et=rd.bind("mcJetEta"), nt=rd.bind("mcJetNConst");
    const int jr=rd.bind(prefix+"recoJetPt"+(Method=="ICS"?"":"Corr"));
    const int er=rd.bind(prefix+"recoJetEta"), nr=rd.bind(prefix+"recoJetNConst");
    const int rr=rd.bind(prefix+"recoJetRho"), gr=rd.bind("gRefMult");
    require(DividedMcDataBy>0,"DividedMcDataBy must be positive");
    // Machine's integer loop condition against a double corresponds to ceil.
    const Long64_t n=Long64_t(ceil(double(rd.tree->GetEntries())/DividedMcDataBy));
    for(Long64_t i=0;i<n;++i) {
        progress(i,n,"MC"); rd.get(i);
        const int cc=code(rd.v(c)); if(cc<0) { ++a.badCent; continue; }
        if(rd.v(maxtrack)>30 || !in(rd.v(dpt),dmin,dmax) || !in(rd.v(rpt),dmin,dmax)) continue;
        const double pz=rd.v(dpt)*sinh(rd.v(deta)), p=rd.v(dpt)*cosh(rd.v(deta));
        const double e=sqrt(p*p+1.864*1.864), y=0.5*log((e+pz)/(e-pz));
        require(isfinite(y),"Non-finite MC D0 rapidity");
        if(abs(y)>0.6 || (abs(rd.v(et))>0.6 && abs(rd.v(er))>0.6) || rd.v(jt)>30) continue;
        if(deleteOneConstituentJets && rd.v(nt)==1 && rd.v(nr)==1) continue;
        require(rd.v(w)>=0,"The smooth-template fit requires nonnegative MC base weights");
        const double ww=rd.v(w)*(fit?(*fit)[cc].raw(rd.v(gr)):1.);
        a.fillNorm(cc,rd.v(gr),ww); // EXACT placement of nCentralityNumbersMC.
        if(abs(rd.v(er))>0.6 || rd.v(nr)<1 || !in(rd.v(jr),jmin,jmax)) continue;
        if(deleteOneConstituentJets && rd.v(nr)==1) continue;
        a.fill(cc,rd.v(gr),rd.v(rr),ww);
    }
}
struct FitGroup {
    int lo=0, hi=0;
    Moment d,m;
    vector<pair<double,double>> cells; // MC-weighted gRefMult mean, sum of weights
};
Moment rangeSum(const map<int,Moment>& cells,int lo,int hi) {
    Moment r;
    for(auto it=cells.lower_bound(lo);it!=cells.end() && it->first<hi;++it) r.add(it->second);
    return r;
}
double quantile(const map<int,Moment>& cells,double fraction) {
    double total=0; for(const auto& kv:cells) total+=kv.second.s;
    double cum=0;
    for(const auto& kv:cells) {
        cum+=kv.second.s;
        if(cum>=fraction*total && kv.second.s>0) return kv.second.r/kv.second.s;
    }
    return cells.empty()?0:cells.rbegin()->first+0.5;
}
GrefFit fitGref(const Sample& data,const Sample& mc,int c,int requested) {
    GrefFit result;
    const auto& d=data.normG[c]; const auto& m=mc.normG[c];
    const double nd=data.norm[c].s, nm=mc.norm[c].s;
    if(d.empty() || m.empty() || !(nd>0) || !(nm>0)) { result.reason="empty_class_unity"; return result; }
    result.neff=nd*nd/data.norm[c].q0;
    const int lo=min(d.begin()->first,m.begin()->first), hi=max(d.rbegin()->first,m.rbegin()->first)+1;
    const double qlo=quantile(m,.01), qhi=quantile(m,.99);
    result.center=(qlo+qhi)/2.; result.scale=max(.5,(qhi-qlo)/2.);
    // Quantile bin boundaries depend only on MC and the total effective sample size,
    // not on local positive/negative data fluctuations. All tails remain included.
    const int nb=max(4,min(12,int(result.neff/15.)));
    vector<int> edges={lo}; double cumulative=0; int target=1;
    for(const auto& kv:m) {
        cumulative+=kv.second.s;
        if(target<nb && cumulative>=nm*target/nb) {
            const int edge=kv.first+1;
            if(edge>edges.back() && edge<hi) edges.push_back(edge);
            while(target<nb && cumulative>=nm*target/nb) ++target;
        }
    }
    if(edges.back()<hi) edges.push_back(hi);
    vector<FitGroup> groups;
    for(size_t j=0;j+1<edges.size();++j) {
        FitGroup g; g.lo=edges[j]; g.hi=edges[j+1];
        g.d=rangeSum(d,g.lo,g.hi); g.m=rangeSum(m,g.lo,g.hi);
        for(auto it=m.lower_bound(g.lo);it!=m.end() && it->first<g.hi;++it)
            if(it->second.s>0) g.cells.push_back({it->second.r/it->second.s,it->second.s});
        groups.push_back(g);
    }
    // A truly empty data cell has no estimated sumw2. Merge it with a neighbour;
    // never delete negative bins or assign a Poisson error to signed contents.
    for(size_t i=0;i<groups.size() && groups.size()>1;) {
        if(groups[i].d.q0>0 && groups[i].m.s>0) { ++i; continue; }
        size_t k=(i+1<groups.size()?i+1:i-1), first=min(i,k), second=max(i,k);
        groups[first].hi=groups[second].hi;
        groups[first].d.add(groups[second].d); groups[first].m.add(groups[second].m);
        groups[first].cells.insert(groups[first].cells.end(),groups[second].cells.begin(),groups[second].cells.end());
        groups.erase(groups.begin()+second); i=first;
    }
    result.edges.push_back(groups.front().lo);
    for(const auto& g:groups) result.edges.push_back(g.hi);
    const double centralFactor=nd/nm;
    auto objective=[&](const double* p) {
        vector<double> prediction(groups.size(),0.); double total=0;
        for(size_t j=0;j<groups.size();++j) for(const auto& cell:groups[j].cells) {
            const double x=max(-1.,min(1.,(cell.first-result.center)/result.scale));
            prediction[j]+=cell.second*exp(p[1]*x+p[2]*x*x);
        }
        for(double v:prediction) total+=v;
        if(!(total>0)) return 1.e100;
        double chi=0;
        for(size_t j=0;j<groups.size();++j) {
            const double variance=groups[j].d.q0+centralFactor*centralFactor*groups[j].m.q0;
            if(!(variance>0)) continue;
            const double delta=groups[j].d.s-nd*exp(p[0])*prediction[j]/total;
            chi+=delta*delta/variance;
        }
        return chi;
    };
    // Baseline chi2 is explicitly evaluated with the preserved nominal yield.
    double zero[3]={0,0,0}; result.chi2Before=objective(zero); result.chi2=result.chi2Before; result.chi2Normalized=result.chi2Before;
    result.ndf=max(0,int(groups.size())-1);
    int degree=requested;
    if(degree==2 && result.neff<150) degree=1;
    if(degree>0 && result.neff<40) degree=0;
    if(requested==0) result.reason="disabled_unity";
    else if(result.neff<40) result.reason="low_statistics_unity";
    else result.reason="fit_failed_unity";
    for(int order=degree;order>=1;--order) {
        if(int(groups.size())<order+3 || qhi-qlo<1.) { result.reason="insufficient_shape_bins_unity"; continue; }
        unique_ptr<ROOT::Math::Minimizer> fitter(ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad"));
        require(bool(fitter),"Cannot create Minuit2 minimizer");
        ROOT::Math::Functor fun(objective,3); fitter->SetFunction(fun);
        fitter->SetPrintLevel(0); fitter->SetMaxFunctionCalls(20000); fitter->SetTolerance(1.e-5); fitter->SetErrorDef(1.);
        fitter->SetLimitedVariable(0,"log_amplitude",0,.02,-2.,2.);
        fitter->SetLimitedVariable(1,"a",0,.02,-4.,4.);
        if(order==2) fitter->SetLimitedVariable(2,"b",0,.02,-2.,2.);
        else fitter->SetFixedVariable(2,"b",0);
        const bool converged=fitter->Minimize(); const bool hesse=fitter->Hesse();
        result.status=fitter->Status(); result.covStatus=fitter->CovMatrixStatus();
        const double* v=fitter->X();
        const bool boundary=abs(v[0])>1.98 || abs(v[1])>3.98 || (order==2 && abs(v[2])>1.98);
        if(!converged || !hesse || result.status!=0 || result.covStatus<3 || boundary || !isfinite(fitter->MinValue())) {
            result.reason=boundary?"boundary_fit_unity":"fit_failed_unity"; continue;
        }
        result.degree=order; result.a=v[1]; result.b=v[2]; result.amplitude=exp(v[0]);
        result.chi2=fitter->MinValue(); result.ndf=int(groups.size())-order-1;
        double normalized[3]={0,result.a,result.b}; result.chi2Normalized=objective(normalized);
        for(int i=0;i<3;++i) for(int j=0;j<3;++j) result.cov[3*i+j]=fitter->CovMatrix(i,j);
        result.reason=order==requested?"fit_ok":"reduced_degree";
        return result;
    }
    return result;
}

vector<Moment> aggregate(const Sample& a,const array<double,NC>& f,int mode) {
    // mode 0=gRefMult (visible bins), 1=nine fine classes, 2=three broad classes.
    vector<Moment> out(mode==0?a.ng:mode==1?9:3);
    for(int c=0;c<NC;++c) {
        if(mode==0) { for(int b=0;b<a.ng;++b) out[b].add(a.ref[c][b+1],f[c]); }
        else out[mode==1?8-c:broad(c)].add(a.selected[c],f[c]);
    }
    return out;
}
TGraphErrors* graph(const vector<Moment>& a,const vector<double>& edges,const TString& name) {
    auto* g=new TGraphErrors; g->SetName(name);
    for(size_t b=0;b<a.size();++b) {
        double y,ey; if(!a[b].mean(y,ey)) continue;
        const int i=g->GetN(); g->SetPoint(i,(edges[b]+edges[b+1])/2,y); g->SetPointError(i,(edges[b+1]-edges[b])/2,ey);
    }
    return g;
}
TGraphErrors* compare(const vector<Moment>& d,const vector<Moment>& m,const vector<double>& edges,const TString& name,bool ratio) {
    auto* g=new TGraphErrors; g->SetName(name);
    for(size_t b=0;b<d.size();++b) {
        double x,ex,y,ey; if(!d[b].mean(x,ex) || !m[b].mean(y,ey) || (ratio && y<=0)) continue;
        const double value=ratio?x/y:x-y;
        const double err=ratio?sqrt(ex*ex/(y*y)+x*x*ey*ey/pow(y,4)):hypot(ex,ey);
        const int i=g->GetN(); g->SetPoint(i,(edges[b]+edges[b+1])/2,value); g->SetPointError(i,(edges[b+1]-edges[b])/2,err);
    }
    return g;
}
void style(TGraphErrors* g,int color,int marker) {
    g->SetLineColor(color); g->SetMarkerColor(color); g->SetMarkerStyle(marker); g->SetMarkerSize(0.65); g->SetLineWidth(1);
}
void panel(const vector<TGraphErrors*>& gs,const vector<const char*>& labels,const char* xt,const char* yt,double xmin,double xmax,double reference=numeric_limits<double>::quiet_NaN()) {
    gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.07); gPad->SetRightMargin(0.025);
    double low=isfinite(reference)?reference:numeric_limits<double>::infinity(), high=isfinite(reference)?reference:-numeric_limits<double>::infinity();
    for(auto* g:gs) for(int i=0;i<g->GetN();++i) { double x,y; g->GetPoint(i,x,y); low=min(low,y-g->GetErrorY(i)); high=max(high,y+g->GetErrorY(i)); }
    if(!isfinite(low) || !isfinite(high)) { low=0; high=1; }
    double span=high-low; if(span<1e-9) span=max(0.1,abs(high)*0.2);
    auto* frame=gPad->DrawFrame(xmin,low-0.12*span,xmax,high+0.50*span,Form(";%s;%s",xt,yt));
    frame->GetXaxis()->SetTitleSize(.048); frame->GetYaxis()->SetTitleSize(.048); frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetXaxis()->CenterTitle(); frame->GetYaxis()->CenterTitle();
    if(isfinite(reference)) { TLine line(xmin,reference,xmax,reference); line.SetLineStyle(2); line.SetLineColor(kGray+1); line.DrawClone(); }
    for(auto* g:gs) g->Draw("P SAME");
    TLegend legend(.20,.76,.97,.92); legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(.029);
    for(size_t i=0;i<gs.size();++i) legend.AddEntry(gs[i],labels[i],gs[i]->GetMarkerStyle()==1?"l":"pe"); legend.DrawClone();
}
void header(TCanvas& cv,const TString& title,const TString& subtitle) {
    cv.cd(); TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(.034); tx.DrawLatex(.035,.966,title);
    tx.SetTextSize(.024); tx.DrawLatex(.035,.928,subtitle);
}
void grid(TCanvas& cv) {
    cv.Clear(); auto* p=new TPad("ri_panels","",0,.04,1,.91); p->Draw(); p->Divide(3,1,.001,.001);
}
void footer(TCanvas& cv) {
    cv.cd(); TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(.017);
    t.DrawLatex(.035,.015,"Signed sWeights retained; errors conditional on weights, independent entries. Stored #rho; occupancy definition unchanged.");
}
void profilePage(TCanvas& cv,const TString& pdf,TFile& out,const Sample& d,const Sample& m,const Sample& fitted,const array<double,NC>& factors,int mode,const TString& subtitle) {
    array<double,NC> one; one.fill(1);
    auto a=aggregate(d,one,mode), b=aggregate(m,one,mode), z=aggregate(m,factors,mode), t=aggregate(fitted,factors,mode);
    vector<double> edges;
    if(mode==0) for(int i=0;i<=d.ng;++i) edges.push_back(i*d.gmax/d.ng);
    else if(mode==1) edges.assign(centEdges,centEdges+10);
    else edges={0,10,40,80};
    const TString tag=mode==0?"gRefMult":mode==1?"centrality9":"centrality3";
    vector<TGraphErrors*> means={graph(a,edges,tag+"_data"),graph(b,edges,tag+"_mc_base"),graph(z,edges,tag+"_mc_reweighted"),graph(t,edges,tag+"_mc_grefFit")};
    vector<TGraphErrors*> ratios={compare(a,b,edges,tag+"_ratio_base",true),compare(a,z,edges,tag+"_ratio_reweighted",true),compare(a,t,edges,tag+"_ratio_grefFit",true)};
    vector<TGraphErrors*> diffs={compare(a,b,edges,tag+"_difference_base",false),compare(a,z,edges,tag+"_difference_reweighted",false),compare(a,t,edges,tag+"_difference_grefFit",false)};
    style(means[0],kBlack,20); style(means[1],kGray+2,24); style(means[2],kBlue+1,25); style(means[3],kGreen+2,26);
    for(auto* v:{&ratios,&diffs}) { style((*v)[0],kGray+2,24); style((*v)[1],kBlue+1,25); style((*v)[2],kGreen+2,26); }
    out.cd(); for(auto* v:{&means,&ratios,&diffs}) for(auto* g:*v) g->Write();
    grid(cv); auto* p=static_cast<TPad*>(cv.GetPrimitive("ri_panels"));
    const char* xt=mode==0?"gRefMult":"Centrality (%)";
    p->cd(1); panel(means,{"D^{0} data","Embedding: base","Embedding: centrality","Embedding: centrality + gRefMult"},xt,"#LT#rho#GT (GeV/c)",edges.front(),edges.back());
    p->cd(2); panel(ratios,{"Data / embedding: base","Data / embedding: centrality","Data / embedding: centrality + gRefMult"},xt,"#LT#rho#GT_{data} / #LT#rho#GT_{emb}",edges.front(),edges.back(),1);
    p->cd(3); panel(diffs,{"Data - embedding: base","Data - embedding: centrality","Data - embedding: centrality + gRefMult"},xt,"#LT#rho#GT_{data} - #LT#rho#GT_{emb} (GeV/c)",edges.front(),edges.back(),0);
    header(cv,mode==0?"Background density vs gRefMult":mode==1?"Background density in nine centrality intervals":"Background density in three centrality intervals",subtitle); footer(cv); cv.Print(pdf);
    cv.Clear(); for(auto* v:{&means,&ratios,&diffs}) for(auto* g:*v) delete g;
}
TGraphErrors* shape(const Sample& a,const array<double,NC>& f,int bc,const TString& name,TFile& out) {
    vector<Moment> bins(a.nr+2); Moment total;
    for(int c=0;c<NC;++c) if(broad(c)==bc) { total.add(a.selected[c],f[c]); for(int i=0;i<a.nr+2;++i) bins[i].add(a.rho[c][i],f[c]); }
    auto* g=new TGraphErrors; g->SetName(name);
    TH1D raw(name+"_raw",";#rho (GeV/c);Weighted entries",a.nr,0,a.rmax); raw.SetDirectory(nullptr); raw.Sumw2();
    for(int i=0;i<a.nr+2;++i) { raw.SetBinContent(i,bins[i].s); raw.SetBinError(i,sqrt(bins[i].q0)); }
    out.cd(); raw.Write();
    if(total.s>0) for(int i=1;i<=a.nr;++i) {
        if(bins[i].q0==0) continue;
        const double fraction=bins[i].s/total.s, width=a.rmax/a.nr;
        // Includes normalization covariance for fixed independent entry weights.
        const double var=(bins[i].q0*(1-2*fraction)+fraction*fraction*total.q0)/(total.s*total.s);
        const int n=g->GetN(); g->SetPoint(n,(i-.5)*width,fraction/width); g->SetPointError(n,.5*width,sqrt(max(0.,var))/width);
    }
    if(bins.front().q0>0 || bins.back().q0>0) cerr << "WARNING: rho under/overflow in " << name << "; included in normalization; increase rhoMax if needed.\n";
    g->Write(); return g;
}
void shapePage(TCanvas& cv,const TString& pdf,TFile& out,const Sample& d,const Sample& m,const Sample& fitted,const array<double,NC>& f,const TString& subtitle) {
    array<double,NC> one; one.fill(1); grid(cv); auto* p=static_cast<TPad*>(cv.GetPrimitive("ri_panels"));
    vector<TGraphErrors*> all;
    for(int c=0;c<3;++c) {
        vector<TGraphErrors*> gs={shape(d,one,c,Form("rho_cent%d_data",c),out),shape(m,one,c,Form("rho_cent%d_mc_base",c),out),shape(m,f,c,Form("rho_cent%d_mc_reweighted",c),out),shape(fitted,f,c,Form("rho_cent%d_mc_grefFit",c),out)};
        style(gs[0],kBlack,20); style(gs[1],kGray+2,24); style(gs[2],kBlue+1,25); style(gs[3],kGreen+2,26);
        double xmax=0;
        for(auto* g:gs) for(int i=0;i<g->GetN();++i) { double x,y; g->GetPoint(i,x,y); xmax=max(xmax,x+g->GetErrorX(i)); }
        xmax=xmax>0?min(d.rmax,xmax*1.08):d.rmax;
        p->cd(c+1); panel(gs,{"D^{0} data","Embedding: base","Embedding: centrality","Embedding: centrality + gRefMult"},"#rho (GeV/c)","(1/#Sigma w) d#Sigma w/d#rho",0,xmax);
        TLatex t; t.SetNDC(); t.SetTextSize(.04); t.DrawLatex(.20,.95,broadLabels[c]); all.insert(all.end(),gs.begin(),gs.end());
    }
    header(cv,"Normalized background-density distributions",subtitle); footer(cv); cv.Print(pdf); cv.Clear(); for(auto* g:all) delete g;
}
void writeFits(TFile& out,const TString& stem,const Sample& d,const Sample& m,const Sample& fitted,
               const GrefFits& fits,const array<double,NC>& factors) {
    ofstream log((stem+"_grefFit.tsv").Data()); require(bool(log),"Cannot create gRefMult fit table");
    log << setprecision(17) << "code\tcent_min\tcent_max\tdegree\tdecision\tNeff_data\tfit_status\tcov_status\tchi2_before\tchi2_fit\tchi2_preserved_yield\tndf\tfit_amplitude\ta\tb\ta_error\tb_error\tg_center\tg_scale\texact_normalization\tmin_extra_weight\tmax_extra_weight\tcentrality_factor\ttarget_yield\tfitted_yield\n";
    for(int c=0;c<NC;++c) {
        const auto& f=fits[c]; int cb=9-c;
        log << c << '\t' << centEdges[cb-1] << '\t' << centEdges[cb] << '\t' << f.degree << '\t' << f.reason << '\t' << f.neff << '\t' << f.status << '\t' << f.covStatus << '\t' << f.chi2Before << '\t' << f.chi2 << '\t' << f.chi2Normalized << '\t' << f.ndf << '\t' << f.amplitude << '\t' << f.a << '\t' << f.b << '\t' << sqrt(max(0.,f.cov[4])) << '\t' << sqrt(max(0.,f.cov[8])) << '\t' << f.center << '\t' << f.scale << '\t' << f.normalization << '\t' << f.minWeight << '\t' << f.maxWeight << '\t' << factors[c] << '\t' << d.norm[c].s << '\t' << factors[c]*fitted.norm[c].s << '\n';
        out.cd();
        const TString xx="TMath::Max(-1.,TMath::Min(1.,(x-[2])/[3]))";
        const TString formula="[4]*exp([0]*("+xx+")+[1]*("+xx+")*("+xx+"))";
        const double xmax=m.normG[c].empty()?1000:m.normG[c].rbegin()->first+1.;
        TF1 extra(Form("gRefExtraWeight_code%d",c),formula,0,xmax);
        extra.SetParameters(f.a,f.b,f.center,f.scale,1./f.normalization);
        extra.SetTitle("Extra shape weight;gRefMult;extra weight (mean 1 in MC normalization sample)"); extra.Write();
        extra.SetName(Form("gRefTotalWeight_code%d",c)); extra.SetParameter(4,factors[c]/f.normalization);
        extra.SetTitle("Centrality times extra shape weight;gRefMult;factor multiplying MC weightCentrality"); extra.Write();
        TMatrixDSym covariance(3); for(int i=0;i<3;++i) for(int j=0;j<3;++j) covariance(i,j)=f.cov[3*i+j];
        covariance.Write(Form("gRefFitCovariance_code%d",c));
        TNamed decision(Form("gRefFitDecision_code%d",c),f.reason.c_str()); decision.Write();
        int bins=1;
        for(const auto* s:{&d,&m,&fitted}) if(!s->normG[c].empty()) bins=max(bins,s->normG[c].rbegin()->first+1);
        auto save=[&](const Sample& sample,const char* label,double mult) {
            TH1D hist(Form("gRefNorm_%s_code%d",label,c),";gRefMult;Weighted entries",bins,0,bins); hist.SetDirectory(nullptr); hist.Sumw2();
            for(const auto& kv:sample.normG[c]) { hist.SetBinContent(kv.first+1,mult*kv.second.s); hist.SetBinError(kv.first+1,abs(mult)*sqrt(kv.second.q0)); }
            hist.Write();
        };
        save(d,"data",1); save(m,"mc_base",1); save(m,"mc_centrality",factors[c]); save(fitted,"mc_grefFit",factors[c]);
        // Cache normalization sumw2 and exact yield preservation separately from diagnostic cuts.
        TH1D totals(Form("gRefNorm_totals_code%d",c),";Sample;Sum of weights",3,0,3); totals.SetDirectory(nullptr); totals.Sumw2();
        totals.GetXaxis()->SetBinLabel(1,"Data"); totals.GetXaxis()->SetBinLabel(2,"MC centrality"); totals.GetXaxis()->SetBinLabel(3,"MC centrality+gRef");
        totals.SetBinContent(1,d.norm[c].s); totals.SetBinContent(2,factors[c]*m.norm[c].s); totals.SetBinContent(3,factors[c]*fitted.norm[c].s);
        totals.SetBinError(1,sqrt(d.norm[c].q0)); totals.SetBinError(2,abs(factors[c])*sqrt(m.norm[c].q0)); totals.SetBinError(3,abs(factors[c])*sqrt(fitted.norm[c].q0)); totals.Write();
    }
}
TGraphErrors* grefDensity(const Sample& s,int c,const vector<double>& edges,const TString& name) {
    auto* gr=new TGraphErrors; gr->SetName(name);
    const double total=s.norm[c].s, qtotal=s.norm[c].q0;
    if(!(total>0)) return gr;
    for(size_t i=0;i+1<edges.size();++i) {
        Moment bin=rangeSum(s.normG[c],int(edges[i]),int(edges[i+1]));
        double p=bin.s/total, width=edges[i+1]-edges[i];
        const double variance=(bin.q0*(1-2*p)+p*p*qtotal)/(total*total);
        int n=gr->GetN(); gr->SetPoint(n,(edges[i]+edges[i+1])/2.,p/width);
        gr->SetPointError(n,width/2.,sqrt(max(0.,variance))/width);
    }
    return gr;
}
TGraphErrors* densityRatio(TGraphErrors* d,TGraphErrors* m,const TString& name) {
    auto* gr=new TGraphErrors; gr->SetName(name);
    if(d->GetN()!=m->GetN()) return gr;
    for(int i=0;i<d->GetN();++i) {
        double x,y,xx,z; d->GetPoint(i,x,y); m->GetPoint(i,xx,z); if(!(z>0)) continue;
        int n=gr->GetN(); gr->SetPoint(n,x,y/z);
        gr->SetPointError(n,d->GetErrorX(i),hypot(d->GetErrorY(i)/z,y*m->GetErrorY(i)/(z*z)));
    }
    return gr;
}
void fitPages(TCanvas& cv,const TString& pdf,TFile& out,const Sample& d,const Sample& m,const Sample& fitted,
              const GrefFits& fits,const array<double,NC>& factors,const TString& subtitle) {
    (void) factors;
    for(int c=8;c>=0;--c) {
        const auto& f=fits[c];
        vector<double> edges=f.edges;
        if(edges.size()<2) edges={0,1};
        auto* gd=grefDensity(d,c,edges,Form("gRefFitQA_data_code%d",c));
        auto* gm=grefDensity(m,c,edges,Form("gRefFitQA_mc_centrality_code%d",c));
        auto* gf=grefDensity(fitted,c,edges,Form("gRefFitQA_mc_grefFit_code%d",c));
        style(gd,kBlack,20); style(gm,kBlue+1,25); style(gf,kGreen+2,26);
        auto* before=densityRatio(gd,gm,Form("gRefFitQA_ratio_before_code%d",c));
        auto* after=densityRatio(gd,gf,Form("gRefFitQA_ratio_after_code%d",c));
        style(before,kBlue+1,25); style(after,kGreen+2,26);
        auto* weights=new TGraphErrors; weights->SetName(Form("gRefExtraWeight_graph_code%d",c));
        for(int i=0;i<=200;++i) { double g=edges.front()+(edges.back()-edges.front())*i/200.; weights->SetPoint(i,g,f.eval(g)); weights->SetPointError(i,0,0); }
        style(weights,kGreen+2,1);
        out.cd(); for(auto* gr:{gd,gm,gf,before,after,weights}) gr->Write();
        grid(cv); auto* pads=static_cast<TPad*>(cv.GetPrimitive("ri_panels"));
        pads->cd(1); panel({gd,gm,gf},{"D^{0} data","Embedding: centrality","Embedding: centrality + gRefMult"},"gRefMult","Normalized density",edges.front(),edges.back());
        pads->cd(2); panel({before,after},{"Data / embedding: centrality","Data / embedding: centrality + gRefMult"},"gRefMult","Data / embedding",edges.front(),edges.back(),1);
        pads->cd(3); panel({weights},{"Extra weight (centrality yield preserved)"},"gRefMult","Extra gRefMult weight",edges.front(),edges.back(),1); weights->Draw("L SAME");
        // Keep fit annotations above the panels, clear of the weight curve.
        TString decision=f.reason.c_str(); decision.ReplaceAll("_"," ");
        header(cv,Form("gRefMult fit: %.0f-%.0f%% (centrality code %d)",centEdges[8-c],centEdges[9-c],c),subtitle);
        cv.cd(); TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(.021);
        t.DrawLatex(.035,.886,Form("degree %d; N_{eff}=%.1f; fit #chi^{2}/ndf = %.1f/%d; %s",f.degree,f.neff,f.chi2,f.ndf,decision.Data()));
        cv.cd(); TLatex foot; foot.SetNDC(); foot.SetTextFont(42); foot.SetTextSize(.017);
        foot.DrawLatex(.035,.015,"Fit QA uses the training sample. Display errors condition on fitted weights; parameter covariance saved in ROOT. #rho was not fitted.");
        cv.Print(pdf); cv.Clear(); for(auto* gr:{gd,gm,gf,before,after,weights}) delete gr;
    }
}

} // namespace RhoInvestDetail

void RhoInvest(double d0Min=-1, double d0Max=-1,
               double jetMin=-1.e9, double jetMax=1.e9,
               const char* sWeightFile="", const char* outputStem="",
               int gRefBins=120, double gRefMax=600,
               int rhoBins=200, double rhoMax=100, int fitDegree=2) {
    using namespace RhoInvestDetail;
    try {
        require(!gSystem->AccessPathName("config.h"),"Run from the project directory containing config.h");
        if(d0Min<0) d0Min=minPtD0Cut; if(d0Max<0) d0Max=maxPtD0Cut;
        require(d0Min<d0Max && jetMin<jetMax,"Invalid pT interval");
        require(gRefBins>0 && rhoBins>0 && gRefMax>0 && rhoMax>0,"Invalid plot binning");
        require(Method=="ICS" || Method=="AREA","Method must be ICS or AREA");
        require(fitDegree>=0 && fitDegree<=2,"fitDegree must be 0, 1 or 2");
        const TString dataFile=TString(sWeightFile).Length()?TString(sWeightFile):"./Output/"+outRoot+"_"+Method+_sys+"2.root";
        const TString mcFile=McJetsFileData;
        const bool jetCut=jetMin>-1.e8 || jetMax<1.e8;
        TString stem=TString(outputStem).Length()?TString(outputStem):Form("RhoInvest/output/RhoInvest_%s_D0_%g_%g",Method.Data(),d0Min,d0Max);
        if(!TString(outputStem).Length()) stem+=Form("_gRefFit%d",fitDegree);
        if(!TString(outputStem).Length() && jetCut) stem+=Form("_Jet_%g_%g",jetMin,jetMax);
        require(gSystem->mkdir(gSystem->DirName(stem),true)==0 || !gSystem->AccessPathName(gSystem->DirName(stem)),"Cannot create output directory");
        cout << "Data: " << dataFile << "\nMC: " << mcFile << "\nD0: [" << d0Min << ", " << d0Max << ") GeV/c\n";
        cout << "MC base = weightCentrality; MC weighted = base * Ndata(c)/Nmc(c). No FONLL weighting.\n";
        Sample data(gRefBins,gRefMax,rhoBins,rhoMax), mc(gRefBins,gRefMax,rhoBins,rhoMax);
        readData(data,dataFile,d0Min,d0Max,jetMin,jetMax); readMC(mc,mcFile,d0Min,d0Max,jetMin,jetMax);
        require(data.accepted>0 && mc.accepted>0,"No selected diagnostic entries in data or MC");
        array<double,NC> factors;
        for(int c=0;c<NC;++c) {
            require(data.norm[c].s>=0,"Negative net data yield in centrality code "+to_string(c)+"; cannot derive a physical centrality weight");
            if(data.norm[c].q0==0 && mc.norm[c].q0==0) { factors[c]=0; continue; }
            require(mc.norm[c].s>0,"Non-positive MC denominator in centrality code "+to_string(c));
            factors[c]=data.norm[c].s/mc.norm[c].s;
            cout << "centrality code " << c << ": data=" << data.norm[c].s << ", mc=" << mc.norm[c].s << ", factor=" << factors[c] << endl;
        }
        GrefFits fits;
        for(int c=0;c<NC;++c) {
            fits[c]=fitGref(data,mc,c,fitDegree);
            cout << "gRefMult fit code " << c << ": " << fits[c].reason
                 << ", degree=" << fits[c].degree << ", Neff=" << fits[c].neff
                 << ", chi2/ndf=" << fits[c].chi2 << "/" << fits[c].ndf << endl;
        }
        Sample fitted(gRefBins,gRefMax,rhoBins,rhoMax);
        // Continuous event-level weights. This second MC pass avoids assuming
        // independence of rho and gRefMult, or quantizing final weights into cells.
        readMC(fitted,mcFile,d0Min,d0Max,jetMin,jetMax,&fits);
        for(int c=0;c<NC;++c) if(mc.norm[c].s>0) {
            fits[c].normalization=fitted.norm[c].s/mc.norm[c].s;
            require(fits[c].normalization>0,"Invalid fitted-weight normalization");
            fitted.scaleClass(c,1./fits[c].normalization);
            require(abs(fitted.norm[c].s-mc.norm[c].s)<1.e-9*mc.norm[c].s,"Centrality yield preservation failed");
            fits[c].minWeight=numeric_limits<double>::infinity(); fits[c].maxWeight=0;
            for(int i=0;i<=1000;++i) {
                double w=fits[c].eval(fits[c].center+fits[c].scale*(2.*i/1000.-1.));
                fits[c].minWeight=min(fits[c].minWeight,w); fits[c].maxWeight=max(fits[c].maxWeight,w);
            }
        }
        TH1::AddDirectory(false); gStyle->SetOptStat(0); gStyle->SetTextFont(42);
        TFile out(stem+".root","RECREATE"); require(!out.IsZombie(),"Cannot create output ROOT file");
        ofstream table((stem+"_centrality.tsv").Data()); require(bool(table),"Cannot create TSV");
        table << setprecision(17) << "code\tcent_min\tcent_max\tNdata\tNmc\tfactor\tselected_data\tselected_mc_base\tselected_mc_weighted\n";
        TH1D hd("centrality_norm_data",";Centrality (%);Weighted entries",9,centEdges), hm("centrality_norm_mc_base",";Centrality (%);Weighted entries",9,centEdges), hw("centrality_factor",";Centrality (%);N_{data}/N_{MC}",9,centEdges), ha("centrality_norm_mc_weighted",";Centrality (%);Weighted entries",9,centEdges);
        for(int c=0;c<NC;++c) {
            int b=9-c; table << c << '\t' << centEdges[b-1] << '\t' << centEdges[b] << '\t' << data.norm[c].s << '\t' << mc.norm[c].s << '\t' << factors[c] << '\t' << data.selected[c].s << '\t' << mc.selected[c].s << '\t' << factors[c]*mc.selected[c].s << '\n';
            hd.SetBinContent(b,data.norm[c].s); hd.SetBinError(b,sqrt(data.norm[c].q0));
            hm.SetBinContent(b,mc.norm[c].s); hm.SetBinError(b,sqrt(mc.norm[c].q0));
            hw.SetBinContent(b,factors[c]); hw.SetBinError(b,0);
            ha.SetBinContent(b,factors[c]*mc.norm[c].s); ha.SetBinError(b,abs(factors[c])*sqrt(mc.norm[c].q0));
        }
        hd.Write(); hm.Write(); hw.Write(); ha.Write();
        writeFits(out,stem,data,mc,fitted,fits,factors);
        TString meta=Form("Source: Machine.C commit 9541f387d53d150b063875c9ee3a01fa646632da\nData=%s\nMC=%s\nMethod=%s\nD0 interval=[%g,%g)\nJet interval=[%g,%g)\nDividedMcDataBy=%d\ndeleteOneConstituentJets=%d\n",dataFile.Data(),mcFile.Data(),Method.Data(),d0Min,d0Max,jetMin,jetMax,int(DividedMcDataBy),int(deleteOneConstituentJets));
        meta+="Data weight=sWeightSignal*weightD0Efficiency*weightCentrality*weightDoubleCount\nMC base=weightCentrality; MC weighted=base*Ndata(code)/Nmc(code). No FONLL.\nNormalization samples follow Machine.LoadDataRealParams. Diagnostics add reco |eta|<=0.6, Nconst>=1, optional jet pT, and optional one-constituent rejection. Data inherits upstream sPlot selections (including stored neutral-fraction cut in reference source); MC inherits Machine truth D0/jet preselection. This is candidate/jet-weighted, not a unique-event MB average.\nSigned weights retained. Means use sum(w*rho)/sum(w); errors use sum(w^2*(rho-mean)^2)/sum(w)^2. Ratios assume independent data/MC at fixed weights. Centrality-map and sPlot-fit uncertainty and repeated-event correlations are not propagated. Non-positive net-weight mean bins omitted.\nRho distributions include flow bins in normalization; means include all finite rho. gRefMult plot excludes flow bins, centrality averages include them. No occupancy recalculation.\n";
        meta+=Form("gRefMult fit requested degree=%d. Model=MC template * exp(a*x+b*x*x), with x clamped to [-1,1] using MC 1%% and 99%% quantiles. Fits use weighted Gaussian chi2 with fixed variance=sumw2_data+(Ndata/Nmc)^2*sumw2_MC. Negative data bins retained; empty-data cells merged. MC template integrated through unit-width cells at their MC-weighted mean gRefMult. Final weights evaluated continuously in a second MC pass, with exact per-class normalization before diagnostic cuts. Fit errors conditional on input sWeights; diagnostics do NOT propagate fit/map uncertainty. No rho information is used in the fit. In-sample fit QA is not independent validation. Requested quadratic needs Neff>=150; linear needs Neff>=40; failed/boundary fits fall back to lower degree, ultimately unity. See _grefFit.tsv for every decision.\n",fitDegree);
        meta+=Form("Selected data=%lld (negative weights=%lld); MC=%lld; skipped invalid MC centrality=%lld\n",data.accepted,data.negative,mc.accepted,mc.badCent);
        TNamed metadata("RhoInvest_metadata",meta.Data()); metadata.Write();
        ofstream notes((stem+"_settings.txt").Data()); require(bool(notes),"Cannot create settings file"); notes << meta.Data();
        TString subtitle=Form("%s; %g #leq p_{T,D^{0}} < %g GeV/c; |#eta_{jet}^{reco}| #leq 0.6",Method.Data(),d0Min,d0Max);
        if(jetCut) subtitle+=Form("; %g #leq p_{T,jet}^{reco} < %g GeV/c",jetMin,jetMax);
        else subtitle+="; no extra jet-p_{T} cut";
        TCanvas cv("ri_canvas","RhoInvest",1500,600); const TString pdf=stem+".pdf"; cv.Print(pdf+"[");
        for(int mode=0;mode<3;++mode) profilePage(cv,pdf,out,data,mc,fitted,factors,mode,subtitle);
        shapePage(cv,pdf,out,data,mc,fitted,factors,subtitle);
        fitPages(cv,pdf,out,data,mc,fitted,fits,factors,subtitle); cv.Print(pdf+"]");
        out.Close(); cout << "Saved " << pdf << ", " << stem << ".root, _centrality.tsv, _settings.txt\n";
    } catch(const exception& e) {
        cerr << "RhoInvest ERROR: " << e.what() << endl;
        if(gROOT->IsBatch()) gSystem->Exit(1);
    }
}
