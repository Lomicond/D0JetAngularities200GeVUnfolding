// Compare two REAL StPicoD0JetAnaMaker productions, without changing either file.
// Put next to RhoInvest.C; run from the project directory (next to config.h):
// root -l -b -q 'RhoInvest/CompareOccupancyReal.C()'
// Default inputs are resolved in ./, then ./Data/ (or pass explicit paths).
// Optional arguments: oldFile, newFile, d0Min, d0Max, jetMin, jetMax, outputStem.
// No config.h / RooUnfold dependency. D0 candidates include combinatorial background:
// weightCentrality ONLY; no sWeights, D0 efficiency, FONLL, or fitted gRef weights.
// rho before the jet-pT cut avoids conditioning exclusively on migrated jets.
// "Equal class mix" means 25% of each NORMALIZED 10%-wide candidate distribution.
// Per-event plots use hCentralityW, filled before D0 selection in the maker.
// Histogram flow bins enter normalization and are reported in the summary TSV.
// Errors are indicative independent-entry errors, with fixed event denominators.
// Repeated candidates/events, cross-production overlap, and denominator uncertainty
// are NOT propagated; ratio error bars must not be used as a significance test.
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TBranch.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TNamed.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TClass.h>
#include <array>
#include <vector>
#include <string>
#include <memory>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace OccupancyReal {
struct Obs { const char *name,*leaf,*title; int bins; double lo,hi; bool ptCut; };
struct Moment {
 double s=0,s2=0,sx=0,s2x=0,s2xx=0; Long64_t n=0;
 void fill(double x,double w){s+=w;s2+=w*w;sx+=w*x;s2x+=w*w*x;s2xx+=w*w*x*x;++n;}
 double mean()const{return s>0?sx/s:0;}
 double error()const{double m=mean();return s>0?std::sqrt(std::max(0.,s2xx-2*m*s2x+m*m*s2))/s:0;}
};
struct Sample {
 std::unique_ptr<TFile> file; TTree* tree=nullptr;
 std::vector<std::array<std::unique_ptr<TH1D>,4> > hist;
 std::vector<std::array<Moment,4> > moments;
 std::array<double,4> events{{0,0,0,0}};
 bool hasEvents=false; Long64_t selected=0;
};
std::string resolve(const char* p){
 if(!gSystem->AccessPathName(p))return p;
 std::string q=std::string("Data/")+p;
 if(!gSystem->AccessPathName(q.c_str()))return q;
 throw std::runtime_error(std::string("Input not found: ")+p+" (also tried "+q+")");
}
TLeaf* leaf(TTree* t,const char* name,bool required=true){
 TLeaf* l=t->GetLeaf(name);
 if(!l){TBranch* b=t->GetBranch(name);if(b&&b->GetListOfLeaves()->GetEntries()==1)l=(TLeaf*)b->GetListOfLeaves()->At(0);}
 if(!l&&required)throw std::runtime_error(std::string("Missing scalar branch: ")+name);
 if(l)t->SetBranchStatus(l->GetBranch()->GetName(),1);
 return l;
}
// Locate a unique event-counter histogram, including maker subdirectories.
void findCounters(TDirectory* d,std::vector<TH1*>& out){
 TIter next(d->GetListOfKeys());TKey* k;
 while((k=(TKey*)next())){
  if(std::string(k->GetName())=="hCentralityW"){
   TH1* h=dynamic_cast<TH1*>(d->Get(k->GetName()));
   if(h&&h->GetDimension()==1&&std::find(out.begin(),out.end(),h)==out.end())out.push_back(h);
  }else if(TClass::GetClass(k->GetClassName())&&TClass::GetClass(k->GetClassName())->InheritsFrom(TDirectory::Class())){
   findCounters((TDirectory*)d->Get(k->GetName()),out);
  }
 }
}
void open(Sample& s,const char* name){
 std::string path=resolve(name);s.file.reset(TFile::Open(path.c_str(),"READ"));
 if(!s.file||s.file->IsZombie())throw std::runtime_error("Cannot open "+path);
 s.tree=dynamic_cast<TTree*>(s.file->Get("jets"));
 if(!s.tree)throw std::runtime_error("Missing jets tree in "+path);
 std::vector<TH1*> counters;findCounters(s.file.get(),counters);
 if(counters.size()==1){
  TH1* h=counters[0];s.hasEvents=true;
  for(int c=0;c<4;++c){int b=h->FindBin(c);s.events[c]=h->GetBinContent(b);
   if(std::abs(h->GetBinCenter(b)-c)>1.e-6||!(s.events[c]>0))s.hasEvents=false;
  }
 }
 std::cout<<path<<": "<<s.tree->GetEntries()<<" tree entries; usable hCentralityW="<<s.hasEvents<<std::endl;
}
void scan(Sample& s,int index,const std::vector<Obs>& obs,double dmin,double dmax,double jmin,double jmax){
 TTree* t=s.tree;t->SetBranchStatus("*",0);
 TLeaf *cent=leaf(t,"centrality"),*w=leaf(t,"weightCentrality"),*dp=leaf(t,"d0Pt"),
       *dy=leaf(t,"d0Rapidity"),*jp=leaf(t,"ICS_jetPt"),*eta=leaf(t,"ICS_jetEta"),*nc=leaf(t,"ICS_nJetConst");
 std::vector<TLeaf*> ls;for(const auto& o:obs)ls.push_back(leaf(t,o.leaf));
 s.hist.resize(obs.size());s.moments.resize(obs.size());
 for(size_t v=0;v<obs.size();++v)for(int c=0;c<4;++c){const auto& o=obs[v];
  s.hist[v][c].reset(new TH1D(Form("raw_%d_%s_c%d",index,o.name,c),o.title,o.bins,o.lo,o.hi));
  s.hist[v][c]->SetDirectory(nullptr);s.hist[v][c]->Sumw2();
 }
 for(Long64_t i=0;i<t->GetEntries();++i){
  if(i%1000000==0)std::cout<<"Sample "<<index<<": "<<i<<"/"<<t->GetEntries()<<std::endl;
  if(t->GetEntry(i)<=0)throw std::runtime_error("Failed to read tree entry");
  int c=int(std::lround(cent->GetValue()));if(c<0||c>3)continue;
  double ww=w->GetValue(),d=dp->GetValue(),y=dy->GetValue(),p=jp->GetValue(),e=eta->GetValue();
  if(!std::isfinite(ww)||ww<0)throw std::runtime_error("Nonfinite/negative weightCentrality");
  if(!std::isfinite(d)||!std::isfinite(y)||!std::isfinite(e)||!std::isfinite(p))continue;
  if(d<dmin||d>=dmax||std::abs(y)>=.6||std::abs(e)>=.6||nc->GetValue()<1)continue;
  ++s.selected;
  for(size_t v=0;v<obs.size();++v){if(obs[v].ptCut&&(p<jmin||p>=jmax))continue;
   double x=ls[v]->GetValue();if(!std::isfinite(x))continue;
   s.hist[v][c]->Fill(x,ww);s.moments[v][c].fill(x,ww);
  }
 }
}
std::unique_ptr<TH1D> clone(const TH1D* h,const TString& name){
 auto* p=(TH1D*)h->Clone(name);p->SetDirectory(nullptr);return std::unique_ptr<TH1D>(p);
}
// Propagate the correlation between a bin and its own total for shape errors.
bool normalize(TH1D& h){
 double sum=h.Integral(0,h.GetNbinsX()+1),q=0;
 for(int b=0;b<=h.GetNbinsX()+1;++b)q+=std::pow(h.GetBinError(b),2);
 if(!(sum>0))return false;
 for(int b=0;b<=h.GetNbinsX()+1;++b){double p=h.GetBinContent(b)/sum, qb=std::pow(h.GetBinError(b),2);
  double var=(qb*(1-2*p)+p*p*q)/(sum*sum);
  h.SetBinContent(b,p);h.SetBinError(b,std::sqrt(std::max(0.,var)));
 }return true;
}
// mode 0: shape; mode 1: per weighted event. panel 0 broad, 1 equal mix, 2..5 fine.
std::unique_ptr<TH1D> distribution(const Sample& s,size_t v,int panel,int mode,const TString& name){
 auto h=clone(s.hist[v][0].get(),name);h->Reset();
 if(panel>=2){int c=5-panel;h->Add(s.hist[v][c].get());
  if(mode)h->Scale(1./s.events[c]);else normalize(*h);
 }else if(panel==0){double n=0;for(int c=0;c<4;++c){h->Add(s.hist[v][c].get());n+=s.events[c];}
  if(mode)h->Scale(1./n);else normalize(*h);
 }else{
  for(int c=0;c<4;++c){auto part=clone(s.hist[v][c].get(),name+Form("_c%d",c));
   if(mode)part->Scale(1./s.events[c]);else if(!normalize(*part))throw std::runtime_error("Empty fine class: cannot build equal-class shape mix");
   h->Add(part.get(),.25);
  }
 }
 h->Scale(1.,"width");return h;
}
void drawPanel(TPad* outer,TH1D& old,TH1D& now,const char* title,const char* xTitle,bool yield){
 outer->cd();auto* up=new TPad(Form("up_%s",old.GetName()),"",0,.27,1,1);up->SetBottomMargin(.02);up->SetLeftMargin(.17);up->Draw();
 auto* down=new TPad(Form("down_%s",old.GetName()),"",0,0,1,.27);down->SetTopMargin(.03);down->SetBottomMargin(.35);down->SetLeftMargin(.17);down->Draw();
 up->cd();old.SetTitle(title);old.SetLineColor(kBlue+1);old.SetMarkerColor(kBlue+1);old.SetMarkerStyle(24);
 now.SetLineColor(kRed+1);now.SetMarkerColor(kRed+1);now.SetMarkerStyle(20);old.SetMarkerSize(.45);now.SetMarkerSize(.45);
 old.GetXaxis()->SetLabelSize(0);old.GetXaxis()->SetTitle("");old.GetYaxis()->SetTitle(yield?"(1/N_{evt}^{w}) dN_{cand}^{w}/dx":"Normalized density");
 old.GetYaxis()->SetTitleSize(.043);old.GetYaxis()->SetTitleOffset(1.65);
 double ymax=0;for(int b=1;b<=old.GetNbinsX();++b)ymax=std::max(ymax,std::max(old.GetBinContent(b)+old.GetBinError(b),now.GetBinContent(b)+now.GetBinError(b)));
 old.SetMinimum(0);old.SetMaximum(ymax>0?1.3*ymax:1);TH1* drawnOld=old.DrawCopy("E1");TH1* drawnNow=now.DrawCopy("E1 SAME");
 TLegend leg(.48,.75,.94,.90);leg.SetBorderSize(0);leg.SetFillStyle(0);leg.SetTextSize(.038);
 leg.AddEntry(drawnOld,"With occupancy C","lep");leg.AddEntry(drawnNow,"Without occupancy C","lep");leg.DrawClone();
 down->cd();TH1D frame("ratioframe","",old.GetNbinsX(),old.GetXaxis()->GetXmin(),old.GetXaxis()->GetXmax());frame.SetDirectory(nullptr);
 frame.GetXaxis()->SetTitle(xTitle);frame.GetYaxis()->SetTitle("No C / C");frame.GetXaxis()->SetTitleSize(.11);frame.GetXaxis()->SetLabelSize(.10);
 frame.GetYaxis()->SetTitleSize(.10);frame.GetYaxis()->SetLabelSize(.09);frame.GetYaxis()->SetTitleOffset(.65);frame.GetYaxis()->SetNdivisions(503);
 TGraphErrors gr;double min=.8,max=1.2;
 for(int b=1;b<=old.GetNbinsX();++b){double a=old.GetBinContent(b),z=now.GetBinContent(b);if(a<=0)continue;
  double ratio=z/a,err=std::hypot(now.GetBinError(b)/a,z*old.GetBinError(b)/(a*a));int n=gr.GetN();
  gr.SetPoint(n,old.GetBinCenter(b),ratio);gr.SetPointError(n,0,err);
  if(a>2*old.GetBinError(b)&&z>2*now.GetBinError(b)){min=std::min(min,ratio-err);max=std::max(max,ratio+err);}
 }
 frame.SetMinimum(std::max(0.,min-.05));frame.SetMaximum(std::max(1.25,max+.05));frame.DrawCopy();
 TLine line(old.GetXaxis()->GetXmin(),1,old.GetXaxis()->GetXmax(),1);line.SetLineStyle(2);line.DrawClone();gr.SetMarkerStyle(20);gr.SetMarkerSize(.4);gr.DrawClone("P SAME");
}
}

void CompareOccupancyReal(const char* oldFile="Output_new_dca1_13102026.root",
                          const char* newFile="Output_noOcu_real_22112026.root",
                          double d0Min=1,double d0Max=10,double jetMin=5,double jetMax=20,
                          const char* outputStem=""){
 using namespace OccupancyReal;
 try{
  if(!(d0Min<d0Max&&jetMin<jetMax))throw std::runtime_error("Invalid pT interval");
  Sample old,now;open(old,oldFile);open(now,newFile);
  std::vector<Obs> obs={
   {"rho_prePt","backgroundDensity","#rho [GeV/c per unit area]",100,0,25,false},
   {"rho_selected","backgroundDensity","#rho [GeV/c per unit area]",100,0,25,true},
   {"jetPt","ICS_jetPt","p_{T,jet}^{ICS} [GeV/c]",60,0,30,false},
   {"LHA","ICS_lambda1_0_5","#lambda_{0.5}^{1}",50,0,1,true},
   {"girth","ICS_lambda1_1","#lambda_{1}^{1}",50,0,.5,true},
   {"gRefMult","gRefMult","gRefMult",80,0,400,false}
  };
  if(leaf(old.tree,"backgroundDensityM",false)&&leaf(now.tree,"backgroundDensityM",false))
   obs.push_back({"rhoM_selected","backgroundDensityM","#rho_{m} [GeV/c per unit area]",60,0,3,true});
  if(leaf(old.tree,"ICS_lambda1_0",false)&&leaf(now.tree,"ICS_lambda1_0",false))
   obs.push_back({"lambda0","ICS_lambda1_0","#lambda_{0}^{1}",60,0,1.2,true});
  scan(old,0,obs,d0Min,d0Max,jetMin,jetMax);scan(now,1,obs,d0Min,d0Max,jetMin,jetMax);
  if(!old.selected||!now.selected)throw std::runtime_error("No selected peripheral entries");
  TString stem=outputStem;if(stem.IsNull())stem=Form("RhoInvest/output/OccupancyReal_D0_%g_%g_Jet_%g_%g",d0Min,d0Max,jetMin,jetMax);
  TString dir=gSystem->DirName(stem);if(gSystem->AccessPathName(dir)&&gSystem->mkdir(dir,true)!=0)throw std::runtime_error("Cannot create output directory");
  TFile out(stem+".root","RECREATE");if(out.IsZombie())throw std::runtime_error("Cannot create ROOT output");
  std::ofstream table((stem+"_summary.tsv").Data());if(!table)throw std::runtime_error("Cannot create summary");
  table<<std::setprecision(12)<<"sample\tobservable\tcent_min\tcent_max\tentries\tsum_w\tweighted_events\tmean\tmean_error_indicative\tunderflow_w\toverflow_w\n";
  for(int k=0;k<2;++k){const Sample& s=k?now:old;for(size_t v=0;v<obs.size();++v)for(int c=3;c>=0;--c){
   const auto& m=s.moments[v][c];TH1D* h=s.hist[v][c].get();out.cd();h->Write();
   table<<(k?"noC":"withC")<<'\t'<<obs[v].name<<'\t'<<70-10*c<<'\t'<<80-10*c<<'\t'<<m.n<<'\t'<<m.s<<'\t';
   if(s.hasEvents)table<<s.events[c];else table<<"nan";
   table<<'\t'<<m.mean()<<'\t'<<m.error()<<'\t'<<h->GetBinContent(0)<<'\t'<<h->GetBinContent(h->GetNbinsX()+1)<<'\n';
  }}
  std::string settings="OLD="+std::string(old.file->GetName())+"\nNEW="+now.file->GetName()+
   Form("\nD0 pT=[%g,%g), |yD0|<0.6; ICS |eta|<0.6, Nconst>=1; selected jet pT=[%g,%g).",d0Min,d0Max,jetMin,jetMax)+
   "\nRaw candidate QA; weightCentrality only. No sWeights or D0 efficiency corrections. No extra mass cut.\nShapes include flows in total normalization. Equal class mix: 25% per fine class.\nPer-event denominator: hCentralityW (all accepted events before D0 selection).\nErrors: independent entries/productions, fixed event denominators; repeated events and sample overlap ignored.\nOnly occupancy must differ for attribution; compare cuts, input run coverage and code versions externally.\n";
  out.cd();TNamed provenance("settings",settings.c_str());provenance.Write();std::ofstream info((stem+"_settings.txt").Data());info<<settings;
  const bool perEvent=old.hasEvents&&now.hasEvents;
  if(!perEvent)std::cerr<<"WARNING: missing/ambiguous/invalid hCentralityW: per-event yield page omitted; shapes still available."<<std::endl;
  gStyle->SetOptStat(0);gStyle->SetTextFont(42);TCanvas cv("occupancyCompare","Occupancy comparison",1350,900);TString pdf=stem+".pdf";cv.Print(pdf+"[");
  const char* labels[]={"40-80%: original class mixture","40-80%: equal class mixture","40-50%","50-60%","60-70%","70-80%"};
  for(size_t v=0;v<obs.size();++v){int modes=(perEvent&&std::string(obs[v].name)=="jetPt")?2:1;
   for(int mode=0;mode<modes;++mode){cv.Clear();cv.cd();
    TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(.024);
    text.DrawLatex(.025,.972,Form("%s | %s | D^{0} p_{T}: [%g,%g) GeV/c%s",obs[v].name,mode?"per weighted event":"shapes",d0Min,d0Max,obs[v].ptCut?Form(" | jet p_{T}: [%g,%g)",jetMin,jetMax):" | before jet-pT cut"));
    text.SetTextSize(.017);text.DrawLatex(.025,.018,"Raw D^{0} candidates, weightCentrality only; indicative errors; flow bins retained in normalization.");
    auto* grid=new TPad("panels","",0,.04,1,.94);grid->Draw();grid->Divide(3,2,.003,.003);
    for(int panel=0;panel<6;++panel){auto a=distribution(old,v,panel,mode,Form("withC_%s_%d_%d",obs[v].name,mode,panel));auto b=distribution(now,v,panel,mode,Form("noC_%s_%d_%d",obs[v].name,mode,panel));
     out.cd();a->Write();b->Write();drawPanel((TPad*)grid->cd(panel+1),*a,*b,labels[panel],obs[v].title,mode);
    }cv.Print(pdf);
   }
  }
  cv.Print(pdf+"]");out.Close();std::cout<<"Wrote "<<pdf<<", "<<stem<<".root, _summary.tsv, _settings.txt"<<std::endl;
 }catch(const std::exception& e){std::cerr<<"CompareOccupancyReal ERROR: "<<e.what()<<std::endl;throw;}
}
