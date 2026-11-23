#ifndef EXPORT_PRIOR_EXPLORER_H
#define EXPORT_PRIOR_EXPLORER_H
// Called in a temporary copy of Machine.C after LoadDataReal(), before unfolding.
// Writes ordinary histograms: the Python GUI does not need RooUnfold.
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TObjString.h"
#include "RooUnfoldResponse.h"
#include <stdexcept>
#include <algorithm>
#include "TDirectory.h"
#include <cmath>
#include <memory>

namespace PEExport {
inline int bins(const TH1* h) {
    return h->GetNbinsX()*(h->GetDimension()==2?h->GetNbinsY():1);
}
inline int global(const TH1* h,int k) {
    return h->GetDimension()==2?h->GetBin(k%h->GetNbinsX()+1,k/h->GetNbinsX()+1):k+1;
}
inline void sameAxis(const TAxis* a,const TAxis* b) {
    if(a->GetNbins()!=b->GetNbins())throw std::runtime_error("PriorExplorer: axis size mismatch");
    for(int i=1;i<=a->GetNbins()+1;++i)
        if(std::abs(a->GetBinLowEdge(i)-b->GetBinLowEdge(i))>1e-9)
            throw std::runtime_error("PriorExplorer: data/response bin edges differ");
}

// Categories include truth underflow, overflow, and unavailable truth values.
// Flat reco indexing always has jet pT as the fastest coordinate.
struct FakeMaps {
    std::unique_ptr<TH2D> pt, joint;
};
inline FakeMaps& maps(int c,int v,const RooUnfoldResponse& r) {
    static FakeMaps store[3][8];
    auto &m=store[c][v];
    if(!m.pt) {
        const auto *t=r.Htruth();int nx=t->GetNbinsX(),ny=t->GetDimension()==2?t->GetNbinsY():0;
        int nr=bins(r.Hmeasured());
        m.pt.reset(new TH2D(Form("pe_fake_pt_%d_%d",c,v),"Fake reco vs truth pT category",nr,0,nr,nx+3,0,nx+3));
        m.pt->SetDirectory(nullptr);m.pt->Sumw2();
        if(ny) {
            int nj=(nx+2)*(ny+2)+1;
            m.joint.reset(new TH2D(Form("pe_fake_joint_%d_%d",c,v),"Fake reco vs extended truth cell",nr,0,nr,nj,0,nj));
            m.joint->SetDirectory(nullptr);m.joint->Sumw2();
        }
    }
    return m;
}
inline void recordFake(int c,int v,const RooUnfoldResponse& r,
                       double rx,double ry,double tx,double ty,double weight) {
    const auto *t=r.Htruth(),*m=r.Hmeasured();
    if(!std::isfinite(weight)||weight<0)throw std::runtime_error("PriorExplorer: invalid fake event weight");
    int ix=m->GetXaxis()->FindFixBin(rx),iy=m->GetDimension()==2?m->GetYaxis()->FindFixBin(ry):1;
    if(ix<1||ix>m->GetNbinsX()||iy<1||iy>(m->GetDimension()==2?m->GetNbinsY():1))
        throw std::runtime_error("PriorExplorer: fake reco outside response axes");
    int i=(iy-1)*m->GetNbinsX()+ix-1,nx=t->GetNbinsX();
    bool validX=std::isfinite(tx)&&tx>0;
    int bx=validX?t->GetXaxis()->FindFixBin(tx):nx+2;
    auto &h=maps(c,v,r);h.pt->Fill(i+.5,bx+.5,weight);
    if(t->GetDimension()==2) {
        int ny=t->GetNbinsY();bool validY=std::isfinite(ty)&&ty>-900;
        int j=(nx+2)*(ny+2);
        if(validX&&validY)j=t->GetYaxis()->FindFixBin(ty)*(nx+2)+bx;
        h.joint->Fill(i+.5,j+.5,weight);
    }
}
inline void writeFakeMaps(int c,int v,const RooUnfoldResponse& r) {
    auto &h=maps(c,v,r);const auto *f=r.Hfakes();
    for(int i=0;i<bins(r.Hmeasured());++i) {
        double expected=f?f->GetBinContent(global(f,i)):0;
        for(const TH2D *a:{h.pt.get(),h.joint.get()}) {
            if(!a)continue;
            double sum=0;for(int j=1;j<=a->GetNbinsY();++j)sum+=a->GetBinContent(i+1,j);
            if(std::abs(sum-expected)>1e-7*std::max(1.,std::abs(expected)))
                throw std::runtime_error("PriorExplorer: fake truth map does not reproduce Hfakes; check instrumentation");
        }
    }
    h.pt->Write("fake_pt");if(h.joint)h.joint->Write("fake_joint");
}
inline void model(TDirectory* d,const RooUnfoldResponse& r,const TH1& data) {
    d->cd(); const TH1* t=r.Htruth(); const TH1* m=r.Hmeasured(); const TH1* f=r.Hfakes();
    const TH2* a=r.Hresponse();
    if(!t||!m||!a)throw std::runtime_error("PriorExplorer: uninitialised response");
    int nt=bins(t),nr=bins(m);
    if(a->GetNbinsX()!=nr||a->GetNbinsY()!=nt)
        throw std::runtime_error("PriorExplorer: overflow/unexpected flattening; require UseOverflow=false");
    if(data.GetDimension()!=m->GetDimension())throw std::runtime_error("PriorExplorer: data dimension mismatch");
    sameAxis(data.GetXaxis(),m->GetXaxis());
    if(m->GetDimension()==2)sameAxis(data.GetYaxis(),m->GetYaxis());
    TH1D truth("truth","Matched plus missed truth",nt,0,nt);
    TH1D reco("reco","Matched plus fake reco",nr,0,nr);
    TH1D fake("fakes","Nominal fake template",nr,0,nr);
    TH1D measured("data","ICS weighted measured counts",nr,0,nr);
    for(int j=0;j<nt;++j) {
        double v=t->GetBinContent(global(t,j));
        if(!std::isfinite(v)||v<0)throw std::runtime_error("PriorExplorer: invalid truth");
        truth.SetBinContent(j+1,v);truth.SetBinError(j+1,t->GetBinError(global(t,j)));
        double col=0;
        for(int i=0;i<nr;++i) {
            double x=a->GetBinContent(i+1,j+1);
            if(!std::isfinite(x)||x<0)throw std::runtime_error("PriorExplorer: invalid migration");
            col+=x;
        }
        if(col>v+1e-7*std::max(1.,v))throw std::runtime_error("PriorExplorer: matched > truth");
    }
    for(int i=0;i<nr;++i) {
        int k=global(m,i);double v=m->GetBinContent(k),fv=f?f->GetBinContent(global(f,i)):0;
        double pred=fv;for(int j=0;j<nt;++j)pred+=a->GetBinContent(i+1,j+1);
        if(!std::isfinite(v)||!std::isfinite(fv)||fv<0||std::abs(pred-v)>1e-7*std::max(1.,std::abs(v)))
            throw std::runtime_error("PriorExplorer: Hresponse projection + Hfakes != Hmeasured");
        reco.SetBinContent(i+1,v);fake.SetBinContent(i+1,fv);
        fake.SetBinError(i+1,f?f->GetBinError(global(f,i)):0);
        measured.SetBinContent(i+1,data.GetBinContent(global(&data,i)));
        measured.SetBinError(i+1,data.GetBinError(global(&data,i)));
    }
    a->Write("migration");t->Write("truth_axes");m->Write("reco_axes");
    truth.Write();reco.Write();fake.Write();measured.Write();
}
inline void save(const char* path,const char* metadata,RooUnfoldResponse* r1,
                 RooUnfoldResponse (*r2)[8],TH1D* d1,TH2D (*d2)[7]) {
    TDirectory::TContext context;
    TFile out(path,"RECREATE");
    if(out.IsZombie())throw std::runtime_error("PriorExplorer: cannot write export");
    TObjString info(metadata);info.Write("metadata");
    const char* names[8]={"jetPt","z","lambda1","lambda1p5","lambda2","lambda3","lambda0p5","pTD"};
    for(int c=0;c<3;++c)for(int v=0;v<8;++v) {
        out.cd();auto d=out.mkdir(Form("c%d_%s",c,names[v]));
        if(v==0) {model(d,r1[c],d1[c]);writeFakeMaps(c,v,r1[c]);}
        else {model(d,r2[c][v-1],d2[c][v-1]);writeFakeMaps(c,v,r2[c][v-1]);}
    }
    out.Close();
}
}
#endif
