#!/usr/bin/env python3
"""Interactive truth-prior / fixed-response forward folding for Machine.C.
Run --help. GUI needs numpy, matplotlib, tkinter; preparing needs ROOT/RooUnfold
and uproot. Bundles include the exact binned response and ICS data from Machine.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import queue
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from dataclasses import dataclass
import numpy as np

KEYS = ['jetPt','z','lambda1','lambda1p5','lambda2','lambda3','lambda0p5','pTD']
LABELS = ['Jet pT (1D)','z','lambda 1','lambda 1.5','lambda 2','lambda 3','lambda 0.5','pT dispersion']
AXLABEL = {'jetPt':r'$p_{T,jet}$ [GeV/c]','z':r'$z$','lambda1':r'$\lambda^1_1$',
 'lambda1p5':r'$\lambda^1_{1.5}$','lambda2':r'$\lambda^1_2$','lambda3':r'$\lambda^1_3$',
 'lambda0p5':r'$\lambda^1_{0.5}$','pTD':r'$p_T^D$'}
CENT = ['0-10%','10-40%','40-80%']

@dataclass
class Model:
    migration: np.ndarray
    truth: np.ndarray
    fakes: np.ndarray
    data: np.ndarray
    variance: np.ndarray
    tx: np.ndarray
    ty: np.ndarray
    rx: np.ndarray
    ry: np.ndarray
    reco: np.ndarray

    def __post_init__(self):
        for name in ('migration','truth','fakes','data','variance','tx','ty','rx','ry','reco'):
            setattr(self,name,np.asarray(getattr(self,name),dtype=float))
        for name in ('tx','ty','rx','ry'):
            a=getattr(self,name)
            if a.size and (a.size<2 or not np.all(np.isfinite(a)) or not np.all(np.diff(a)>0)):
                raise ValueError('Invalid '+name+' edges')
        nt=(len(self.tx)-1)*max(1,len(self.ty)-1); nr=(len(self.rx)-1)*max(1,len(self.ry)-1)
        if self.migration.shape!=(nr,nt) or self.truth.shape!=(nt,) or any(getattr(self,n).shape!=(nr,) for n in ('fakes','data','variance','reco')):
            raise ValueError('Incompatible matrix/vector/axis dimensions')
        for n in ('migration','truth','fakes','data','variance','reco'):
            a=getattr(self,n)
            if not np.all(np.isfinite(a)) or (n!='data' and np.any(a<0)):
                raise ValueError('Invalid values in '+n)
        if np.any(self.migration.sum(axis=0)>self.truth+1e-7*np.maximum(1,self.truth)):
            raise ValueError('Matched counts exceed truth (check misses/flattening)')
        if not np.allclose(self.migration.sum(axis=1)+self.fakes,self.reco,rtol=1e-7,atol=1e-8):
            raise ValueError('Migration + fakes does not reproduce nominal reco')
        if self.truth.sum()<=0 or self.reco.sum()<=0: raise ValueError('Empty response')
        self.response=np.divide(self.migration,self.truth[None,:],out=np.zeros_like(self.migration),where=self.truth[None,:]>0)
        self.response.setflags(write=False)
        for n in ('migration','truth','fakes','data','variance','reco'):
            getattr(self,n).setflags(write=False)

    @property
    def two_d(self): return len(self.ty)>0

    def coordinates(self):
        x=(self.tx[1:]+self.tx[:-1])/2
        y=(self.ty[1:]+self.ty[:-1])/2 if self.two_d else np.array([0.])
        # RooUnfold uses x (jet pT) as the fastest flat coordinate.
        return np.tile(x,len(y)),np.repeat(y,len(x))

    def vary(self,mode='MC tilt',a=0.,b=0.,c=0.,d=0.,p0=3.,n=5.):
        x,y=self.coordinates(); logx=np.log(np.maximum(x,1e-6)/8.)
        yn=(2*(y-self.ty[0])/(self.ty[-1]-self.ty[0])-1) if self.two_d else np.zeros_like(x)
        if mode=='Power law':
            if p0<=0 or n<=0: raise ValueError('p0 and n must be positive')
            # Integrate target density over each pT bin; retain MC conditional y|pT.
            nodes,weights=np.polynomial.legendre.leggauss(12)
            lo,hi=self.tx[:-1],self.tx[1:]
            pts=(lo[:,None]+hi[:,None])/2+(hi-lo)[:,None]*nodes/2
            target=((1+pts/p0)**(-n)*weights).sum(axis=1)*(hi-lo)/2
            marginal=self.truth.reshape((-1,len(lo))).sum(axis=0)
            wt=np.divide(target,marginal,out=np.zeros_like(target),where=marginal>0)
            base=np.tile(wt,max(1,len(self.ty)-1))
            exponent=c*yn+d*logx*yn
        else:
            base=np.ones_like(x);exponent=a*logx+b*logx**2+c*yn+d*logx*yn
        if not np.all(np.isfinite(exponent)) or np.max(np.abs(exponent))>500:
            raise ValueError('Weight exponent too large; reduce deformation')
        w=base*np.exp(exponent);new=self.truth*w
        if new.sum()<=0: raise ValueError('No supported truth bins')
        w*=self.truth.sum()/new.sum(); new=self.truth*w
        return new,w

    def predict(self,settings):
        truth,weights=self.vary(**{k:settings[k] for k in ('mode','a','b','c','d','p0','n')})
        fk=settings['fake']; raw=self.response@truth+fk*self.fakes
        total=self.data.sum()
        if total<=0: raise ValueError('Nonpositive total data yield; cannot use this normalization')
        scale=total/(raw.sum() if settings['normalize'] else self.reco.sum())*settings['amplitude']
        nominal_scale=total/self.reco.sum()
        return dict(truth=truth*scale,weights=weights,reco=raw*scale,
                    fake=self.fakes*fk*scale,nominal_truth=self.truth*nominal_scale,
                    nominal_reco=self.reco*nominal_scale,scale=scale)

    def project(self,values,level,axis):
        x,y=(self.tx,self.ty) if level=='truth' else (self.rx,self.ry)
        a=np.asarray(values).reshape((max(1,len(y)-1),len(x)-1))
        return a.sum(axis=0 if axis==0 else 1)

@dataclass
class Bundle:
    path: Path
    metadata: dict
    models: dict

    @property
    def label(self):
        return f"D0 {self.metadata['d0_min']:g}-{self.metadata['d0_max']:g} | {self.metadata.get('prior_order','?')} | {self.path.name}"


def save_bundle(path,metadata,models):
    arr={'metadata':np.array(json.dumps(metadata,ensure_ascii=False))}
    for key,m in models.items():
        for field in Model.__dataclass_fields__: arr[key+'__'+field]=getattr(m,field)
    np.savez_compressed(path,**arr)


def load_bundle(path):
    path=Path(path)
    with np.load(path,allow_pickle=False) as z:
        metadata=json.loads(str(z['metadata']));models={}
        if metadata.get('format')!=1: raise ValueError('Unsupported bundle format')
        for c in range(3):
            for v in KEYS:
                key=f'c{c}_{v}'
                if key+'__truth' in z:
                    models[key]=Model(**{n:z[key+'__'+n] for n in Model.__dataclass_fields__})
    if not models:raise ValueError('No response models found')
    return Bundle(path,metadata,models)


def convert_export(source,destination):
    import uproot
    with uproot.open(source) as f:
        meta=json.loads(str(f['metadata']));models={}
        for c in range(3):
            for v in KEYS:
                key=f'c{c}_{v}';g=f[key]; ta=g['truth_axes'];ra=g['reco_axes']
                var=g['data'].variances()
                if var is None: raise ValueError('Missing measured Sumw2')
                models[key]=Model(g['migration'].values(),g['truth'].values(),g['fakes'].values(),
                    g['data'].values(),var,ta.axes[0].edges(),ta.axes[1].edges() if len(ta.axes)==2 else [],
                    ra.axes[0].edges(),ra.axes[1].edges() if len(ra.axes)==2 else [],g['reco'].values())
    save_bundle(destination,meta,models)
    return load_bundle(destination)


def cpp(s): return json.dumps(str(s),ensure_ascii=True)


def patched_source(machine,header,metadata,export_path,before):
    source=machine.read_text();matches=list(re.finditer(r'\bvoid\s+(Machine2|Machine)\s*\(',source))
    if len(matches)!=1: raise ValueError('Expected one void Machine(...) or Machine2(...) entry point')
    main=matches[0]; name=main.group(1)
    # Require the known loading sequence. Refuse an unfamiliar version instead of guessing.
    tail=source[main.start():]
    pattern=r'(?m)^\s*if\s*\(\s*!ClosureTest\s*\)\s*LoadDataReal\s*\(\s*\)\s*;'
    hooks=list(re.finditer(pattern,tail))
    if len(hooks)!=1: raise ValueError('Cannot locate unique real-data load in Machine entry point; send the current macro for adaptation')
    hook=hooks[0];pos=main.start()+hook.end()
    injection='\n    PEExport::save('+cpp(export_path)+','+cpp(json.dumps(metadata))+',rurResponse,rurResponse2D,hRealData,hRealData2D);\n    return; // PriorExplorer: stop before unfolding and plotting final results\n'
    source=source[:pos]+injection+source[pos:]
    # Headers in the temporary copy still resolve relative to the original macro.
    def absolute_include(m):
        file=machine.parent/m.group(1)
        return '#include '+cpp(file.resolve()) if file.exists() else m.group(0)
    source=re.sub(r'#include\s+"([^"\n]+)"',absolute_include,source)
    source='#include '+cpp(header)+'\n'+source
    setup=['ClosureTest=false;','UseCachedRM=false;','FillStandardRM=true;','FillCacheRM=false;',
           'Unfold2D=true;','MissingJets=true;','FakeJets=true;','UseOverflow=false;','Method="ICS";']
    optional={'crossClosure':'false','ForwardFoldCrossClosureQA':'false','ResponseWeightComparisonQA':'false'}
    # Set flags at entry; the driver passes no optional binning override.
    body_start=source.index('{',source.index('void '+name+'('))
    setup_code='\n'+''.join('\n    '+s for s in setup)+'\n'
    for flag,value in optional.items():
        setup_code+=f'    if(gROOT->GetGlobal({cpp(flag)})) gROOT->ProcessLine({cpp(flag+"="+value+";")});\n'
    setup_code+='    if(!gROOT->GetGlobal("cutBefore")) throw std::runtime_error("PriorExplorer requires a Machine/config version with cutBefore");\n'
    setup_code+='    gROOT->ProcessLine('+cpp('cutBefore='+('true' if before else 'false')+';')+');\n'
    source=source[:body_start+1]+setup_code+source[body_start+1:]
    source+='\nvoid PriorExplorerBuild() {\n  '+name+'(true,false,-30,4,nullptr,"PriorExplorer",'+str(metadata['d0_min'])+','+str(metadata['d0_max'])+');\n}\n'
    return source


def prepare(machine,dlo,dhi,outdir,before=True,root_cmd='root',log=print):
    machine=Path(machine).expanduser().resolve();outdir=Path(outdir).expanduser().resolve()
    header=Path(__file__).with_name('ExportPriorExplorer.h')
    if not machine.is_file() or not header.is_file():raise ValueError('Missing Machine.C or ExportPriorExplorer.h')
    if not 1<=dlo<dhi<=10:raise ValueError('D0 range must be inside 1-10 GeV/c')
    repo=machine.parent.parent
    if not (repo/'config.h').is_file():raise ValueError('Place Machine.C in its usual Unfolding/ directory under config.h')
    outdir.mkdir(parents=True,exist_ok=True)
    job=Path(tempfile.mkdtemp(prefix='prepare_',dir=outdir))
    meta=dict(format=1,d0_min=float(dlo),d0_max=float(dhi),method='ICS',
              prior_order='cutBefore' if before else 'cutAfter',machine=str(machine),
              machine_sha256=hashlib.sha256(machine.read_bytes()).hexdigest(),
              config_sha256=hashlib.sha256((repo/'config.h').read_bytes()).hexdigest(),
              data='Output/Output_ICS2.root, loaded by LoadDataReal()',
              created=time.strftime('%Y-%m-%dT%H:%M:%S'),overflow=False,
              note='Exact Machine histograms before unfolding; all MC; closure off; no negative-bin clipping; no extra prior tilt')
    # Isolate conventional Machine diagnostic outputs from the user's results.
    for d in ('Output','OutputPdf','Tests','Test'): (job/d).mkdir()
    if (repo/'Data').exists():(job/'Data').symlink_to(repo/'Data',target_is_directory=True)
    for name in ('FONLL_D0Pt_0_31.root','FONLL_D0Pt_0_10.root','Output_ICS2.root','Output_D02.root'):
        src=repo/'Output'/name
        if not src.is_file():raise ValueError(f'Required Machine input missing: {src}')
        (job/'Output'/name).symlink_to(src)
    # Preserve local ROOT logon setup (e.g. RooUnfold loading) for the isolated job.
    if (repo/'.rootlogon.C').is_file():
        (job/'.rootlogon.C').write_text('{ gROOT->Macro('+cpp(repo/'.rootlogon.C')+'); }\n')
    export=job/'export.root';macro=job/'PriorExplorerBuild.C'
    macro.write_text(patched_source(machine,header.resolve(),meta,export,before))
    log('Preparing from '+str(machine)+'; log: '+str(job/'prepare.log'))
    env=os.environ.copy(); env['ROOT_INCLUDE_PATH']=str(repo)+':'+str(machine.parent)+':'+env.get('ROOT_INCLUDE_PATH','')
    with (job/'prepare.log').open('w') as stream:
        process=subprocess.Popen([root_cmd,'-l','-b','-q',str(macro)],cwd=job,env=env,stdout=stream,stderr=subprocess.STDOUT)
        code=process.wait()
    if code or not export.is_file():raise RuntimeError(f'ROOT preparation failed (exit {code}). Read {job / "prepare.log"}')
    stem=f'PriorExplorer_D0_{dlo:g}_{dhi:g}_{meta["prior_order"]}_{job.name[8:]}'
    destination=outdir/(stem+'.npz')
    result=convert_export(export,destination)
    log('Ready: '+str(destination))
    return result


def defaults():
    return dict(mode='MC tilt',a=0.,b=0.,c=0.,d=0.,p0=3.,n=5.,fake=1.,amplitude=1.,normalize=True)


def demo_bundle(path):
    """Synthetic detector + deliberately different truth; never analysis data."""
    rng=np.random.default_rng(125);models={}
    for c in range(3):
        for v in KEYS:
            tx=np.array([1,2,3,4,5,7,9,11,13,15,20.]);rx=np.array([-5,-2,0,1,2,3,4,5,6,7,8,10,13,16,20,25.])
            ty=np.array([0,.05,.1,.2,.4,.7,1.]) if v!='jetPt' else np.array([]);ry=ty.copy()
            x=np.tile((tx[1:]+tx[:-1])/2,max(1,len(ty)-1));y=np.repeat((ty[1:]+ty[:-1])/2,len(tx)-1) if len(ty) else np.zeros_like(x)
            xr=np.tile((rx[1:]+rx[:-1])/2,max(1,len(ry)-1));yr=np.repeat((ry[1:]+ry[:-1])/2,len(rx)-1) if len(ry) else np.zeros_like(xr)
            density=np.exp(-.5*((xr[:,None]-x[None,:])/(3-.5*c))**2)
            if len(ty):density*=np.exp(-.5*((yr[:,None]-y[None,:])/.13)**2)
            widths=np.tile(np.diff(rx),max(1,len(ry)-1))
            if len(ry):widths*=np.repeat(np.diff(ry),len(rx)-1)
            density*=widths[:,None];density/=density.sum(axis=0);density*=.75+.15*x[None,:]/20
            t=2e4*(1+x/3)**-4*np.exp(-4*y)
            f=5*np.exp(-xr**2/20);reco=density@t+f
            alternate=t*np.exp(.5*np.log(x/8)+.6*y);alternate*=t.sum()/alternate.sum()
            truth_reco=density@alternate+f;var=np.maximum(truth_reco,1.)
            data=truth_reco+rng.normal(0,np.sqrt(var))
            models[f'c{c}_{v}']=Model(density*t,t,f,data,var,tx,ty,rx,ry,reco)
    save_bundle(path,dict(format=1,d0_min=1.,d0_max=5.,prior_order='SYNTHETIC DEMO',method='ICS-like toy',demo=True),models)
    return load_bundle(path)


def draw(fig,model,key,settings,title):
    result=model.predict(settings);fig.clear();axes=fig.subplots(2,3)
    for row,(axis,variable) in enumerate(((0,'jetPt'),(1,key))):
        if row and not model.two_d:
            for ax in axes[row]:ax.axis('off')
            axes[row,1].text(.5,.5,'1D response\nSelect z or an angularity for a 2D view',ha='center',va='center',transform=axes[row,1].transAxes)
            continue
        tx=model.tx if axis==0 else model.ty;rx=model.rx if axis==0 else model.ry
        t=model.project(result['truth'],'truth',axis);t0=model.project(result['nominal_truth'],'truth',axis)
        rec=model.project(result['reco'],'reco',axis);r0=model.project(result['nominal_reco'],'reco',axis)
        fk=model.project(result['fake'],'reco',axis);data=model.project(model.data,'reco',axis)
        err=np.sqrt(model.project(model.variance,'reco',axis));xc=(rx[:-1]+rx[1:])/2
        at,ar,aq=axes[row]
        at.stairs(t0/np.diff(tx),tx,color='0.55',linestyle='--',label='Nominal MC')
        at.stairs(t/np.diff(tx),tx,color='#da6b16',linewidth=2,label='Trial prior')
        ar.stairs(r0/np.diff(rx),rx,color='0.55',linestyle='--',label='Nominal reco')
        ar.stairs(rec/np.diff(rx),rx,color='#1976b5',linewidth=2,label='Folded trial + fakes')
        ar.stairs(fk/np.diff(rx),rx,color='#9362a3',linestyle=':',label='Fake component')
        ar.errorbar(xc,data/np.diff(rx),yerr=err/np.diff(rx),fmt='o',color='black',markersize=3,label='ICS data')
        good=err>0;res=(rec[good]-data[good])/err[good]
        aq.axhline(0,color='0.5',lw=1);aq.axhspan(-1,1,color='#e8eef2')
        aq.plot(xc[good],res,'o-',color='#1976b5',markersize=3)
        aq.set_ylabel('(prediction - data) / data error');aq.set_title('Residual (diagonal errors only)',fontsize=10)
        at.set_title('Truth prior');ar.set_title('Reco comparison')
        for a in (at,ar):
            a.set_ylabel('Weighted yield / bin width');a.legend(fontsize=8);a.grid(alpha=.18)
            if settings.get('log',False):
                a.set_yscale('symlog',linthresh=max(1e-3,float(np.max(np.abs(data/np.diff(rx))))*.001))
        for a in (at,ar,aq):a.set_xlabel(AXLABEL[variable])
    fig.suptitle(title,fontsize=12)
    fig.text(.5,.013,'R fixed; bin-constant truth weights. Full acceptance projections. Fake shape fixed. No unfolding / no fit probability.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.035,1,.94),pad=1.8)
    return result


class App:
    def __init__(self,root,bundles):
        import tkinter as tk
        from tkinter import ttk
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
        self.tk=tk;self.ttk=ttk;self.root=root;self.bundles=[];self.pending=None;self.events=queue.Queue();self.busy=False
        root.title('Prior Explorer | fixed response + ICS data');root.geometry('1460x960')
        top=ttk.Frame(root,padding=6);top.pack(fill='x')
        ttk.Button(top,text='Open prepared bundle(s)',command=self.open_files).pack(side='left')
        ttk.Button(top,text='Prepare from Machine.C',command=self.prepare_dialog).pack(side='left',padx=6)
        ttk.Button(top,text='Save current trial',command=self.save).pack(side='left')
        ttk.Button(top,text='Load trial settings',command=self.load_settings).pack(side='left',padx=6)
        self.selection=tk.StringVar();self.combo=ttk.Combobox(top,textvariable=self.selection,state='readonly',width=60);self.combo.pack(side='left',padx=8)
        self.combo.bind('<<ComboboxSelected>>',lambda e:self.reset())
        select=ttk.Frame(root,padding=5);select.pack(fill='x')
        self.cent=tk.StringVar(value=CENT[0]);self.obs=tk.StringVar(value=LABELS[0]);self.mode=tk.StringVar(value='MC tilt')
        for label,var,choices in [('Centrality',self.cent,CENT),('Response',self.obs,LABELS),('Prior family',self.mode,['MC tilt','Power law'])]:
            ttk.Label(select,text=label).pack(side='left',padx=(8,3));box=ttk.Combobox(select,textvariable=var,values=choices,state='readonly',width=18);box.pack(side='left');box.bind('<<ComboboxSelected>>',lambda e:self.schedule())
        self.normalize=tk.BooleanVar(value=True);self.log=tk.BooleanVar(value=False)
        ttk.Checkbutton(select,text='Match total reco yield',variable=self.normalize,command=self.schedule).pack(side='left',padx=8)
        ttk.Checkbutton(select,text='Log / signed log',variable=self.log,command=self.schedule).pack(side='left')
        ttk.Button(select,text='Reset prior',command=self.reset).pack(side='left',padx=8)
        body=ttk.Frame(root);body.pack(fill='both',expand=True)
        controls=ttk.Frame(body,padding=10,width=225);controls.pack(side='left',fill='y');controls.pack_propagate(False)
        self.values={};self.outputs={}
        sliders=[('a','pT tilt a',-4,4,0),('b','pT curvature b',-3,3,0),('p0','Power law p0 [GeV/c]',.3,20,3),('n','Power law n',1,15,5),('c','Observable tilt c',-4,4,0),('d','pT-observable correlation d',-3,3,0),('fake','Fake normalization',0,2,1),('amplitude','Overall yield factor',.5,1.5,1)]
        for key,label,lo,hi,val in sliders:
            row=ttk.Frame(controls);row.pack(fill='x',pady=(8,0));ttk.Label(row,text=label).pack(side='left')
            out=ttk.Label(row,text=f'{val:.3f}');out.pack(side='right');self.outputs[key]=out
            var=tk.DoubleVar(value=val);self.values[key]=var
            ttk.Scale(controls,from_=lo,to=hi,variable=var,command=lambda v:self.schedule()).pack(fill='x')
        ttk.Label(controls,text='MC tilt:\nw = exp[a L + b L² + c u + d L u]\nL = log(pT / 8 GeV/c)\nu = observable scaled to [-1,1]\n\nPower law replaces the pT marginal:\nf(pT) ~ (1 + pT/p0)^(-n)\nRetains nominal observable | pT\nbefore c and d are applied.\n\nWeights are evaluated once per\ntruth bin and held constant within it.\n\na,b: MC tilt only\np0,n: power law only\nc,d: 2D models only',wraplength=205,justify='left').pack(pady=16)
        plot=ttk.Frame(body);plot.pack(side='left',fill='both',expand=True)
        self.fig=Figure(figsize=(12,7.8),dpi=100);self.canvas=FigureCanvasTkAgg(self.fig,master=plot);self.canvas.get_tk_widget().pack(fill='both',expand=True)
        NavigationToolbar2Tk(self.canvas,plot)
        self.status=tk.StringVar(value='Open a bundle or prepare one from Machine.C.');ttk.Label(root,textvariable=self.status,padding=6,wraplength=1400).pack(fill='x')
        for b in bundles:self.add(b)
        root.after(200,self.poll);self.schedule()

    def add(self,bundle):
        self.bundles.append(bundle);self.combo.configure(values=[b.label for b in self.bundles]);self.combo.current(len(self.bundles)-1);self.schedule()
    def current(self):
        i=self.combo.current()
        if i<0:raise ValueError('Open a prepared bundle first')
        b=self.bundles[i];key=f'c{CENT.index(self.cent.get())}_{KEYS[LABELS.index(self.obs.get())]}'
        if key not in b.models:raise ValueError('This model is not in the bundle')
        return b,b.models[key],key
    def settings(self):
        return dict(mode=self.mode.get(),**{k:v.get() for k,v in self.values.items()},normalize=self.normalize.get(),log=self.log.get())
    def schedule(self):
        if self.pending:self.root.after_cancel(self.pending)
        self.pending=self.root.after(100,self.update)
    def update(self):
        self.pending=None
        for k,v in self.values.items():self.outputs[k].configure(text=f'{v.get():.3f}')
        if not self.bundles:return
        try:
            b,m,key=self.current();s=self.settings();v=key.split('_',1)[1]
            title=f'D0 {b.metadata["d0_min"]:g}-{b.metadata["d0_max"]:g} | {b.metadata.get("prior_order","?")} | {self.cent.get()} | {self.obs.get()} | '+('reco yield matched' if s['normalize'] else 'fixed normalization')
            r=draw(self.fig,m,v,s,title);self.canvas.draw_idle()
            valid=m.variance>0;res=(r['reco'][valid]-m.data[valid])/np.sqrt(m.variance[valid]);w=r['weights'][m.truth>0]
            self.status.set(f'R is fixed. Truth weights {w.min():.3g} to {w.max():.3g}; reco/data total {r["reco"].sum()/m.data.sum():.4f}; diagonal residual RMS {np.sqrt(np.mean(res**2)):.3g} (descriptive only). '+('SYNTHETIC DEMO.' if b.metadata.get('demo') else 'MC/data-fit covariance not included.'))
        except Exception as e:self.status.set(str(e))
    def reset(self):
        self.mode.set("MC tilt")
        s=defaults()
        for k,v in self.values.items():v.set(s[k])
        self.schedule()
    def open_files(self):
        from tkinter import filedialog,messagebox
        for p in filedialog.askopenfilenames(filetypes=[('Prior bundles','*.npz')]):
            try:self.add(load_bundle(p))
            except Exception as e:messagebox.showerror('Cannot open bundle',str(e))
    def save(self):
        from tkinter import filedialog,messagebox
        try:
            b,m,key=self.current();s=self.settings();r=m.predict(s)
            path=filedialog.asksaveasfilename(defaultextension='.json',initialfile='prior_trial.json',filetypes=[('Trial settings','*.json')])
            if not path:return
            p=Path(path);payload=dict(bundle=str(b.path.resolve()),metadata=b.metadata,model=key,settings=s,
                note='Exploratory data-dependent prior, not a validated unfolding result; weights constant in truth bins; fixed fake shape')
            p.write_text(json.dumps(payload,indent=2,ensure_ascii=False))
            np.savez_compressed(p.with_suffix('.npz'),**r,tx=m.tx,ty=m.ty,rx=m.rx,ry=m.ry,data=m.data,data_variance=m.variance)
            self.update()  # Flush any pending slider redraw before saving the figure.
            self.fig.savefig(p.with_suffix('.png'),dpi=170)
            self.status.set('Saved settings, numerical trial and PNG: '+str(p.with_suffix('')))
        except Exception as e:messagebox.showerror('Save failed',str(e))
    def load_settings(self):
        from tkinter import filedialog,messagebox
        path=filedialog.askopenfilename(filetypes=[('Trial settings','*.json')])
        if not path:return
        try:
            payload=json.loads(Path(path).read_text());target=Path(payload['bundle']).resolve()
            found=next((i for i,b in enumerate(self.bundles) if b.path.resolve()==target),None)
            if found is None:self.add(load_bundle(target))
            else:self.combo.current(found)
            c,key=payload['model'].split('_',1);self.cent.set(CENT[int(c[1:])]);self.obs.set(LABELS[KEYS.index(key)])
            s=payload['settings'];self.mode.set(s['mode']);self.normalize.set(s['normalize']);self.log.set(s.get('log',False))
            for k,v in self.values.items():v.set(s[k])
            self.schedule()
        except Exception as e:messagebox.showerror('Cannot load settings',str(e))
    def prepare_dialog(self):
        from tkinter import filedialog,messagebox
        if self.busy:messagebox.showinfo('Preparing','A preparation job is already running.');return
        tk,ttk=self.tk,self.ttk;win=tk.Toplevel(self.root);win.title('Prepare exact Machine response + ICS data');win.geometry('690x360')
        vals={k:tk.StringVar(value=v) for k,v in [('machine',str(Path.cwd()/'Unfolding/Machine.C')),('out',str(Path.cwd()/'PriorExplorerBundles')),('lo','1'),('hi','5'),('root','root')]}
        for row,(key,title) in enumerate([('machine','Machine.C / Machine2.C'),('out','Output directory'),('lo','D0 pT minimum'),('hi','D0 pT maximum'),('root','ROOT executable')]):
            ttk.Label(win,text=title).grid(row=row,column=0,padx=8,pady=8,sticky='w');ttk.Entry(win,textvariable=vals[key],width=62).grid(row=row,column=1,padx=5)
            if key=='machine':ttk.Button(win,text='...',command=lambda:vals['machine'].set(filedialog.askopenfilename(filetypes=[('ROOT macro','*.C')]) or vals['machine'].get())).grid(row=row,column=2)
        before=tk.BooleanVar(value=True);ttk.Checkbutton(win,text='cutBefore=true (D0 cut, then reshape)',variable=before).grid(row=5,column=0,columnspan=3,pady=8)
        ttk.Label(win,text='Uses your own loading functions. Builds once, stops before unfolding.\nOriginal source is unchanged; diagnostics/logs go to an isolated prepare directory.').grid(row=6,column=0,columnspan=3,pady=8)
        def start():
            try:args=(vals['machine'].get(),float(vals['lo'].get()),float(vals['hi'].get()),vals['out'].get(),before.get(),vals['root'].get())
            except ValueError as e:messagebox.showerror('Invalid range',str(e));return
            self.busy=True;win.destroy();self.status.set('Building response; this can take as long as the MC loading in Machine.C.')
            def worker():
                try:self.events.put(('bundle',prepare(*args,log=lambda s:self.events.put(('status',s)))))
                except Exception as e:self.events.put(('error',str(e)))
                finally:self.events.put(('done',None))
            threading.Thread(target=worker,daemon=True).start()
        ttk.Button(win,text='Prepare and open',command=start).grid(row=7,column=1,pady=8)
    def poll(self):
        from tkinter import messagebox
        try:
            while True:
                kind,value=self.events.get_nowait()
                if kind=='bundle':self.add(value)
                elif kind=='error':messagebox.showerror('Preparation failed',value);self.status.set(value)
                elif kind=='done':self.busy=False
                else:self.status.set(value)
        except queue.Empty:pass
        self.root.after(200,self.poll)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundles',nargs='*',help='Prepared .npz bundles')
    p.add_argument('--prepare',metavar='MACHINE',help='Prepare a bundle using the local Machine.C, then exit')
    p.add_argument('--d0',nargs=2,type=float,default=[1,5],metavar=('MIN','MAX'))
    p.add_argument('--cut-after',action='store_true',help='Use cutBefore=false for preparation; default true')
    p.add_argument('--outdir',default='PriorExplorerBundles');p.add_argument('--root',default='root')
    p.add_argument('--demo',action='store_true',help='Open a clearly labelled synthetic example')
    p.add_argument('--render',metavar='PNG',help='Headless rendering for checking a prepared bundle or demo')
    a=p.parse_args()
    if a.prepare:
        prepare(a.prepare,*a.d0,a.outdir,not a.cut_after,a.root);return
    bundles=[load_bundle(x) for x in a.bundles]
    if a.demo:
        tmp=Path(tempfile.mkdtemp(prefix='prior_explorer_demo_'));bundles.append(demo_bundle(tmp/'demo.npz'))
    if a.render:
        if not bundles:p.error('--render needs a bundle or --demo')
        import matplotlib;matplotlib.use('Agg')
        from matplotlib.figure import Figure
        b=bundles[-1];fig=Figure(figsize=(14,8),dpi=120);draw(fig,b.models['c0_lambda1'],'lambda1',defaults(),b.label+' | 0-10%');fig.savefig(a.render);return
    import tkinter as tk
    root=tk.Tk();App(root,bundles);root.mainloop()

if __name__=='__main__':main()
