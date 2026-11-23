"""Finite grid search of fixed-response priors; no optimizer or GUI dependencies."""
import itertools
import heapq
import math
import time
from contextlib import nullcontext
import numpy as np
import fake_weights

PARAMETERS = ('a', 'b', 'p0', 'n', 'c', 'd', 'fake', 'amplitude')
DEFAULT_GRID = dict(a=(-4,4,.2), b=(-1,1,.1), p0=(.5,10,.5), n=(2,10,.5),
                    c=(-2,2,.5), d=(-1,1,.5), fake=(.5,1.5,.1), amplitude=(.8,1.2,.05))

def grid_values(lo, hi, step):
    lo,hi,step=map(float,(lo,hi,step))
    if not all(map(math.isfinite,(lo,hi,step))) or hi<lo or step<=0:
        raise ValueError('Grid needs finite min <= max and step > 0')
    count=int(math.floor((hi-lo)/step+1e-10))+1
    if count>2_000_000: raise ValueError('Too many points on one grid axis')
    # The last point is <= max; max is included only when reached by the step.
    return lo+step*np.arange(count,dtype=float)

def plans(settings, two_d, grid, families):
    unknown=set(grid)-set(PARAMETERS)
    if unknown:raise ValueError('Unknown grid parameters: '+', '.join(sorted(unknown)))
    result=[]
    if not families or len(set(families))!=len(families): raise ValueError('Choose distinct prior families')
    for mode in families:
        if mode not in ('MC tilt','Power law'): raise ValueError('Unknown prior family')
        keys=list(('a','b') if mode=='MC tilt' else ('p0','n'))
        if two_d: keys+=['c','d']
        keys+=['fake','amplitude']
        axes={}
        for key in keys:
            values=grid_values(*grid[key]) if key in grid else np.array([float(settings[key])])
            if key in ('p0','n','amplitude') and np.any(values<=0): raise ValueError(key+' must be positive')
            if key=='fake' and np.any(values<0): raise ValueError('Fake normalization must be nonnegative')
            axes[key]=values
        total=math.prod(len(v) for v in axes.values())
        result.append((mode,axes,total))
    if sum(r[2] for r in result)>2_000_000: raise ValueError('More than 2,000,000 combinations; increase steps or fix parameters')
    return result

def predictions(model, base, mode, keys, points):
    """Vectorized equivalent of Model.predict, one row per candidate."""
    count=len(points)
    params={k:np.full(count,float(base[k])) for k in PARAMETERS}
    params.update({k:points[:,i] for i,k in enumerate(keys)})
    x,y=model.coordinates(); logx=np.log(np.maximum(x,1e-6)/8.)
    u=(2*(y-model.ty[0])/(model.ty[-1]-model.ty[0])-1) if model.two_d else np.zeros_like(x)
    exponent=params['c'][:,None]*u+params['d'][:,None]*(logx*u)
    if mode=='MC tilt':
        exponent+=params['a'][:,None]*logx+params['b'][:,None]*logx**2
        shape=np.ones((count,len(x)))
    else:
        nodes,quad=np.polynomial.legendre.leggauss(12)
        lo,hi=model.tx[:-1],model.tx[1:]
        pts=(lo[:,None]+hi[:,None])/2+(hi-lo)[:,None]*nodes/2
        integral=np.sum((1+pts[None,:,:]/params['p0'][:,None,None])**(-params['n'][:,None,None])*quad,axis=2)*(hi-lo)/2
        marginal=model.truth.reshape((-1,len(lo))).sum(axis=0)
        ratio=np.divide(integral,marginal[None,:],out=np.zeros_like(integral),where=marginal[None,:]>0)
        shape=np.tile(ratio,(1,max(1,len(model.ty)-1)))
    valid=np.all(np.isfinite(exponent)&(np.abs(exponent)<=500),axis=1)
    with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
        weights=shape*np.exp(np.clip(exponent,-500,500))
        truth=weights*model.truth[None,:]
        norm=model.truth.sum()/truth.sum(axis=1)[:,None]
        weights*=norm;truth*=norm
        fk=fake_weights.predict(model,weights,dict(base,mode=mode),c_values=params['c'],d_values=params['d'])
        raw=truth@model.response.T+params['fake'][:,None]*fk
        denominator=raw.sum(axis=1) if base['normalize'] else np.full(count,model.reco.sum())
        reco=raw*(model.data.sum()/denominator*params['amplitude'])[:,None]
    valid &= np.all(np.isfinite(reco),axis=1)&(denominator>0)
    return reco,valid

def scan(model,settings,grid,families,progress=None,cancel=None,chunk_size=128):
    if chunk_size<1: raise ValueError('Positive chunk size required')
    if model.data.sum()<=0: raise ValueError('Total data yield must be positive')
    plan=plans(settings,model.two_d,grid,families)
    good=model.variance>0
    if not np.any(good): raise ValueError('No reco bins with positive data errors')
    baseline=model.predict(settings)['reco']
    baseline_score=float(np.sum((baseline[good]-model.data[good])**2/model.variance[good]))
    total=sum(p[2] for p in plan);done=0;valid_count=0;best=None;top=[];start=time.monotonic();stopped=False
    try:
        from threadpoolctl import threadpool_limits
        limits=threadpool_limits(limits=1,user_api='blas')
    except ImportError: limits=nullcontext()
    with limits:
        for mode,axes,_ in plan:
            keys=list(axes);iterator=itertools.product(*axes.values())
            while True:
                if cancel is not None and cancel.is_set():stopped=True;break
                block=list(itertools.islice(iterator,chunk_size))
                if not block:break
                points=np.asarray(block);pred,valid=predictions(model,settings,mode,keys,points)
                with np.errstate(over='ignore',invalid='ignore'):
                    scores=np.sum((pred[:,good]-model.data[None,good])**2/model.variance[None,good],axis=1)
                valid &= np.isfinite(scores);valid_count+=int(valid.sum())
                for i in np.flatnonzero(valid):
                    score=float(scores[i]);index=done+int(i)
                    if len(top)<25 or score < -top[0][0]:
                        trial=dict(settings,mode=mode,**{k:float(points[i,j]) for j,k in enumerate(keys)})
                        record=dict(score=score,index=index,settings=trial)
                        if len(top)<25:heapq.heappush(top,(-score,-index,record))
                        else:heapq.heapreplace(top,(-score,-index,record))
                        if best is None or score<best['score']:best=record
                done+=len(block)
                if progress:progress(done,total,None if best is None else best['score'])
            if stopped:break
    boundary=[]
    if best:
        axes=next(a for mode,a,_ in plan if mode==best['settings']['mode'])
        boundary=[k for k,v in axes.items() if len(v)>1 and (np.isclose(best['settings'][k],v[0]) or np.isclose(best['settings'][k],v[-1]))]
    return dict(best=best,top=[r[2] for r in sorted(top,key=lambda r:(-r[0],-r[1]))],
                completed=not stopped,planned=total,evaluated=done,invalid=done-valid_count,
                seconds=time.monotonic()-start,baseline_score=baseline_score,
                bins_used=int(good.sum()),bins_excluded=int((~good).sum()),boundary_parameters=boundary,
                initial_settings=dict(settings),grid={k:list(map(float,v)) for k,v in grid.items()},families=list(families),
                actual_axes={mode:{k:v.tolist() for k,v in axes.items()} for mode,axes,_ in plan},
                criterion='sum((prediction-data)^2 / data_variance); full reco vector; diagonal errors only, no fit probability')
