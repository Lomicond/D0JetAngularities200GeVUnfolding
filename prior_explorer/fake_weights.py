"""Fake templates conditional on actual truth coordinates, including overflow.
The truth-bin weights already contain the normalization used for matched/missed.
"""
import numpy as np
POLICIES=('Truth weights / keep outside','Truth weights / edge extension','Fixed')

def predict(model,weights,settings,c_values=None,d_values=None):
    weights=np.atleast_2d(weights);batch=len(weights)
    policy=settings.get('fake_policy','Fixed')
    if policy not in POLICIES:raise ValueError('Unknown fake policy')
    if policy=='Fixed':return np.broadcast_to(model.fakes,(batch,len(model.fakes))).copy()
    if not model.fake_pt.size:
        raise ValueError('This old bundle has no fake truth maps. Prepare a new v3 bundle or explicitly select Fixed.')
    nx=len(model.tx)-1;ny=max(1,len(model.ty)-1)
    cube=weights.reshape(batch,ny,nx)
    c=np.broadcast_to(settings.get('c',0) if c_values is None else c_values,(batch,))
    d=np.broadcast_to(settings.get('d',0) if d_values is None else d_values,(batch,))
    # Factorization w(x,y)=C*g(x)*h(x,y). Outside the y range, keep h=1
    # (or extend its edge). Never discard the valid jet-pT part just because y
    # failed acceptance. This is continuous as c,d approach zero.
    xweights=cube[:,0,:].copy()
    if model.two_d:
        u0=2*((model.ty[0]+model.ty[1])/2-model.ty[0])/(model.ty[-1]-model.ty[0])-1
        lx=np.log(np.maximum((model.tx[:-1]+model.tx[1:])/2,1e-6)/8.)
        xweights*=np.exp(-u0*(c[:,None]+d[:,None]*lx))
    marginal=model.truth.reshape(ny,nx).sum(axis=0)
    if settings['mode']=='Power law':xweights[:,marginal==0]=1.
    wx=np.pad(xweights,((0,0),(1,1)),mode='edge')
    if policy=='Truth weights / keep outside':wx[:,[0,-1]]=1.
    wx=np.concatenate((wx,np.ones((batch,1))),axis=1)
    out=wx@model.fake_pt.T
    if model.two_d and np.any((c!=0)|(d!=0)):
        if not model.fake_joint.size:raise ValueError('Missing joint fake truth map')
        full=cube.copy()
        if settings['mode']=='Power law':full[:,:,marginal==0]=1.
        extended=np.pad(full,((0,0),(1,1),(1,1)),mode='edge')
        delta=extended-wx[:,None,:-1]
        if policy=='Truth weights / keep outside':
            delta[:,[0,-1],:]=0.;delta[:,:,[0,-1]]=0.
        djoint=np.concatenate((delta.reshape(batch,-1),np.zeros((batch,1))),axis=1)
        out+=djoint@model.fake_joint.T
    # Roundoff cancellation can produce tiny negative numbers at zero.
    out=np.maximum(out,0.)
    return out


def coverage(model):
    """Nominal fractions for interpretation, not statistical estimates."""
    if not model.fake_pt.size:return None
    total=model.fakes.sum()
    if total<=0:return dict(pt_outside=0.,pt_unknown=0.,joint_outside=0.,joint_unknown=0.)
    pt=model.fake_pt.sum(axis=0);result=dict(pt_outside=float((pt[0]+pt[-2])/total),pt_unknown=float(pt[-1]/total))
    if model.two_d:
        joint=model.fake_joint.sum(axis=0);core=joint[:-1].reshape(len(model.ty)+1,len(model.tx)+1)
        result.update(joint_outside=float((core.sum()-core[1:-1,1:-1].sum())/total),joint_unknown=float(joint[-1]/total))
    return result
