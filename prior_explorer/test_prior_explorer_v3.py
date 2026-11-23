"""Numerical regression checks; no GUI, ROOT or analysis input required."""
import unittest
from pathlib import Path
import tempfile
import numpy as np
import prior_explorer_v3 as pe
import fake_weights as fw
import scan_engine as se
import threading

class ForwardFoldTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.bundle=pe.demo_bundle(Path(self.tmp.name)/'toy.npz')
        self.m=self.bundle.models['c0_lambda1']
    def tearDown(self): self.tmp.cleanup()
    def test_identity_and_efficiency(self):
        m=self.m
        np.testing.assert_allclose(m.response@m.truth+m.fakes,m.reco,rtol=1e-12)
        self.assertTrue(np.all(m.response.sum(axis=0)<=1))
        result=m.predict(pe.defaults())
        np.testing.assert_allclose(result['reco'],result['nominal_reco'])
        np.testing.assert_allclose(result['weights'],1)
    def test_fixed_response_and_weighted_counts(self):
        m=self.m;old=m.response.copy();s=pe.defaults();s.update(a=.8,b=-.2,c=.4,d=.1)
        t,w=m.vary(**{k:s[k] for k in ('mode','a','b','c','d','p0','n')})
        np.testing.assert_allclose((m.migration*w[None,:])/t[None,:],m.response)
        np.testing.assert_array_equal(m.response,old)
        np.testing.assert_allclose(t.sum(),m.truth.sum())
    def test_x_fastest_and_projection(self):
        m=self.m; nx=len(m.tx)-1;ny=len(m.ty)-1
        x,y=m.coordinates()
        np.testing.assert_allclose(x[:nx],(m.tx[:-1]+m.tx[1:])/2)
        self.assertTrue(np.all(y[:nx]==y[0]))
        a=np.arange(nx*ny,dtype=float)
        np.testing.assert_equal(m.project(a,'truth',0),[sum(a[j*nx+i] for j in range(ny)) for i in range(nx)])
        self.assertAlmostEqual(m.project(a,'truth',1).sum(),a.sum())
    def test_power_law_bin_integrals(self):
        m=self.m;p0=3;n=5;t,w=m.vary(mode='Power law',p0=p0,n=n)
        target=p0/(1-n)*((1+m.tx[1:]/p0)**(1-n)-(1+m.tx[:-1]/p0)**(1-n))
        target*=m.truth.sum()/target.sum()
        np.testing.assert_allclose(m.project(t,'truth',0),target,rtol=1e-9)
    def test_signed_data_and_fake_template(self):
        m=self.m; s=pe.defaults();s['normalize']=False
        s['fake']=1;r1=m.predict(s);s['fake']=0;r0=m.predict(s)
        np.testing.assert_allclose(r1['reco']-r0['reco'],r1['fake'],rtol=1e-12,atol=1e-12)
        data=m.data.copy();data[0]=-3
        kw={n:getattr(m,n) for n in pe.Model.__dataclass_fields__};kw['data']=data
        signed=pe.Model(**kw);self.assertEqual(signed.data[0],-3)
        self.assertTrue(np.all(np.isfinite(signed.predict(pe.defaults())['reco'])))
    def test_invalid_response_rejected(self):
        m=self.m;kw={n:getattr(m,n) for n in pe.Model.__dataclass_fields__};kw['truth']=m.truth*.2
        with self.assertRaises(ValueError):pe.Model(**kw)


class ScanTests(ForwardFoldTests):
    def test_vectorized_matches_individual(self):
        for mode in ['MC tilt','Power law']:
            for normalize in [False,True]:
                s=pe.defaults();s.update(normalize=normalize,mode=mode)
                keys=['a','b','p0','n','c','d','fake','amplitude']
                points=np.array([[.4,-.3,3.,5.,-.7,.2,1.,.9],[-1.3,.2,2.1,4.,.6,-.4,.7,1.2]])
                reco,valid=se.predictions(self.m,s,mode,keys,points)
                self.assertTrue(np.all(valid))
                for i,pt in enumerate(points):
                    trial=dict(s,**dict(zip(keys,pt)))
                    np.testing.assert_allclose(reco[i],self.m.predict(trial)['reco'],rtol=1e-12,atol=1e-10)
    def test_recover_known_grid_point_with_fakes(self):
        m=self.m;s=pe.defaults();target=dict(s,a=.4,b=-.2,c=.5,d=-.5)
        expected=m.predict(target)['reco']
        kw={n:getattr(m,n) for n in pe.Model.__dataclass_fields__};kw['data']=expected
        toy=pe.Model(**kw)
        grid=dict(a=(0,.8,.2),b=(-.4,0,.2),c=(0,1,.5),d=(-1,0,.5))
        report=se.scan(toy,s,grid,['MC tilt'],chunk_size=17)
        self.assertEqual(report['evaluated'],135);self.assertTrue(report['completed'])
        self.assertLess(report['best']['score'],1e-18)
        for k in grid:self.assertAlmostEqual(report['best']['settings'][k],target[k])
        self.assertEqual(report['boundary_parameters'],[])
        scores=[r['score'] for r in report['top']];self.assertEqual(scores,sorted(scores))
    def test_signed_data_zero_errors_and_bruteforce(self):
        m=self.m;kw={n:getattr(m,n).copy() for n in pe.Model.__dataclass_fields__}
        kw['data'][0]=-2;kw['variance'][1]=0;m=pe.Model(**kw);s=pe.defaults()
        grid=dict(a=(-.5,.5,.5),b=(-.2,.2,.2));report=se.scan(m,s,grid,['MC tilt'])
        good=m.variance>0;scores=[]
        for a in [-.5,0,.5]:
            for b in [-.2,0,.2]:
                r=m.predict(dict(s,a=a,b=b))['reco'];scores.append(np.sum((r[good]-m.data[good])**2/m.variance[good]))
        self.assertAlmostEqual(report['best']['score'],min(scores),places=8)
        self.assertEqual(report['bins_excluded'],1)
    def test_grid_limits_and_cancel(self):
        np.testing.assert_allclose(se.grid_values(0,1,.3),[0,.3,.6,.9])
        with self.assertRaises(ValueError):se.grid_values(1,0,.1)
        with self.assertRaises(ValueError):se.plans(pe.defaults(),True,{'a':(0,10,.000001)},['MC tilt'])
        stop=threading.Event();stop.set()
        report=se.scan(self.m,pe.defaults(),{},['MC tilt'],cancel=stop)
        self.assertFalse(report['completed']);self.assertIsNone(report['best']);self.assertEqual(report['evaluated'],0)
    def test_two_families_and_1d_inactive_parameters(self):
        m=self.bundle.models['c0_jetPt'];s=pe.defaults()
        grid=dict(a=(0,1,1),b=(0,1,1),p0=(2,3,1),n=(4,5,1),c=(-4,4,.1),d=(-3,3,.1))
        r=se.scan(m,s,grid,['MC tilt','Power law'])
        self.assertEqual(r['evaluated'],8);self.assertNotIn('c',r['actual_axes']['MC tilt'])


class FakeWeightTests(ForwardFoldTests):
    def test_batch_matches_single_for_all_policies(self):
        for model in [self.m,self.bundle.models['c0_jetPt']]:
            for policy in fw.POLICIES:
                for mode in ['MC tilt','Power law']:
                    s=pe.defaults();s.update(mode=mode,fake_policy=policy)
                    keys=['a','b','p0','n','c','d','fake','amplitude']
                    points=np.array([[.4,-.3,3.,5.,-.7,.2,1.,.9],[-1.3,.2,2.1,4.,0.,0.,.7,1.2]])
                    reco,valid=se.predictions(model,s,mode,keys,points)
                    self.assertTrue(np.all(valid))
                    for i,pt in enumerate(points):
                        np.testing.assert_allclose(reco[i],model.predict(dict(s,**dict(zip(keys,pt))))['reco'],rtol=1e-11,atol=1e-9)
    def test_exact_pt_known_outside_and_unknown(self):
        m=self.bundle.models['c0_jetPt'];s=pe.defaults();s.update(a=-1,b=.2)
        t,w=m.vary(a=-1,b=.2);nx=len(m.tx)-1
        expected=m.fake_pt[:,1:nx+1]@w+m.fake_pt[:,0]+m.fake_pt[:,-2]+m.fake_pt[:,-1]
        np.testing.assert_allclose(fw.predict(m,w,s)[0],expected)
        s['fake_policy']='Truth weights / edge extension'
        expected=m.fake_pt[:,1:nx+1]@w+m.fake_pt[:,0]*w[0]+m.fake_pt[:,-2]*w[-1]+m.fake_pt[:,-1]
        np.testing.assert_allclose(fw.predict(m,w,s)[0],expected)
    def test_pt_only_reweights_fakes_outside_second_axis(self):
        m=self.m;s=pe.defaults();s.update(a=-1)
        _,w=m.vary(a=-1);nx=len(m.tx)-1
        expected=m.fake_pt[:,1:nx+1]@w[:nx]+m.fake_pt[:,0]+m.fake_pt[:,-2]+m.fake_pt[:,-1]
        np.testing.assert_allclose(fw.predict(m,w,s)[0],expected)
        # No discontinuous fake-policy change at c=0.
        _,wnear=m.vary(a=-1,c=1e-10)
        np.testing.assert_allclose(fw.predict(m,wnear,dict(s,c=1e-10))[0],expected,rtol=1e-8,atol=1e-8)
    def test_legacy_requires_explicit_fixed_policy(self):
        kw={n:getattr(self.m,n) for n in pe.Model.__dataclass_fields__};kw.update(fake_pt=None,fake_joint=None)
        m=pe.Model(**kw)
        with self.assertRaises(ValueError):m.predict(pe.defaults())
        r=m.predict(dict(pe.defaults(),fake_policy='Fixed'))
        np.testing.assert_allclose(r['reco'],r['nominal_reco'])
    def test_corrupt_fake_map_rejected(self):
        kw={n:getattr(self.m,n) for n in pe.Model.__dataclass_fields__};kw['fake_pt']=self.m.fake_pt*2
        with self.assertRaises(ValueError):pe.Model(**kw)
    def test_fixed_response_but_changed_fake_fraction(self):
        m=self.m;s=pe.defaults();r=m.predict(s);old=m.response.copy()
        s['a']=-2;trial=m.predict(s)
        np.testing.assert_array_equal(old,m.response)
        self.assertFalse(np.allclose(r['fake'],trial['fake']))
    def test_joint_weights_against_explicit_category_loop(self):
        m=self.m;nx=len(m.tx)-1;ny=len(m.ty)-1
        for policy in fw.POLICIES[:2]:
            s=dict(pe.defaults(),a=-.6,b=.2,c=.7,d=-.3,fake_policy=policy)
            _,w=m.vary(**{k:s[k] for k in ('mode','a','b','c','d','p0','n')})
            cube=w.reshape(ny,nx);x=(m.tx[:-1]+m.tx[1:])/2;y=(m.ty[:-1]+m.ty[1:])/2
            u=2*(y-m.ty[0])/(m.ty[-1]-m.ty[0])-1
            gx=cube[0]/np.exp(u[0]*(s['c']+s['d']*np.log(x/8)))
            expected=m.fake_pt[:,-1].copy()
            for by in range(ny+2):
                for bx in range(nx+2):
                    if policy.endswith('keep outside') and bx in (0,nx+1):v=1.
                    elif policy.endswith('keep outside') and by in (0,ny+1):v=gx[bx-1]
                    else:v=cube[min(max(by-1,0),ny-1),min(max(bx-1,0),nx-1)]
                    expected+=m.fake_joint[:,by*(nx+2)+bx]*v
            np.testing.assert_allclose(fw.predict(m,w,s)[0],expected,rtol=1e-12,atol=1e-10)

if __name__=='__main__':unittest.main()
