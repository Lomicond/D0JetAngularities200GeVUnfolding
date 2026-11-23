"""Numerical regression checks; no GUI, ROOT or analysis input required."""
import unittest
from pathlib import Path
import tempfile
import numpy as np
import prior_explorer as pe

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

if __name__=='__main__':unittest.main()
