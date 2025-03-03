import pytest
import numpy as np

from damask import Crystal

class TestCrystal:

    def test_double_to_lattice(self):
        c = Crystal(lattice='cF')
        with pytest.raises(KeyError):
            c.to_lattice(direction=np.ones(3),plane=np.ones(3))

    def test_double_to_frame(self):
        c = Crystal(lattice='cF')
        with pytest.raises(KeyError):
            c.to_frame(uvw=np.ones(3),hkl=np.ones(3))

    @pytest.mark.parametrize('lattice,a,b,c,alpha,beta,gamma',
                            [
                             ('aP',0.5,2.0,3.0,0.8,0.5,1.2),
                             ('mP',1.0,2.0,3.0,np.pi/2,0.5,np.pi/2),
                             ('oI',0.5,1.5,3.0,np.pi/2,np.pi/2,np.pi/2),
                             ('tP',0.5,0.5,3.0,np.pi/2,np.pi/2,np.pi/2),
                             ('hP',1.0,None,1.6,np.pi/2,np.pi/2,2*np.pi/3),
                             ('cF',1.0,1.0,None,np.pi/2,np.pi/2,np.pi/2),
                            ])
    def test_bases_contraction(self,lattice,a,b,c,alpha,beta,gamma):
        c = Crystal(lattice=lattice,
                    a=a,b=b,c=c,
                    alpha=alpha,beta=beta,gamma=gamma)
        assert np.allclose(np.eye(3),np.einsum('ik,jk',c.basis_real,c.basis_reciprocal))


    @pytest.mark.parametrize('keyFrame,keyLattice',[('uvw','direction'),('hkl','plane'),])
    @pytest.mark.parametrize('vector',np.array([
                                                [1.,1.,1.],
                                                [-2.,3.,0.5],
                                                [0.,0.,1.],
                                                [1.,1.,1.],
                                                [2.,2.,2.],
                                                [0.,1.,1.],
                                               ]))
    @pytest.mark.parametrize('lattice,a,b,c,alpha,beta,gamma',
                            [
                             ('aP',0.5,2.0,3.0,0.8,0.5,1.2),
                             ('mP',1.0,2.0,3.0,np.pi/2,0.5,np.pi/2),
                             ('oI',0.5,1.5,3.0,np.pi/2,np.pi/2,np.pi/2),
                             ('tP',0.5,0.5,3.0,np.pi/2,np.pi/2,np.pi/2),
                             ('hP',1.0,1.0,1.6,np.pi/2,np.pi/2,2*np.pi/3),
                             ('cF',1.0,1.0,1.0,np.pi/2,np.pi/2,np.pi/2),
                            ])
    def test_to_frame_to_lattice(self,lattice,a,b,c,alpha,beta,gamma,vector,keyFrame,keyLattice):
        c = Crystal(lattice=lattice,
                    a=a,b=b,c=c,
                    alpha=alpha,beta=beta,gamma=gamma)
        assert np.allclose(vector,
                           c.to_frame(**{keyFrame:c.to_lattice(**{keyLattice:vector})}))
