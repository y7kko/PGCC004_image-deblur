import numpy as np
import matplotlib.pyplot as plt
from .regu import kronecker_split,nagy_solve,psf_split,vec
class Nagy():
    def __init__(self,A=None,B=None):
        self.A = A
        self.B = B
        self.svd_is_computed=False
        pass

    def from_psf(self,psf,mid_idx=None,force_positive=False):
        if mid_idx==None:
            mid_idx = tuple(_len//2 for _len in psf.shape)
        self.A, self.B = kronecker_split(psf,mid_idx=mid_idx,force_positive=force_positive)
        return self

    def get_gcv(self,G,grid_size=100,min_reg=-7,max_reg=0):
        self._compute_svd(self.A,self.B,G)
        reg_grid = np.logspace(min_reg,max_reg,num=grid_size)
        V = np.zeros(len(reg_grid))
        for l_idx,lam in enumerate(reg_grid):
            V[l_idx] = ((self.gh/(self.S**2+lam**2))**2).sum()
            V[l_idx] /= (((np.arange(len(self.S))+1)/(self.S**2+lam**2))**2).sum()
        return reg_grid,V
    
    def fit(self,G,reg):
        self._compute_svd(self.A,self.B,G)
        self.sol= self.Vb @ (self.Gh*self.Ss/(self.Ss**2 + reg**2)) @ self.Va.T
        return self

    def fit_and_return(self,G,reg):
        self.fit(G,reg)
        return self.sol

    def _compute_svd(self,A,B,G):
        if self.svd_is_computed:
            return 
        
        self.Ua,self.Sa,self.Va = np.linalg.svd(A,)
        self.Va = self.Va.T
        self.Ub,self.Sb,self.Vb = np.linalg.svd(B,)
        self.Vb = self.Vb.T

        self.Gh = ((self.Ub.T @ G) @ self.Ua)
        self.gh = vec(self.Gh)

        self.Ss = np.outer(self.Sb,self.Sa)
        self.S = np.kron(self.Sa,self.Sb)

        self.svd_is_computed=True

            
