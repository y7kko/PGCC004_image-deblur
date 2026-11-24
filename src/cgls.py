import numpy as np
from numpy import fft
from numpy.linalg import norm
import matplotlib.pyplot as plt

class CGLS():
    def __init__(self,max_iter=50):
        self.max_iter = max_iter
        self.res_norm = np.zeros(max_iter)
        self.sol_norm = np.zeros(max_iter)
        pass

    def _foward(self,H,X): #matvec
        return fft.ifft2(fft.fft2(H)*fft.fft2(X)).real

    def _residue(self,H,X,X_true):
        return X_true - self._foward(H,X)

    def _adjoint(self,H,res): #rmatvec
        H_w = (fft.fft2(H)).conj()
        adj_w = H_w * fft.fft2(res)
        out = fft.ifft2(adj_w).real
        return out

    def fit(self,H,G): 
        """Executes the CGLS regularization method

        Args:
            H (ndarray): blur kernel
            G (ndarray): blurred image

        Returns:
            self
        """
        X0 = np.zeros(G.shape) #x0
        res = self._residue(H,X0,G) # residual
        s = self._adjoint(H,res) # Adjoint

        sol = X0.copy()
        p = s.copy()
        self.cgls_sols = np.zeros([self.max_iter,*list(sol.shape)])
        for iter in range(self.max_iter):
            q = self._foward(H,p)

            alpha = (norm(s) / norm(q))**2
            sol = sol + (alpha * p)
            res = res - (alpha * q)

            self.cgls_sols[iter] = np.roll(sol,tuple(dim//2 for dim in sol.shape),axis=[0,1])

            s_old = s
            s = self._adjoint(H,res)
            beta = (norm(s) / norm(s_old))**2
            p = s + beta*p

            self.res_norm[iter] = norm(res)
            self.sol_norm[iter] = norm(sol)
        return self

    def plot_lcurve(self,):
        plt.plot(self.res_norm,self.sol_norm,marker='x')
        plt.ylabel('$\\|X\\|_2$')
        plt.xlabel('$\\|F_\\mathrm{true}- H*X\\|_2$')

    def __getitem__(self, key):
        return self.cgls_sols[key]
