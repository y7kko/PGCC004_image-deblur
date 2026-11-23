import numpy as np
import warnings

def toep(arr,diag_idx):
    """
    Transforma array unidimensional em matriz de toeplitz, conforme convenção de
    Kamm, Nagy
    """
    if arr.ndim > 1:
        raise ValueError('`arr` há de ser unidimensional')
    arr_len = len(arr)
    output = np.zeros([arr_len,arr_len])
    for idx in range(output.shape[1]):
        idx_roll = idx-diag_idx
        arr_roll = np.roll(arr,idx_roll)
        clip_idx = np.max([0,idx_roll])
        arr_roll[:clip_idx] = 0
        output[:,idx] = arr_roll
    return output

def blur_img(img,A,B):
    return ( B @ img ) @ A.T

def gen_delta(shape,coord=None):
    img_delta = np.zeros(shape) #Dk é delta de kronecker
    if coord is None:
        coord = ((np.asarray(img_delta.shape)-1)/2).astype(int)

    img_delta[coord,coord] = 1
    return img_delta

def vec(x:np.ndarray):
    x = np.atleast_2d(x)
    return x.reshape(np.prod(x.shape),1,order='F').squeeze()

def kronecker_split(P,mid_idx):
    rank = np.linalg.matrix_rank(P)
    if rank != 1:
        warnings.warn(f"rank(P) = {rank}, this function only works with rank = 1")

    
    u,s,v = np.linalg.svd(P)
    v = v.T
    
    if np.sum(u) < 0:
        u1 = -u
    if np.sum(v) < 0:
        v1 = -v

    a1 = np.sqrt(s[0])*u1[:,0]
    b1 = np.sqrt(s[0])*v1[:,0]

    A = toep(a1,mid_idx[0])
    B = toep(b1,mid_idx[0])
    return A,B

