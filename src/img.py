from PIL import Image
import numpy as np

def read_img(path):
    return np.array(Image.open(path))

def gaussian_psf(N=256,radius=7.0,split_xy=False):
    x = np.arange(N) - N // 2
    y = np.arange(N) - N // 2
    X, Y = np.meshgrid(x, y)

    if split_xy:
        Py_a = np.exp(-y**2 / (2*radius**2))
        Px_a = np.exp(-x**2 / (2*radius**2))
        return Px_a,Py_a

    PSF = np.exp(-((X)**2 / (2*radius**2) + (Y)**2 / (2*radius**2)))
    return PSF
