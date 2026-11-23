from PIL import Image
import numpy as np

def read_img(path):
    return np.array(Image.open(path))

