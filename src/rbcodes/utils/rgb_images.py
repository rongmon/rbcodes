import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def stiff_rgb(r, g, b,
              min_percent=0.01, max_percent=99.99,
              gamma=2.2,
              color_balance=True,
              savefile=None, dpi=300):
    """
    Emulate STIFF's RGB image generation using gamma + percentile stretch.

    Parameters
    ----------
    r, g, b : 2D numpy arrays
        Red, green, and blue channel images (same shape).
    min_percent : float
        Lower percentile for contrast stretch (e.g., 0.01).
    max_percent : float
        Upper percentile for contrast stretch (e.g., 99.99).
    gamma : float
        Gamma correction value (e.g., 2.2).
    color_balance : bool
        Whether to normalize each channel by its median.
    savefile : str
        If provided, save the RGB image to this path.
    dpi : int
        Resolution for saved image.

    Returns
    -------
    rgb : 3D numpy array
        RGB image array in shape (height, width, 3).
    """
    # Stack to handle percentile clipping together
    all_data = np.stack([r, g, b], axis=0)
    vmin = np.percentile(all_data, min_percent)
    vmax = np.percentile(all_data, max_percent)

    def stretch(img):
        img = np.clip(img, vmin, vmax)
        img = (img - vmin) / (vmax - vmin)
        img = img ** (1 / gamma)
        return img

    r_s = stretch(r)
    g_s = stretch(g)
    b_s = stretch(b)

    # Optional color balance (normalize each channel)
    if color_balance:
        r_s /= np.median(r_s)
        g_s /= np.median(g_s)
        b_s /= np.median(b_s)

    # Re-normalize to [0, 1] after color balance
    rgb_stack = np.stack([r_s, g_s, b_s], axis=-1)
    rgb_stack /= np.max(rgb_stack)

    if savefile:
        plt.imsave(savefile, rgb_stack, dpi=dpi)

    return rgb_stack
