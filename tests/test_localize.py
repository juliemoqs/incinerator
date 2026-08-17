"""Test the Localization Method"""

import numpy as np
import pandas as pd
from astropy.utils.data import get_pkg_data_filename

import incinerator.localize as loc
from incinerator.utils import coords_to_pixels

file_name = get_pkg_data_filename("data/test.fits.gz")


def test_localize():
    """
    Test that the Localize class correctly loads a Kepler TPF and performs transit localization.
    """

    tces = pd.DataFrame({'star_id':6922244.01,
                             'period':3.5224,
                             'mes':np.nan,
                             't0':121.1194228+2454833,
                             'tdur':.13},index=[0])
    

    loc_obj = loc.Localize.from_tpf_info(file_name,tces,'6922244',mission='kepler')

    #testing dimensions of flux, time, and flux_err
    assert (loc_obj.flux.ndim == 3 and loc_obj.flux_err.ndim == 3)
    assert (loc_obj.time.size == loc_obj.pix.shape[0] and loc_obj.time.size == loc_obj.pix_err.shape[0])
    
    #testing dimensions of design matrix
    loc_obj.build_design_matrix()
    assert(loc_obj.design_matrix.shape[0] == 1626 and loc_obj.design_matrix.shape[1] == 92)

    loc_obj.solve_transit_weights()

    #testing that the coordinate to pixel transformation is being done correctly
    pix_col, pix_row = coords_to_pixels(loc_obj.wcs,loc_obj.ra_targ,loc_obj.dec_targ)
    assert np.round(pix_col,decimals=2) == 2.76
    assert np.round(pix_row,decimals=2) == 2.19

    fit = loc_obj.fit_to_heatmap(method='prf',which_tce=0)

    #extract fit metrics
    x0 = fit[1]['centerx']
    y0 = fit[1]['centery']
    x0_err = fit[1]['sigmax']
    y0_err = fit[1]['sigmay']
    #checking to make sure the pos and unc are positive
    assert x0 > 0 and y0 > 0 and x0_err > 0 and y0_err > 0

    #calculating differences
    dx = x0 - pix_col
    dy = y0 - pix_row

    #finding how many sigma away the center of the fit is from the actual location
    nsigma = np.sqrt((dx/x0_err)**2+(dy/y0_err)**2)
    #checking nsigma
    assert np.round(nsigma,decimals=0) == 5.0



