"""Test the Localization Method"""

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

import incinerator.localize as loc
from incinerator.localize import Localize
from incinerator.utils import coords_to_pixels

file_name = get_pkg_data_filename("data/test.fits.gz")


def test_load_tpf_info():
    """
    Making sure that header is being properly read
    """
    info = Localize._get_kepler_tpf_info(file_name)

    #testing info is accurately taken from hdulist
    assert info['ra_targ'] == 281.28812
    assert info['dec_targ'] == 42.45108
    assert info['origin_row'] == 190
    assert info['origin_col'] == 674
    assert info['channel'] == 3



def test_localize():
    """
    Check that the Localize class methods work
    """
    hdulist = fits.open(file_name)

    #extract time, flux cube, and flux cube error
    time = hdulist[1].data['TIME']
    flux = hdulist[1].data['FLUX']
    flux_err = hdulist[1].data['FLUX_ERR']
    
    #testing dimensions of flux, time, and flux_err
    assert (flux.ndim == 3 and flux_err.ndim == 3)
    assert (time.size == flux.shape[0] and time.size == flux_err.shape[0])

    tces = pd.DataFrame({'star_id':6922244.01,
                             'period':3.5224,
                             'mes':np.nan,
                             't0':121.1194228,
                             'tdur':.13},index=[0])
    

    loc_obj = loc.Localize(time,flux,flux_err,tces.to_numpy(),'6922244',
                       mission='kepler', file_name = file_name)
    
    #testing dimensions of design matrix
    loc_obj.build_design_matrix()
    assert(loc_obj.design_matrix.shape[0] == 1626 and loc_obj.design_matrix.shape[1] == 92)

    loc_obj.solve_transit_weights()

    #testing that the coordinate to pixel transformation is being done correctly
    pix_col, pix_row = coords_to_pixels(loc_obj.wcs,loc_obj.ra_targ,loc_obj.dec_targ)
    assert np.round(pix_col,decimals=2) == 2.75
    assert np.round(pix_row,decimals=2) == 2.18

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



