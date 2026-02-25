#do i need to import os? 

import astropy.units as u
import numpy as np
from astropy.coordinates import Angle, SkyCoord
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
from lamatrix import Spline


def get_p_tdur_t0(tce):
    return tce[1], tce[4], tce[3]



def create_spline_design_matrix(time,tdur,order=3):

    #normalizing the time range
    t = (time - time.mean())/(time.max() - time.min())

    #setting the spline knots
    knot_spacing = (3.*tdur) / (time.max() - time.min())
    knots = np.arange(np.min(t) - knot_spacing * (order - 1), np.max(t) + (knot_spacing * (order - 1)) + knot_spacing, knot_spacing)

    #creating spline design matrix
    model = Spline('x', knots=knots, order=order)
    dM = model.design_matrix(x=t)

    return dM



#THIS NEEDS TO BE ADDED IN
def create_CBV_design_matrix(time,mask):
    #will be added in the future, can be used instead of 

    return 1.



def create_polynomial_design_matrix(time):

    #normalizing the time range
    t = (time - time.mean())/(time.max() - time.min())
    
    #compute spacing between points
    dt = np.diff(t)/np.median(np.diff(t))

    #identify where spacing is abnormally large
    breakpoints = np.where(dt > 4)[0] + 1

    #polynomial trends
    dM = np.asarray([t**idx for idx in range(3)]).T


    #create masks for each segment defined by breakpoints
    seg_masks = np.asarray([np.isin(np.arange(len(dM)), i) for i in np.array_split(np.arange(len(dM)), breakpoints)])
        
    #create piece wise design matrix
    dM2 = np.hstack(np.asarray([dM * m[:, None] for m in seg_masks]))

    return dM2



def coords_to_pixels(wcs,ra,dec):

    #initializing skycoord object
    coord = SkyCoord(ra, dec, unit='deg')

    #converting coordinates to pixel using wcs
    pix = wcs.world_to_pixel(coord)
    pix_col = pix[0]
    pix_row = pix[1]

    return pix_col, pix_row



#ADD A CACHE -- OKAY MAYBE NO CACHE USE VIZIER INSTEAD
def query_gaia_background(wcs,ra_targ,dec_targ,radius=25):

    #initializing skycoord object using the Kepler target coords
    coord = SkyCoord(ra=ra_targ, dec=dec_targ, unit=(u.degree, u.degree), frame='icrs')

    #doing a cone search to find all the stars in a cone radius around the Kepler target star
    search = Gaia.cone_search_async(coord, radius=u.Quantity(radius, u.arcsecond))
    #extracting the results 
    res = search.get_results()

    #initializing new skycoord object with all the star coordinates
    coord_bkgd = SkyCoord(ra=res['ra'], dec=res['dec'], unit=(u.degree, u.degree), frame='icrs')
    pix_bkgd = wcs.all_world2pix(coord_bkgd.ra,coord_bkgd.dec, 0)

    return coord_bkgd, pix_bkgd



def query_vizier_background(wcs,ra_targ,dec_targ,tpf_shape):

    #initializing skycoord object using the kepler target coords
    coord = SkyCoord(ra=ra_targ, dec=dec_targ, unit=(u.degree, u.degree), frame='icrs')

    #using TPF shape to get the radius
    rad =(np.nanmax(tpf_shape) * 0.5 ) * 4
    #rad = 0.5 * np.sqrt(tpf_shape[0]**2 + tpf_shape[1]**2) * 3.98


    #getting the catalog of stars using Vizier
    #Catalog = Vizier(columns=['ID','Tmag','pmRA','pmDEC','RAJ2000','DEJ2000','RAOdeg','DEOdeg'],row_limit=-1) -- for TESS v8.2
    Catalog = Vizier(columns=['ID','Tmag','pmRA','pmDEC','RA_ICRS','DE_ICRS'],row_limit=-1) # -- for GAIA DR3
    
    result = Catalog.query_region(coord,radius=Angle(rad, unit="arcsec"),catalog='I/355/gaiadr3') #'IV/39/tic82'

    if len(result) == 0:
        return None
    
    res = result[0]

    #initializing new skycoord object with all the star coordinates -- EVENTUALLY TAKE INTO ACCOUNT THE PM
    coord_bkgd = SkyCoord(ra=res['RA_ICRS'], dec=res['DE_ICRS'], unit=(u.degree, u.degree), frame='icrs')
    pix_bkgd = wcs.all_world2pix(coord_bkgd.ra,coord_bkgd.dec, 0)

    return coord_bkgd, pix_bkgd



def prf_residual(params,prf,data,data_err,origin,shape):#(params, prf, data, data_err):

    #extracting parameter values
    amp = params['amplitude'].value
    row_pix = params['centery'].value
    col_pix = params['centerx'].value

    #getting the ccd location of the target
    row_ccd = origin[0] + row_pix
    col_ccd = origin[1] + col_pix

    #evaluating the prf model
    model = prf.evaluate(targets=[(row_ccd,col_ccd)],origin=origin,shape=shape)[0,:,:]

    #calcualting the residual
    res = ((amp*model - data)/data_err).ravel()#/data_err[tce_num]).ravel()

    #masking out NaNs
    mask = np.isfinite(res)
    
    return res[mask]



def all_prf_residual(params,prf,data,data_err,origin,shape,tces):#(params, prf, data, data_err, origin, shape, tces):
    
    #initializing residuals list
    res = []

    #extracting parameter values
    row_pix = params['centery'].value
    col_pix = params['centerx'].value

    #getting the ccd location of the target
    row_ccd = origin[0] + row_pix
    col_ccd = origin[1] + col_pix

    #evaluating the prf model
    model = prf.evaluate(targets=[(row_ccd,col_ccd)],origin=origin,shape=shape)[0,:,:]
    
    #looping through tces
    for i in range(len(tces)):
        #extracting the amplitude param value for tce
        amp = params[f'amplitude_{i}'].value
        #calculating the residual for tce
        res_i = ((amp*model - data[i])/data_err[i]).ravel()
        res.append(res_i)

    #concatenating all the residuals
    res = np.concatenate(res)

    #masking out any NaNs
    mask = np.isfinite(res)

    return res[mask]



