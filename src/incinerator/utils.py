import astropy.units as u
import diskcache as dc
import numpy as np
from astropy.coordinates import Angle, SkyCoord
from astropy.time import Time
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
from lamatrix import Spline

#creating cache for vizier query
cache = dc.Cache('gaia_cache')


def get_p_tdur_t0(tce):
    """
    Extract the orbital period, transit duration, and transit epoch (t0)
    from a TCE NumPy array.

    Parameters
    ----------
    tce : numpy.ndarray
        One-dimensional array containing TCE parameters in a fixed order.
        This function assumes:
            tce[1] = orbital period
            tce[4] = transit duration
            tce[3] = transit epoch (t0)

    Returns
    -------
    tuple : 
        A tuple (period, tdur, t0).
    """
    return tce[1], tce[4], tce[3]



def create_spline_design_matrix(time,tdur,order=3,spacing_mult = 3.):
    """
    Construct a piecewise spline design matrix.

    Parameters
    ----------
    time : numpy.ndarray
        One-dimensional array of time values.
    tdur : float
        Transit duration in the same units as `time`. Used to set
        the knot spacing.
    order : int, optional
        Order of the spline (default is 3 for cubic splines).
    spacing_mult : int, optional
        Multiplier applied to the transit duration (tdur) to set spline knot spacing (default 3.).

    Returns
    -------
    numpy.ndarray :
        The spline design matrix evaluated at the normalized time values.
    """

    #normalizing the time range
    t = (time - time.mean())/(time.max() - time.min())

    #compute spacing between points
    dt = np.diff(t)/np.median(np.diff(t))

    #setting the spline knots
    knot_spacing = (spacing_mult*tdur) / (time.max() - time.min())

    #identify where spacing is larger then knot spacing
    breakpoints = np.where(dt > 4)[0] + 1

    #split into segments based on breakpoints
    segments = np.array_split(np.arange(len(t)), breakpoints)

    #list to store each segment's spline design matrix
    dM_list = []

    #loop over each segment
    for seg in segments:
        t_seg = t[seg]

        #skip tiny segments that are too short for the spline order
        #if len(t_seg) < order + 1:
            #continue

        #setting up the knots    
        knots = np.arange(np.min(t_seg) - knot_spacing * (order - 1), np.max(t_seg) + (knot_spacing * (order - 1)) + knot_spacing, knot_spacing)

        #creating spline design matrix
        model = Spline('x', knots=knots, order=order)
        dM_seg = model.design_matrix(x=t_seg)

        #expand the segment design matrix to full length
        dM_full = np.zeros((len(t),dM_seg.shape[1]))
        dM_full[seg,:] = dM_seg

        #append matrix to list
        dM_list.append(dM_full)

    #create piece wise design matrix
    dM2 = np.hstack(dM_list)

    return dM2



#THIS NEEDS TO BE ADDED IN
def create_CBV_design_matrix(time,mask):
    #will be added in the future, can be used instead in addition to spline and polynomial 

    return 1.



def create_polynomial_design_matrix(time):
    """
    Construct a piecewise quadratic polynomial design matrix.

    Parameters
    ----------
    time : numpy.ndarray
        One-dimensional array of time values.

    Returns
    -------
    numpy.ndarray : 
        Piecewise polynomial design matrix with separate quadratic
        trends for each continuous time segment.
    """

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
    """
    Convert sky coordinates (RA, Dec) to detector pixel coordinates.

    Parameters
    ----------
    wcs : astropy.wcs.WCS
        World Coordinate System transformation object.
    ra : float 
        Right ascension in degrees.
    dec : float
        Declination in degrees.

    Returns
    -------
    tuple : 
        (pix_col, pix_row) pixel coordinates corresponding to the
        input sky coordinates.
    """

    #initializing skycoord object
    coord = SkyCoord(ra, dec, unit='deg')

    #converting coordinates to pixel using wcs
    pix = wcs.world_to_pixel(coord)
    pix_col = pix[0]
    pix_row = pix[1]

    return pix_col, pix_row



#better to use query_vizier_background (which is the default) as this has not been updated to account for proper motion
def query_gaia_background(wcs,ra_targ,dec_targ,radius=25):
    """
    Query Gaia sources around a target coordinate and convert them to pixels.

    Parameters
    ----------
    wcs : astropy.wcs.WCS
        World Coordinate System transformation object.
    ra_targ : float
        Target right ascension in degrees.
    dec_targ : float
        Target declination in degrees.
    radius : float, optional
        Cone search radius in arcseconds (default 25).

    Returns
    -------
    coord_bkgd : astropy.coordinates.SkyCoord
        Sky coordinates of Gaia sources within the search radius.
    pix_bkgd : tuple
        Pixel coordinates of the sources as returned by WCS
        (typically arrays of x and y positions).
    """

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



#@cache.memoize(expire=2592000) #adding a cache for the vizier query
def query_vizier_background(ra_targ,dec_targ,tpf_shape,mission,radius=60):
    """
    Query background sources from Vizier

    Parameters
    ----------
    ra_targ : float
        Target right ascension in degrees.
    dec_targ : float
        Target declination in degrees.
    tpf_shape : tuple 
        Shape of the target pixel file (e.g., (n_rows, n_cols)),
        used to estimate the search radius.
    radius : float
        Search radius in arcseconds. Used when mission is not specified 
        and pixel size is unknown (default 60).
    mission : string
        Mission name. If mission is 'tess' or 'kepler', 
        search radius is based on tpf_shape rather than radius.

    Returns
    -------
    query_result : astropy.table.Table or None
        Table of sources from Gaia DR3 within the search radius, including
        'ID', 'Gmag', 'pmra', 'pmdec', 'RA_ICRS', 'DE_ICRS'. Returns None if no sources found.
    
    """

    #initializing skycoord object using the kepler target coords
    coord = SkyCoord(ra=ra_targ, dec=dec_targ, unit=(u.degree, u.degree), frame='icrs')

    if mission == 'kepler':
        #using TPF shape to get the radius
        rad =((np.nanmax(tpf_shape) + 4) * 0.5) * 4

    elif mission == 'tess':
        #using TPF shape to get the radius
        rad =((np.nanmax(tpf_shape) + 4) * 0.5) * 21

    else:
        #since pixel size is unknown, using a set radius
        rad = radius

    #getting the catalog of stars using Vizier
    Catalog = Vizier(columns=['DR3Name','Gmag','pmRA','pmDE','RA_ICRS','DE_ICRS'],
                     column_filters={'Gmag':'<=25'},row_limit=-1) # -- for GAIA DR3
    
    #querying
    query = Catalog.query_region(coord,radius=Angle(rad, unit="arcsec"),catalog='I/355/gaiadr3')

    if len(query) == 0:
        return None
    
    res = query[0]

    return res



def propagate_query(query_result, wcs, epoch):
    """
    Propagate Gaia query coordinates to a given epoch and convert to pixel coordinates.

    Parameters
    ----------
    query_result : dict or Table
        Source data containing at least:
        - 'RA_ICRS', 'DE_ICRS' in degrees
        - 'pmRA', 'pmDE' in mas/yr
    wcs : astropy.wcs.WCS
        World Coordinate System transformation object.
    epoch : astropy.time.Time
        Observation epoch to propagate coordinates to.

    Returns
    -------
    coord_bkgd : astropy.coordinates.SkyCoord or None
        Sky coordinates of catalog sources within the search region,
        or None if no sources are found.
    pix_bkgd : tuple
        Pixel coordinates of the sources as returned by WCS
        (arrays of x and y positions).
    """
    #initializing new skycoord object with all the star coordinates and proper motions
    coord_gaia = SkyCoord(ra=query_result['RA_ICRS'], dec=query_result['DE_ICRS'], 
                          pm_ra_cosdec=query_result['pmRA'], pm_dec=query_result['pmDE'], 
                          frame='icrs', obstime=Time('J2016'))
    
    #propogating to observation time
    coord_bkgd = coord_gaia.apply_space_motion(new_obstime=epoch)

    #converting to pixel space
    pix_bkgd = wcs.all_world2pix(coord_bkgd.ra,coord_bkgd.dec, 0)

    return coord_bkgd, pix_bkgd



def propagate_coord(ra_targ, dec_targ, query_result, epoch, match_radius=2.0*u.arcsec):
    """
    Propagate coordinates to a given epoch and convert to pixel coordinates.

    Parameters
    ----------
    ra_targ : float
        Target right ascension in degrees.
    dec_targ : float
        Target declination in degrees.
    query_result : dict or Table
        Source data containing at least:
        - 'RA_ICRS', 'DE_ICRS' in degrees
        - 'pmRA', 'pmDEC' in mas/yr
    end_epoch : astropy.time.Time
        Epoch to propagate to.
    start_epoch : astropy.time.Time
        Epoch of the input coordinate. Default Time("J2000")
    match_radius : Quantity
        Maximum separation allowed for Gaia match.

    Returns
    -------
    ra_new, dec_new : float
        Propagated coordinates in degrees.
    """
    
    #initializing target SkyCoord obj
    target = SkyCoord(ra=ra_targ*u.deg, dec=dec_targ*u.deg)

    #initializing query table SkyCoord obj
    gaia_coords = SkyCoord(ra=query_result["RA_ICRS"],dec=query_result["DE_ICRS"],
                           pm_ra_cosdec=query_result["pmRA"],pm_dec=query_result["pmDE"], 
                            frame="icrs",obstime=Time("J2016"))

    #finding closest matches
    idx, sep, _ = target.match_to_catalog_sky(gaia_coords)

    #if sep > match_radius:
        #raise ValueError("No Gaia match within match_radius")

    #getting matched row from the query 
    row = query_result[idx]

    #creating new SkyCoord object with the target coord and matched proper motion
    coord = SkyCoord(ra=row["RA_ICRS"]*u.deg,dec=row["DE_ICRS"]*u.deg,
                     pm_ra_cosdec=row["pmRA"]*u.mas/u.yr,pm_dec=row["pmDE"]*u.mas/u.yr, 
                     frame="icrs",obstime=Time("J2016"))

    #propagating coordinate
    coord_new = coord.apply_space_motion(new_obstime=epoch)

    return coord_new.ra.deg, coord_new.dec.deg



def prf_residual(params,prf,data,data_err,origin,shape):
    """
    Compute the flattened, uncertainty-weighted residual for PRF fitting.

    Parameters
    ----------
    params : lmfit.Parameters
        Fit parameters containing 'amplitude', 'centerx', and 'centery'.
    prf : object
        PRF model object with an `evaluate` method.
    data : numpy.ndarray
        2D array of observed pixel values.
    data_err : numpy.ndarray
        2D array of per-pixel uncertainties.
    origin : tuple
        (row_origin, col_origin) CCD reference position.
    shape : tuple
        Shape of the image stamp used for PRF evaluation.

    Returns
    -------
    numpy.ndarray
        Flattened array of finite, uncertainty-weighted residuals.
    """

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
    res = ((amp*model - data)/data_err).ravel()

    #masking out NaNs
    mask = np.isfinite(res)

    res = res[mask]

    if len(res) == 0:
        raise ValueError('PRF residual failed: no valid pixels for this TCE.')
    
    return res



def all_prf_residual(params,prf,data,data_err,origin,shape,tces):
    """
    Compute concatenated, uncertainty-weighted residuals for multiple TCEs.

    Parameters
    ----------
    params : lmfit.Parameters
        Fit parameters containing 'centerx', 'centery', and
        'amplitude_i' for each TCE index i.
    prf : object
        PRF model object with an `evaluate` method.
    data : numpy.ndarray
        Array of 2D pixel data for each TCE (shape: n_tces × rows × cols).
    data_err : numpy.ndarray
        Array of per-pixel uncertainties for each TCE.
    origin : tuple
        (row_origin, col_origin) CCD reference position.
    shape : tuple
        Shape of the image stamp used for PRF evaluation.
    tces : sequence
        Collection of TCEs used to determine the number of amplitudes.

    Returns
    -------
    numpy.ndarray
        Flattened array of finite, concatenated residuals across all TCEs.
    """

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

        #masking out nans
        mask_i = np.isfinite(res_i)
        res_i = res_i[mask_i]

        #making sure that its not fully null
        if res_i.size == 0:
            continue

        res.append(res_i)

    if len(res) == 0:
        raise ValueError("All TCE PRF residuals failed: no valid pixels in any TCE")

    #concatenating all the residuals and returning
    return np.concatenate(res)



def quarters_prf_residual(params,prf_list,data_list,data_err_list,wcs_list,origin_list,shape_list,tce_mapping_list,valid_idxs,which_tce):
    
    #intializing residuals list
    res = []

    #extracting parameter values
    ra = params['centerra'].value
    dec = params['centerdec'].value

    #iterating through the valid quarters for this tce
    for i in valid_idxs:
        which_col = tce_mapping_list[i].get(which_tce)
        if which_col is None:
            continue

        #turning ra and dec to pixel coords for the specific quarter
        col_pix, row_pix = coords_to_pixels(wcs_list[i],ra,dec)

        #getting the ccd location for the target
        origin = origin_list[i]
        row_ccd = origin[0] + row_pix
        print(row_ccd)
        col_ccd = origin[1] + col_pix

        #evaluating the prf model
        model = prf_list[i].evaluate(targets=[(row_ccd, col_ccd)],origin=origin,shape=shape_list[i])[0, :, :]

        #extracting the amplitude param value for quarter
        amp = params[f'amplitude_{i}'].value

        #calculating the residual for quarter
        res_i = ((amp * model - data_list[i][which_col]) / data_err_list[i][which_col]).ravel()

        #masking out nans
        mask_i = np.isfinite(res_i)
        res_i = res_i[mask_i]

        #making sure that its not fully null
        if res_i.size == 0:
            continue

        res.append(res_i)

    if len(res) == 0:
        raise ValueError("Multiple Quarter PRF residual failed: no valid pixels in any quarter")

    #concatenating all the residuals and returning
    return np.concatenate(res)





