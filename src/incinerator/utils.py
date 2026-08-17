import astropy.units as u
import diskcache as dc
import numpy as np
import pandas as pd
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table
from astropy.time import Time
from astropy.units import Quantity
from astropy.wcs import WCS
from astroquery.vizier import Vizier
from lamatrix import Spline
from lmfit import Parameters

#creating cache for vizier query
cache = dc.Cache('vizier_cache')


def get_p_tdur_t0(tce: pd.Series) -> tuple[float, float, float]:
    """
    Extract the orbital period, transit duration, and transit epoch (t0)
    from a TCE NumPy array.

    Parameters
    ----------
    tce : pandas.Series
        Row from a TCE DataFrame containing the columns ``period``,
        ``tdur``, and ``t0``. ``period`` and ``tdur`` are in days,
        and ``t0`` is in BJD.

    Returns
    -------
    tuple[float, float, float]
        A tuple containing the period, transit duration, and transit
        epoch (t0), in that order.
    """
    return tce['period'], tce['tdur'], tce['t0']



def create_spline_design_matrix(time: np.ndarray,tdur: float,order: int = 3,spacing_mult: float = 3.0) -> np.ndarray:
    """
    Construct a piecewise spline design matrix.

    Parameters
    ----------
    time : numpy.ndarray
        One-dimensional array of time values.
    tdur : float
        Transit duration in the same units as ``time``. Used to set
        the spline knot spacing.
    order : int, optional
        Order of the spline. Default is 3 for cubic splines.
    spacing_mult : float, optional
        Multiplier applied to the transit duration to set the spline
        knot spacing. Default is 3.0.

    Returns
    -------
    numpy.ndarray
        Spline design matrix evaluated at the normalized time values.
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



def create_polynomial_design_matrix(time: np.ndarray) -> np.ndarray:
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
    dM2 = np.hstack([dM * m[:, None] for m in seg_masks])

    return dM2



def coords_to_pixels(wcs: WCS,ra: float,dec: float) -> tuple[float, float]:
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
    tuple[float, float]
        Pixel column and row coordinates corresponding to the input
        sky coordinates.
    """

    #initializing skycoord object
    coord = SkyCoord(ra, dec, unit='deg')

    #converting coordinates to pixel using wcs
    pix = wcs.world_to_pixel(coord)
    pix_col = pix[0]
    pix_row = pix[1]

    return pix_col, pix_row



@cache.memoize(expire=2592000) #adding a cache for the vizier query
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



def propagate_query(query_result: Table, wcs: WCS,epoch: Time) -> tuple[SkyCoord, tuple[np.ndarray, np.ndarray]]:
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
    coord_bkgd : astropy.coordinates.SkyCoord
        Sky coordinates of catalog sources propagated to the observation epoch.
    pix_bkgd : tuple[numpy.ndarray, numpy.ndarray]
        Pixel coordinates of the sources as arrays of column and row positions.
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



def propagate_coord(ra_targ: float, dec_targ: float, query_result: Table, epoch: Time, match_radius: Quantity = 2.0 * u.arcsec) -> tuple[float, float]:
    """
    Propagate a target coordinate to a given epoch using matched Gaia proper motions.

    Parameters
    ----------
    ra_targ : float
        Target right ascension in degrees.
    dec_targ : float
        Target declination in degrees.
    query_result : astropy.table.Table
        Gaia source data containing at least ``RA_ICRS``, ``DE_ICRS``,
        ``pmRA``, and ``pmDE``.
    epoch : astropy.time.Time
        Epoch to which the coordinate is propagated.
    match_radius : astropy.units.Quantity, optional
        Maximum separation allowed for the Gaia match.
        Default is 2 arcseconds.

    Returns
    -------
    tuple[float, float]
        Propagated right ascension and declination in degrees.

    Raises
    ------
    ValueError
        If no Gaia source is found within ``match_radius`` of the target.
    """
    
    #initializing target SkyCoord obj
    target = SkyCoord(ra=ra_targ*u.deg, dec=dec_targ*u.deg)

    #initializing query table SkyCoord obj
    gaia_coords = SkyCoord(ra=query_result["RA_ICRS"],dec=query_result["DE_ICRS"],
                           pm_ra_cosdec=query_result["pmRA"],pm_dec=query_result["pmDE"], 
                            frame="icrs",obstime=Time("J2016"))

    #finding closest matches
    idx, sep, _ = target.match_to_catalog_sky(gaia_coords)

    if sep > match_radius:
        raise ValueError("No Gaia match within match_radius")

    #getting matched row from the query 
    row = query_result[idx]

    #creating new SkyCoord object with the target coord and matched proper motion
    coord = SkyCoord(ra=row["RA_ICRS"]*u.deg,dec=row["DE_ICRS"]*u.deg,
                     pm_ra_cosdec=row["pmRA"]*u.mas/u.yr,pm_dec=row["pmDE"]*u.mas/u.yr, 
                     frame="icrs",obstime=Time("J2016"))

    #propagating coordinate
    coord_new = coord.apply_space_motion(new_obstime=epoch)

    return coord_new.ra.deg, coord_new.dec.deg



def prf_residual(params: Parameters, prf, data: np.ndarray, data_err: np.ndarray, origin: tuple[float, float], shape: tuple[int, int]) -> np.ndarray:
    """
    Compute the flattened, uncertainty-weighted residual for PRF fitting.

    Parameters
    ----------
    params : lmfit.Parameters
        Fit parameters containing 'amplitude', 'centerx', and 'centery'.
    prf : lkprf PRF model
        PRF model object used to evaluate the expected spatial flux
        distribution. For example, a ``KeplerPRF`` or ``TESSPRF`` object.
    data : numpy.ndarray
        2D array of observed pixel values.
    data_err : numpy.ndarray
        2D array of per-pixel uncertainties.
    origin : tuple[float, float]
        CCD reference position as ``(row_origin, col_origin)``.
    shape : tuple[int, int]
        Shape of the image stamp as ``(n_rows, n_cols)``.

    Returns
    -------
    numpy.ndarray
        Flattened array of finite, uncertainty-weighted residuals.

    Raises
    ------
    ValueError
        If no finite residuals can be calculated from the input pixels.
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



def all_prf_residual(params: Parameters, prf, data: np.ndarray, data_err: np.ndarray, origin: tuple[float, float], shape: tuple[int, int], valid_tces: np.ndarray, tce_mapping: dict[int, int]) -> np.ndarray:
    """
    Compute concatenated, uncertainty-weighted residuals for multiple TCEs.

    Parameters
    ----------
    params : lmfit.Parameters
        Fit parameters containing ``centerx``, ``centery``, and
        ``amplitude_i`` for each TCE index ``i``.
    prf : lkprf PRF model
        PRF model object used to evaluate the expected spatial flux
        distribution, such as a ``KeplerPRF`` or ``TESSPRF`` object.
    data : numpy.ndarray
        Array of 2D pixel data for each TCE with shape
        ``(n_tces, n_rows, n_cols)``.
    data_err : numpy.ndarray
        Array of per-pixel uncertainties with the same shape as ``data``.
    origin : tuple[float, float]
        CCD reference position as ``(row_origin, col_origin)``.
    shape : tuple[int, int]
        Shape of the image stamp as ``(n_rows, n_cols)``.
    valid_tces : numpy.ndarray
        Indices of TCEs that produce non-zero transit models and are
        included in the fit.
    tce_mapping : dict[int, int]
        Mapping from the original TCE index to its column index in
        ``data`` and ``data_err``.

    Returns
    -------
    numpy.ndarray
        Flattened array of finite, concatenated residuals across all valid TCEs.

    Raises
    ------
    ValueError
        If no valid residuals can be calculated for any TCE.
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
    for i in valid_tces:
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



def quarters_prf_residual(params: Parameters, prf_list: list, data_list: list[np.ndarray], data_err_list: list[np.ndarray], wcs_list: list[WCS], origin_list: list[tuple[float, float]], shape_list: list[tuple[int, int]], tce_mapping_list: list[dict[int, int]], valid_idxs: list[int], which_tce: int) -> np.ndarray:
    """
    Compute concatenated, uncertainty-weighted residuals for a TCE
    across multiple quarters or sectors.

    Parameters
    ----------
    params : lmfit.Parameters
        Fit parameters containing ``centerra``, ``centerdec``, and
        ``amplitude_i`` for each valid quarter or sector.
    prf_list : list
        List of lkprf PRF model objects, one for each observation.
    data_list : list[numpy.ndarray]
        List of transit depth maps for each observation.
    data_err_list : list[numpy.ndarray]
        List of per-pixel uncertainties corresponding to ``data_list``.
    wcs_list : list[astropy.wcs.WCS]
        List of WCS objects used to convert sky coordinates to pixel
        coordinates for each observation.
    origin_list : list[tuple[float, float]]
        List of CCD reference positions as
        ``(row_origin, col_origin)`` for each observation.
    shape_list : list[tuple[int, int]]
        List of image stamp shapes as ``(n_rows, n_cols)``.
    tce_mapping_list : list[dict[int, int]]
        List of mappings from TCE indices to their corresponding
        columns in the transit depth maps.
    valid_idxs : list[int]
        Indices of observations containing the TCE being fitted.
    which_tce : int
        Index of the TCE being fitted.

    Returns
    -------
    numpy.ndarray
        Flattened array of finite, concatenated residuals across all
        valid observations.

    Raises
    ------
    ValueError
        If no valid residuals can be calculated for any observation.
    """
    
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