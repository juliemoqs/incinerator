import batman
import lkprf
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from lmfit import Parameters, minimize
from lmfit.models import Gaussian2dModel
from matplotlib.gridspec import GridSpec
from scipy.spatial.distance import mahalanobis

from .utils import (
    all_prf_residual,
    coords_to_pixels,
    create_polynomial_design_matrix,
    create_spline_design_matrix,
    get_p_tdur_t0,
    prf_residual,
    propagate_coord,
    propagate_query,
    quarters_prf_residual,
    query_vizier_background,
)


#initializing Localize Class -- single quarter localization
class Localize(object):

    #initiliazing -- data parser and loader
    def __init__(self, time: np.ndarray, flux: np.ndarray,
             flux_err: np.ndarray, tces: pd.DataFrame, id: str,
             ra_targ: float,
             dec_targ: float,
             wcs: WCS, epoch: Time,
             mission: str | None = None,
             ra_bonus: float | None = None,
             dec_bonus: float | None = None,
             prf: object | None = None,
             headers: list[fits.Header] | None = None,
             **mission_keywords) -> None:
        """
        Initialize the Localize object.

        Parameters
        ----------
        time : numpy.ndarray
            1D array of time values in BJD.
        flux : numpy.ndarray
            3D flux cube with shape (n_time, n_rows, n_cols).
        flux_err : numpy.ndarray
            3D array of per-pixel uncertainties matching `flux`.
        tces : pandas.DataFrame
            DataFrame containing Threshold Crossing Event (TCE) information.
            Must contain the columns `period`, `t0`, and `tdur`. The `period`
            and `tdur` values must be in days, and `t0` must be in BJD.
        id : str
            Target identifier. Do not include the mission prefix (TIC/KIC/etc.).
        ra_targ : float
            J2000 right ascension of the target, in degrees.
        dec_targ : float
            J2000 declination of the target, in degrees.
        wcs : astropy.wcs.WCS
            WCS solution for pixel-to-sky transformations.
        epoch : astropy.time.Time
            Epoch of the input coordinates.
        mission : str, optional
            Mission name (e.g., 'kepler', 'tess').
        ra_bonus : float, optional
            J2000 right ascension of a source different from the nominal
            mission target, in degrees.
        dec_bonus : float, optional
            J2000 declination of a source different from the nominal
            mission target, in degrees.
        prf : object, optional
            Precomputed PRF model.
        headers : list of astropy.io.fits.Header, optional
            FITS headers from the original FITS file. Each entry corresponds
            to the header of a single HDU.
        **mission_keywords
            Additional mission-specific keyword arguments stored directly in
            the instance dictionary. Common entries include sector/quarter,
            channel, CCD, camera, and other instrument metadata required for
            PRF evaluation.

        Raises
        ------
        ValueError
            If `flux` is not 3D, if the length of `time` does not match the
            first axis of `flux`, or if `tces` does not contain the required
            columns.
        TypeError
            If `tces` is not a pandas DataFrame.
        """
        
        self.time = time
        self.flux = flux
        self.flux_err = flux_err
        self.tces = tces
        self.wcs = wcs
        self.mission = mission
        self.ra_bonus = ra_bonus
        self.dec_bonus = dec_bonus
        self.ra_targ = ra_targ
        self.dec_targ = dec_targ
        self.prf = prf
        self.id = id
        self.headers = headers
        self.epoch = epoch

        self.origin_row: float | None = None
        self.origin_col: float | None = None
        self.quarter: int | None = None
        self.sector: int | None = None
        self.weights: np.ndarray | None = None
        self.weights_err: np.ndarray | None = None

        self.__dict__.update(mission_keywords)

        #checking data types and shapes
        if not isinstance(self.tces, pd.DataFrame):
            raise TypeError("tces must be a pandas DataFrame")

        if self.flux.ndim != 3:
            raise ValueError("Flux must be a 3D array: (time, y, x)")

        if self.flux.shape[0] != self.time.size:
            raise ValueError("Time length must match flux first axis")

        required_tce_columns = {"period", "t0", "tdur"}
        if not required_tce_columns.issubset(self.tces.columns):
            raise ValueError(f"tces must contain the columns: {required_tce_columns}")

        self.shape = self.flux.shape[1:]

        #search for nearby stars and propagate coordinates
        self._search_and_propagate()

        #clean the data
        self._clean_data()


    
    @staticmethod
    def from_tpf_info(filename: str, tces: pd.DataFrame, id: str,
                  mission: str, ra_bonus: float | None = None,
                  dec_bonus: float | None = None) -> "Localize":
        """
        Create a Localize object from a Kepler or TESS Target Pixel File (TPF).
    
        Parameters
        ----------
        filename : str
            Path to the TPF FITS file.
        tces : pandas.DataFrame
            DataFrame containing Threshold Crossing Event (TCE) information
            associated with the target.
        id : str
            Identifier for the target.
        mission : str
            Name of the mission ('kepler' or 'tess').
        ra_bonus : float, optional
            J2000 right ascension of a source different from the nominal
            mission target, in degrees.
        dec_bonus : float, optional
            J2000 declination of a source different from the nominal
            mission target, in degrees.
    
        Returns
        -------
        Localize
            Initialized Localize object containing the data and metadata
            extracted from the TPF.
        """
        
        mission = mission.lower()
        if mission not in {"kepler", "tess"}:
            raise ValueError("mission must be 'kepler' or 'tess'")

        #reads in the file 
        with fits.open(filename) as hdulist:

            time = hdulist[1].data['time']
            flux = hdulist[1].data['flux']
            flux_err = hdulist[1].data['flux_err']

            headers = [hdu.header for hdu in hdulist]

            #gets all the info we need whoop whoop
            ra_targ= hdulist[0].header['ra_obj']
            dec_targ = hdulist[0].header['dec_obj']
            wcs = WCS(hdulist[2].header)
            origin_row = hdulist[1].header['2CRV4P']
            origin_col = hdulist[1].header['1CRV4P']

            if mission == 'kepler':
                channel = hdulist[0].header['channel']
                quarter = hdulist[0].header['quarter']

                time_offset = 2454833 

                #initializing prf for the specific channel
                prf = lkprf.KeplerPRF(channel = channel)

                mission_keywords = {'channel':channel,
                                        'origin_row':origin_row,
                                        'origin_col':origin_col,
                                        'quarter':quarter}

            elif mission == 'tess':
                camera = hdulist[0].header['camera']
                ccd = hdulist[0].header['ccd']
                sector = hdulist[0].header['sector']

                time_offset = 2457000

                #initializing prf for the specific camera and ccd chip 
                prf = lkprf.TESSPRF(camera = camera, ccd = ccd)

                mission_keywords = {'camera':camera,'ccd':ccd,
                                    'origin_row':origin_row,
                                    'origin_col':origin_col,
                                    'sector':sector}


        #getting the epoch of the observation
        epoch = Time(np.nanmin(time) + time_offset, format='jd')

        #converting time to bjd
        time = time + time_offset 

        return Localize(time, flux, flux_err, tces, id=id, wcs=wcs, mission=mission,
                        ra_bonus=ra_bonus, dec_bonus=dec_bonus, ra_targ=ra_targ, dec_targ=dec_targ,
                        prf=prf, epoch=epoch, headers=headers,**mission_keywords)



    def _clean_data(self) -> None:
            """
            Remove invalid cadences and pixels from the flux cube.
    
            Attributes
            ----------
            good_pix_mask : numpy.ndarray
                2D boolean mask identifying pixels with nonzero total flux.
            good_cad_mask : numpy.ndarray
                1D boolean mask identifying cadences with finite time, flux,
                and flux uncertainties.
            pix : numpy.ndarray
                2D cleaned flux array with shape
                (n_good_cadences, n_good_pixels).
            pix_err : numpy.ndarray
                2D cleaned uncertainty array matching `pix`.
            time : numpy.ndarray
                Filtered time array containing only good cadences.
            """
            
            flux = self.flux
            time = self.time
            flux_err = self.flux_err
    
            #good pixels: exclude NaNs and zeros
            good_pix = np.nansum(flux, axis=0) != 0
            #good cadences: time, flux, and flux_err is finite
            good_cad = np.isfinite(time) & np.isfinite(flux[:, good_pix]).all(axis=(1)) & np.isfinite(flux_err[:, good_pix]).all(axis=(1))
    
            #save masks for later use
            self.good_pix_mask = good_pix
            self.good_cad_mask = good_cad
    
            #extract the good cadences and pixels
            self.pix = flux[good_cad][:, good_pix]
            self.pix_err = flux_err[good_cad][:, good_pix]
            self.time = time[good_cad]



    def _search_and_propagate(self) -> None:
        """
        Query Gaia around the target and propagate coordinates to the
        observation epoch.
    
        Attributes
        ----------
        gaia_sources : astropy.table.Table
            Gaia sources returned from the cone search.
        gaia_coord : astropy.coordinates.SkyCoord
            Gaia source coordinates propagated to the observation epoch.
        gaia_pix : numpy.ndarray
            Pixel coordinates of the Gaia sources.
        ra_targ, dec_targ : float
            Target coordinates in degrees propagated to the observation
            epoch.
        ra_bonus, dec_bonus : float
            Bonus source coordinates in degrees propagated to the
            observation epoch, if provided.
        """

        #doing a cone search by querying gaia dr3 with vizier -- search is cached
        self.gaia_sources = query_vizier_background(self.ra_targ,self.dec_targ,self.shape,mission=self.mission)

        #propagating cone search coordinates to observation epoch and converting into pixels
        self.gaia_coord, self.gaia_pix = propagate_query(self.gaia_sources, self.wcs, self.epoch)

        #propagating target coordinates that were either extracted from header or given as inputs in initialization
        self.ra_targ, self.dec_targ = propagate_coord(self.ra_targ, self.dec_targ, self.gaia_sources, self.epoch)

        if self.ra_bonus is not None and self.dec_bonus is not None:
            #propagating bonus coordinates if given as inputs in initialization
            self.ra_bonus, self.dec_bonus = propagate_coord(self.ra_bonus, self.dec_bonus, self.gaia_sources, self.epoch)



    def _get_transit(self, which_tce: int, method: str = "box", tdur_frac: float = 3.0, **batman_kws) -> np.ndarray:
        """
        Generate a transit model mask for a specified TCE.

        Parameters
        ----------
        which_tce : int
            Index of the TCE in `self.tces`.
        method : str, optional
            Method used to generate the transit model. Options are
            ``"box"`` or ``"batman"``. Default is ``"box"``.
        tdur_frac : float, optional
            Factor used to determine the width of the box-shaped transit
            mask. The transit mask extends to ``tdur / tdur_frac`` from
            the transit center. Default is 3.0.
        **batman_kws
            Additional parameters required to generate the ``batman`` transit
            model. Required keywords are ``rp``, ``a``, ``inc``, ``ecc``,
            ``w``, ``u``, and ``limb_dark``.

        Returns
        -------
        numpy.ndarray
            1D transit model evaluated at the object's time array. For
            the box model, values are -1 during transit and 0 elsewhere.
            For the ``batman`` model, the transit flux is shifted so that the
            out-of-transit baseline is approximately 0.

        Raises
        ------
        ValueError
            If `method` is not ``"box"`` or ``"batman"``.
        """
        
        time = self.time
        period,tdur,t0 = get_p_tdur_t0(self.tces.iloc[which_tce])

        #fold lc
        ph = ((time - t0)) % period

        #shift phase to make it go from -P/2 to P/2 instead of 0 to P
        ph[ph > (0.5 * period)] -= period
        
        if method == 'box':

            #creating transit mask
            transit =  - (np.abs(ph) < (tdur/tdur_frac)).astype(float) #0 out of transit, -1 in transit

        elif method == 'batman':

            #create object to store parameters
            params = batman.TransitParams()

            #set up parameters
            params.t0 = t0 
            params.per = period
            params.rp = batman_kws['rp']
            params.a = batman_kws['a']
            params.inc = batman_kws['inc']
            params.ecc = batman_kws['ecc']
            params.w = batman_kws['w']
            params.u = batman_kws['u']
            params.limb_dark = batman_kws['limb_dark']
        
            #build model
            m = batman.TransitModel(params,time,supersample_factor=5,exp_time=np.nanmedian(np.diff(time)))

            flux = m.light_curve(params)

            #shift out-of-transit baseline to 0 for the design matrix
            transit = flux - np.median(flux)

        else:
            raise ValueError("method must be 'box' or 'batman'")

        return transit



    def build_design_matrix(self, order: int = 3, spacing_mult: float = 3.0,
                        tdur_frac: float = 3.0, transit_method: str = "box",
                        batman_kws: list[dict] | None = None) -> None:
        """
        Construct the full design matrix.

        Parameters
        ----------
        order : int, optional
            Order of the spline used to model short-term variability.
            Default is 3.
        spacing_mult : float, optional
            Multiplier applied to the maximum transit duration to set the
            spline knot spacing. Default is 3.0.
        tdur_frac : float, optional
            Factor used to determine the width of the box-shaped transit
            model. The transit mask extends to ``tdur / tdur_frac`` from
            the transit center. Default is 3.0.
        transit_method : str, optional
            Method used to generate the transit component. Options are
            ``"box"`` and ``"batman"``. Default is ``"box"``.
        batman_kws : list of dict, optional
            List of dictionaries containing additional arguments passed to
            the ``batman`` transit model. One dictionary is required for
            each TCE when ``transit_method="batman"``.

        Raises
        ------
        ValueError
            If `transit_method` is not ``"box"`` or ``"batman"``, or if
            `transit_method="batman"` and `batman_kws` is not provided.

        Attributes
        ----------
        valid_tces : numpy.ndarray
            Indices of TCEs in `self.tces` that produce non-zero transit
            models and are included as columns in the design matrix.
        tce_mapping : dict
            Mapping from TCE indices in the `self.tces` DataFrame to their
            corresponding column indices in the design matrix.
        design_matrix : numpy.ndarray or None
            Final design matrix combining transit, polynomial, and spline
            components. Set to None if no valid transits are found.
    
        """
        
        time = self.time

        tdur_list = []

        transit_cols = []
        valid_tce_indices = []
        
        #iterating through all the tces
        for i in range(len(self.tces)):
            period,tdur,t0 = get_p_tdur_t0(self.tces.iloc[i])
            tdur_list.append(tdur)
            if transit_method == 'batman':
                if batman_kws is None:
                    raise ValueError("batman_kws must be provided when transit_method='batman'")
                else:
                    batman_kws_i = batman_kws[i]
                
                #get transit
                transit = self._get_transit(i,method=transit_method,**batman_kws_i)

            else:
                #get transit
                transit = self._get_transit(i,method=transit_method,tdur_frac=tdur_frac)

            #making sure the tce is actually transiting during this quarter
            if not np.allclose(transit, 0, atol=1e-10):
                #appending the transit if its not all 0s
                transit_cols.append(transit[:, None])
                #appending the tce index to keep track of valid tces for each quarter
                valid_tce_indices.append(i)

        if not transit_cols:
            self.valid_tces = np.array([], dtype=int)
            self.tce_mapping = {}
            self.design_matrix = None
            return

        transits = np.hstack(transit_cols)

        tdur_long = np.nanmax(tdur_list)

        #get the piecewise spline portion -- short-term variability
        spline = create_spline_design_matrix(time,tdur_long,order,spacing_mult)

        #get piecewise polynomial portion -- long-term variability
        polynomial = create_polynomial_design_matrix(time)

        #combine piecewise trends, spline, and transit signal into a final design matrix
        dM = np.hstack([transits, polynomial, spline])

        self.valid_tces = np.array(valid_tce_indices, dtype=int)
        self.tce_mapping = {tce_idx: col_idx for col_idx, tce_idx in enumerate(self.valid_tces)}
        self.design_matrix = dM
    

    
    def _solve_weights(self) -> None:
        """
        Solve for per-pixel linear model weights and uncertainties.

        Attributes
        ----------
        weights : numpy.ndarray or None
            Array of fitted weights with shape (n_pix, n_params).
            Set to None if no design matrix is available.
        weights_err : numpy.ndarray or None
            Array of 1-sigma uncertainties on the fitted weights with
            shape (n_pix, n_params). Set to None if no design matrix is
            available.
        """
        
        pix = self.pix
        pix_err= self.pix_err
        dM = self.design_matrix

        if dM is None:
            self.weights = None
            self.weights_err = None
            return

        n_pix = pix.shape[1]
        n_params = dM.shape[1]
    
        self.weights = np.zeros((n_pix, n_params))
        self.weights_err = np.zeros((n_pix, n_params))

        #solving for the weights
        for idx, y, e in zip(range(n_pix), pix.T, pix_err.T):
            sigma_w_inv = dM.T.dot(dM / e[:, None]**2)
            B = dM.T.dot(y / e**2)
            self.weights[idx] = np.linalg.solve(sigma_w_inv, B)
            self.weights_err[idx] = np.sqrt(np.diag(np.linalg.inv(sigma_w_inv)))        



    def _solve_transit(self, which_tce: int) -> tuple[np.ndarray, np.ndarray]:
        """
        Extract the spatial transit depth solution for a given TCE.
    
        Parameters
        ----------
        which_tce : int
            Index of the TCE corresponding to the transit component in
            the design matrix.
    
        Returns
        -------
        transit_weight : numpy.ndarray
            2D array of fitted transit depths per pixel.
        transit_weight_err : numpy.ndarray
            2D array of 1-sigma uncertainties on the transit depths.
        """

        self._solve_weights()
        good_pix = self.good_pix_mask

        if self.weights is None or self.weights_err is None:
            raise RuntimeError("Weights were not initialized.")

        #getting the weights for the transit depth
        transit_weight, transit_weight_err = np.zeros((2, *self.flux.shape[1:])) * np.nan
        transit_weight[good_pix] = self.weights[:, which_tce]
        transit_weight_err[good_pix] = self.weights_err[:, which_tce]

        return transit_weight, transit_weight_err
    


    def solve_transit_weights(self) -> None:
        """
        Solve for the spatial transit depth maps for all TCEs.

        Attributes
        ----
        self.transit_weights : list of numpy.ndarray
            List of 2D transit depth maps, one per TCE.
        self.transit_weights_err : list of numpy.ndarray
            List of 2D 1-sigma uncertainty maps corresponding to
            each transit depth map.
        """

        amps = []
        amps_err = []
        for i in self.valid_tces:
            tce_col = self.tce_mapping[i]
            amp_i, amp_err_i = self._solve_transit(tce_col)
            amps.append(amp_i)
            amps_err.append(amp_err_i)

        self.transit_weights = amps
        self.transit_weights_err = amps_err



    def fit_to_heatmap(self,model_func=None,which_tce: int | None = None,method: str = "prf",initial_loc: tuple[float, float] | None = None,**fit_kw) -> tuple:
        """
        Fit a model to the per-pixel transit depth heatmap.

        Parameters
        ----------
        model_func : callable, optional
            Custom residual function used when ``method="custom"``.
        which_tce : int, optional
            Index of the TCE to fit. For ``method="prf"``, if None, a joint
            PRF fit is performed across all valid TCEs.
        method : str, optional
            Fitting method. Must be ``"prf"``, ``"2dgaussian"``, or
            ``"custom"``.
        initial_loc : tuple of float, optional
            Initial guess for the centroid position as ``(col, row)`` in pixel
            coordinates. If None, the position is initialized using ``ra_bonus``
            and ``dec_bonus`` if provided, otherwise ``ra_targ`` and ``dec_targ``.
        **fit_kw
            Additional keyword arguments passed to the custom model.

        Returns
        -------
        result : lmfit.model.ModelResult or lmfit.minimizer.MinimizerResult
            Full fit result object.
        fit_metrics : dict
            Dictionary containing the best-fit centroid position and
            corresponding 1-sigma uncertainties. Keys are ``"centerx"``,
            ``"centery"``, ``"sigmax"``, and ``"sigmay"``.

        Raises
        ------
        ValueError
            If an unsupported fitting method is provided, if ``which_tce``
            is required but not provided, or if a required argument for
            the selected fitting method is missing.
        KeyError
            If ``which_tce`` does not correspond to a valid TCE for the
            current quarter or sector.
        RuntimeError
            If an unexpected error occurs during fitting.
        """
        
        amp = self.transit_weights.copy()
        amp_err = self.transit_weights_err.copy()

        #initializing metrics dictionary
        fit_metrics = {}

        #getting pix coord
        if initial_loc is not None:
            pix_col, pix_row = initial_loc
        elif self.ra_bonus is not None and self.dec_bonus is not None:
            pix_col, pix_row = coords_to_pixels(self.wcs, self.ra_bonus, self.dec_bonus)
        else:
            pix_col, pix_row = coords_to_pixels(self.wcs, self.ra_targ, self.dec_targ)

        try:

            if method == 'prf':

                params = Parameters()
                #constraining params
                params.add('centery', value=pix_row, min=0, max=self.shape[0]-1)
                params.add('centerx', value=pix_col, min=0, max=self.shape[1]-1)

                #doing the fit for only one tce, if a tce number is given 
                if which_tce is not None:

                    which_col = self.tce_mapping.get(which_tce)
                    if which_col is None:
                        raise KeyError('This TCE does not transit in this quarter/sector.')

                    #constraining the amplitude param
                    params.add('amplitude',value=np.nanmax(amp[which_col]),min=0)

                    if self.origin_row is None or self.origin_col is None:
                        raise ValueError("origin_row and origin_col must be provided for a PRF fit.")

                    #inputs for prf_residual function
                    prf_fit_kw = {'prf':self.prf,'data':amp[which_col],'data_err':amp_err[which_col],
                                  'origin':(self.origin_row,self.origin_col),'shape':self.shape}

                    #minimizing to get the best fit results using prf_residual
                    result = minimize(prf_residual,params,kws=prf_fit_kw)

                else:
                    #constraining the amplitude param for all tces
                    for i in self.valid_tces:
                        col_i = self.tce_mapping[i]
                        params.add(f'amplitude_{i}',value=np.nanmax(amp[col_i]),min=0)

                    #inputs for all_prf_residual function
                    prf_fit_kw = {'prf':self.prf,'data':amp,'data_err':amp_err,
                                  'origin':(self.origin_row,self.origin_col),'shape':self.shape,
                                  'valid_tces':self.valid_tces, 'tce_mapping': self.tce_mapping}

                    #minimizing to get the best fit results using prf_residual
                    result = minimize(all_prf_residual,params,kws=prf_fit_kw)


                #extracing most important fit metrics 
                fit_metrics['centerx'] = result.params['centerx'].value
                fit_metrics['sigmax'] = result.params['centerx'].stderr
                fit_metrics['centery'] = result.params['centery'].value
                fit_metrics['sigmay'] = result.params['centery'].stderr


            elif method == 'custom':
                print('i dont know if this will actually work so good luck :)')

                if which_tce is None:
                    raise ValueError("which_tce must be provided.")

                which_col = self.tce_mapping.get(which_tce)
                
                if which_col is None:
                    raise KeyError("This TCE does not transit in this quarter/sector.")

                params = Parameters()
                #constraining params
                params.add('centery', value=pix_row, min=0, max=self.shape[0]-1)
                params.add('centerx', value=pix_col, min=0, max=self.shape[1]-1)
                params.add('amplitude',value=np.nanmax(amp[which_col]),min=0)

                #minimizing to get the best fit results
                result = minimize(model_func,params,kws=fit_kw)#,grid=False)

                #extracting most important fit metrics
                fit_metrics['centerx'] = result.params['centerx'].value
                fit_metrics['sigmax'] = result.params['centerx'].stderr
                fit_metrics['centery'] = result.params['centery'].value
                fit_metrics['sigmay'] = result.params['centery'].stderr


            elif method == '2dgaussian':
                if which_tce is None:
                    raise ValueError("which_tce must be provided for a 2D Gaussian fit.")

                which_col = self.tce_mapping.get(which_tce)

                if which_col is None:
                    raise KeyError("This TCE does not transit in this quarter/sector.")

                #need to make sure there are no NaNs
                mask = np.isfinite(amp[which_col])

                #create rows, cols, and amplitude
                rows,cols = np.indices(amp[which_col].shape)
                x = cols[mask].ravel()
                y = rows[mask].ravel()
                z = amp[which_col][mask].ravel()

                #using lmfit's 2d Gaussian model -- this is hardcoded
                model_func = Gaussian2dModel()
                params = model_func.guess(z,x=x,y=y)

                #constraining params
                params['amplitude'].set(min=None, max=None)
                params['sigmax'].set(min=0.5)
                params['sigmay'].set(min=.5)

                #fitting model to get results
                result = model_func.fit(z,params,x=x,y=y)

                #extracting most important fit metrics 
                fit_metrics['centerx'] = result.params['centerx'].value
                fit_metrics['sigmax'] = result.params['centerx'].stderr
                fit_metrics['centery'] = result.params['centery'].value
                fit_metrics['sigmay'] = result.params['centery'].stderr

            else:
                raise ValueError("method should be 'prf', '2dgaussian', or 'custom'.")
            
        except ValueError as e:
            if which_tce is not None:
                raise ValueError(f"Fit failed for TCE {which_tce} using method '{method}': {e}") from e
            else:
                raise ValueError(f"Fit failed for all TCEs using method '{method}': {e}") from e
        except Exception as e:
            if which_tce is not None:
                raise RuntimeError(f"Unexpected error while fitting TCE {which_tce} using method '{method}': {e}") from e
            else:
                raise RuntimeError(f"Unexpected error while fitting all TCEs using method '{method}': {e}") from e

        #returns the whole report and the most important fit metrics
        return result, fit_metrics



    def get_offset(self, fit_metrics: dict) -> tuple[float, float, float, float]:
        """
        Compute the centroid offset in units of sigma.

        Parameters
        ----------
        fit_metrics : dict
            Dictionary containing the fitted centroid positions and their
            uncertainties. Must contain ``"centerx"``, ``"centery"``,
            ``"sigmax"``, and ``"sigmay"``.
        Returns
        -------
        nsigma : float
            Radial offset between the fitted centroid and the expected source
            position, expressed in units of the propagated positional
            uncertainty.
        dist : float
            Euclidean distance between the fitted centroid and the expected
            source position in pixels.
        sigma_dist : float
            Propagated 1-sigma uncertainty on the radial distance in pixels.
        mah_dist : float
            Mahalanobis distance between the fitted centroid and the expected
            source position.

        Raises
        ------
        KeyError
            If a required key is missing from ``fit_metrics``.
        """
        #extract fit metrics
        x0 = fit_metrics['centerx']
        y0 = fit_metrics['centery']
        x0_err = fit_metrics['sigmax']
        y0_err = fit_metrics['sigmay']
        
        #getting the pixel info for the given ra and dec
        if self.ra_bonus is not None and self.dec_bonus is not None:
            pix_col, pix_row = coords_to_pixels(self.wcs,self.ra_bonus,self.dec_bonus)
        else:
            pix_col, pix_row = coords_to_pixels(self.wcs,self.ra_targ,self.dec_targ)

        #calculating differences
        dx = x0 - pix_col
        dy = y0 - pix_row

        dist = (dx**2+dy**2)**(1/2)
        sigma_dist = 1/dist * (dx**2 * x0_err**2 + dy**2 * y0_err**2)**(1/2)

        cov = np.array([[x0_err**2,0],[0,y0_err**2]])
        inv_cov = np.linalg.inv(cov)
        x = np.array([x0,y0])
        y = np.array([pix_col,pix_row])
        mah_dist = mahalanobis(x,y,inv_cov)

        #finding how many sigma away the center of the fit is from the actual location
        nsigma = dist/sigma_dist

        return nsigma, dist, sigma_dist, mah_dist



    def plot_heatmap(self,metrics: dict,which_tce: int,savefig: bool = False,save_directory: str = ".") -> None:
        """
        Generate a localization report plot for a given TCE.

        Parameters
        ----------
        metrics : dict
            Dictionary containing the fitted centroid positions and their
            uncertainties. Must contain ``"centerx"``, ``"centery"``,
            ``"sigmax"``, and ``"sigmay"``.
        which_tce : int
            Index of the TCE in ``self.tces`` to visualize.
        savefig : bool, optional
            If True, save the figure to disk. Default is False.
        save_directory : str, optional
            Directory where the figure will be saved. Default is the
            current directory.

        Raises
        ------
        KeyError
            If ``which_tce`` does not correspond to a valid TCE for the
            current quarter or sector.
        """

        which_col = self.tce_mapping.get(which_tce)
        if which_col is None:
            raise KeyError('This TCE does not transit in this quarter/sector.')

        #extracting metrics
        x0 = metrics['centerx']
        x0_err = metrics['sigmax']
        y0 = metrics['centery']
        y0_err = metrics['sigmay']

        #getting the pixel info for the given ra and dec
        if self.ra_bonus is not None and self.dec_bonus is not None:
            pix_col, pix_row = coords_to_pixels(self.wcs,self.ra_bonus,self.dec_bonus)
        else:
            pix_col, pix_row = coords_to_pixels(self.wcs,self.ra_targ,self.dec_targ)

        #getting offset
        nsigma,dist,sigma_dist,mah_dist = self.get_offset(metrics)

        fig, axes = plt.subplots(1,3, figsize=(16,4), layout='compressed', 
                         gridspec_kw={'width_ratios': [1.2, 1.2, 0.8]})

        axes[1].plot(pix_col,pix_row, 'rx')

        axes[1].errorbar(x0,y0,xerr=x0_err,yerr=y0_err,fmt='o',c='k',markersize=4,capsize=1.5)

        cax=axes[1].imshow(self.transit_weights[which_col]/self.transit_weights_err[which_col], cmap='inferno',origin='lower')
        axes[1].scatter(self.gaia_pix[0],self.gaia_pix[1],c='w',edgecolors='k',s=60)
        axes[1].tick_params(axis='both', labelsize=13)
        axes[1].set_xlabel('Pixel Number', fontsize=14)
        axes[1].set_ylabel('Pixel Number', fontsize=14)
        axes[1].set_xlim(-0.5,self.shape[1]-0.5)
        axes[1].set_ylim(-0.5,self.shape[0]-0.5)
        
        cbar1 = fig.colorbar(cax, ax=axes[1], fraction=0.046, pad=0.04)
        cbar1.ax.tick_params(labelsize=13)
        cbar1.set_label('Depth / Depth Error', fontsize=16)


        cax1 = axes[0].imshow(np.nanmean(self.flux, axis=0),origin='lower')
        axes[0].plot(pix_col,pix_row, 'rx')

        axes[0].scatter(self.gaia_pix[0],self.gaia_pix[1],c='w',edgecolors='k',s=60)
        axes[0].tick_params(axis='both', labelsize=13)
        axes[0].set_xlabel('Pixel Number', fontsize=14)
        axes[0].set_ylabel('Pixel Number', fontsize=14)
        axes[0].set_xlim(-0.5,self.shape[1]-0.5)
        axes[0].set_ylim(-0.5,self.shape[0]-0.5)
        
        if self.shape[0] < 6 or self.shape[1] < 6:
            axes[1].set_xticks(np.arange(0,self.shape[1],1))
            axes[1].set_yticks(np.arange(0,self.shape[0],1))
            axes[0].set_xticks(np.arange(0,self.shape[1],1))
            axes[0].set_yticks(np.arange(0,self.shape[0],1))
        else:
            axes[1].set_xticks(np.arange(0,self.shape[1],2))
            axes[1].set_yticks(np.arange(0,self.shape[0],2))
            axes[0].set_xticks(np.arange(0,self.shape[1],2))
            axes[0].set_yticks(np.arange(0,self.shape[0],2))
        
        cbar2 = fig.colorbar(cax1, ax=axes[0], fraction=0.046, pad=0.04)
        cbar2.ax.tick_params(labelsize=13)
        cbar2.set_label('Flux [e-/sec]', fontsize=16)

        #side panel with offset info
        axes[2].axis('off')

        report_text = (f'Euclidean Dist:\n {dist:.4g} pixels\n\n'
                       f'Euclidean Dist Error:\n {sigma_dist:.4g} pixels\n\n'
                       f'Offset:\n {nsigma:.4g}'+r'$\sigma$'+f'\n\n'
                       f'Mahalanobis Dist:\n {mah_dist:.4g}'+r'$\sigma$')
        
        axes[2].text(0, 0.5, report_text,fontsize=14,verticalalignment='center')

        
        if self.mission=='kepler':
            cat_name='KIC'
            obs_label=f'Quarter {self.quarter}'
            obs_prefix=f"{obs_label.replace(' ', '').lower()}_"
        elif self.mission=='tess':
            cat_name='TIC'
            obs_label=f'Sector {self.sector}'
            obs_prefix=f"{obs_label.replace(' ', '').lower()}_"
        else:
            cat_name='CANDIDATE'
            obs_label=''
            obs_prefix=''

        
        plt.suptitle(f'{obs_label} InCINERATOr Report for {cat_name} {self.id}.0{which_tce+1}', fontsize=18)

        if savefig:
            savefile = f'{obs_prefix}incinerator_report_{cat_name}{self.id}.0{which_tce+1}.png'

            plt.savefig(f'{save_directory}/{savefile}', dpi=200, pad_inches=0.5, bbox_inches='tight')


    
#initializing MultiLocalize Class -- multi-quarter joint fit built from Localize objects
class MultiLocalize(object):
    
    #initializing -- loading in Localize objects
    def __init__(self, loc_object_list: list[Localize]) -> None:
        """
        Initialize a collection of Localize objects.

        Parameters
        ----------
        loc_object_list : list[Localize]
            List of Localize objects corresponding to different quarters
            or sectors. Objects with no valid TCEs are excluded.

        Raises
        ------
        ValueError
            If no Localize objects contain valid TCEs.
        """

        #filtering out quarters that have no transits at all
        valid_loc_list = [loc for loc in loc_object_list if len(loc.valid_tces) > 0]
        if len(valid_loc_list) == 0:
            raise ValueError("No Localize objects contain valid TCEs.")

        #extracting information from each Localize object given 
        self.loc_object_list = valid_loc_list
        self.n_observations = len(valid_loc_list)
        self.tces = valid_loc_list[0].tces
        self.ra_targ = valid_loc_list[0].ra_targ
        self.dec_targ = valid_loc_list[0].dec_targ
        self.ra_bonus = valid_loc_list[0].ra_bonus
        self.dec_bonus = valid_loc_list[0].dec_bonus
        self.gaia_coord = valid_loc_list[0].gaia_coord
        self.mission = valid_loc_list[0].mission
        self.id = valid_loc_list[0].id
        self.shape_list = [loc.shape for loc in valid_loc_list]
        self.origin_list = [(loc.origin_row, loc.origin_col) for loc in valid_loc_list]
        self.wcs_list = [loc.wcs for loc in valid_loc_list]
        self.transit_weights_list = [loc.transit_weights for loc in valid_loc_list]
        self.transit_weights_err_list = [loc.transit_weights_err for loc in valid_loc_list]
        self.prf_list = [loc.prf for loc in valid_loc_list]
        self.valid_tces_list = [loc.valid_tces for loc in valid_loc_list]
        self.tce_mapping_list = [loc.tce_mapping for loc in valid_loc_list]

        if self.mission == 'kepler':
            self.quarter_list = [loc.quarter for loc in valid_loc_list]

        if self.mission == 'tess':
            self.sector_list = [loc.sector for loc in valid_loc_list]

    

    def fit_to_quarters(self, which_tce: int | None = None, initial_location: tuple[float, float] | None = None) -> tuple:
        """
        Fit a joint PRF model to the transit depth maps across observations.

        Parameters
        ----------
        which_tce : int
            Index of the TCE in ``self.tces`` to fit.
        initial_location : tuple of float, optional
            Initial guess for the sky position as ``(ra, dec)``. If None, the
            position is initialized using ``ra_bonus`` and ``dec_bonus`` if
            provided, otherwise ``ra_targ`` and ``dec_targ``.

        Returns
        -------
        result : lmfit.minimizer.MinimizerResult
            Full fit result object.
        fit_metrics : dict
            Dictionary containing the fitted sky position and its
            1-sigma uncertainties. Contains ``"centerra"``,
            ``"sigmara"``, ``"centerdec"``, and ``"sigmadec"``.

        Raises
        ------
        KeyError
            If the requested TCE is not present in any observation.
        ValueError
            If the PRF fit fails or the requested TCE cannot be fit.
        RuntimeError
            If an unexpected error occurs during the PRF fit.
        """

        amp_list = self.transit_weights_list.copy()
        amp_err_list = self.transit_weights_err_list.copy()

        #initializing metrics dictionary
        fit_metrics = {}

        params = Parameters()
        #constraining params
        if initial_location is not None:
            ra, dec = initial_location
        elif self.ra_bonus is not None and self.dec_bonus is not None:
            ra, dec = self.ra_bonus, self.dec_bonus
        else:
            ra, dec = self.ra_targ, self.dec_targ

        params.add('centerra', value=ra)
        params.add('centerdec', value=dec)
        
        #doing the fit for only one tce, if a tce number is given 
        if which_tce is not None:
            #create a subset of all information we need for minimizer based on which quarters/sectors contain this tce
            valid_idxs = [i for i in range(self.n_observations) if which_tce in self.valid_tces_list[i]]
            if not valid_idxs:
                raise KeyError(f"TCE {which_tce} is not present in any observation.")
            
            #constraining the amplitude param for all quarters
            for i in valid_idxs:
                which_col = self.tce_mapping_list[i].get(which_tce)
                if which_col is None:
                    continue
                params.add(f'amplitude_{i}',value=np.nanmax(amp_list[i][which_col]),min=0)
            
            #inputs for all_prf_residual function
            prf_fit_kw = {'prf_list':self.prf_list,'data_list':amp_list,'data_err_list':amp_err_list,
                          'origin_list':self.origin_list,'wcs_list':self.wcs_list,'shape_list':self.shape_list,
                          'tce_mapping_list':self.tce_mapping_list,'valid_idxs':valid_idxs,'which_tce':which_tce}
            
            try:
                #minimizing to get the best fit results using prf_residual
                result = minimize(quarters_prf_residual,params,kws=prf_fit_kw)#,grid=False) 
                
            except ValueError as e:
                raise ValueError(f"PRF fit failed for TCE {which_tce}: {e}") from e

            except Exception as e:
                raise RuntimeError(f"Unexpected error during PRF fit for TCE {which_tce}: {e}") from e
            
        else:
            raise ValueError("Only supports one tce at a time for the moment. You need to specify a tce.")
        
        #extracing most important fit metrics 
        fit_metrics['centerra'] = result.params['centerra'].value
        fit_metrics['sigmara'] = result.params['centerra'].stderr
        fit_metrics['centerdec'] = result.params['centerdec'].value
        fit_metrics['sigmadec'] = result.params['centerdec'].stderr

        #returns the who report and the most important fit metrics
        return result, fit_metrics
    


    def get_quarters_offset(self,fit_metrics: dict[str, float]) -> tuple[float, float, float, float]:
        """
        Calculate the angular offset between the fitted and reference positions.

        Parameters
        ----------
        fit_metrics : dict
            Dictionary containing the fitted sky position and its
            1-sigma uncertainties. Required keys are ``'centerra'``,
            ``'centerdec'``, ``'sigmara'``, and ``'sigmadec'``.

        Returns
        -------
        nsigma : float
            Offset between the fitted and reference positions in units
            of the propagated 1-sigma uncertainty.
        dist : float
            Angular separation between the fitted and reference positions
            in degrees.
        sigma_dist : float
            1-sigma uncertainty on the angular separation in degrees.
        mah_dist : float
            Mahalanobis distance between the fitted and reference positions.

        Raises
        ------
        KeyError
            If a required key is missing from ``fit_metrics``.
        """
        
        #extract fit metrics
        ra0 = fit_metrics['centerra'] 
        dec0 = fit_metrics['centerdec'] 
        ra0_err = fit_metrics['sigmara'] 
        dec0_err = fit_metrics['sigmadec']

        #calculating differences
        if self.ra_bonus is not None and self.dec_bonus is not None:
            reference_ra = self.ra_bonus
            reference_dec = self.dec_bonus
        else:
            reference_ra = self.ra_targ
            reference_dec = self.dec_targ

        cos_dec = np.cos(np.radians(reference_dec))
        dra = (ra0 - reference_ra) * cos_dec
        ddec = dec0 - reference_dec  

        #calculating distance and its propagated uncertainty
        dist = np.sqrt(dra**2 + ddec**2)
        sigma_dist = (np.sqrt(dra**2 * ra0_err**2+ ddec**2 * dec0_err**2) / dist)

        #calculating mahalanobis distance
        cov = np.array([[(ra0_err*cos_dec)**2, 0],[0, dec0_err**2]])
        inv_cov = np.linalg.inv(cov)
        x = np.array([ra0*cos_dec, dec0])
        y = np.array([reference_ra*cos_dec, reference_dec])
        mah_dist = mahalanobis(x, y, inv_cov)

        #finding how many sigma away the center of the fit is from the actual location
        nsigma = dist/sigma_dist

        return nsigma, dist, sigma_dist, mah_dist



    def plot_multi_quarter(self,fit_metrics: dict[str, float],which_tce: int,savefig: bool = False,save_directory: str = '.') -> None:
        """
        Generate a localization plot combining fits from multiple observations.

        Parameters
        ----------
        fit_metrics : dict
            Dictionary containing the joint fitted sky position and its
            1-sigma uncertainties. Required keys are ``'centerra'``,
            ``'centerdec'``, ``'sigmara'``, and ``'sigmadec'``.
        which_tce : int
            Index of the TCE to visualize.
        savefig : bool, optional
            If True, save the figure to disk (default False).
        save_directory : str, optional
            Directory where the figure will be saved (default '.').

        Raises
        ------
        KeyError
            If ``which_tce`` is not present in an observation.
        ValueError
            If a per-observation PRF fit fails.
        RuntimeError
            If an unexpected error occurs during a per-observation fit.
        """

        ra0 = fit_metrics['centerra'] 
        dec0 = fit_metrics['centerdec'] 
        ra0err = fit_metrics['sigmara'] 
        dec0err = fit_metrics['sigmadec'] 

        #calculating offset, distance, and uncertainity of distance
        nsigma, dist, sigma_dist, mah_dist = self.get_quarters_offset(fit_metrics)

        #get catalog position
        if self.ra_bonus is not None and self.dec_bonus is not None:
            catalog_ra = self.ra_bonus 
            catalog_dec = self.dec_bonus 
        else:
            catalog_ra = self.ra_targ 
            catalog_dec = self.dec_targ 

        per_quart_ra_fit = []
        per_quart_dec_fit = []
        for i in range(self.n_observations):
            fit = self.loc_object_list[i].fit_to_heatmap(method='prf',which_tce=which_tce)
            x0 = fit[1]['centerx']
            y0 = fit[1]['centery']
            sky = self.wcs_list[i].pixel_to_world(x0,y0)
            per_quart_ra_fit.append(sky.ra.degree)
            per_quart_dec_fit.append(sky.dec.degree)

        #use catalog dec for projection
        cos_dec = np.cos(np.radians(catalog_dec))

        #calculate differences
        ra_fit_diff = (ra0 - catalog_ra) * cos_dec
        dec_fit_diff = (dec0 - catalog_dec)

        per_quarter_ra_diff = (np.array(per_quart_ra_fit) - np.array(catalog_ra)) * cos_dec
        per_quarter_dec_diff = (np.array(per_quart_dec_fit) - np.array(catalog_dec))

        ra_neighbors_diff = (self.gaia_coord.ra.deg - catalog_ra) * cos_dec
        dec_neighbors_diff = (self.gaia_coord.dec.deg - catalog_dec)
        
        #convert to arcseconds
        ra_fit_diff *= 3600
        dec_fit_diff *= 3600
        ra0err *= 3600
        dec0err *= 3600
        ra_neighbors_diff *= 3600
        dec_neighbors_diff *= 3600
        per_quarter_ra_diff *= 3600
        per_quarter_dec_diff *= 3600

        fig = plt.figure(figsize=(8,6))
        gs = GridSpec(1,2, width_ratios=[3,1], figure=fig, wspace=0.05)

        #main plot
        ax = fig.add_subplot(gs[0])

        #plotting per quarter fits as well
        ax.scatter(per_quarter_ra_diff, per_quarter_dec_diff,color='blue', s=30, alpha=0.7, label='Per-quarter fits')

        #fit location in terms of arcsec offset from catalog position 
        ax.errorbar(ra_fit_diff, dec_fit_diff, xerr= ra0err, yerr=dec0err,fmt='o',c='k',markersize=4,capsize=1.5)
        
        #ax.scatter(ra_neighbors_diff,dec_neighbors_diff, color='blue', alpha=0.8)

        #catalog at origin
        ax.scatter(0, 0, color='orange', marker='*', s=400)

        ax.set_xlabel(r'$\Delta$ RA [arcsec]', fontsize=14) 
        ax.set_ylabel(r'$\Delta$ Dec [arcsec]',fontsize=14)
        ax.set_aspect('equal', 'box')

        #side panel with offset info
        ax_text = fig.add_subplot(gs[1])
        ax_text.axis('off')

        report_text = (f'Euclidean Dist:\n {dist*3600:.4g} arcsec\n\n'
                       f'Eucl Dist Error:\n {sigma_dist*3600:.4g} arcsec\n\n'
                       f'Offset:\n {nsigma:.4g}'+r'$\sigma$'+f'\n\n'
                       f'Mahalanobis Dist:\n {mah_dist:.4g}'+r'$\sigma$')
        
        ax_text.text(0, 0.5, report_text,fontsize=14,verticalalignment='center')

        if self.mission=='kepler':
            cat_name='KIC'
        elif self.mission=='tess':
            cat_name='TIC'
        else:
            cat_name='CANDIDATE'
        
        plt.suptitle('Multi Quarter InCINERATOr Report for '+cat_name+' '+self.id+'.0'+str(which_tce+1), fontsize=18)

        if savefig:
            plt.savefig(save_directory+'/multi_quarter_incinerator_report_'+cat_name+self.id+'.0'+str(which_tce+1)+'.png', dpi=200, pad_inches=0.5, bbox_inches='tight')