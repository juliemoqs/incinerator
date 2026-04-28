import batman
import lkprf
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from lmfit import Parameters, minimize
from lmfit.models import Gaussian2dModel
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Ellipse
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

    # ADD DATA TYPES AND TCES AS DATAFRAME NOT NUMPY ARRAY 
    #initiliazing -- data parser and loader
    def __init__(self, time, flux, flux_err, tces, id, 
                 wcs=None, mission=None, file_name=None, 
                 ra_bonus=None, dec_bonus=None,
                 ra_targ=None , dec_targ=None, 
                 prf=None, epoch=None):
        """
        Initialize the Localize object.

        Parameters
        ----------
        time : numpy.ndarray
            1D array of time values.
        flux : numpy.ndarray
            3D flux cube with shape (n_time, n_rows, n_cols).
        flux_err : numpy.ndarray
            3D array of per-pixel uncertainties matching `flux`.
        tces : numpy.ndarray
            Array of TCE parameters. Typically derived from a pandas
            DataFrame and converted to a NumPy array (one row per TCE).
        id : str
            Target identifier. Do not include mission prefix (TIC/KIC/etc.)
        wcs : astropy.wcs.WCS
            WCS solution for pixel-to-sky transformations. If not provided, one is generated from file header for Kepler and TESS.
        mission : str, optional
            Mission name (e.g., 'kepler', 'tess').
        file_name : str, optional
            Path to a TPF or mission-specific file. Not optional for Kepler and TESS.
        ra_bonus, dec_bonus : float, optional
            Optional J2000 sky coordinate of a source different from the nominal mission target. 
        ra_targ, dec_targ : float
            J2000 target sky coordinates in degrees. If not provided, one is generated from file header for Kepler and TESS.
        prf : object, optional
            Precomputed PRF model. If not provided, one is generated for Kepler and TESS.
        epoch : astropy.time.Time, optional
            Epoch of the input coordinate. If not provided, one is generated for Kepler and TESS.

        Raises
        ------
        ValueError
            If flux is not 3D or time length does not match the flux cube.
            Also raised if `file_name` is provided without `mission`.
        """
        
        self.time = time
        self.flux = flux
        self.flux_err = flux_err
        self.tces = tces
        self.wcs = wcs
        self.mission = mission
        self.ra_bonus = ra_bonus
        self.dec_bonus = dec_bonus
        self.file_name = file_name
        self.ra_targ = ra_targ
        self.dec_targ = dec_targ
        self.prf = prf
        self.id = id
        self.shape = self.flux.shape[1:]

        #making sure data cube is 3D
        if self.flux.ndim != 3:
            raise ValueError("Flux must be a 3D array: (time,y,x)")
        if self.flux.shape[0] != self.time.size:
            raise ValueError("Time length must match flux first axis")
        
        #getting information for Kepler tpf or TESS tesscube
        if file_name is not None:
            if mission is None:
                raise ValueError("mission must be specified when providing file_name")

            self.load_tpf_info()

        #getting the epoch of the observation
        if mission == "kepler":
            self.epoch = Time(np.nanmedian(self.time) + 2454833, format='jd')
        
        elif mission == "tess":
            self.epoch = Time(np.nanmedian(self.time) + 2457000, format='jd')
        
        else:
            if epoch is None:
                raise ValueError("Please provide an astropy Time object for epoch.")
            else:
                self.epoch = epoch

        #need to get the prf model
        self.prf = self.get_prf_model()

        #search for nearby stars and propagate all coordinates
        self._search_and_propagate()

        #cleaning the data
        self._clean_data()

    
 
    def _clean_data(self):
        """
        Remove invalid cadences and pixels from the flux cube.

        Attributes
        --------------
        good_pix_mask : numpy.ndarray
            Boolean mask identifying pixels with nonzero total flux.
        good_cad_mask : numpy.ndarray
            Boolean mask identifying cadences with finite time, flux,
            and flux uncertainties.
        pix : numpy.ndarray
            Cleaned flux array containing only good cadences and pixels.
        pix_err : numpy.ndarray
            Cleaned uncertainty array matching `pix`.
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

        #save masks for later
        self.good_pix_mask = good_pix
        self.good_cad_mask = good_cad

        #extract the flux and flux_err arrays only for good pixels and finite times
        self.pix = flux[good_cad][:, good_pix]
        self.pix_err = flux_err[good_cad][:, good_pix]
        
        self.time = time[good_cad]
        
    
    # REMOVE THIS!
    def load_tpf_info(self):
        """
        Load mission-specific metadata from a target pixel file. The returned values are added directly to the object's
        attributes.

        Raises
        ------
        ValueError
            If the mission is not supported.
        """

        if self.mission == "kepler":
            info = self._get_kepler_tpf_info(self.file_name)
        elif self.mission == "tess":
            info = self._get_tesscut_info(self.file_name)
        else:
            raise ValueError(f"Unsupported mission: {self.mission}")

        self.__dict__.update(info)


    #MAKE THIS GENERAL AND ASK FOR MISSION -- RETURN THE LOCALIZE OBJECT
    @staticmethod
    def from_kepler_tpf_info(file_name):
        """
        Extract relevant metadata from a Kepler target pixel file (TPF).

        Parameters
        ----------
        file_name : str
            Path to the Kepler TPF FITS file.

        Returns
        -------
        dict
            Dictionary containing:
            - 'ra_targ' : target right ascension (deg)
            - 'dec_targ' : target declination (deg)
            - 'wcs' : WCS solution
            - 'channel' : CCD channel number
            - 'origin_row' : detector row origin
            - 'origin_col' : detector column origin
            - 'quarter' : observation quarter number
            - 'headers' : a list of headers, one per HDU in the FITS file. Includes the primary header and all extension headers.
        """
        
        #reads in the file 
        hdulist = fits.open(file_name)

        #gets all the info we need whoop whoop
        info = {'ra_targ': hdulist[0].header['ra_obj'],
                'dec_targ': hdulist[0].header['dec_obj'],
                'wcs': WCS(hdulist[2].header),
                'channel': hdulist[0].header['channel'],
                'origin_row': hdulist[1].header['2CRV4P'],
                'origin_col': hdulist[1].header['1CRV4P'],
                'quarter': hdulist[0].header['quarter'],
                'headers':[hdu.header for hdu in hdulist]}

        hdulist.close()

        return info
    


    @staticmethod
    def _get_tesscut_info(file_name):
        """
        Extract relevant metadata from a TESSCUT file.

        Parameters
        ----------
        file_name : str
            Path to the TESSCUT FITS file.

        Returns
        -------
        dict
            Dictionary containing:
            - 'ra_targ' : target right ascension (deg)
            - 'dec_targ' : target declination (deg)
            - 'wcs' : WCS solution
            - 'camera' : Camera number
            - 'ccd' : CCD Chip number
            - 'origin_row' : detector row origin
            - 'origin_col' : detector column origin
            - 'sector' : observation sector number
            - 'headers' : a list of headers, one per HDU in the FITS file. Includes the primary header and all extension headers.
        """

        #reads in the file 
        hdulist = fits.open(file_name)

        #gets all the info we need whoop whoop
        info = {'ra_targ': hdulist[0].header['ra_obj'],
                'dec_targ': hdulist[0].header['dec_obj'],
                'wcs': WCS(hdulist[2].header),
                'camera': hdulist[0].header['camera'],
                'ccd': hdulist[0].header['ccd'],
                'origin_row': hdulist[1].header['2CRV4P'],
                'origin_col': hdulist[1].header['1CRV4P'],
                'sector' : hdulist[0].header['sector'],
                'headers':[hdu.header for hdu in hdulist]}

        hdulist.close()

        return info
    
    

    def get_prf_model(self):
        """
        Initialize and return the mission-specific PRF model.

        Returns
        -------
        object or None
            PRF model instance if supported, otherwise None.
        """
        
        if self.mission == 'kepler':
            #run get_kepler_prf
            prf = self._get_kepler_prf()

        elif self.mission == 'tess':
            #run get_tess_prf
            prf = self._get_tess_prf()

        else:
            return

        return prf



    def _get_kepler_prf(self):
        """
        Construct the Kepler PRF model for the target channel.

        Returns
        -------
        lkprf.KeplerPRF
            Kepler PRF model for the specified channel.
        """
        
        #initializing prf for the specific channel
        prf = lkprf.KeplerPRF(channel = self.channel)

        return prf



    def _get_tess_prf(self):
        """
        Construct the TESS PRF model for the target camera and ccd chip.

        Returns
        -------
        lkprf.TESSPRF
            TESS PRF model for the specified camera and ccd chip.
        """

        #initializing prf for the specific camera and ccd chip 
        prf = lkprf.TESSPRF(camera = self.camera, ccd = self.ccd)

        return prf
    


    def _search_and_propagate(self):
        """
        Query Gaia around the target and propagate coordinates to the
        observation epoch.

        Attributes
        -------
        gaia_sources : astropy.table.Table
            Gaia sources returned from the cone search.
        gaia_coord : astropy.coordinates.SkyCoord
            Gaia source coordinates at the observation epoch.
        gaia_pix : ndarray
            Pixel coordinates of the Gaia sources.
        ra_targ, dec_targ : float
            Target coordinates propagated to the observation epoch.
        ra_bonus, dec_bonus : float
            Bonus coordinates propagated to the observation epoch if provided.
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



    def _get_transit(self,which_tce,method='box',tdur_frac=3.,**batman_kws):
        """
        Generate a transit model mask for a specified TCE.

        Parameters
        ----------
        which_tce : int
            Index of the TCE in `self.tces`.

        Returns
        -------
        numpy.ndarray
            1D array with value -1 during transit and 0 elsewhere,
            evaluated at the object's time array.
        """
        
        time = self.time
        #good_cad = self.good_cad_mask
        period,tdur,t0 = get_p_tdur_t0(self.tces[which_tce])

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

            #creating transit mask for design-matrix
            transit = flux - np.median(flux) #out of transit baseline needs to be at 0


        else:
            raise ValueError('Method should be box or batman')


        return transit



    def build_design_matrix(self, order=3, spacing_mult = 3., tdur_frac = 3., transit_method='box',batman_kws=None,cbv=False):
        """
        Construct the full design matrix.

        Parameters
        ----------
        order : int, optional
            Order of the spline used when `method='spline'` (default 3).
        spacing_mult : int, optional
            Multiplier applied to the transit duration (tdur) to set spline knot spacing when `method='spline'` (default 3.).
        transit_method : str, optional
            Method used for the transit component. Supports 'box' and 'batman'. Default is 'box'.
        batman_kws : list, optional
            List of dictionaries for each tce of additional arguments passed to the `batman` transit model when `'transit_method='batman'`.

        Raises
        ------
        ValueError
            If an unsupported `transit_method` is provided.

        Attributes
        ----
        self.valid_tces : numpy.ndarray
            Indices of self.tces that produced non-zero transit models and are included as columns in the design matrix.
        self.tce_mapping : dict
            Mapping of tce index in the tces dataframe to the column index in the design matrix. 
        self.design_matrix : numpy.ndarray
            Final design matrix combining transit, polynomial, and spline components.
        """
        
        time = self.time

        tdur_list = []

        transit_cols = []
        valid_tce_indices = []
        #iterating through all the tces
        for i in range(len(self.tces)):
            period,tdur,t0 = get_p_tdur_t0(self.tces[i])
            tdur_list.append(tdur)
            if transit_method == 'batman':
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
            
            if len(transit_cols) > 0:
                transits = np.hstack(transit_cols)

            else:
                #no transits in this quarter → skip quarter entirely
                transits = None
                self.valid_tces = []
                self.design_matrix = None
                return

        tdur_long = np.nanmax(tdur_list)

        #get the piecewise spline portion -- short-term variability
        spline = create_spline_design_matrix(time,tdur_long,order,spacing_mult)
        
        if cbv is True:
            #get cbv portion
            raise ValueError('CBVs will be added in the future')

        #get piecewise polynomial portion -- long-term variability
        polynomial = create_polynomial_design_matrix(time)

        #combine piecewise trends, spline, and transit signal into a final design matrix
        dM = np.hstack([transits, polynomial, spline])

        self.valid_tces = np.array(valid_tce_indices, dtype=int)
        self.tce_mapping = {tce_idx: col_idx for col_idx, tce_idx in enumerate(self.valid_tces)}
        self.design_matrix = dM
    

    
    def _solve_weights(self):
        """
        Solve for per-pixel linear model weights and uncertainties.

        Attributes
        ----
        self.weights : numpy.ndarray
            Array of fitted weights with shape (n_pix, n_params).
        self.weights_err : numpy.ndarray
            Array of 1-sigma uncertainties on the fitted weights
            with shape (n_pix, n_params).
        """
        
        pix = self.pix
        pix_err= self.pix_err
        dM = self.design_matrix

        if self.design_matrix is None:
            self.weights = None
            self.weights_err = None
            return

        n_pix = pix.shape[1]
        n_params = dM.shape[1]
    
        self.weights, self.weights_err = np.zeros((2, n_pix, n_params))

        #solving for the weights
        for idx, y, e in zip(range(n_pix), pix.T, pix_err.T):
            sigma_w_inv = dM.T.dot(dM / e[:, None]**2)
            B = dM.T.dot(y / e**2)
            self.weights[idx] = np.linalg.solve(sigma_w_inv, B)
            self.weights_err[idx] = np.sqrt(np.diag(np.linalg.inv(sigma_w_inv)))        



    def _solve_transit(self,which_tce):
        """
        Extract the spatial transit depth solution for a given TCE.

        Parameters
        ----------
        which_tce : int
            Index of the TCE corresponding to the transit component in the design matrix.

        Returns
        -------
        transit_weight : numpy.ndarray
            2D array of fitted transit depths per pixel.
        transit_weight_err : numpy.ndarray
            2D array of 1-sigma uncertainties on the transit depths.
        """

        self._solve_weights()
        good_pix = self.good_pix_mask

        #getting the weights for the transit depth
        transit_weight, transit_weight_err = np.zeros((2, *self.flux.shape[1:])) * np.nan
        transit_weight[good_pix] = self.weights[:, which_tce]
        transit_weight_err[good_pix] = self.weights_err[:, which_tce]

        return transit_weight, transit_weight_err
    


    def solve_transit_weights(self):
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
        for i in self.valid_tces:#range(len(self.tces)):
            col_i = self.tce_mapping[i]
            amp_i, amp_err_i = self._solve_transit(col_i)
            amps.append(amp_i)
            amps_err.append(amp_err_i)

        self.transit_weights = amps
        self.transit_weights_err = amps_err



    def fit_to_heatmap(self, model_func=None, which_tce = None, method=None,**fit_kw):
        """
        Fit a model to the per-pixel transit depth heatmap.

        The spatial transit depth map(s) are fit using one of the following:
        - 'prf'        : Fit a mission PRF model (single TCE or joint fit).
        - '2dgaussian' : Fit an unconstrained 2D Gaussian model.
        - 'custom'     : Fit a user-supplied residual function.

        Parameters
        ----------
        model_func : callable, optional
            Custom residual function (used only if `method='custom'`).
        which_tce : int, optional
            Index of the TCE to fit. If None and `method='prf'`, a joint PRF fit is performed across all TCEs.
        method : str
            Fitting method: 'prf', '2dgaussian', or 'custom'.
        **fit_kw
            Additional keyword arguments passed to the custom model.

        Returns
        -------
        result : lmfit.ModelResult or lmfit.MinimizerResult
            Full fit result object.
        fit_metrics : dict
            Dictionary containing key localization metrics:
            - 'centerx', 'centery' : best-fit centroid position
            - 'sigmax', 'sigmay'   : 1-sigma uncertainties on centroid
        Raises
        ------
        ValueError
            If an unsupported fitting method is provided.
        """
        
        amp = self.transit_weights.copy()
        amp_err = self.transit_weights_err.copy()

        #initializing metrics dictionary
        fit_metrics = {}

        #getting pix coord for target star
        pix_col, pix_row = coords_to_pixels(self.wcs,self.ra_targ,self.dec_targ)

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

                    #inputs for prf_residual function
                    prf_fit_kw = {'prf':self.prf,'data':amp[which_col],'data_err':amp_err[which_col],
                                  'origin':(self.origin_row,self.origin_col),'shape':self.shape}

                    #minimizing to get the best fit results using prf_residual
                    result = minimize(prf_residual,params,kws=prf_fit_kw) 

                else:
                    #constraining the amplitude param for all tces
                    for i in self.valid_tces:#range(len(self.tces)):
                        col_i = self.tce_mapping[i]
                        params.add(f'amplitude_{i}',value=np.nanmax(amp[col_i]),min=0)

                    #inputs for all_prf_residual function
                    prf_fit_kw = {'prf':self.prf,'data':amp,'data_err':amp_err,
                                  'origin':(self.origin_row,self.origin_col),'shape':self.shape,'tces':self.tces}

                    #minimizing to get the best fit results using prf_residual
                    result = minimize(all_prf_residual,params,kws=prf_fit_kw) 


                #extracing most important fit metrics 
                fit_metrics['centerx'] = result.params['centerx'].value
                fit_metrics['sigmax'] = result.params['centerx'].stderr
                fit_metrics['centery'] = result.params['centery'].value
                fit_metrics['sigmay'] = result.params['centery'].stderr


            elif method == 'custom':
                print('i dont know if this will actually work so good luck :)')

                params = Parameters()
                #constraining params
                params.add('centery', value=pix_row, min=0, max=self.shape[0]-1)
                params.add('centerx', value=pix_col, min=0, max=self.shape[1]-1)
                params.add('amplitude',value=np.nanmax(amp[which_tce]),min=0)

                #minimizing to get the best fit results
                result = minimize(model_func,params,kws=fit_kw)

                #extracting most important fit metrics
                fit_metrics['centerx'] = result.params['centerx'].value
                fit_metrics['sigmax'] = result.params['centerx'].stderr
                fit_metrics['centery'] = result.params['centery'].value
                fit_metrics['sigmay'] = result.params['centery'].stderr


            elif method == '2dgaussian':
                #need to make sure there are no NaNs
                mask = np.isfinite(amp[which_tce])

                #create rows, cols, and amplitude
                rows,cols = np.indices(amp[which_tce].shape)
                x = cols[mask].ravel()
                y = rows[mask].ravel()
                z = amp[which_tce][mask].ravel()

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
                raise ValueError(f"Fit failed for TCE {which_tce} with Method {method}. {e}") from e
            else:
                raise ValueError(f"Fit failed for all TCEs with Method {method}. {e}") from e
        except Exception as e:
            if which_tce is not None:
                raise RuntimeError(f"Unexpected Fit failure for all TCEs with Method {method}: {e}") from e


        #returns the whole report and the most important fit metrics
        return result, fit_metrics



    def get_offset(self,fit_metrics):
        """
        Compute the centroid offset in units of sigma.

        Parameters
        ----------
        fit_metrics : dict
            Dictionary containing fitted centroid values and uncertainties
            ('centerx', 'centery', 'sigmax', 'sigmay').

        Attributes
        ----
        self.offset : float
            Radial offset significance in units of sigma.
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
    


    def _get_all_metrics(self,pix_col,pix_row,model_func=None, method='prf',**fit_kw):
        """
        Run heatmap fitting and return localization metrics.

        Parameters
        ----------
        pix_col : float
            Reference pixel column position.
        pix_row : float
            Reference pixel row position.
        model_func : callable, optional
            Custom residual function (used if `method='custom'`).
        method : str, optional
            Fitting method ('prf', '2dgaussian', or 'custom').
        **fit_kw
            Additional keyword arguments passed to the fitting routine.

        Returns
        -------
        dict
            Dictionary containing centroid positions, uncertainties,
            and the radial offset significance ('offset').
        """
        
        #getting the fit metrics
        _, fit_metrics = self.fit_to_heatmap(model_func=model_func, method=method,**fit_kw)

        all_metrics = fit_metrics.copy()

        #calculating the amount of sigma away
        offset, _, _, _ = self.get_offset(fit_metrics,pix_col,pix_row)

        #adding offset to metrics dictionary
        all_metrics['offset'] = offset

        return all_metrics



    def plot_heatmap(self,metrics,which_tce,savefig=False,save_directory='.'):
        """
        Generate a localization report plot for a given TCE.

        Parameters
        ----------
        metrics : dict
            Dictionary containing fitted centroid values and uncertainties
            ('centerx', 'centery', 'sigmax', 'sigmay').
        which_tce : int
            Index of the TCE in `self.tce` to visualize.
        savefig : bool, optional
            If True, save the figure to disk (default False).
        save_directory : str, optional
            Directory where the figure will be saved.
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

        fig, axes = plt.subplots(1,3, figsize=(16,4))

        axes[1].plot(pix_col,pix_row, 'rx')

        axes[1].errorbar(x0,y0,xerr=x0_err,yerr=y0_err,fmt='o',c='k',markersize=4,capsize=1.5)

        cax=axes[1].imshow(self.transit_weights[which_col]/self.transit_weights_err[which_col], cmap='inferno',origin='lower')#,vmin=0)
        axes[1].scatter(self.gaia_pix[0],self.gaia_pix[1],c='w',edgecolors='k',s=60)
        axes[1].tick_params(axis='both', labelsize=13)
        axes[1].set_xlabel('Pixel Number', fontsize=14)
        axes[1].set_ylabel('Pixel Number', fontsize=14)
        axes[1].set_xlim(-0.5,self.shape[1]-0.5)
        axes[1].set_ylim(-0.5,self.shape[0]-0.5)
        
        cbar1 = fig.colorbar(cax, ax=axes[1])
        cbar1.ax.tick_params(labelsize=13)
        cbar1.set_label('Depth / Depth Error', fontsize=16)


        cax1 = axes[0].imshow(np.nanmean(self.flux, axis=0),origin='lower')#, norm=colors.LogNorm()) 
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
        
        cbar2 = fig.colorbar(cax1, ax=axes[0])
        cbar2.ax.tick_params(labelsize=13)
        cbar2.set_label('Flux [e-/sec]', fontsize=16)

        #side panel with offset info
        axes[2].axis('off')

        report_text = (f'Offset:\n {nsigma:.4g}'+r'$\sigma$'+f'\n\n'
                       f'Euclidean Dist:\n {dist:.4g} pixels\n\n'
                       f'Euclidean Dist Error:\n {sigma_dist:.4g} pixels\n\n'
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
    def __init__(self, loc_object_list):

        #filtering out quarters that have no transits at all
        valid_loc_list = [loc for loc in loc_object_list if len(loc.valid_tces) > 0]

        #extracting information from each Localize object given 
        self.loc_object_list = valid_loc_list
        self.n_quarters = len(valid_loc_list)
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

    

    def fit_to_quarters(self, model_func=None, which_tce=None, method=None,**fit_kw):

        amp_list = self.transit_weights_list.copy()
        amp_err_list = self.transit_weights_err_list.copy()

        #initializing metrics dictionary
        fit_metrics = {}

        if method == 'prf':

            params = Parameters()
            #constraining params
            if self.ra_bonus is not None and self.dec_bonus is not None:
                params.add('centerra', value=self.ra_bonus)
                params.add('centerdec', value=self.dec_bonus)

            else:
                params.add('centerra', value=self.ra_targ)
                params.add('centerdec', value=self.dec_targ)
            

            #doing the fit for only one tce, if a tce number is given 
            if which_tce is not None:
                #create a subset of all information we need for minimizer based on which quarters/sectors contain this tce
                valid_idxs = [i for i in range(self.n_quarters) if which_tce in self.valid_tces_list[i]]
                print(valid_idxs)

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
                    result = minimize(quarters_prf_residual,params,kws=prf_fit_kw) 
                    
                except ValueError as e:
                    raise ValueError(f"PRF fit failed for TCE {which_tce}. {e}") from e

                except Exception as e:
                    raise RuntimeError(f"Unexpected PRF failure for TCE {which_tce}. {e}") from e
                
            else:
                raise ValueError("Only supports one tce at a time for the moment.")

        else:
            raise ValueError("Method only supports 'prf' for the moment.")
        
        #extracing most important fit metrics 
        fit_metrics['centerra'] = result.params['centerra'].value
        fit_metrics['sigmara'] = result.params['centerra'].stderr
        fit_metrics['centerdec'] = result.params['centerdec'].value
        fit_metrics['sigmadec'] = result.params['centerdec'].stderr

        
        #returns the who report and the most important fit metrics
        return result, fit_metrics
    


    def get_quarters_offset(self,fit_metrics):

        #extract fit metrics
        ra0 = fit_metrics['centerra'] 
        dec0 = fit_metrics['centerdec'] 
        ra0_err = fit_metrics['sigmara'] 
        dec0_err = fit_metrics['sigmadec']

        cov = np.array([[ra0_err**2,0],[0,dec0_err**2]])
        inv_cov = np.linalg.inv(cov)

        #calculating differences
        if self.ra_bonus is not None and self.dec_bonus is not None:
            dra = (ra0 - self.ra_bonus)
            ddec = (dec0 - self.dec_bonus)

            dist = (dra**2 * np.cos(np.radians(self.dec_bonus))**2 + ddec**2)**(1/2)
            sigma_dist = 1/dist * (dra**2 * np.cos(np.radians(self.dec_bonus))**4 * ra0_err**2 + ddec**2 * dec0_err**2)**(1/2)

            x = np.array([ra0,dec0])
            y = np.array([self.ra_bonus,self.dec_bonus])
            mah_dist = mahalanobis(x,y,inv_cov)

        else:
            dra = (ra0 - self.ra_targ)*np.cos(np.radians(self.dec_targ))
            ddec = (dec0 - self.dec_targ)

            dist = (dra**2 * np.cos(np.radians(self.dec_targ))**2 + ddec**2)**(1/2)
            sigma_dist = 1/dist * (dra**2 * np.cos(np.radians(self.dec_targ))**4 * ra0_err**2 + ddec**2 * dec0_err**2)**(1/2)

            x = np.array([ra0,dec0])
            y = np.array([self.ra_targ,self.dec_targ])
            mah_dist = mahalanobis(x,y,inv_cov)


        #finding how many sigma away the center of the fit is from the actual location
        nsigma = dist/sigma_dist

        return nsigma, dist, sigma_dist, mah_dist



    def plot_multi_quarter(self,fit_metrics,which_tce,savefig=False,save_directory='.'):
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
        for i in range(self.n_quarters):
            fit = self.loc_object_list[i].fit_to_heatmap(method='prf',which_tce=which_tce)
            x0 = fit[1]['centerx']
            y0 = fit[1]['centery']
            sky = self.wcs_list[i].pixel_to_world(x0,y0)
            print(sky)
            per_quart_ra_fit.append(sky.ra.degree)
            per_quart_dec_fit.append(sky.dec.degree)

        #use catalog dec for projection
        dec0_rad = np.deg2rad(catalog_dec)

        #calculate differences
        ra_fit_diff = (ra0 - catalog_ra) * np.cos(dec0_rad)
        dec_fit_diff = (dec0 - catalog_dec)

        per_quarter_ra_diff = (per_quart_ra_fit - catalog_ra) * np.cos(dec0_rad)
        per_quarter_dec_diff = (per_quart_dec_fit - catalog_dec)

        ra_neighbors_diff = (self.gaia_coord.ra.deg - catalog_ra) * np.cos(dec0_rad)
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
        ax.scatter(0, 0, color='orange', marker='*', s=400)#, label='Catalog position')

        #3 sigma ellipse around fit 
        ax.add_patch(Ellipse((ra_fit_diff, dec_fit_diff),width=3*mah_dist,height=3*mah_dist,fill=False,color='c'))

        ax.set_xlabel(r'$\Delta$ RA [arcsec]', fontsize=14) 
        ax.set_ylabel(r'$\Delta$ Dec [arcsec]',fontsize=14)
        ax.set_aspect('equal', 'box')

        #side panel with offset info
        ax_text = fig.add_subplot(gs[1])
        ax_text.axis('off')

        report_text = (f'Offset:\n {nsigma:.4g}'+r'$\sigma$'+f'\n\n'
                       f'Euclidean Dist:\n {dist*3600:.4g} arcsec\n\n'
                       f'Eucl Dist Error:\n {sigma_dist*3600:.4g} arcsec\n\n'
                       f'Mahalanobis Dist:\n {mah_dist:.4g}'+r'$\sigma$')
        
        ax_text.text(0, 0.5, report_text,fontsize=14,verticalalignment='center')#,family='monospace')

        #r = np.hypot(ra_fit_diff, dec_fit_diff)
        #ax.set_xlim(-r-ra0err,r+ra0err)
        #ax.set_ylim(-r-dec0err-1,r+dec0err+1)
        

        if self.mission=='kepler':
            cat_name='KIC'
        elif self.mission=='tess':
            cat_name='TIC'
        else:
            cat_name='CANDIDATE'
        
        plt.suptitle('Multi Quarter InCINERATOr Report for '+cat_name+' '+self.id+'.0'+str(which_tce+1), fontsize=18)

        if savefig:
            plt.savefig(save_directory+'/multi_quarter_incinerator_report_'+cat_name+self.id+'.0'+str(which_tce+1)+'.png', dpi=200, pad_inches=0.5, bbox_inches='tight')

        