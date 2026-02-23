#do i need to import os? 

import numpy as np
import matplotlib.pyplot as plt
from lmfit.models import Gaussian2dModel
from lmfit import Parameters, minimize
from astropy.io import fits
from astropy.wcs import WCS
import lkprf
from matplotlib import colors

from .utils import *

#initializing Class
class Localize(object):

    #initiliazing -- data parser and loader
    def __init__(self, time, flux, flux_err, tces, id, 
                 wcs=None, mission=None, file_name=None, 
                 ra_bonus=None, dec_bonus=None,
                 ra_targ=None , dec_targ=None, prf=None):
        
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

        #need to get the prf model
        self.prf = self.get_prf_model()

        #cleaning the data
        self._clean_data()

    
 
    def _clean_data(self):
        
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
        
    

    def load_tpf_info(self):
        #mission = mission.lower()

        if self.mission == "kepler":
            info = self._get_kepler_tpf_info(self.file_name)
        elif self.mission == "tess":
            info = self._get_tess_tpf_info(self.file_name)
        else:
            raise ValueError(f"Unsupported mission: {self.mission}")

        self.__dict__.update(info)



    @staticmethod
    def _get_kepler_tpf_info(file_name):
        
        #reads in the file 
        hdulist = fits.open(file_name)

        #gets all the info we need whoop whoop
        info = {'ra_targ': hdulist[0].header['ra_obj'],
                'dec_targ': hdulist[0].header['dec_obj'],
                'wcs': WCS(hdulist[2].header),
                #'shape': hdulist[1].data['Flux'].shape[1:],
                'channel': hdulist[0].header['channel'],
                'origin_row': hdulist[1].header['2CRV4P'],
                'origin_col': hdulist[1].header['1CRV4P']}

        hdulist.close()

        return info
    


    #THIS NEEDS TO BE FINISHED
    @staticmethod
    def _get_tesscube_info():

        return print('will eventually be added in :)')
    
    

    def get_prf_model(self):
        if self.mission == 'kepler':
            #run get_kepler_prf
            prf = self._get_kepler_prf()
            #print('YAY KEPLER')

        #THIS IS NOT FINISHED
        elif self.mission == 'tess':
            #run get_tess_prf
            print('YAY TESS, this will be added in the future!')
            prf = None

        else:
            return

        return prf



    def _get_kepler_prf(self):

        #initializing prf for the specific channel
        prf = lkprf.KeplerPRF(channel = self.channel)

        return prf



    #THIS NEEDS TO BE FINISHED
    def _get_tess_prf():
    #this will be added in, in the future
        print('YAY TESS, this will be added later!')
        return 1.



    def _get_transit(self,which_tce):
        
        time = self.time
        #good_cad = self.good_cad_mask
        period,tdur,t0 = get_p_tdur_t0(self.tces[which_tce])
        
        #fold lc
        #ph = ((time[good_cad] - t0)) % period
        ph = ((time - t0)) % period

        #shift phase to make it go from -P/2 to P/2 instead of 0 to P
        ph[ph > (0.5 * period)] -= period

        #creating transit mask
        transit =  - (np.abs(ph) < (tdur/2.)).astype(float) 

        return transit


#this is save in case we don't want to fit all transits at once
    #def build_design_matrix(self, tce_num, method='spline', order=3):
    #
    #    period,tdur,t0 = get_p_tdur_t0(self.tces[tce_num])
    #
    #    time = self.time
    #    good_cad = self.good_cad_mask
#
    #    #get transit
    #    transit = self._get_transit(tce_num)
#
    #    #get the systematics portion
    #    if method == 'spline':
    #       #get spline portion
    #        systematics = create_spline_design_matrix(time,good_cad,tdur,order=3)
    #
    #    elif method == 'CBV':
    #        #get cbv portion
    #        raise ValueError('Please use spline for the time being.')
#
    #    else:
    #        raise ValueError('Method should be spline or CBV.')
#
#
    #    #get piecewise polynomial portion -- generally stellar variability
    #    polynomial = create_polynomial_design_matrix(time,good_cad)
#
    #    #combine piecewise trends, spline, and transit signal into a final design matrix
    #    dM = np.hstack([transit[:, None], polynomial, systematics])
#
    #    self.design_matrix = dM



    def build_design_matrix(self, method='spline', order=3):
        
        time = self.time
        #good_cad = self.good_cad_mask

        tdur_list = []
        #iterating through all the tces
        for i in range(len(self.tces)):
            period,tdur,t0 = get_p_tdur_t0(self.tces[i])
            tdur_list.append(tdur)

            #get transit
            transit = self._get_transit(i)
            
            if i == 0:
                transits = transit[:,None]
            
            else:
                transits = np.hstack([transits, transit[:,None]])

        tdur_long = np.nanmax(tdur_list)

        #get the systematics portion
        if method == 'spline':
            #get spline portion
            systematics = create_spline_design_matrix(time,tdur_long,order)
        
        elif method == 'CBV':
            #get cbv portion
            raise ValueError('Please use spline for the time being.')

        else:
            raise ValueError('Method should be spline or CBV.')


        #get piecewise polynomial portion -- generally stellar variability
        polynomial = create_polynomial_design_matrix(time)

        #combine piecewise trends, spline, and transit signal into a final design matrix
        dM = np.hstack([transits, polynomial, systematics])

        self.design_matrix = dM
    


    #def _solve_weights(self):
    #    #idk if i want this because what if i wanna try a different design matrix? 
    #    #if hasattr(self, "weights"):
    #        #return  # already solved
#
    #    pix = self.pix
    #    pix_err= self.pix_err
    #    dM = self.design_matrix
#
    #    n_pix = pix.shape[1]
    #
    #    self.weights, self.weights_err = np.zeros((2, pix.shape[1], dM.shape[1]))
#
    #    print("Design matrix shape:", dM.shape)
#
    #    #solving for the weights
    #    for idx, y, e in zip(range(n_pix), pix.T, pix_err.T):
    #        sigma_w_inv = dM.T.dot(dM / e[:, None]**2)
    #        print("Rank:", np.linalg.matrix_rank(sigma_w_inv))
    #        print("Condition number:", np.linalg.cond(sigma_w_inv))
    #        B = dM.T.dot(y / e**2)
    #        self.weights[idx] = np.linalg.solve(sigma_w_inv, B)
    #        self.weights_err[idx] = np.sqrt(np.diag(np.linalg.inv(sigma_w_inv)))        

    
    def _solve_weights(self):
        #idk if i want this because what if i wanna try a different design matrix? 
        #if hasattr(self, "weights"):
            #return  # already solved

        pix = self.pix
        pix_err= self.pix_err
        dM = self.design_matrix

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

        self._solve_weights()
        good_pix = self.good_pix_mask

        #getting the weights for the transit depth
        transit_weight, transit_weight_err = np.zeros((2, *self.flux.shape[1:])) * np.nan
        transit_weight[good_pix] = self.weights[:, which_tce]
        transit_weight_err[good_pix] = self.weights_err[:, which_tce]

        return transit_weight, transit_weight_err
    


    def solve_transit_weights(self):

        amps = []
        amps_err = []
        for i in range(len(self.tces)):
            amp_i, amp_err_i = self._solve_transit(i)
            amps.append(amp_i)
            amps_err.append(amp_err_i)

        self.transit_weights = amps
        self.transit_weights_err = amps_err



    def fit_to_heatmap(self, model_func=None, which_tce = None, method=None,**fit_kw):
        #MODEL_FUNC ONLY FOR CUSTOM 
        #IF METHOD IS PRF NEED TO SET ALL TCES TO TRUE OR FALSE -- ADAPT CODE TO READ THIS TO CALL PRF OR ALL_PRF
        
        amp = self.transit_weights.copy()
        amp_err = self.transit_weights_err.copy()

        #initializing metrics dictionary
        fit_metrics = {}

        #getting pix coord for target star
        pix_col, pix_row = coords_to_pixels(self.wcs,self.ra_targ,self.dec_targ)

        if method == 'prf':

            params = Parameters()
            #constraining params
            params.add('centery', value=pix_row, min=0, max=self.shape[0]-1)
            params.add('centerx', value=pix_col, min=0, max=self.shape[1]-1)

            #doing the fit for only one tce, if a tce number is given 
            if which_tce is not None:

                #constraining the amplitude param
                params.add(f'amplitude',value=np.nanmax(amp[which_tce]),min=0)
                
                #inputs for prf_residual function
                prf_fit_kw = {'prf':self.prf,'data':amp[which_tce],'data_err':amp_err[which_tce],
                              'origin':(self.origin_row,self.origin_col),'shape':self.shape}

                #minimizing to get the best fit results using prf_residual
                result = minimize(prf_residual,params,kws=prf_fit_kw) 
            
            else:
                #constraining the amplitude param for all tces
                for i in range(len(self.tces)):
                    params.add(f'amplitude_{i}',value=np.nanmax(amp[i]),min=0)
                
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
            params.add(f'amplitude',value=np.nanmax(amp[which_tce]),min=0)

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

        #returns the who report and the most important fit metrics
        return result, fit_metrics



    #some type of report for the metrics, most specifically the sigma offset
    def get_offset(self,fit_metrics,pix_col,pix_row):

        #extract fit metrics
        x0 = fit_metrics['centerx']
        y0 = fit_metrics['centery']
        x0_err = fit_metrics['sigmax']
        y0_err = fit_metrics['sigmay']

        #calculating differences
        dx = x0 - pix_col
        dy = y0 - pix_row

        #finding how many sigma away the center of the fit is from the actual location
        nsigma = np.sqrt((dx/x0_err)**2+(dy/y0_err)**2)

        self.offset = nsigma
    


    def _get_all_metrics(self,pix_col,pix_row,model_func=None, method='prf',**fit_kw):
        #this can eventually be modified to add chi2 and whatever other metrics we want : )
        
        #getting the fit metrics
        _, fit_metrics = self.fit_to_heatmap(model_func=model_func, method=method,**fit_kw)

        all_metrics = fit_metrics.copy()

        #calculating the amount of sigma away
        self.get_offset(fit_metrics,pix_col,pix_row)

        #adding offset to metrics dictionary
        all_metrics['offset'] = self.offset

        return all_metrics



    def plot_heatmap(self,metrics,which_tce,savefig=False,save_directory='.'):

        #extracting metrics
        x0 = metrics['centerx']
        x0_err = metrics['sigmax']
        y0 = metrics['centery']
        y0_err = metrics['sigmay']

        if self.ra_bonus is not None and self.dec_bonus is not None:
            #getting the pixel info for the given ra and dec
            pix_col, pix_row = coords_to_pixels(self.wcs,self.ra_bonus,self.dec_bonus)
        
        else:
            pix_col, pix_row = coords_to_pixels(self.wcs,self.ra_targ,self.dec_targ)

        #getting offset
        self.get_offset(metrics,pix_col,pix_row)

        #search for background stars using a cone search
        _, pix_bkgd = query_vizier_background(self.wcs,self.ra_targ,self.dec_targ,self.shape)

        fig, axes = plt.subplots(1,2, figsize=(11,4))

        axes[0].plot(pix_col,pix_row, 'rx')

        axes[0].errorbar(x0,y0,xerr=x0_err,yerr=y0_err,fmt='o',c='k',markersize=2,capsize=1.5)

        cax=axes[0].imshow(self.transit_weights[which_tce]/self.transit_weights_err[which_tce], cmap='inferno')#,vmin=0)
        axes[0].scatter(pix_bkgd[0],pix_bkgd[1],c='w',edgecolors='k',s=60)
        axes[0].tick_params(axis='both', labelsize=13)
        axes[0].set_xlabel('Pixel Number', fontsize=14)
        axes[0].set_ylabel('Pixel Number', fontsize=14)
        #axes[0].set_title(f'Offset: {self.offset:.3f}$\sigma$',fontsize=14)
        
        cbar1 = fig.colorbar(cax, ax=axes[0])
        cbar1.ax.tick_params(labelsize=13)
        cbar1.set_label('Depth / Depth Error', fontsize=16)


        cax1 = axes[1].imshow(np.nanmean(self.flux, axis=0))#, norm=colors.LogNorm()) 
        axes[1].plot(pix_col,pix_row, 'rx')

        axes[1].scatter(pix_bkgd[0],pix_bkgd[1],c='w',edgecolors='k',s=60)
        axes[1].tick_params(axis='both', labelsize=13)
        axes[1].set_xlabel('Pixel Number', fontsize=14)
        axes[1].set_ylabel('Pixel Number', fontsize=14)
        
        cbar2 = fig.colorbar(cax1, ax=axes[1])
        cbar2.ax.tick_params(labelsize=13)
        cbar2.set_label('Flux [e-/sec]', fontsize=16)

        fig.text(0.5, -0.1, f'Offset = {self.offset:.3f}$\sigma$', ha='center',fontsize=14)
        
        if self.mission=='kepler':
            cat_name='KIC'
        elif self.mission=='tess':
            cat_name='TIC'
        else:
            cat_name='CANDIDATE'
        
        plt.suptitle('INCINERATOR Report for '+cat_name+' '+self.id+'.0'+str(which_tce+1), fontsize=18)

        if savefig == True:
            plt.savefig(save_directory+'/incinerator_report_'+cat_name+self.id+'.0'+str(which_tce+1)+'.png', dpi=200, pad_inches=0.5, bbox_inches='tight')


    
    #THIS NEEDS TO BE ADDED
    #need some type of visualization functionnnnn
    def get_localization_report(self):

        return 1.