# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import astropy as ap
import astropy.units
from astropy.table import Table
import astropy.constants as apc
import time
import sys
import scipy
import matplotlib.ticker as ti
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import sqlalchemy
import sqlalchemy.ext.declarative
import sqlalchemy.orm
import sqlalchemy.orm.query
import matplotlib
from sqlalchemy import and_, or_
from astropy import units as u
from astropy.coordinates import Angle
#from evtk.hl import gridToVTK

seconds_per_day = 60*60*24
batch_size = 1e8

# ==============================================================================
def calculate_FWHM(X,Y):
    """Calculate FWHM from arrays"""
    # Taken from http://stackoverflow.com/questions/10582795/finding-the-full-width-half-maximum-of-a-peak
    # Create 'difference' array by subtracting half maximum
    d = Y - (np.amax(Y) / 2)
    # Find the points where the difference is positive
    indexes = np.where(d>0)[0]
    # The first and last positive points are the edges of the peak
    return abs(X[indexes[-1]] - X[indexes[0]]) 

def calculate_centroid(X,Y, bounds=None):
    """Returns the centroid position, with optional flux interval"""
    bins = X
    vals = Y
    centroid_total = np.sum(vals)
    centroid_position = np.sum(np.multiply(bins,vals))/centroid_total

    if bounds is not None:
        bound_width = bounds/2
        bound_min = -1
        bound_max = -1
        value_total = 0
        for index, value in enumerate(vals):
            value_total += value
            if value_total/centroid_total >= 0.5+bound_width:
                bound_max = bins[index]
                break
        value_total = centroid_total
        for index,value in enumerate(vals[::-1]):
            value_total -= value
            if value_total/centroid_total <= 0.5-bound_width:
                bound_min = bins[len(bins)-1-index]
                break
        return centroid_position, bound_min, bound_max
    else:
        return centroid_position

def calculate_modal_value(X,Y):
    """Find the modal delay"""
    return X[np.argmax(Y)]

def calculate_midpoints(X):
    X_midp = np.zeros(shape=len(X)-1)
    for i in range(0,len(X)-1):
        X_midp[i] = (X[i] + X[i+1])/ 2
    return X_midp

# ==============================================================================       
def calculate_delay(angle, phase, radius, timescale):
    """Delay relative to continuum observed at angle for emission at radius"""
    # Calculate delay compared to continuum photons
    # Draw plane at r_rad_min out. Find x projection of disk position.
    # Calculate distance travelled to that plane from the current disk position
    # Delay relative to continuum is thus (distance from centre to plane)
    # + distance from centre to point
    # vr_disk     = np.array([r_rad.value*np.cos(r_phase), 0.0]) * u.m
    vr_disk     = np.array([radius*np.cos(phase), 0.0])
    vr_normal   = np.array([np.cos(angle), np.sin(angle)])
    vr_plane    = radius * vr_normal
    # return (np.dot((vr_plane - vr_disk), vr_normal) * u.m / apc.c.value).to(timescale)
    return (np.dot((vr_plane - vr_disk), vr_normal) / apc.c.value) / seconds_per_day
    

def keplerian_velocity(mass, radius):
    """Calculates Keplerian velocity at radius"""
    return np.sqrt(ap.constants.G.value * mass / radius)

def path_to_delay(path):
    """Converts path to time delay"""
    return path / apc.c.value

def doppler_shift_wave(line, vel):
    """Converts passed line and velocity into red/blueshifted wavelength"""
    return line * apc.c.value / (apc.c.value - vel)

def doppler_shift_vel(line, wave):
    """Converts passed red/blueshifted wave into velocity"""
    if wave > line:
        return -1*apc.c.value * (1 - (line / wave))
    else:
        return apc.c.value * ((line / wave) - 1)

def clamp(minimum, x, maximum):
    return max(minimum, min(x, maximum))

# ==============================================================================       
class TransferFunction:    
    # FAKE INIT
    def __init__(self, database, filename, luminosity=None, dimensions=None, template=None, template_different_line=False, template_different_spectrum=False):
        """Initialises the TF, taking the query it is to execute, the dimensions to bin by, and the delay range bins"""
        assert dimensions is not None or template is not None,\
            "Must provide either dimensions or another TF to copy them from!"
        # self._query = database.query(Photon.Wavelength, Photon.Delay, Photon.Weight, Photon.X, Photon.Z)
        start = time.clock()

        self._database = database
        Session = sqlalchemy.orm.sessionmaker(bind=self._database)
        self._session = Session()
        self._query = self._session.query(Photon.Wavelength, Photon.Delay, Photon.Weight)

        self._delay_dynamic_range = None
        self._velocity = None
        self._line_list = None
        self._line_wave = None
        self._line_num = None
        self._delay_range = None
        self._luminosity=None
        self._filename = filename
        self._dimensions = dimensions
        self._bins_vel = None
        self._bins_wave = None
        self._bins_delay = None
        self._flux = None
        self._flux_w = None
        self._count = None
        self._response_map = None
        self._wave_range = None
        self._spectrum = None

        if template is not None:
            print("Templating '{}' off of '{}'...".format(self._filename, template._filename))
            self._dimensions = template._dimensions
            self._delay_dynamic_range = template._delay_dynamic_range
            self._bins_vel = template._bins_vel
            self._bins_delay = template._bins_delay
            self._luminosity= template._luminosity

            # If we are templating off of the same line, we want the same wavelength bins
            if template_different_line is False:
                self._bins_wave = template._bins_wave
            if template._line_wave is not None and template_different_line is False:
                self.line(template._line_num, template._line_wave)
            if template._velocity is not None:
                self.velocities(template._velocity)
            if template._wave_range is not None and template_different_line is False:
                self.wavelengths(self._wave_range[0], self._wave_range[1])
            if template._line_list is not None and template_different_line is False:
                self.lines(template._line_list)
            if template._spectrum is not None and template_different_spectrum is False:
                self.spectrum(template._spectrum)

        if luminosity is not None:
            self._luminosity=luminosity
        # print("'{}' successfully created ({:.1f}s)".format(self._filename, time.clock()-start))

    def close_query(self):
        self._session.close()
        del(self._query)
        del(self._session)
        Session = sqlalchemy.orm.sessionmaker(bind=self._database)
        self._session = Session()
        self._query = self._session.query(Photon.Wavelength, Photon.Delay, Photon.Weight)
        return(self)

    def spectrum(self, number):
        self._spectrum = number
        self._query = self._query.filter(Photon.Spectrum == number)
        return self
    def line(self, number, wavelength):
        self._line_wave = wavelength
        self._line_num = number
        self._query = self._query.filter(Photon.Resonance == number)
        return self
    def velocities(self, velocity):
        assert self._line_wave is not None,\
            "Cannot limit doppler shift around a line without specifying a line!"
        self._velocity = velocity
        self._query = self._query.filter(Photon.Wavelength >= doppler_shift_wave(self._line_wave, -velocity), 
                                         Photon.Wavelength <= doppler_shift_wave(self._line_wave,  velocity))
        return self
    def wavelengths(self, wave_min, wave_max):
        assert wave_min < wave_max,\
            "Minimum wavelength must be lower than maximum wavelength!"
        self._wave_range = [wave_min, wave_max]
        self._query = self._query.filter(Photon.Wavelength >= wave_min, Photon.Wavelength <= wave_max)
        return self
    def wavelengths(self, wave_range):
        assert len(wave_range) < 3,\
            "When providing an array, it must be of more than 2 entries! Use wavelength(min, max)."
        self._bins_wave = wave_range
        self._wavelengths(wave_range[0], wave_range[-1])
        return self

    def lines(self, line_list):
        assert len(lines) > 1,\
            "For a single line, use the 'line()' filter rather than 'lines()'!"
        self._line_list = lines
        self._query = self._query.filter(Photon.Resonance.in_(line_list))
        return self
    def delays(self, delay_min, delay_max, unit='d'):
        assert delay_min < delay_max,\
            "Minimum delay must be below maximum delay!"

        if unit in ['d','D','day','Day','days','Days']:
            self._delay_range = [delay_min * seconds_per_day, delay_max * seconds_per_day]
        else:
            self._delay_range = [delay_min, delay_max]
        self._query=self._query.filter(Photon.Delay > self._delay_range[0], Photon.Delay < self._delay_range[1])
        return self
    def cont_scatters(self, scat_min, scat_max=None):
        if scat_max is not None:
            assert scat_min < scat_max,\
                "Minimum continuum scatters must be below maximum scatters!"
        assert scat_min >= 0,\
            "Must select a positive number of continuum scatters" 

        if scat_max is not None:
            self._query = self._query.filter(Photon.ContinuumScatters >= scat_min, Photon.ContinuumScatters <= scat_max)
        else:
            self._query = self._query.filter(Photon.ContinuumScatters == scat_min)
        return self
    def res_scatters(self, scat_min, scat_max=None):
        if scat_max is not None:
            assert scat_min < scat_max,\
                "Minimum resonant scatters must be below maximum scatters!"
        assert scat_min >= 0,\
            "Must select a positive number of resonant scatters"    

        if scat_max is not None:
            self._query = self._query.filter(Photon.ResonantScatters >= scat_min, Photon.ResonantScatters <= scat_max)
        else:
            self._query = self._query.filter(Photon.ResonantScatters == scat_min)
        return self

    def filter(self, *args):
        self._query=self._query.filter(args)
        return self

    def response_map_by_tf(self, low_state, high_state, plot=False, cf_min=1, cf_max=1):
        """Creates a response map from two other transfer functions, to be applied during plotting"""
        # The other two TFs ***must*** have identical bins and both provide ionising luminosity information
        assert self._flux is not None,\
            "You must run the TF query with '.run()' before response mapping it!"
        assert low_state._flux is not None and high_state._flux is not None,\
            "You must run the low and high state TF queries with '.run()' before response mapping using them!"
        assert np.array_equal(self._bins_wave, low_state._bins_wave) and np.array_equal(self._bins_delay, low_state._bins_delay),\
            "Low state TF is binned differently to target TF! Cannot rescale using it."
        assert np.array_equal(self._bins_wave, high_state._bins_wave) and np.array_equal(self._bins_delay, high_state._bins_delay),\
            "High state TF is binned differently to target TF! Cannot rescale using it."
        assert self._luminosity != None,\
            "TF missing continuum luminosity information!"
        assert low_state._luminosity != None,\
            "Low state TF missing continuum luminosity information!"
        assert high_state._luminosity != None,\
            "High state TF missing continuum luminosity information!"
        assert low_state._luminosity <= self._luminosity,\
            "Low state ionising luminosity greater than target TF ionising luminosity!"
        assert high_state._luminosity >= self._luminosity,\
            "High state ionising luminosity lower than target TF ionising luminosity!"

        luminosity_difference = high_state._luminosity - low_state._luminosity

        # If that is true, the map is trivial to construct. We divide the difference in TFs by the luminosity difference
        response_map = ((high_state._flux*high_state._luminosity*cf_max) - (low_state._flux*low_state._luminosity*cf_min)) / luminosity_difference

        #print("Response mapping from TFs with total values {:.2e}, {:.2e}, {:.2e}".format(total_low, total_mid, total_high))

        self._response_map = response_map
        return self


    def FWHM(self, response=False, velocity=True):
        """Calculates the full width half maximum of the TF"""
           
        if velocity:
            midpoints = calculate_midpoints(self._bins_vel)
        else:
            midpoints = calculate_midpoints(self._bins_wave)

        if response:
            return calculate_FWHM(midpoints, np.sum(self._response_map, 0)) 
        else:
            return calculate_FWHM(midpoints, np.sum(self._flux, 0)) 

    def delay(self, response=False, threshold=0, bounds=None):
        """Calculates the delay for the current data"""
        assert threshold < 1 or threshold >= 0,\
            "Threshold is a multiplier to the peak flux! It must be between 0 and 1"

        data = None
        if response:
            data = np.sum(self._response_map,1)
        else:
            data = np.sum(self._flux ,1)
        value_threshold = np.amax(data) * threshold
        delay_midp = calculate_midpoints(self._bins_delay)
        
        delay_weighted = 0
        value_total = 0
        for value, delay in zip(data, delay_midp):
            if value >= value_threshold:
                delay_weighted += value * delay
                value_total    += value
                
        return delay_weighted/value_total

    def run(self, response_map=None, line=None, scaling_factor=1.0, delay_dynamic_range=None, limit=None):
        """Performs a query on the photon DB and bins it"""
        assert response_map is None or line is not None,\
            "Passing a response map but no line information for it!"
        assert response_map is None or self._flux_w is None,\
            "A response map has already been built!"
        assert response_map is None,\
            "Response mapping by location not yet implemented!"
        assert self._flux is None,\
            "TF has already been run!"
        start = time.clock()

        # If we're not in limit mode, fetch all
        data = None
        if limit is None:
            data = np.asarray(self._query.all())
        else:
            print("Limited to {} results...".format(limit))
            data = np.asarray(self._query.limit(limit).all())

        assert len(data) > 0,\
            "No records found!"

        if self._spectrum is not None:
            print("Fetched {} records from '{}' for spectrum {}...".format(len(data), self._filename, self._spectrum))
        else:
            print("Fetched {} records from '{}'...".format(len(data), self._filename))


        # Check if we've already got delay bins from another TF
        if self._bins_delay is None:
            # Data returned as Wavelength, Delay, Weight. Find min and max delays
            if delay_dynamic_range is not None:
                self._delay_dynamic_range = delay_dynamic_range
                range_delay = [0,np.percentile(data[:,1],(1 - (10**(-delay_dynamic_range)))*100)]
                print("Delays up to the {} percentile value, {}d".format((1 - (10**(-delay_dynamic_range)))*100, range_delay[1]/seconds_per_day))
            else:
                range_delay = [0,np.amax(data[:,1])]
            self._bins_delay = np.linspace(range_delay[0], range_delay[1],
                                           self._dimensions[1]+1, endpoint=True, dtype=np.float64)

        # Check if we've already got wavelength bins from another TF
        if self._bins_wave is None:
            # If we have no velocity bins, this is a factory-fresh TF
            if self._bins_vel is None:
                # Data returned as Wavelength, Delay, Weight. Find min and max delays and wavelengths
                range_wave = [np.amin(data[:,0]), np.amax(data[:,0])] 
            # If we do have velocity bins, this was templated off a different line and we need to copy the velocities (but bins are in km! not m!)
            else:
                range_wave = [doppler_shift_wave(self._line_wave, self._bins_vel[0]*1000), doppler_shift_wave(self._line_wave, self._bins_vel[-1]*1000)]
                print("Creating new wavelength bins from template, velocities from {:.2e}-{:.2e} to waves: {:.2f}-{:.2f}".format(self._bins_vel[0], self._bins_vel[-1], range_wave[0], range_wave[1]))

            # Now create the bins for each dimension
            self._bins_wave  = np.linspace(range_wave[0], range_wave[1],
                                           self._dimensions[0]+1, endpoint=True, dtype=np.float64)                     

        # Check if we've already got velocity bins from another TF
        if self._bins_vel is None:
            # If this is a line-based TF, we can set velocity bins up
            if self._line_wave is not None:
                self._bins_vel  = np.linspace(doppler_shift_vel(self._line_wave, range_wave[1]), 
                                    doppler_shift_vel(self._line_wave, range_wave[0]),
                                    self._dimensions[0]+1, endpoint=True, dtype=np.float64)
                # Convert speed from m/s to km/s
                self._bins_vel = np.true_divide(self._bins_vel, 1000.0)


        # Now we bin the photons, weighting them by their photon weights for the luminosity
        self._flux, junk, junk = np.histogram2d(data[:,1], data[:,0], weights=data[:,2], 
                                                    bins=[self._bins_delay, self._bins_wave]) 
        # Keep an unweighted photon count for statistical error purposes                  
        self._count, junk, junk = np.histogram2d(data[:,1], data[:,0],
                                                    bins=[self._bins_delay, self._bins_wave])


        # Scaling factor! Each spectral cycle outputs L photons. If we do 50 cycles, we want a factor of 1/50
        self._flux *= scaling_factor
        # Scale to continuum luminosity
        self._flux /= self._luminosity

        print("'{}' successfully run ({:.1f}s)".format(self._filename,time.clock()-start))
        # Make absolutely sure this data is wiped as it's *HUGE*
        del(data)
        return self
        
    def flux(self, wave, delay):
        """Returns the photon luminosity in this bin"""
        return(self._flux[np.searchsorted(self._bins_delay, delay),
                          np.searchsorted(self._bins_wave, wave)])
    def count(self, wave, delay):
        """Returns the photon count in this bin"""
        return(self._count[np.searchsorted(self._bins_delay, delay),
                           np.searchsorted(self._bins_wave, wave)])

    def response_function_1d(self, days=True):
        """Returns a 1d response function"""
        if days:
            return np.column_stack((calculate_midpoints(self._bins_delay/seconds_per_day), np.sum(self._response_map, 1)))
        else:
            return np.column_stack((calculate_midpoints(self._bins_delay), np.sum(self._response_map, 1)))

    def plot(self, log=False, normalised=False, rescaled=False, velocity=False, name=None, days=True, 
            response_map=False, keplerian=None, dynamic_range=None):
        """Takes the data gathered by calling 'run' and outputs a plot"""
        assert response_map is False or self._response_map is not None,\
            "No data available for response map!"
        assert log is False or response_map is False,\
            "Cannot plot a logarithmic response map!"
        assert normalised is False or rescaled is False,\
            "Cannot be both normalised and rescaled!"
        assert self._bins_wave is not None,\
            "You must run the TF query with '.run()' before plotting it!"

        #matplotlib.rcParams["text.usetex"] = "True" 
        matplotlib.rcParams.update({'font.size': 14})

        start = time.clock()
        if name is not None:
            print("Plotting to file '"+self._filename+"_"+name+".eps'...")
        else:
            print("Plotting to file '"+self._filename+".eps'...")

        if dynamic_range is not None:
            log_range = dynamic_range
        elif self._delay_dynamic_range is not None:
            log_range = self._delay_dynamic_range
        else:
            log_range = 3

        fig = None
        ax_spec = None
        ax_tf = None
        ax_resp = None
        ax_rms = None
        # Set up the multiplot figure and axis
        fig, ((ax_spec, ax_none), (ax_tf, ax_resp)) = plt.subplots(2,2,sharex='col', sharey='row',
            gridspec_kw={'width_ratios':[3,1], 'height_ratios':[1,3]})
        ax_none.axis('off')
        ax_resp.invert_xaxis()
        fig.subplots_adjust(hspace=0, wspace=0)


        if response_map:
            ratio = np.sum(self._response_map)/np.sum(self._flux)
            ratio_exp = np.floor(np.log10(ratio))
            ratio_text = '\n'
            
            if ratio_exp < -1 or ratio_exp > 1:
                ratio_text_exp = r"{}{:.0f}{}".format("{",ratio_exp,"}")
                ratio_text += r"${:.2f}\times 10^{}$".format(ratio/(10**ratio_exp), ratio_text_exp)
            else:
                ratio_text += r"${:.3g}$".format(ratio)

            ax_tf.text(0.05, 0.95, r"$\frac{\Delta L}{L}/\frac{\Delta C}{C}=$"+ratio_text,
                transform=ax_tf.transAxes, fontsize=18, verticalalignment='top', horizontalalignment='left')

        # Set the properties that depend on log and wave/velocity status
        cb_label = None
        cb_label_vars = r""
        cb_label_units = r""
        cb_label_scale= r""
        cb_map = "afmhot_r"

        # Copy the data for later modification.
        data_plot = None
        if response_map:
            data_plot = np.copy(self._response_map)
            print("Total response: {:.3e}".format(np.sum(data_plot)))
            psi_label = r"$\Psi_{R}$"
        else:
            data_plot = np.copy(self._flux)
            print("Total line: {:.3e}".format(np.sum(data_plot)))
            psi_label = r"$\Psi_{T}$"
        cb_label = psi_label
    

        # Set the xlabel and colour bar label - these differ if velocity or not
        x_bin_mult = 1
        bins_x = np.zeros(shape=self._dimensions[0])
        if velocity:
            # We're rescaling the axis to e.g. 10^3 km/s but the colorbar is still in km/s
            # So when we scale by bin width, we need a multiplier on the bin widths
            oom = np.log10(np.amax(self._bins_vel))
            oom = oom - oom%3
            bins_x = self._bins_vel/(10**oom)
            x_bin_mult = 10**oom
            ax_tf.set_xlabel(r'Velocity ($10^{:.0f}$ km s$^{}$)'.format(oom, '{-1}'))
            cb_label_vars = r"($v, \tau$)"
            cb_label_units = r"/ km s$^{-1}$"
        else:
            bins_x = self._bins_wave
            ax_tf.set_xlabel(r'Wavelength $\lambda$ ($\AA$)')
            cb_label_vars += r"($\lambda, \tau$)"
            cb_label_units = r"$\AA$"

        bins_x_midp = np.zeros(shape=self._dimensions[0])
        for i in range(0,self._dimensions[0]):
            bins_x_midp[i] = (bins_x[i] + bins_x[i+1])/ 2

       # Set the ylabel and y bins for whether it's in days or seconds
        if days:
            bins_y = np.true_divide(self._bins_delay, float(seconds_per_day))
            data_plot *= seconds_per_day
            ax_tf.set_ylabel(r'Delay $\tau$ (days)')
            cb_label_units += r' d'
        else:
            bins_y = self._bins_delay
            ax_tf.set_ylabel(r'Delay $\tau$ (seconds)')
            cb_label_units += r' s'

        bins_y_midp = np.zeros(shape=self._dimensions[1])
        for i in range(0,self._dimensions[1]):
            bins_y_midp[i] = (bins_y[i] + bins_y[i+1])/ 2

        # Rescale the values to be luminosity/km s^-1 d or /A d
        for i in range(0, self._dimensions[0]):
            width_x = bins_x[i+1] - bins_x[i]
            for j in range(0, self._dimensions[1]):
                width_y = bins_y[j+1] - bins_y[j]
                data_plot[i][j] /= (width_x * x_bin_mult * width_y)

        # Plot the spectrum and light curve, normalised
        data_plot_spec = np.sum(data_plot, 0)
        data_plot_resp = np.sum(data_plot, 1)
        exponent_spec = np.floor(np.log10(np.amax(data_plot_spec)))
        exponent_resp = np.floor(np.log10(np.amax(data_plot_resp)))
        exponent_resp_text = "{}{:.0f}{}".format("{",exponent_resp,"}")
        exponent_spec_text = "{}{:.0f}{}".format("{",exponent_spec,"}")

 

        resp = ax_resp.plot(data_plot_resp/(10**exponent_resp), bins_y_midp, c='m')
        ax_spec.set_ylabel(r'{}(v) $10^{}/$'.format(psi_label, exponent_spec_text)+"\n"+r'km s$^{-1}$')
        if days:
            ax_resp.set_xlabel(r'{}($\tau$) $10^{}$ / d'.format(psi_label, exponent_resp_text))
        else:
            ax_resp.set_xlabel(r'{}($\tau$) $10^{}$ / s'.format(psi_label, exponent_resp_text))

        if response_map:
            ax_spec.axhline(0, color='grey')
            ax_resp.axvline(0, color='grey')
 
            data_plot_rms = np.sqrt(np.sum(np.square(data_plot), 0) / self._dimensions[0])           
            exponent_rms = np.floor(np.log10(np.amax(data_plot_rms)))
            exponent_rms_text = "{}{:.0f}{}".format("{",exponent_rms,"}")
            maximum_spec = np.amax(data_plot_spec)/np.power(10, exponent_spec)
            maximum_rms = np.amax(data_plot_rms)/np.power(10, exponent_rms)
            data_plot_rms /= np.amax(data_plot_rms)
            data_plot_spec /= np.amax(data_plot_spec)

            rms = ax_spec.plot(bins_x_midp, data_plot_rms, c='c', label=r'RMS {}(v)/{:.2f}$x10^{}$'.format(psi_label, maximum_rms, exponent_rms_text))
            spec = ax_spec.plot(bins_x_midp, data_plot_spec, c='m', label=r'{}(v)/{:.2f}$x10^{}$'.format(psi_label, maximum_spec, exponent_spec_text))
            lg_orig = ax_spec.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
            ax_spec.set_ylabel(r'{}(v)'.format(psi_label)+"\n"+r'km s$^{-1}$')
        else:
            spec = ax_spec.plot(bins_x_midp, data_plot_spec/(10**exponent_spec), c='m')
            ax_spec.set_ylabel(r'{}(v) $10^{}/$'.format(psi_label, exponent_spec_text)+"\n"+r'km s$^{-1}$')

        # If this is a log plot, take log and correct label and limits
        if log:
            cb_max = np.log10(np.amax(data_plot))
            cb_min = cb_max-log_range
            cb_label = r"Log "+cb_label
            data_plot = np.ma.log10(data_plot)
         # Else just scale the data
        else:
            maxval = np.floor(np.log10(np.amax(data_plot)))
            data_plot /= np.power(10,maxval)
            cb_max = np.amax(data_plot)
            cb_min = np.amin(data_plot)
            dummy = "{}{:.0f}{}".format("{",maxval,"}")
            cb_label_scale = r" 10$^{}$".format(dummy)        

        # If this is a response map, it may have a negative component and need a different plot
        if response_map:
            cb_max = np.amax([cb_max, np.abs(cb_min)])
            cb_min = -cb_max
            cb_map = 'RdBu_r' #'seismic'

        # Normalise or rescale the data. If doing neither, put units on cb.
        if normalised:
            data_plot /= np.sum(data_plot)     
            cb_label_units = r""
            cb_label_scale = r""
        elif rescaled:
            data_plot /= np.amax(data_plot)
            cb_label_units = r""
            cb_label_scale = r""

        # Plot the main colourplot for the transfer function
        tf = ax_tf.pcolor(bins_x, bins_y, data_plot,
                          vmin=cb_min, vmax=cb_max, cmap=cb_map)
        ax_tf.set_ylim(bottom=bins_y[0], top=bins_y[-1])
        ax_tf.set_xlim(left=bins_x[0], right=bins_x[-1])
        ax_tf.set_aspect('auto')


        # Add lines for keplerian rotational outflows
        if keplerian is not None:
            resolution  = 1000
            r_angle     = np.radians(keplerian["angle"])
            r_mass_bh   = keplerian["mass"] * apc.M_sun.value
            r_rad_grav  = (6 * apc.G.value * r_mass_bh / np.power(apc.c.value, 2))
            ar_wave     = np.zeros(resolution) # * u.angstrom
            ar_delay    = np.zeros(resolution) # * u.s
            ar_phase    = np.linspace(0, np.pi*2, resolution)
            ar_rad      = np.linspace(keplerian["radius"][0]*r_rad_grav, keplerian["radius"][1]*r_rad_grav,resolution)
            ar_vel      = np.zeros(resolution) 
            r_rad_min   = r_rad_grav * keplerian["radius"][0]
            r_vel_max   = keplerian_velocity( r_mass_bh, r_rad_min)
    
           
    
            # ITERATE OVER INNER EDGE
            for r_phase, r_wave, r_delay, r_vel in np.nditer([ar_phase, ar_wave, ar_delay, ar_vel], op_flags=['readwrite']):
	            r_vel[...]  	= r_vel_max *  np.sin(r_phase) * np.sin(r_angle) / (1e3 * x_bin_mult)
	            r_wave[...]		= doppler_shift_wave( self._line_wave, r_vel )
	            r_delay[...]	= calculate_delay( r_angle, r_phase, r_rad_min, u.day)
            ax_tf.plot(ar_vel, ar_delay, '-', c='m')

            # ITERATE OVER BLUE BOUND
            for r_rad, r_wave, r_delay, r_vel in np.nditer([ar_rad, ar_wave, ar_delay, ar_vel], op_flags=['readwrite']):
                r_rad           = r_rad # * u.m
                r_vel[...]      = keplerian_velocity( r_mass_bh, r_rad ) * np.sin(r_angle) / (1e3 * x_bin_mult)
                r_wave[...]     = doppler_shift_wave( self._line_wave, r_vel )
                r_delay[...]    = calculate_delay( r_angle, np.pi/2, r_rad, u.day )
            ax_tf.plot(ar_vel, ar_delay, '-', c='m')

            # ITERATE OVER RED BOUND
            for r_rad, r_wave, r_delay, r_vel in np.nditer([ar_rad, ar_wave, ar_delay, ar_vel], op_flags=['readwrite']):
                r_rad           = r_rad # * u.m
                r_vel[...]      = -keplerian_velocity( r_mass_bh, r_rad ) * np.sin(r_angle) / (1e3 * x_bin_mult)
                r_wave[...]     = doppler_shift_wave( self._line_wave, r_vel )
                r_delay[...]    = calculate_delay( r_angle, np.pi/2, r_rad, u.day )
            ax_tf.plot(ar_vel, ar_delay, '-', c='m')

        cbar = plt.colorbar(tf, orientation="vertical")
        cbar.set_label(cb_label+cb_label_vars+cb_label_scale+cb_label_units)

        if name is None:
            plt.savefig("{}.eps".format(self._filename),bbox_inches='tight')
        else:
            plt.savefig("{}_{}.eps".format(self._filename, name),bbox_inches='tight')
        print("Successfully plotted ({:.1f}s)".format(time.clock()-start))
        plt.close(fig)
        return self

# ==============================================================================
def open_database(s_file, s_user, s_password):
    ### TRY OPENING THE DATABASE ###
    print ("Opening database '{}'...".format(s_file))

    db_engine = None
    try:
        db_engine = sqlalchemy.create_engine("sqlite:///{}.db".format(s_file))
    except sqlalchemy.exc.SQLAlchemyError as e:
        print(e)
        sys.exit(1)

    ### DOES IT ALREADY EXIST? ###
    # print ("Searching for table 'Photons'...")
    Session = sqlalchemy.orm.sessionmaker(bind=db_engine)
    session = Session()

    start = time.clock()

    try:
        session.query(Photon.Weight).first()
        # If so, we go with what we've found.
        print("Found existing filled photon database '{}'".format(s_file))
    except sqlalchemy.exc.SQLAlchemyError as e:
        # If not, we populate from the delay dump file. This bit is legacy!
        print("No existing filled photon database, reading from file '{}.delay_dump'".format(s_file))
        Base.metadata.create_all(db_engine)

        added = 0
        delay_dump = open("{}.delay_dump".format(s_file), 'r')
        for line in delay_dump:
            if line.startswith('#'): 
                continue
            
            try:
                values = [float(i) for i in line.split()]
            except:
                print("Malformed line: '{}'".format(line))
                continue

            if len(values) is not 13:
                print("Malformed line: '{}'".format(line))
                continue

            session.add(Photon(Wavelength=values[1], Weight=values[2], X=values[3], Y=values[4], Z=values[5],
                            ContinuumScatters=values[6]-values[7], ResonantScatters=values[7], Delay=values[8],
                            Spectrum=values[10], Origin=(values[11]%10), Resonance=values[12], Origin_matom = (values[11]>9)))
            added += 1
            if added > 25000:
                added = 0
                session.commit()

        session.commit()
        session.close()
        del(session)
        print("Successfully read in ({:.1f}s)".format(time.clock()-start))

    return db_engine
# ==============================================================================       

s_user = "root"
s_password = "password"

Base = sqlalchemy.ext.declarative.declarative_base()
class Spectrum(Base):
    __tablename__ = "Spectra"
    id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True)
    angle = sqlalchemy.Column(sqlalchemy.Float)

class Origin(Base):
    __tablename__ = "Origins"
    id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True)
    name = sqlalchemy.Column(sqlalchemy.String)

class Photon(Base):
    __tablename__ = "Photons"
    id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True, autoincrement=True)
    Wavelength = sqlalchemy.Column(sqlalchemy.Float)
    Weight = sqlalchemy.Column(sqlalchemy.Float)
    X = sqlalchemy.Column(sqlalchemy.Float)
    Y = sqlalchemy.Column(sqlalchemy.Float)
    Z = sqlalchemy.Column(sqlalchemy.Float)
    ContinuumScatters = sqlalchemy.Column(sqlalchemy.Integer)
    ResonantScatters = sqlalchemy.Column(sqlalchemy.Integer)
    Delay = sqlalchemy.Column(sqlalchemy.Float)
    Spectrum = sqlalchemy.Column(sqlalchemy.Integer)
    Origin = sqlalchemy.Column(sqlalchemy.Integer)
    Resonance = sqlalchemy.Column(sqlalchemy.Integer)
    Origin_matom = sqlalchemy.Column(sqlalchemy.Boolean)
    #__table_args__ = (sqlalchemy.Index('spec_res', "Spectrum", "Resonance"),)


dims = [50, 50]
kep_sey = {"angle":40, "mass":1e7, "radius":[50,2000]}
kep_qso = {"angle":40, "mass":1e9, "radius":[50,20000]}

def do_tf_plots(tf_list_inp, dynamic_range=None, keplerian=None, name=None, file=None):
    tf_delay = []
    for tf_inp in tf_list_inp:  
        tf_inp.plot(velocity=True, keplerian=keplerian, log=False, name=name)
        tf_inp.plot(velocity=True, keplerian=keplerian, log=True,  name=('log' if name is None else name+"_log"), dynamic_range=dynamic_range)
        tf_delay.append(tf_inp.delay(threshold=0.8))

    if file is not None:
        print("Saving TF plots to file: {}".format(file+"_tf_delay.txt"))
        np.savetxt(file+"_tf_delay.txt", np.array(tf_delay,dtype='float'), header="Delay")

def do_rf_plots(tf_min, tf_mid, tf_max, keplerian=None, name=None, file=None):
    rf_delay = []
    if name is not None:
        name += '_'
    else:
        name = ''
    
    total_min  = np.sum(tf_min._flux).item()
    total_mid  = np.sum(tf_mid._flux).item()
    total_max = np.sum(tf_max._flux).item()

    calibration_factor = total_mid / ((total_min + total_max) / 2)

    tf_mid.response_map_by_tf(tf_min, tf_max,cf_min=1, cf_max=1).plot(velocity=True, response_map=True, keplerian=keplerian, name=name+"resp_mid")
    rf_mid = tf_mid.delay(response=True, threshold=0.8)
    
    tf_mid.response_map_by_tf(tf_min, tf_mid,cf_min=calibration_factor, cf_max=1).plot(velocity=True, response_map=True, keplerian=keplerian, name=name+"resp_low")
    rf_min = tf_mid.delay(response=True, threshold=0.8)
    tf_mid.response_map_by_tf(tf_mid, tf_max,cf_min=1, cf_max=calibration_factor).plot(velocity=True, response_map=True, keplerian=keplerian, name=name+"resp_high")
    rf_max = tf_mid.delay(response=True, threshold=0.8)

    if file is not None:
        print("Saving RF plots to file: {}".format(file+"_rf_delay.txt"))
        np.savetxt(file+"_rf_delay.txt", np.array([rf_min, rf_mid, rf_max],dtype='float'), header="Delay") 

# # ==============================================================================
# # RUN FOR SEYFERT
# # ==============================================================================
# sey100_db = open_database("/home/swm1n12/python_runs/paper1_5548_resp/sey_100", "root", "password")
# sey090_db = open_database("/home/swm1n12/python_runs/paper1_5548_resp/sey_090", "root", "password")
# sey110_db = open_database("/home/swm1n12/python_runs/paper1_5548_resp/sey_110", "root", "password")

lim_sey = 9999999

scale_sey_100 = (1/20)/(20/40)
scale_sey_110 = (1/20)
scale_sey_090 = (1/20)

wave_ha = 6562.81
wave_ha = 6562
wave_hb = 4861
wave_c4 = 1548.18

# ==============================================================================
# RUN FOR SHORT TIMESCALES
# ==============================================================================
# RUN FOR C4
# ------------------------------------------------------------------------------

#sey100_tf_c4 = TransferFunction(sey100_db, "sey100_c4", luminosity=1.043e44, dimensions=dims)
#sey100_tf_c4.line(416, wave_c4).spectrum(0).run(scaling_factor=scale_sey_100, delay_dynamic_range=1.25, limit=lim_sey)
#sey090_tf_c4 = TransferFunction(sey090_db, "sey090_c4", luminosity=0.9391e44, template=sey100_tf_c4, template_different_spectrum=True).spectrum(0).run(scaling_factor=scale_sey_090, limit=lim_sey)
#sey110_tf_c4 = TransferFunction(sey110_db, "sey110_c4", luminosity=1.148e44,  template=sey100_tf_c4, template_different_spectrum=True).spectrum(0).run(scaling_factor=scale_sey_110, limit=lim_sey)

#do_tf_plots([sey090_tf_c4, sey100_tf_c4, sey110_tf_c4], keplerian=kep_sey, file='sey_c4', dynamic_range=2)
#do_rf_plots( sey090_tf_c4, sey100_tf_c4, sey110_tf_c4,  keplerian=kep_sey, file='sey_c4')


# # ------------------------------------------------------------------------------
# # RUN FOR Ha
# # ------------------------------------------------------------------------------
#sey100_tf_ha = TransferFunction(sey100_db, "sey100_ha", luminosity=1.043e44,  template=sey100_tf_c4, template_different_line=True)
#sey100_tf_ha.line(28, wave_ha).spectrum(0).run(scaling_factor=scale_sey_100, limit=lim_sey)
#sey090_tf_ha = TransferFunction(sey090_db, "sey090_ha", luminosity=0.9391e44, template=sey100_tf_ha, template_different_spectrum=True).spectrum(0).run(scaling_factor=scale_sey_090, limit=lim_sey)
#sey110_tf_ha = TransferFunction(sey110_db, "sey110_ha", luminosity=1.148e44,  template=sey100_tf_ha, template_different_spectrum=True).spectrum(0).run(scaling_factor=scale_sey_110, limit=lim_sey)

#do_tf_plots([sey090_tf_c4, sey100_tf_c4, sey110_tf_c4], keplerian=kep_sey, file='sey_ha', dynamic_range=2)
#do_rf_plots( sey090_tf_c4, sey100_tf_c4, sey110_tf_c4,  keplerian=kep_sey, file='sey_ha')

# ------------------------------------------------------------------------------
# RUN FOR Hb
# ------------------------------------------------------------------------------
#sey100_tf_hb = TransferFunction(sey100_db, "sey100_hb", luminosity=1.043e44,  template=sey100_tf_c4, template_different_line=True)
#sey100_tf_hb.line(31, wave_hb).spectrum(0).run(scaling_factor=scale_sey_100, limit=lim_sey)
#sey090_tf_hb = TransferFunction(sey090_db, "sey090_hb", luminosity=0.9391e44, template=sey100_tf_hb, template_different_spectrum=True).spectrum(0).run(scaling_factor=scale_sey_090, limit=lim_sey)
#sey110_tf_hb = TransferFunction(sey110_db, "sey110_hb", luminosity=1.148e44,  template=sey100_tf_hb, template_different_spectrum=True).spectrum(0).run(scaling_factor=scale_sey_110, limit=lim_sey)

#do_tf_plots([sey090_tf_hb, sey100_tf_hb, sey110_tf_hb], keplerian=kep_sey, file='sey_hb', dynamic_range=2)
#do_rf_plots( sey090_tf_hb, sey100_tf_hb, sey110_tf_hb,  keplerian=kep_sey, file='sey_hb')

# ==============================================================================
# RUN FOR LONG TIMESCALES
# ==============================================================================
# RUN FOR C4
# ------------------------------------------------------------------------------
# sey100_tf_c4 = TransferFunction(sey100_db, "sey100_c4", luminosity=1.043e44,  dimensions=dims)
# sey100_tf_c4.line(416, wave_c4).spectrum(0).run(scaling_factor=scale_sey_100, delay_dynamic_range=1.5, limit=lim_sey)
# sey090_tf_c4 = TransferFunction(sey090_db, "sey090_c4", luminosity=0.9391e44, template=sey100_tf_c4, template_different_spectrum=True).spectrum(0).run(scaling_factor=scale_sey_090, limit=lim_sey)
# sey110_tf_c4 = TransferFunction(sey110_db, "sey110_c4", luminosity=1.148e44,  template=sey100_tf_c4, template_different_spectrum=True).spectrum(0).run(scaling_factor=scale_sey_110, limit=lim_sey)

# do_tf_plots([sey090_tf_c4, sey100_tf_c4, sey110_tf_c4], keplerian=kep_sey, name="long", file='sey_c4_long', dynamic_range=2)
# do_rf_plots( sey090_tf_c4, sey100_tf_c4, sey110_tf_c4,  keplerian=kep_sey, name="long", file='sey_c4_long')

# ------------------------------------------------------------------------------
# RUN FOR Ha
# ------------------------------------------------------------------------------
# sey100_tf_ha = TransferFunction(sey100_db, "sey100_ha", luminosity=1.043e44,  template=sey100_tf_c4, template_different_line=True)
# sey100_tf_ha.line(28, wave_ha).spectrum(0).run(scaling_factor=scale_sey_100, limit=lim_sey)
# sey090_tf_ha = TransferFunction(sey090_db, "sey090_ha", luminosity=0.9391e44, template=sey100_tf_ha, template_different_spectrum=True).spectrum(0).run(scaling_factor=scale_sey_090, limit=lim_sey)
# sey110_tf_ha = TransferFunction(sey110_db, "sey110_ha", luminosity=1.148e44,  template=sey100_tf_ha, template_different_spectrum=True).spectrum(0).run(scaling_factor=scale_sey_110, limit=lim_sey)

# do_tf_plots([sey090_tf_ha, sey100_tf_ha, sey110_tf_ha], keplerian=kep_sey, name="long", file='sey_ha_long', dynamic_range=2)
# do_rf_plots( sey090_tf_ha, sey100_tf_ha, sey110_tf_ha,  keplerian=kep_sey, name="long", file='sey_ha_long')

# ------------------------------------------------------------------------------
# RUN FOR Hb
# ------------------------------------------------------------------------------
# sey100_tf_hb = TransferFunction(sey100_db, "sey100_hb", luminosity=1.043e44,  template=sey100_tf_c4, template_different_line=True)
# sey100_tf_hb.line(31, wave_hb).spectrum(0).run(scaling_factor=scale_sey_100, limit=lim_sey)
# sey090_tf_hb = TransferFunction(sey090_db, "sey090_hb", luminosity=0.9391e44, template=sey100_tf_hb, template_different_spectrum=True).spectrum(0).run(scaling_factor=scale_sey_090, limit=lim_sey)
# sey110_tf_hb = TransferFunction(sey110_db, "sey110_hb", luminosity=1.148e44,  template=sey100_tf_hb, template_different_spectrum=True).spectrum(0).run(scaling_factor=scale_sey_110, limit=lim_sey)

# do_tf_plots([sey090_tf_hb, sey100_tf_hb, sey110_tf_hb], keplerian=kep_sey, name="long", file='sey_hb_long', dynamic_range=2)
# do_rf_plots( sey090_tf_hb, sey100_tf_hb, sey110_tf_hb,  keplerian=kep_sey, name="long", file='sey_hb_long')

# ==============================================================================
# RUN FOR 3 DAYS
# ==============================================================================
# RUN FOR C4
# ------------------------------------------------------------------------------
# sey100_tf_c4 = TransferFunction(sey100_db, "sey100_c4", luminosity=1.043e44, dimensions=dims).delays(0,3.2)
# sey100_tf_c4.line(416, wave_c4).spectrum(0).run(scaling_factor=scale_sey_100, delay_dynamic_range=5, limit=lim_sey)
# sey090_tf_c4 = TransferFunction(sey090_db, "sey090_c4", luminosity=0.9391e44, template=sey100_tf_c4, template_different_spectrum=True).spectrum(0).run(scaling_factor=scale_sey_090, limit=lim_sey)
# sey110_tf_c4 = TransferFunction(sey110_db, "sey110_c4", luminosity=1.148e44,  template=sey100_tf_c4, template_different_spectrum=True).spectrum(0).run(scaling_factor=scale_sey_110, limit=lim_sey)

# do_tf_plots([sey090_tf_c4, sey100_tf_c4, sey110_tf_c4], keplerian=kep_sey, name='300', file='sey_c4_300', dynamic_range=2)
# do_rf_plots( sey090_tf_c4, sey100_tf_c4, sey110_tf_c4,  keplerian=kep_sey, name='300', file='sey_c4_300')

# # ------------------------------------------------------------------------------
# # RUN FOR Ha
# # ------------------------------------------------------------------------------
# sey100_tf_ha = TransferFunction(sey100_db, "sey100_ha", luminosity=1.043e44,  template=sey100_tf_c4, template_different_line=True)
# sey100_tf_ha.line(28, wave_ha).spectrum(0).run(scaling_factor=scale_sey_100, limit=lim_sey)
# sey090_tf_ha = TransferFunction(sey090_db, "sey090_ha", luminosity=0.9391e44, template=sey100_tf_ha, template_different_spectrum=True).spectrum(0).run(scaling_factor=scale_sey_090, limit=lim_sey)
# sey110_tf_ha = TransferFunction(sey110_db, "sey110_ha", luminosity=1.148e44,  template=sey100_tf_ha, template_different_spectrum=True).spectrum(0).run(scaling_factor=scale_sey_110, limit=lim_sey)

# do_tf_plots([sey090_tf_ha, sey100_tf_ha, sey110_tf_ha], keplerian=kep_sey, name='300', file='sey_ha_300', dynamic_range=2)
# do_rf_plots( sey090_tf_ha, sey100_tf_ha, sey110_tf_ha,  keplerian=kep_sey, name='300', file='sey_ha_300')

# ==============================================================================
# RUN FOR QSO
# ==============================================================================
qso100_db = open_database("/media/data/Endjinn-Extra/swm1n12/broken_angles/agn_100", "root", "password")
qso090_db = open_database("/media/data/Endjinn-Extra/swm1n12/broken_angles/agn_090", "root", "password")
qso110_db = open_database("/media/data/Endjinn-Extra/swm1n12/broken_angles/agn_110", "root", "password")

lim_qso = 9999999999

scale_qso_100 = 1/20
scale_qso_110 = 1/20
scale_qso_090 = 1/20

# ==============================================================================
# RUN FOR SHORT TIMESCALES
# ==============================================================================
# RUN FOR C4
# ------------------------------------------------------------------------------
# qso100_tf_c4 = TransferFunction(qso100_db, "qso100_c4", luminosity=1.043e46, dimensions=dims)
# qso100_tf_c4.line(445, wave_c4).spectrum(2).run(scaling_factor=scale_qso_100, delay_dynamic_range=1, limit=lim_qso)
# qso110_tf_c4 = TransferFunction(qso110_db, "qso110_c4", luminosity=1.148e46,  template=qso100_tf_c4).run(scaling_factor=scale_qso_110, limit=lim_qso)
# qso090_tf_c4 = TransferFunction(qso090_db, "qso090_c4", luminosity=0.9391e46, template=qso100_tf_c4).run(scaling_factor=scale_qso_090, limit=lim_qso)

# do_tf_plots([qso090_tf_c4, qso100_tf_c4, qso110_tf_c4], keplerian=kep_qso, file='qso_c4', dynamic_range=2)
# do_rf_plots(qso090_tf_c4,  qso100_tf_c4, qso110_tf_c4,  keplerian=kep_qso, file='qso_c4')

# # ------------------------------------------------------------------------------
# # RUN FOR Ha
# # ------------------------------------------------------------------------------
# qso100_tf_ha = TransferFunction(qso100_db, "qso100_ha", luminosity=1.043e46,  template=qso100_tf_c4, template_different_line=True)
# qso100_tf_ha.line(44, wave_ha).run(scaling_factor=scale_qso_100, limit=lim_qso)
# qso090_tf_ha = TransferFunction(qso090_db, "qso090_ha", luminosity=0.9391e46, template=qso100_tf_ha).run(scaling_factor=scale_qso_090, limit=lim_qso)
# qso110_tf_ha = TransferFunction(qso110_db, "qso110_ha", luminosity=1.148e46,  template=qso100_tf_ha).run(scaling_factor=scale_qso_110, limit=lim_qso)

# do_tf_plots([qso090_tf_ha, qso100_tf_ha, qso110_tf_ha], keplerian=kep_qso, file='qso_ha', dynamic_range=2)
# do_rf_plots( qso090_tf_ha, qso100_tf_ha, qso110_tf_ha,  keplerian=kep_qso, file='qso_ha')

# # ------------------------------------------------------------------------------
# # RUN FOR Hb
# # ------------------------------------------------------------------------------
# qso100_tf_hb = TransferFunction(qso100_db, "qso100_hb", luminosity=1.043e46,  template=qso100_tf_c4, template_different_line=True)
# qso100_tf_hb.line(49, wave_hb).run(scaling_factor=scale_qso_100, limit=lim_qso)
# qso090_tf_hb = TransferFunction(qso090_db, "qso090_hb", luminosity=0.9391e46, template=qso100_tf_hb).run(scaling_factor=scale_qso_090, limit=lim_qso)
# qso110_tf_hb = TransferFunction(qso110_db, "qso110_hb", luminosity=1.148e46,  template=qso100_tf_hb).run(scaling_factor=scale_qso_110, limit=lim_qso)

# do_tf_plots([qso090_tf_hb, qso100_tf_hb, qso110_tf_hb], keplerian=kep_qso, file='qso_hb', dynamic_range=2)
# do_rf_plots( qso090_tf_hb, qso100_tf_hb, qso110_tf_hb,  keplerian=kep_qso, file='qso_hb')

# ==============================================================================
# RUN FOR LONG TIMESCALES
# ==============================================================================
# RUN FOR C4
# ------------------------------------------------------------------------------
#qso100_tf_c4 = TransferFunction(qso100_db, "qso100_c4", luminosity=1.043e46,  dimensions=dims)
#qso100_tf_c4.line(445, wave_c4).spectrum(2).run(scaling_factor=scale_qso_100, delay_dynamic_range=2, limit=lim_qso)
#qso090_tf_c4 = TransferFunction(qso090_db, "qso090_c4", luminosity=0.9391e46, template=qso100_tf_c4).run(scaling_factor=scale_qso_090, limit=lim_qso)
#qso110_tf_c4 = TransferFunction(qso110_db, "qso110_c4", luminosity=1.148e46,  template=qso100_tf_c4).run(scaling_factor=scale_qso_110, limit=lim_qso)

#do_tf_plots([qso090_tf_c4, qso100_tf_c4, qso110_tf_c4], keplerian=kep_qso, name="long", file='qso_c4_long', dynamic_range=3)
#do_rf_plots( qso090_tf_c4, qso100_tf_c4, qso110_tf_c4,  keplerian=kep_qso, name="long", file='qso_c4_long')

# ------------------------------------------------------------------------------
# RUN FOR Ha
# ------------------------------------------------------------------------------
#qso100_tf_ha = TransferFunction(qso100_db, "qso100_ha", luminosity=1.043e46,  template=qso100_tf_c4, template_different_line=True)
#qso100_tf_ha.line(44, wave_ha).run(scaling_factor=scale_qso_100, limit=lim_qso)
#qso090_tf_ha = TransferFunction(qso090_db, "qso090_ha", luminosity=0.9391e46, template=qso100_tf_ha).run(scaling_factor=scale_qso_090, limit=lim_qso)
#qso110_tf_ha = TransferFunction(qso110_db, "qso110_ha", luminosity=1.148e46,  template=qso100_tf_ha).run(scaling_factor=scale_qso_110, limit=lim_qso)

#do_tf_plots([qso090_tf_ha, qso100_tf_ha, qso110_tf_ha], keplerian=kep_qso, name="long", file='qso_ha_long', dynamic_range=3)
#do_rf_plots( qso090_tf_ha, qso100_tf_ha, qso110_tf_ha,  keplerian=kep_qso, name="long", file='qso_ha_long')

# ------------------------------------------------------------------------------
# RUN FOR Hb
# ------------------------------------------------------------------------------
#qso100_tf_hb = TransferFunction(qso100_db, "qso100_hb", luminosity=1.043e46,  template=qso100_tf_c4, template_different_line=True)
#qso100_tf_hb.line(49, wave_hb).run(scaling_factor=scale_qso_100, limit=lim_qso)
#qso090_tf_hb = TransferFunction(qso090_db, "qso090_hb", luminosity=0.9391e46, template=qso100_tf_hb).run(scaling_factor=scale_qso_090, limit=lim_qso)
#qso110_tf_hb = TransferFunction(qso110_db, "qso110_hb", luminosity=1.148e46,  template=qso100_tf_hb).run(scaling_factor=scale_qso_110, limit=lim_qso)#, limit=lim_qso)

#do_tf_plots([qso090_tf_hb, qso100_tf_hb, qso110_tf_hb], keplerian=kep_qso, name="long", file='qso_hb_long', dynamic_range=3)
#do_rf_plots( qso090_tf_hb, qso100_tf_hb, qso110_tf_hb,  keplerian=kep_qso, name="long", file='qso_hb_long')

# ==============================================================================
# RUN FOR 300 DAYS
# ==============================================================================
# RUN FOR C4
# ------------------------------------------------------------------------------
# qso100_tf_c4 = TransferFunction(qso100_db, "qso100_c4", luminosity=1.043e46,  dimensions=dims).delays(0,320)
# qso100_tf_c4.line(445, wave_c4).spectrum(2).run(scaling_factor=scale_qso_100, delay_dynamic_range=5, limit=lim_qso)
# qso090_tf_c4 = TransferFunction(qso090_db, "qso090_c4", luminosity=0.9391e46, template=qso100_tf_c4).run(scaling_factor=scale_qso_090, limit=lim_qso)
# qso110_tf_c4 = TransferFunction(qso110_db, "qso110_c4", luminosity=1.148e46,  template=qso100_tf_c4).run(scaling_factor=scale_qso_110, limit=lim_qso)

# do_tf_plots([qso090_tf_c4, qso100_tf_c4, qso110_tf_c4], keplerian=kep_qso, name="300", file='qso_c4_300', dynamic_range=3)
# do_rf_plots( qso090_tf_c4, qso100_tf_c4, qso110_tf_c4,  keplerian=kep_qso, name="300", file='qso_c4_300')

# tf_mid = qso100_tf_c4.delay(threshold=0.8)
# tfc_mid = qso100_tf_c4.delay(threshold=0.0)
# tf_min = qso090_tf_c4.delay(threshold=0.8)
# tfc_min = qso090_tf_c4.delay(threshold=0.0)
# tf_max = qso110_tf_c4.delay(threshold=0.8)
# tfc_max = qso110_tf_c4.delay(threshold=0.0)

# np.savetxt("centroid_delay.txt", np.array([[tf_min, tf_mid, tf_max], [tfc_min, tfc_mid, tfc_max]],dtype='float'), header="Delay_norm Delay_full") 

# ------------------------------------------------------------------------------
# RUN FOR Ha
# ------------------------------------------------------------------------------
# qso100_tf_ha = TransferFunction(qso100_db, "qso100_ha", luminosity=1.043e46,  template=qso100_tf_c4, template_different_line=True)
# qso100_tf_ha.line(44, wave_ha).run(scaling_factor=scale_qso_100, limit=lim_qso)
# qso090_tf_ha = TransferFunction(qso090_db, "qso090_ha", luminosity=0.9391e46, template=qso100_tf_ha).run(scaling_factor=scale_qso_090, limit=lim_qso)
# qso110_tf_ha = TransferFunction(qso110_db, "qso110_ha", luminosity=1.148e46,  template=qso100_tf_ha).run(scaling_factor=scale_qso_110, limit=lim_qso)

# do_tf_plots([qso090_tf_ha, qso100_tf_ha, qso110_tf_ha], keplerian=kep_qso, name="300", file='qso_ha_300', dynamic_range=3)
# do_rf_plots( qso090_tf_ha, qso100_tf_ha, qso110_tf_ha,  keplerian=kep_qso, name="300", file='qso_ha_300')

# outputs = qso100_tf_ha.response_function_1d()
# np.savetxt("qso100_1drf.txt",outputs, header='Days Response')
# print(qso100_tf_ha._response_map[qso100_tf_ha._response_map < 0])

# ==============================================================================
# RUN FOR f-vs-Angle DAYS
# ==============================================================================
qso100_db = open_database("/media/data/Endjinn-Extra/swm1n12/agn_angle_hb_100", "root", "password")
qso090_db = open_database("/media/data/Endjinn-Extra/swm1n12/agn_angle_hb_090", "root", "password")
qso110_db = open_database("/media/data/Endjinn-Extra/swm1n12/agn_angle_hb_110", "root", "password")


angles = [20, 30, 40, 50, 60, 75]
lim_qso = 99999999

scale_qso_090 = 1/10
scale_qso_100 = 1/10
scale_qso_110 = 1/10

lines = [
    {"line": 49,  "name":'hb',  "wave": wave_hb},
 #   {"line": 44,  "namaaae":'ha',  "wave": wave_ha},
 #   {"line": 445,  "name":'c4',  "wave": wave_c4},
 #   {"line": 131, "name":'mg2', "wave": 2796.35},
 #   {"line": 383, "name":'he2', "wave": 1640.47}
]

for line in lines:
    delay_tf = np.zeros(len(angles))
    delay_rf = np.zeros(len(angles))
    factor_rf = np.zeros(len(angles))
    factor_tf = np.zeros(len(angles))

    for angle in range(0,len(angles)):
        tf_mid = TransferFunction(qso100_db, "qso_hb_angle_{}".format(angles[angle]), luminosity=1.043e46,  dimensions=dims).delays(0,320)
        tf_mid.line(line["line"], line["wave"]).spectrum(angle).run(scaling_factor=scale_qso_100, delay_dynamic_range=5, limit=lim_qso).plot(velocity=True)
        tf_min = TransferFunction(qso090_db, "qso090", luminosity=0.9391e46, template=tf_mid).run(scaling_factor=scale_qso_090, limit=lim_qso).plot(velocity=True)
        tf_max = TransferFunction(qso110_db, "qso110", luminosity=1.148e46,  template=tf_mid).run(scaling_factor=scale_qso_110, limit=lim_qso).plot(velocity=True)
                                       
        M_actual = 1e9 * apc.M_sun.value

        delay_rf[angle] = tf_mid.response_map_by_tf(tf_min, tf_max,cf_min=1, cf_max=1).plot(velocity=True, response_map=True, name='resp').delay(response=True, threshold=0.8)        
        M_calc = apc.c.value * delay_rf[angle] * np.power(1e3 * tf_mid.FWHM(response=True),2) / apc.G.value
        factor_rf[angle] = (M_actual / M_calc)
        print("R - Delay: {} FWHM: {} M_calc: {} f: {}".format(delay_rf[angle]/seconds_per_day, tf_mid.FWHM(response=True), M_calc, factor_rf[angle]))

        delay_tf[angle] = tf_mid.delay(threshold=0.8)
        M_calc = apc.c.value * delay_tf[angle] * np.power(1e3 * tf_mid.FWHM(response=False),2) / apc.G.value
        factor_tf[angle] = (M_actual / M_calc)      
        print("T - Delay: {} FWHM: {} M_calc: {} f: {}".format(delay_tf[angle]/seconds_per_day, tf_mid.FWHM(response=False), M_calc, factor_tf[angle]))

    np.savetxt("qso_angle_"+line["name"]+".txt", np.column_stack((angles, delay_rf/seconds_per_day, delay_tf/seconds_per_day, factor_rf, factor_tf)), header="Angle d_RF d_TF f_RF f_TF")

sys.exit(1)

def cart2cyl(x,y):
    r = np.sqrt(x*x + y*y)
    theta = np.arctan2(y, x)
    return r, theta

# r, theta, z
print("Querying")
p110 = np.asarray(sey110_db.connect().query(Photon.X, Photon.Y, Photon.Z, Photon.Weight).filter(Photon.Delay < 1*seconds_per_day, Photon.Resonance == 416, Photon.Wavelength > 1550, \
    Photon.Wavelength < doppler_shift_wave(1550, 5e6)).all())
p090 = np.asarray(sey090_db.connect().query(Photon.X, Photon.Y, Photon.Z, Photon.Weight).filter(Photon.Delay < 1*seconds_per_day, Photon.Resonance == 416, Photon.Wavelength > 1550, \
    Photon.Wavelength < doppler_shift_wave(1550, 5e6)).all())
p100 = np.asarray(sey100_db.connect().query(Photon.X, Photon.Y, Photon.Z, Photon.Weight).filter(Photon.Delay < 1*seconds_per_day, Photon.Resonance == 416, Photon.Wavelength > 1550, \
    Photon.Wavelength < doppler_shift_wave(1550, 5e6)).all())

a090 = np.delete(p090, np.s_[-1], axis=1)
a110 = np.delete(p110, np.s_[-1], axis=1)
h090, edges = np.histogramdd(a090, weights=p090[:,3], bins=100)
h110, edges = np.histogramdd(a110, weights=p110[:,3], bins=edges)
hDiff = (h110 - h090) / (50*(1.096e44 - 0.9912e44))
hMask = np.zeros(hDiff.shape)
MaskLim = np.amax(np.abs(hDiff))/1e7 
hMask[hDiff >  MaskLim] = 1
hMask[hDiff < -MaskLim] = 1
gridToVTK("./difference", edges[0], edges[1], edges[2], cellData={"Diff":hDiff, "Mask":hMask})

for photon in p110:
    photon[0], photon[1] = cart2cyl(photon[0], photon[1])
for photon in p090:
    photon[0], photon[1] = cart2cyl(photon[0], photon[1])
for photon in p100:
    photon[0], photon[1] = cart2cyl(photon[0], photon[1])

d100, xedges, yedges = np.histogram2d(p100[:,1], p100[:,0], bins=50, weights=p100[:,3])
ax = plt.subplot(111, projection='polar')
pc = ax.pcolormesh(xedges, yedges, d100, cmap='viridis')
ax.set_rmax(3e15)
ax.set_rlabel("Radius (cm)")
ax.set_rticks([1e15, 2e15, 3e15])  # less radial ticks
ax.grid(True)
cbar = plt.colorbar(pc, orientation="vertical")
cbar.set_label(r"{$\Delta$ L$_{line}$/$\Delta$ C")
plt.savefig("cylindrical.eps", bbox_inches='tight')

d090, xedges, yedges = np.histogram2d(p090[:,1], p090[:,0], bins=50, weights=p090[:,3])
d110, xedges, yedges = np.histogram2d(p110[:,1], p110[:,0], bins=50, weights=p110[:,3])
dDiff = (d110 - d090) / (50*(1.096e44 - 0.9912e44))
clims = np.amax(np.abs(dDiff))
pc = ax.pcolormesh(xedges, yedges, dDiff, cmap='RdBu_r', vmin=(-clims), vmax=clims)
ax.set_rmax(3e15)
ax.set_rlabel("Radius (cm)")
ax.set_rticks([1e15, 2e15, 3e15])  # less radial ticks
ax.grid(True)
cbar = plt.colorbar(pc, orientation="vertical")
cbar.set_label(r"{$\Delta L$_{line}$/$\Delta$ C")
plt.savefig("cylindrical_diff.eps", bbox_inches='tight')
