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

seconds_per_day = 60*60*24

# ==============================================================================       
def calculate_delay(r_angle, r_phase, r_rad, timescale):
    """Delay relative to continuum observed at angle for emission at radius"""
    # Calculate delay compared to continuum photons
    # Draw plane at r_rad_min out. Find x projection of disk position.
    # Calculate distance travelled to that plane from the current disk position
    # Delay relative to continuum is thus (distance from centre to plane)
    # + distance from centre to point
    # vr_disk     = np.array([r_rad.value*np.cos(r_phase), 0.0]) * u.m
    vr_disk     = np.array([r_rad*np.cos(r_phase), 0.0])
    vr_normal   = np.array([np.cos(r_angle), np.sin(r_angle)])
    vr_plane    = r_rad * vr_normal
    # return (np.dot((vr_plane - vr_disk), vr_normal) * u.m / apc.c.value).to(timescale)
    return (np.dot((vr_plane - vr_disk), vr_normal) / apc.c.value) / seconds_per_day
    

def keplerian_velocity(r_mass, r_rad):
    """Calculates Keplerian velocity at radius"""
    return np.sqrt(ap.constants.G.value * r_mass / r_rad)

def path_to_delay(r_path):
    """Converts path to time delay"""
    return r_path / apc.c.value

def doppler_shift_wave(line, vel):
    """Converts passed line and velocity into red/blueshifted wavelength"""
    return line * apc.c.value / (apc.c.value + vel)
    
def doppler_shift_vel(line, wave):
    """Converts passed red/blueshifted wave into velocity"""
    if wave > line:
        return -1*apc.c.value * (1 - (line / wave))
    else:
        return apc.c.value * ((line / wave) - 1)

def clamp(minimum, x, maximum):
    return max(minimum, min(x, maximum))

# ==============================================================================       
class ResponseMap:
    def __init__(self, low_state, low_state_lx, high_state, high_state_lx, line_list, low_scaling_factor=1.0, high_scaling_factor=1.0):
        """Initialises the response map, taking the low and high states and all the lines in them"""
        assert low_state_lx < high_state_lx,\
            "Low and high states are in the wrong order!"
        self._checked_positions = 0
        self._bins_x = np.loadtxt("wind_grid_x.txt")
        self._bins_z = np.loadtxt("wind_grid_z.txt")
        # The last bin is inexplicably 0, so up it
        self._bins_x[0] = self._bins_x[1]-1
        self._bins_z[0] = self._bins_z[1]-1       
       	self._bins_x[-1] = self._bins_x[-2]+1
       	self._bins_z[-1] = self._bins_z[-2]+1
        self._name_lo = low_state
        self._name_hi = high_state

        self._map = {}
        self._dimensions = [len(self._bins_x)-1, len(self._bins_z)-1]
        lx_diff = high_state_lx - low_state_lx

        # For each requested line, read in and reshape arrays then take scaled difference
        for line in line_list:
            temp = Table.read(low_state+line[1], format='ascii')
            self._wind = np.reshape(temp['inwind'], self._dimensions)
            line_lo = np.reshape(temp['var'], self._dimensions) * low_scaling_factor
            temp = Table.read(high_state+line[1], format='ascii')
            line_hi = np.reshape(temp['var'], self._dimensions) * high_scaling_factor
            self._map['L_'+line[0]] = (line_hi - line_lo)/lx_diff  

    def index(self, x_pos, z_pos):
        """Given an position in the wind, returns the position index"""
        x_index = np.searchsorted(self._bins_x,np.abs(x_pos))-1
        z_index = np.searchsorted(self._bins_z,z_pos)-1
        # First and last bins are dummies so clip to real bins
        x_index = clamp(1, x_index, len(self._bins_x)-2)
        z_index = clamp(1, z_index, len(self._bins_z)-2)
        return [x_index, z_index]

    def weight_for(self, line, x_pos, z_pos):
        x_index, z_index = self.index(x_pos, z_pos)
        return self._map["L_"+line][x_index][z_index]

    def plot(self, line, x_range=[None, None], z_range=[None, None]):
        fig, ax1 = plt.subplots(1)
        ax1.set_xlabel('Log X (cm)')
        ax1.set_ylabel('Log Z (cm)')
        if x_range[0] is None:
            x_range[0] = np.log10(self._bins_x[0])
        if x_range[1] is None:
            x_range[1] = np.log10(self._bins_x[-1])
        if z_range[0] is None:
            z_range[0] = np.log10(self._bins_z[0])
        if z_range[1] is None:
            z_range[1] = np.log10(self._bins_z[-1])
        cb_range = np.amax([np.amax(self._map['L_'+line]),np.abs(np.amin(self._map['L_'+line]))])

        ax1.set_xlim(x_range[0], x_range[1])
        ax1.set_ylim(z_range[0], z_range[1])

        im1 = ax1.pcolor(np.log10(self._bins_x), np.log10(self._bins_z), self._map['L_'+line].T, cmap='seismic', vmin=-cb_range, vmax=cb_range)
        cb1 = fig.colorbar(im1, ax=ax1)
        cb1.set_clim(-cb_range, cb_range)
        cb1.set_label(r"${\Delta } L_{line} / {\Delta} L_{ion}$")
        fig.savefig("{}-{}_{}_map.eps".format(self._name_lo, self._name_hi, line))
        plt.close(fig)

# ==============================================================================       
class TransferFunction:    
	# # FAKE INIT
	# def __init(database, ):
	# 	self._query = database.query(Photon.Wavelength, Photon.Delay, Photon.Weight, Photon.X, Photon.Z)

	# def line(wavelength, number):
	# 	self._line_wave = wavelength
	# 	self._query = self._query.filter(Photon.Resonance == number)
	# 	return self
	# def velocities(velocity):
	# 	assert self._line_wave is not None,\
	# 		"Cannot limit doppler shift around a line without specifying a line!"
	# 	self._query = self._query.filter(Photon.Wavelength >= doppler_shift_wave(self._line_wave, -velocity), 
	# 									 Photon.Wavelength <= doppler_shift_wave(self._line_wave,  velocity))
	# 	return self
	# def wavelengths(wave_min, wave_max):
	# 	assert wave_min < wave_max,\
	# 		"Minimum wavelength must be lower than maximum wavelength!"
	# 	self._query = self._query.filter(Photon.Wavelength >= wave_min, Photon.Wavelength <= wave_max)
	# 	return self
	# def lines(line_list):
	# 	assert len(lines) > 1,\
	# 		"For a single line, use the 'line()' filter rather than 'lines()'!"
	# 	self._query = self._query.filter(Photon.Resonance.in_(line_list))
	# 	return self
	# def delays(delay_min, delay_max, unit='d'):
	# 	assert delay_min > delay_max,\
	# 		"Minimum delay must be below maximum delay!"
	# 	if unit is in ['d','D','day','Day','days','Days']:
	# 		self._query=self._query.filter(Photon.Delay > delay_min*seconds_per_day, Photon.Delay < delay_max*seconds_per_day)
	# 	else:
	# 		self._query=self._query.filter(Photon.Delay > delay_min, Photon.Delay < delay_max)
	# 	return self

    def __init__(self, query, filename, continuum_luminosity=None, dimensions=None, template=None):
        """Initialises the TF, taking the query it is to execute, the dimensions to bin by, and the delay range bins"""
        assert dimensions is not None or template is not None,\
            "Must provide either dimensions or another TF to copy them from!"
        
        self._query = query
        self._line_wave = None
        self._delay_range = None
        self.continuum_luminosity=continuum_luminosity
        self._filename = filename
        self._dimensions = dimensions
        self._bins_vel = None
        self._bins_wave = None
        self._bins_delay = None
        self._flux = None
        self._flux_w = None
        self._count = None
        self._response_map = None

        if template is not None:
            self._dimensions = template._dimensions
            self._bins_vel = template._bins_vel
            self._bins_wave = template._bins_wave
            self._bins_delay = template._bins_delay
            if template._line_wave is not None:
                self._line_wave = template._line_wave

    def line(self, line_wave):
        """Flags this TF as being for a single line, and records its wavelength"""
        self._line_wave = line_wave
        return self

    def response_map_by_tf(self, low_state, high_state, plot=False):
        """Creates a response map from two other transfer functions, to be applied during plotting"""
        # The other two TFs ***must*** have identical bins and both provide ionising luminosity information
        assert self._flux is not None,\
            "You must run the TF query with '.run()' before response mapping it!"
        assert low_state._flux is not None and high_state._flux is not None,\
            "You must run the low and high state TF queries with '.run()' before response mapping using them!"
        assert self._flux_w is None,\
            "A response map has already been built!"
        assert np.array_equal(self._bins_wave, low_state._bins_wave) and np.array_equal(self._bins_delay, low_state._bins_delay),\
            "Low state TF is binned differently to target TF! Cannot rescale using it."
        assert np.array_equal(self._bins_wave, high_state._bins_wave) and np.array_equal(self._bins_delay, high_state._bins_delay),\
            "High state TF is binned differently to target TF! Cannot rescale using it."
        assert self.continuum_luminosity != None,\
            "TF missing ionising luminosity information!"
        assert low_state.continuum_luminosity != None,\
            "Low state TF missing ionising luminosity information!"
        assert high_state.continuum_luminosity != None,\
            "High state TF missing ionising luminosity information!"
        assert low_state.continuum_luminosity < self.continuum_luminosity,\
            "Low state ionising luminosity greater than target TF ionising luminosity!"
        assert high_state.continuum_luminosity > self.continuum_luminosity,\
            "High state ionising luminosity lower than target TF ionising luminosity!"

        # If that is true, the map is trivial to construct. We divide the difference in TFs by the luminosity difference
        luminosity_difference = high_state.continuum_luminosity - low_state.continuum_luminosity
        response_map = (high_state._flux - low_state._flux) / luminosity_difference
        self._response_map = response_map
        return self

    def run(self, response_map=None, line=None, scaling_factor=1.0, delay_dynamic_range=2.0):
        """Performs a query on the photon DB and bins it"""
        assert response_map is None or line is not None,\
            "Passing a response map but no line information for it!"
        assert response_map is None or self._flux_w is None,\
            "A response map has already been built!"
        assert response_map is not None,\
            "Response mapping by location not yet implemented!"
        data = np.asarray(self._query.all())
        assert len(data) > 0,\
            "No records found!"

        # Check if we've already got bins from another TF
        if self._bins_wave is None:
            # Data returned as Wavelength, Delay, Weight. Find min and max delays and wavelengths
            range_delay = [0,np.percentile(data[:,1],1 - (10**(-delay_dynamic_range))*100)]
            range_wave = [np.amin(data[:,0]), np.amax(data[:,0])] 

            # Now create the bins for each dimension
            self._bins_vel  = np.linspace(doppler_shift_vel(self._line_wave, range_wave[1]), 
                                          doppler_shift_vel(self._line_wave, range_wave[0]),
                                          self._dimensions[0]+1, endpoint=True, dtype=np.float64)
            self._bins_wave  = np.linspace(range_wave[0], range_wave[1],
                                           self._dimensions[0]+1, endpoint=True, dtype=np.float64)                     
            self._bins_delay = np.linspace(range_delay[0], range_delay[1],
                                           self._dimensions[1]+1, endpoint=True, dtype=np.float64)

            # Convert delay from seconds to days, speed from m/s to km/s
            self._bins_vel = np.true_divide(self._bins_vel, 1000.0)

        # Now we bin the photons, weighting them by their photon weights for the luminosity
        self._flux, junk, junk = np.histogram2d(data[:,1], data[:,0], 
                                                    weights=data[:,2], 
                                                    bins=[self._bins_delay, self._bins_wave]) 
        # Keep an unweighted photon count for statistical error purposes                  
        self._count, junk, junk = np.histogram2d(data[:,1], data[:,0], 
                                                    bins=[self._bins_delay, self._bins_wave])
        # Scaling factor! Each spectral cycle outputs L photons. If we do 50 cycles, we want a factor of 1/50
        self._flux *= scaling_factor
        # Scale to continuum luminosity
        self._flux /= self._continuum_luminosity

        return self
        
    def flux(self, wave, delay):
        """Returns the photon luminosity in this bin"""
        return(self._flux[np.searchsorted(self._bins_delay, delay),
                          np.searchsorted(self._bins_wave, wave)])
    def count(self, wave, delay):
        """Returns the photon count in this bin"""
        return(self._count[np.searchsorted(self._bins_delay, delay),
                           np.searchsorted(self._bins_wave, wave)])

    def plot(self, log=False, dynamic_range=3, normalised=False, rescaled=False, velocity=False, name=None, days=True, 
            response=False, response_map=False, keplerian=None):
        """Takes the data gathered by calling 'run' and outputs a plot"""
        assert response_map_only is False or response is False,\
            "Cannot plot a response function *and* the response map!"
        assert response is False or self._flux_w is not None,\
            "No data available for response function!"
        assert response_map_only is False or self._flux_w is not None,\
            "No data available for response map!"
        assert response_map_only is False or self._response_map is not None,\
            "Response mapping by local luminosity does not produce a plottable map!"
        assert log is False or (response is False and response_map_only is False),\
            "Cannot plot a logarithmic response function or map!"
        assert normalised is False or rescaled is False,\
            "Cannot be both normalised and rescaled!"
        assert self._bins_wave is not None,\
            "You must run the TF query with '.run()' before plotting it!"

        if name is not None:
            print("Plotting to file "+self._filename+"_"+name+".eps")
        else:
            print("Plotting to file "+self._filename+".eps")

        # Set up the multiplot figure and axis
        fig, ((ax_spec, ax_none), (ax_tf, ax_resp)) = plt.subplots(2,2,sharex='col', sharey='row',
            gridspec_kw={'width_ratios':[3,1], 'height_ratios':[1,3]})
        fig.subplots_adjust(hspace=0, wspace=0)

        # Set the properties that depend on log and wave/velocity status
        cb_label = r"$\psi_{em}$"
        cb_label_vars = r""
        cb_label_units = r""
        cb_label_scale= r""
        cb_map = "afmhot_r"

        # Copy the data for later modification.
        data_plot = np.copy(self._flux)
        data_plot_w = None
        if response:
            # If this is a response plot, we need the weighted data too
            data_plot_w = np.copy(self._flux_w)
            cb_label = r"$\psi$"
        elif response_map_only:
            data_plot = np.copy(self._response_map)
            cb_label = r"$\Delta L_{line}$ / $\Delta L_{ion}$"

        # Set the xlabel and colour bar label - these differ if velocity or not
        x_bin_mult = 1
        if velocity:
            oom = np.log10(np.amax(self._bins_vel))
            oom = oom - oom%3
            # We're rescaling the axis to e.g. 10^3 km/s but the colorbar is still in km/s
            # So when we scale by bin width, we need a multiplier on the bin widths
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
                if response:
                    data_plot_w[i][j] /= (width_x * x_bin_mult * width_y)

        # If this is a log plot, take log and correct label and limits
        if log:
            cb_max = np.ma.log10(cb_max)
            cb_min = cb_max-dynamic_range
            cb_label = r"Log "+cb_label
            data_plot = np.log10(data_plot)
        # If this is a response plot, scale the weighted data
        elif response:
            maxval = np.floor(np.log10(np.amax(data_plot_w)))
            data_plot_w /= np.power(10,maxval)
            cb_max = np.amax(data_plot_w)
            cb_min = np.amin(data_plot_w)
            dummy = "{}{:.0f}{}".format("{",maxval,"}")
            cb_label_scale = r" 10$^{}$".format(dummy)
        # If this is not a response plot, scale the unweighted data
        else:
            maxval = np.floor(np.log10(np.amax(data_plot)))
            data_plot /= np.power(10,maxval)
            cb_max = np.amax(data_plot)
            cb_min = np.amin(data_plot)
            dummy = "{}{:.0f}{}".format("{",maxval,"}")
            cb_label_scale = r" 10$^{}$".format(dummy)        

        # If this is a response function, it may have a negative component and need a different plot
        if response or response_map:
            cb_max = np.amax([cb_max, np.abs(cb_min)])
            cb_min = -cb_max
            cb_map = 'seismic'

        # Normalise or rescale the data. If doing neither, put units on cb.
        if normalised:
            data_plot /= np.sum(data_plot)
            if response:
                data_plot_w /= np.sum(data_plot_w)      
            cb_label_units = r""
            cb_label_scale = r""
        elif rescaled:
            data_plot /= np.amax(data_plot)
            if response:
                data_plot_w /= np.amax(data_plot_w)
            cb_label_units = r""
            cb_label_scale = r""
        elif response_map_only:
            cb_label_scale = r""
            cb_label_units = r""

        # Plot the main colourplot for the transfer function
        tf = None
        if response:
            tf = ax_tf.pcolor(bins_x, bins_y, data_plot_w,
                             vmin=cb_min, vmax=cb_max, cmap=cb_map)
        else:
            tf = ax_tf.pcolor(bins_x, bins_y, data_plot,
                             vmin=cb_min, vmax=cb_max, cmap=cb_map)
        ax_tf.set_ylim(top=bins_y[-1])
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
            ar_rad      = np.linspace(keplerian["rad_min"]*r_rad_grav, keplerian["rad_max"]*r_rad_grav,resolution)
            ar_vel      = np.zeros(resolution) 

            # ITERATE OVER BLUE BOUND
            for r_rad, r_wave, r_delay, r_vel in np.nditer([ar_rad, ar_wave, ar_delay, ar_vel], op_flags=['readwrite']):
                r_rad           = r_rad # * u.m
                r_vel[...]      = keplerian_velocity( r_mass_bh, r_rad ) * np.sin(r_angle) / (1e6)
                r_wave[...]     = doppler_shift_wave( self._line_wave, r_vel )
                r_delay[...]    = calculate_delay( r_angle, np.pi/2, r_rad, u.day )
            ax_tf.plot(ar_vel, ar_delay, '-', c='C0')

            # ITERATE OVER RED BOUND
            for r_rad, r_wave, r_delay, r_vel in np.nditer([ar_rad, ar_wave, ar_delay, ar_vel], op_flags=['readwrite']):
                r_rad           = r_rad # * u.m
                r_vel[...]      = -keplerian_velocity( r_mass_bh, r_rad ) * np.sin(r_angle) / (1e6)
                r_wave[...]     = doppler_shift_wave( self._line_wave, r_vel )
                r_delay[...]    = calculate_delay( r_angle, np.pi/2, r_rad, u.day )
            ax_tf.plot(ar_vel, ar_delay, '-', c='C0')

        # Plot the spectrum and light curve, normalised
        data_plot_spec = np.sum(data_plot, 0) / np.sum(data_plot)
        data_plot_resp = np.sum(data_plot, 1) / np.sum(data_plot)
        spec = ax_spec.plot(bins_x_midp, data_plot_spec, c='C0', label='Emissivity')
        resp = ax_resp.plot(data_plot_resp, bins_y_midp, c='C0', label='_nolegend_')

        # If this is a response plot, include the emissivity spectra and light curve for comparison
        if response:
            data_plot_spec_w = np.sum(data_plot_w, 0) / np.sum(data_plot_w)
            data_plot_resp_w = np.sum(data_plot_w, 1) / np.sum(data_plot_w)
            spec_w = ax_spec.plot(bins_x_midp, data_plot_spec_w, c='C1', label='Response')
            resp_w = ax_resp.plot(data_plot_resp_w, bins_y_midp, c='C1', label='_nolegend_')
            lg_orig = ax_spec.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)

        # Clear the ticks as they're a mess
        ax_spec.tick_params(labelbottom='off', labelleft='off', labeltop='off', labelright='off', left='off', bottom='off')
        ax_spec.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter(''))
        ax_resp.tick_params(labelbottom='off', labelleft='off', labeltop='off', labelright='off', bottom='off', left='off')
        ax_resp.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter(''))
        ax_none.axis('off')
        cbar = plt.colorbar(tf, orientation="vertical")
        cbar.set_label(cb_label+cb_label_vars+cb_label_scale+cb_label_units)

        if name is None:
            plt.savefig("{}.eps".format(self._filename),bbox_inches='tight')
        else:
            plt.savefig("{}_{}.eps".format(self._filename, name),bbox_inches='tight')
        plt.close(fig)

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
    print ("Searching for table 'Photons'...")
    Session = sqlalchemy.orm.sessionmaker(bind=db_engine)
    dbc = Session()

    start = time.clock()

    try:
        dbc.query(Photon.Weight).first()
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
            values = line.split()
            if len(values) < 12:
                continue
            
            for index,value in enumerate(values):
                values[index] = float(value)
            del values[0]
            del values[8]
            values[5] = (values[5] - values[6])
            matom_bool = False
            if(values[9]>= 10):
                values[9] = values[9] - 10
                matom_bool = True

            dbc.add(Photon(Wavelength=values[0], Weight=values[1], X=values[2], Y=values[3], Z=values[4],
                            ContinuumScatters=values[5], ResonantScatters=values[6], Delay=values[7],
                            Spectrum=values[8], Origin=values[9], Resonance=values[10], Origin_matom = matom_bool))
            added += 1
            if added > 10000:
                added = 0
                dbc.commit()
        dbc.commit()
        print("Input took {:.2f} seconds".format(time.clock()-start))

    return db_engine, dbc
# ==============================================================================       

#s_root = "ngc5548_1e00_obs_2"
s_user = "root"
s_password = "password"
#s_file = "test_large"

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


sey100_engine, sey100_dbc = open_database("/Users/swm1n12/python_runs/paper1_5548_resp/sey_100", "root", "password")
sey_lo_engine, sey_lo_dbc = open_database("/Users/swm1n12/python_runs/paper1_5548_resp/sey_090", "root", "password")
sey_hi_engine, sey_hi_dbc = open_database("/Users/swm1n12/python_runs/paper1_5548_resp/sey_110", "root", "password")

# sey100_engine, sey100_dbc = open_database("/Users/swm1n12/python_runs/paper1_5548_resp/sey_100", "root", "password")
# sey095_engine, sey095_dbc = open_database("/Users/swm1n12/python_runs/paper1_5548_resp/sey_095", "root", "password")
# sey105_engine, sey105_dbc = open_database("/Users/swm1n12/python_runs/paper1_5548_resp/sey_105", "root", "password")

# sey100_engine, sey100_dbc = open_database("sey_100", "root", "password")
# sey095_engine, sey095_dbc = open_database("sey_095", "root", "password")
# sey105_engine, sey105_dbc = open_database("sey_105", "root", "password")

# should be around 1475 - 1625 
c4_wave_lim = [doppler_shift_wave(1550, 12e6), doppler_shift_wave(1550, -12e6)]
ha_wave_lim = [doppler_shift_wave(6563, 12e6), doppler_shift_wave(6563, -12e6)]
dims = [50, 50]

# ------------------------------------------------------------------------------
sey100_query_ha = sey100_dbc.query(Photon.Wavelength, Photon.Delay, Photon.Weight, Photon.X, Photon.Z)
sey100_query_ha = sey100_query_ha.filter(and_(Photon.Wavelength > ha_wave_lim[0], Photon.Wavelength < ha_wave_lim[1]))
sey100_query_ha = sey100_query_ha.filter(and_(Photon.Delay < (20 * seconds_per_day), Photon.Resonance == 28))
sey100_tf_ha = TransferFunction(sey100_query_ha, "sey_100", dimensions=dims, continuum_luminosity=1.043e44).line(6563).run(scaling_factor=0.2)

sey_lo_query_ha = sey_lo_dbc.query(Photon.Wavelength, Photon.Delay, Photon.Weight)
sey_lo_query_ha = sey_lo_query_ha.filter(and_(Photon.Wavelength > ha_wave_lim[0], Photon.Wavelength < ha_wave_lim[1]))
sey_lo_query_ha = sey_lo_query_ha.filter(and_(Photon.Delay < (20 * seconds_per_day), Photon.Resonance == 28))
sey_lo_tf_ha = TransferFunction(sey_lo_query_ha, "sey_090", continuum_luminosity=0.9391e44, template=sey100_tf_ha).run()
sey_lo_tf_ha.plot(velocity=True, name="ha_vel")

sey_hi_query_ha = sey_hi_dbc.query(Photon.Wavelength, Photon.Delay, Photon.Weight)
sey_hi_query_ha = sey_hi_query_ha.filter(and_(Photon.Wavelength > ha_wave_lim[0], Photon.Wavelength < ha_wave_lim[1]))
sey_hi_query_ha = sey_hi_query_ha.filter(and_(Photon.Delay < (20 * seconds_per_day), Photon.Resonance == 28))
sey_hi_tf_ha = TransferFunction(sey_hi_query_ha, "sey_110", continuum_luminosity=1.148e44, template=sey100_tf_ha).run()

sey100_tf_ha.response_map_by_tf(sey_lo_tf_ha, sey_hi_tf_ha)
sey100_tf_ha.plot(velocity=True, name="ha_vel")
sey100_tf_ha.plot(velocity=True, name="ha_vel_respmap", response_map_only=True)

sey100_tf_ha = TransferFunction(sey100_query_ha, "sey_100", dimensions=dims, continuum_luminosity=1.043e44).line(6563).run(scaling_factor=1/20, delay_dynamic_range=3)
sey100_tf_ha.plot(velocity=True, name="ha_vel_long")
sey_lo_tf_ha = TransferFunction(sey_lo_query_ha, "sey_090", continuum_luminosity=9.391e43, template=sey100_tf_ha).run(scaling_factor=1/20)
sey_hi_tf_ha = TransferFunction(sey_hi_query_ha, "sey_110", continuum_luminosity=1.148e44, template=sey100_tf_ha).run(scaling_factor=1/20)
sey100_tf_ha.response_map_by_tf(sey_lo_tf_ha, sey_hi_tf_ha)
sey100_tf_ha.plot(velocity=True, name="ha_vel_long_respmap", response_map_only=True)

# ------------------------------------------------------------------------------
sey100_query_c4 = sey100_dbc.query(Photon.Wavelength, Photon.Delay, Photon.Weight, Photon.X, Photon.Z)
sey100_query_c4 = sey100_query_c4.filter(and_(Photon.Wavelength > c4_wave_lim[0], Photon.Wavelength < c4_wave_lim[1]))
sey100_query_c4 = sey100_query_c4.filter(and_(Photon.Delay < (20 * seconds_per_day), Photon.Resonance == 416))
sey100_tf_c4 = TransferFunction(sey100_query_c4, "sey_100", dimensions=dims, continuum_luminosity=1.043e44).line(1550).run(scaling_factor=1/20)
sey100_tf_c4.plot(velocity=True, name="c4_vel")

sey_lo_query_c4 = sey_lo_dbc.query(Photon.Wavelength, Photon.Delay, Photon.Weight)
sey_lo_query_c4 = sey_lo_query_c4.filter(and_(Photon.Wavelength > c4_wave_lim[0], Photon.Wavelength < c4_wave_lim[1]))
sey_lo_query_c4 = sey_lo_query_c4.filter(and_(Photon.Delay < (20 * seconds_per_day), Photon.Resonance == 416))
sey_lo_tf_c4 = TransferFunction(sey_lo_query_c4, "sey_090", continuum_luminosity=0.9391e44, template=sey100_tf_c4).run(scaling_factor=1/20)
sey_lo_tf_c4.plot(velocity=True, name="vel")

sey_hi_query_c4 = sey_hi_dbc.query(Photon.Wavelength, Photon.Delay, Photon.Weight)
sey_hi_query_c4 = sey_hi_query_c4.filter(and_(Photon.Wavelength > c4_wave_lim[0], Photon.Wavelength < c4_wave_lim[1]))
sey_hi_query_c4 = sey_hi_query_c4.filter(and_(Photon.Delay < (20 * seconds_per_day), Photon.Resonance == 416))
sey_hi_tf_c4 = TransferFunction(sey_hi_query_c4, "sey_110", continuum_luminosity=1.148e44, template=sey100_tf_c4).run(scaling_factor=1/20)
sey_hi_tf_c4.plot(velocity=True, name="vel")

sey100_tf_c4.response_map_by_tf(sey_lo_tf_c4, sey_hi_tf_c4)
sey100_tf_c4.plot(velocity=True, name="c4_vel_respmap", response_map_only=True)

sey100_tf_c4 = TransferFunction(sey100_query_c4, "sey_100", dimensions=dims, continuum_luminosity=1.043e44).line(1550).run(scaling_factor=1/20, delay_dynamic_range=3)
sey100_tf_c4.plot(velocity=True, name="c4_vel_long")
sey_lo_tf_c4 = TransferFunction(sey_lo_query_c4, "sey_090", continuum_luminosity=9.391e43, template=sey100_tf_c4).run(scaling_factor=1/20)
sey_hi_tf_c4 = TransferFunction(sey_hi_query_c4, "sey_110", continuum_luminosity=1.148e44, template=sey100_tf_c4).run(scaling_factor=1/20)
sey100_tf_c4.response_map_by_tf(sey_lo_tf_c4, sey_hi_tf_c4)
sey100_tf_c4.plot(velocity=True, name="c4_vel_long_respmap", response_map_only=True)

# line_list = [   ["Ha", ".lineH1.3-2.dat"],
#                 ["Hb", ".lineH1.4-2.dat"],
# #                ["Lb", ".lineH1.3-1.dat"],
# #                ["Lg", ".lineH1.4-1.dat"],
#                 ["C4", ".lineC4.dat"]   ]
# map_95_to_105 = ResponseMap("sey_095", 0.95e43, "sey_105", 1.05e43, line_list)
# map_95_to_105.plot('C4', x_range=[None, 15.7], z_range=[None, 15.2])
# map_95_to_105.plot('Ha', x_range=[None, 15.7], z_range=[None, 15.2])
# map_95_to_105.plot('Hb', x_range=[None, 15.7], z_range=[None, 15.2])
# sey100_tf2_c4 = TransferFunction(sey100_query_c4, "sey_100", dimensions=dims, ionising_luminosity=1.00e43).line(1550)
# sey100_tf2_c4.run(map_95_to_105, line="C4", scaling_factor=0.2)
# keplerian_info = { "mass":1e7, "rad_min":50, "rad_max":1000, "angle":40}
# sey100_tf2_c4.plot(velocity=True, name="c4_vel")
# sey100_tf2_c4.plot(velocity=True, name="c4_vel_keplerian", keplerian=keplerian_info)
# sey100_tf2_c4.plot(velocity=True, name="c4_vel_resp_local", response=True)
