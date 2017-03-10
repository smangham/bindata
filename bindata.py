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


seconds_per_day = 60*60*24

# ==============================================================================       
def doppler_shift_wave(wave, vel):
    """Doppler shifts the passed wavelength"""
    return r_wave * (1 + (vel / apc.c.value))

def doppler_shift_vel(wave, line):
    """Converts passed wavelength into velocity"""
    return apc.c.value * ((wave/line) - 1)

# ==============================================================================       
class ResponseMap:
    def __init__(self, low_state, low_state_lx, high_state, high_state_lx, line_list, low_scaling_factor=1.0, high_scaling_factor=1.0):
        """Initialises the response map, taking the low and high states and all the lines in them"""
        assert low_state_lx < high_state_lx,\
            "Low and high states are in the wrong order!"

        self._bins_x = np.loadtxt("wind_grid_x.txt")
        self._bins_z = np.loadtxt("wind_grid_z.txt")
        self._map = {}
        self._dimensions = [len(self._bins_x), len(self._bins_z)]
        lx_diff = high_state_lx - low_state_lx

        # For each requested line, read in and reshape arrays then take scaled difference
        for line in line_list:
            temp = Table.read(low_state+line[1], format='ascii')
            line_lo = np.reshape(temp['var'], self._dimensions) * low_scaling_factor
            temp = Table.read(high_state+line[1], format='ascii')
            line_hi = np.reshape(temp['var'], self._dimensions) * high_scaling_factor
            self._map['L_'+line[0]] = (line_hi - line_lo)/lx_diff  
    
    def index(self, x_pos, z_pos):
        """Given an position in the wind, returns the position index"""
        x_index = np.searchsorted(self._bins_x,x)-1
        z_index = np.searchsorted(self._bins_z,z)-1
        if x_index < 0: 
            x_index = 0
        if z_index < 0: 
            z_index = 0
        return [x_index, z_index]

# ==============================================================================       
class TransferFunction:    
    def __init__(self, query, filename, delay_dynamic_range=2.0, ionising_luminosity=None, dimensions=None, bins_from=None):
        """Initialises the TF, taking the query it is to execute, the dimensions to bin by, and the delay range bins"""
        assert dimensions is not None or bins_from is not None,\
            "Must provide either dimensions or another TF to copy them from!"
        
        self._query = query
        self._line_wave = None
        self._delay_range = 1 - (10**(-delay_dynamic_range))

        self._ionising_luminosity=ionising_luminosity
        self._filename = filename
        if bins_from is not None:
            self._dimensions = bins_from._dimensions
            self._bins_vel = bins_from._bins_vel
            self._bins_wave = bins_from._bins_wave
            self._bins_delay = bins_from._bins_delay
        else:
            self._dimensions = dimensions
            self._bins_vel = None
            self._bins_wave = None
            self._bins_delay = None
        self._flux = None
        self._flux_w = None
        self._count = None
        self._response_map = None

    def line(self, line_wave):
        """Flags this TF as being for a single line, and records its wavelength"""
        self._line_wave = line_wave
        return self

    def response_map_by_tf(self, low_state, high_state, plot=False):
        """Creates a responsivity map from two other transfer functions, to be applied during plotting"""
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
        assert self._ionising_luminosity != None,\
            "TF missing ionising luminosity information!"
        assert low_state._ionising_luminosity != None,\
            "Low state TF missing ionising luminosity information!"
        assert high_state._ionising_luminosity != None,\
            "High state TF missing ionising luminosity information!"
        assert low_state._ionising_luminosity < self._ionising_luminosity,\
            "Low state ionising luminosity greater than target TF ionising luminosity!"
        assert high_state._ionising_luminosity > self._ionising_luminosity,\
            "High state ionising luminosity lower than target TF ionising luminosity!"

        # If that is true, the map is trivial to construct. We divide the difference in TFs by the luminosity difference
        luminosity_difference = high_state._ionising_luminosity - low_state._ionising_luminosity
        response_map = (high_state._flux - low_state._flux) / luminosity_difference
   
        # Multiply by the response map
        self._flux_w = np.multiply(self._flux, response_map)
        self._response_map = response_map

        return self

    def run(self, response_map=None, line=None, scaling_factor=1.0):
        """Performs a query on the photon DB and bins it"""
        data = np.asarray(self._query.all())
        update_bins = False
        assert len(data) > 0,\
            "No records found!"
        assert response_map is None or line is not None,\
            "Passing a response map but no line information for it!"
        assert response_map is None or self._flux_w is None,\
            "A response map has already been built!"

        # Check if we've already got bins from another TF
        if self._bins_wave is None:
            # Data returned as Wavelength, Delay, Weight. Find min and max delays and wavelengths
            range_delay = [0,np.percentile(data[:,1],self._delay_range*100)]
            range_wave = [np.amin(data[:,0]), np.amax(data[:,0])] 

            # Now create the bins for each dimension
            self._bins_vel  = np.linspace(doppler_shift_vel(range_wave[0], self._line_wave), 
                                          doppler_shift_vel(range_wave[1], self._line_wave),
                                          self._dimensions[0]+1, endpoint=True, dtype=np.float64)
            self._bins_wave  = np.linspace(range_wave[0], range_wave[1],
                                           self._dimensions[0]+1, endpoint=True, dtype=np.float64)                     
            self._bins_delay = np.linspace(range_delay[0], range_delay[1],
                                           self._dimensions[1]+1, endpoint=True, dtype=np.float64)

            # Convert delay from seconds to days, speed from m/s to km/s
            self._bins_vel = np.true_divide(self._bins_vel, 1000.0)
            self._bins_delay_d = np.true_divide(self._bins_delay, float(seconds_per_day))

        # Now we bin the photons, weighting them by their photon weights for the luminosity
        self._flux, junk, junk = np.histogram2d(data[:,1], data[:,0], 
                                                    weights=data[:,2], 
                                                    bins=[self._bins_delay, self._bins_wave]) 
        # Keep an unweighted photon count for statistical error purposes                  
        self._count, junk, junk = np.histogram2d(data[:,1], data[:,0], 
                                                    bins=[self._bins_delay, self._bins_wave])
        # If there are any points with no response, set a dummy value, then copy to velocity version
        self._flux[self._flux == 0] = 1e-34
        # Scaling factor! Each spectral cycle outputs L photons. If we do 50 cycles, we want a factor of 1/50
        self._flux *= scaling_factor

        # If we've been given a response map, keep a weighted array too
        if response_map is not None:
            response_weights = np.ones(len(data))
            for index, photon in enumerate(data):
                x, z = response_map.index(photon[3], photon[4])
                response_weights[index] = response_map._map['L_'+line][x][z]

            self._flux_w = np.zeros(shape=self._dimensions, dtype=np.float64)
            self._flux_w, junk, junk = np.histogram2d(data[:,1], data[:,0], 
                                                        weights=data[:,2]*response_weights, 
                                                        bins=[self._bins_delay, self._bins_wave])    
            self._flux_w[self._flux_w == 0] = 1e-34   
            self._flux_w *= scaling_factor

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
            response=False, response_map_only=False):
        """Takes the data gathered by calling 'run' and outputs a plot"""
        assert response_map_only is False or response is False,\
            "Cannot plot a response function *and* the response map!"
        assert response is False or self._flux_w is not None,\
            "No data available for response function!"
        assert response_map_only is False or self._flux_w is not None,\
            "No data available for response map!"
        assert log is False or (response is False and response_map_only is False),\
            "Cannot plot a logarithmic response function or map!"
        assert normalised is False or rescaled is False,\
            "Cannot be both normalised and rescaled!"
        assert self._bins_wave is not None,\
            "You must run the TF query with '.run()' before plotting it!"

        # Set up the multiplot figure and axis
        fig, ((ax_spec, ax_none), (ax_tf, ax_resp)) = plt.subplots(2,2,sharex='col', sharey='row',
            gridspec_kw={'width_ratios':[3,1], 'height_ratios':[1,3]})
        fig.subplots_adjust(hspace=0, wspace=0)

        # Set the properties that depend on log and wave/velocity status
        cb_label = r"$\psi_{em}$"
        cb_label_units = ""
        cb_label_scale= ""
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
            cb_label = r"$\delta L_{line}$ / $\delta L_{ion}$"

        # Set the xlabel and colour bar label - these differ if velocity or not
        if velocity:
            oom = np.log10(np.amax(self._bins_vel))
            oom = oom - oom%3
            bins_x = self._bins_vel/(10**oom)
            ax_tf.set_xlabel(r'Velocity ($10^{:.0f}$ km s$^{}$)'.format(oom, '{-1}'))
            cb_label += r"($v, \tau$)"
            cb_label_units = r" erg s$^{-1}$ / km s$^{-1}$"
        else:
            bins_x = self._bins_wave
            ax_tf.set_xlabel(r'Wavelength $\lambda$ ($\AA$)')
            cb_label += r"($\lambda, \tau$)"
            cb_label_units = r" erg s$^{-1}$ / $\AA$"

        bins_x_midp = np.zeros(shape=self._dimensions[0])
        for i in range(0,self._dimensions[0]):
            bins_x_midp[i] = (bins_x[i] + bins_x[i+1])/ 2

       # Set the ylabel and y bins
        if days:
            bins_y = np.true_divide(self._bins_delay, float(seconds_per_day))
            ax_tf.set_ylabel(r'Delay $\tau$ (days)')
            cb_label_units += r' day'
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
                data_plot[i][j] /= (width_x * width_y)
                if response:
                    data_plot_w[i][j] /= (width_x * width_y)

        # If this is a log plot, take log and correct label and limits
        if log:
            cb_max = np.log10(cb_max)
            cb_min = cb_max-dynamic_range
            cb_label = r"Log "+cb_label
            data_plot = np.log10(data_plot)
        elif response:
            maxval = np.floor(np.log10(np.amax(data_plot_w)))
            data_plot_w /= np.power(10,maxval)
            cb_max = np.amax(data_plot_w)
            cb_min = np.amin(data_plot_w)
            dummy = "{}{:.0f}{}".format("{",maxval,"}")
            cb_label_scale = r" 10$^{}$".format(dummy)
        else:
            maxval = np.floor(np.log10(np.amax(data_plot)))
            data_plot /= np.power(10,maxval)
            cb_max = np.amax(data_plot)
            cb_min = np.amin(data_plot)
            dummy = "{}{:.0f}{}".format("{",maxval,"}")
            cb_label_scale = r" 10$^{}$".format(dummy)        

        # If this is a response function, it may have a negative component and need a different plot
        if response or response_map_only:
            cb_max = np.amax([cb_max, np.abs(cb_min)])
            cb_min = -cb_max
            # Default RdBu colourmap goes via grey; we want it to go via white
            cdict = {'red':   ((0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)),
                     'green': ((0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0)),
                     'blue':  ((0.0, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0))}
            cb_map = mcolors.LinearSegmentedColormap('CustomMap', cdict)
            if response:
                cb_label_units += r" erg"

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
        cbar.set_label(cb_label+cb_label_scale+cb_label_units)

        if name is None:
            plt.savefig("{}.eps".format(self._filename))
        else:
            plt.savefig("{}_{}.eps".format(self._filename, name))

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
            #print("Reading line: {}".format(line))
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
sey095_engine, sey095_dbc = open_database("/Users/swm1n12/python_runs/paper1_5548_resp/sey_095", "root", "password")
sey105_engine, sey105_dbc = open_database("/Users/swm1n12/python_runs/paper1_5548_resp/sey_105", "root", "password")

sey100_query = sey100_dbc.query(Photon.Wavelength, Photon.Delay, Photon.Weight)
sey100_query = sey100_query.filter(and_(Photon.Wavelength > 1475, Photon.Wavelength < 1625))
sey100_query = sey100_query.filter(and_(Photon.Delay < (3 * seconds_per_day)))
sey100_tf = TransferFunction(sey100_query, "sey_100", dimensions=[100,100], ionising_luminosity=1.00e43).line(1550).run(scaling_factor=0.2)

sey095_query = sey095_dbc.query(Photon.Wavelength, Photon.Delay, Photon.Weight)
sey095_query = sey095_query.filter(and_(Photon.Wavelength > 1475, Photon.Wavelength < 1625))
sey095_query = sey095_query.filter(and_(Photon.Delay < (3 * seconds_per_day)))
sey095_tf = TransferFunction(sey095_query, "sey_095", ionising_luminosity=0.95e43, bins_from=sey100_tf).line(1550).run()
sey095_tf.plot(velocity=True, name="vel")

sey105_query = sey105_dbc.query(Photon.Wavelength, Photon.Delay, Photon.Weight)
sey105_query = sey105_query.filter(and_(Photon.Wavelength > 1475, Photon.Wavelength < 1625))
sey105_query = sey105_query.filter(and_(Photon.Delay < (3 * seconds_per_day)))
sey105_tf = TransferFunction(sey105_query, "sey_105", ionising_luminosity=1.05e43, bins_from=sey100_tf).line(1550).run()
sey105_tf.plot(velocity=True, name="vel")

sey100_tf.response_map_by_tf(sey095_tf, sey105_tf)
sey100_tf.plot(velocity=True, name="vel")
sey100_tf.plot(velocity=True, name="vel_resp_tf", response=True)
sey100_tf.plot(velocity=True, name="vel_respmap", response_map_only=True)


line_list = [   ["Hα", ".lineH1.3-2.dat"],
                ["Hβ", ".lineH1.4-2.dat"],
                ["Lβ", ".lineH1.3-1.dat"],
                ["Lγ", ".lineH1.4-1.dat"],
                ["C4", ".lineC4.dat"]   ]
map_95_to_105 = ResponseMap("sey_095", 0.95e43, "sey_105", 1.05e43, line_list)
sey100_tf2 = TransferFunction(sey100_query, "sey_100", dimensions=[100,100], ionising_luminosity=1.00e43, scaling_factor=0.2).line(1550)
sey100_tf2.run(map_95_to_105, "C4")
sey100_tf2.plot(velocity=True, name="vel_resp_local")

#query = query.filter(and_(or_(Photon.ResonantScatters > 0, Photon.Origin == 3, Photon.Origin_matom == True)))
#query = query.filter(and_(Photon.Resonance.in_([380, 381, 382, 416])))

