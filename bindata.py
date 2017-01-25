# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import HDFql
import numpy as np
import astropy as ap
import astropy.units
import astropy.constants as apc
import time
import sys
import scipy
import matplotlib.pyplot as plt
import sqlalchemy
import sqlalchemy.ext.declarative
import sqlalchemy.orm
import sqlalchemy.orm.query
import matplotlib
from sqlalchemy import and_, or_

seconds_per_day = 60*60*24

def doppler_shift_wave(wave, vel):
    """Doppler shifts the passed wavelength"""
    return r_wave * (1 + (vel / apc.c.value))

def doppler_shift_vel(wave, line):
    """Converts passed wavelength into velocity"""
    return apc.c.value * ((wave/line) - 1)
       
class TransferFunction:    
    def __init__(self, query, dimensions, delay_dynamic_range=2.0):
        self._query = query
        self._dimensions = dimensions
        self._flux = np.zeros(shape=dimensions, dtype=np.float64)
        self._count = np.zeros(shape=dimensions, dtype=np.intc)
        self._bins_vel = None
        self._bins_vel_midp = None
        self._bins_wave = None
        self._bins_wave_midp = None
        self._bins_delay = None
        self._bins_delay_midp = None
        self._range_velocity = [None, None]
        self._line_wave = None
        self._delay_range = 1 - (10**(-delay_dynamic_range))
        self._as_velocity = False
        self._flux_spec = None
        self._flux_resp = None

    def line(self, line_wave):
        self._line_wave = line_wave
        return self

    def as_velocity(self, vel_min=None, vel_max=None):
        self._range_velocity = [vel_min, vel_max]
        self._as_velocity = True
        return self

    def run(self):
        data = np.asarray(self._query.all())
        if len(data) == 0:
            print("Error: No records found!")
            return

        # Data returned as Wavelength, Delay, Weight
        range_delay = [0,np.percentile(data[:,1],self._delay_range*100)]
        range_wave = [np.amin(data[:,0]), np.amax(data[:,0])] 
        #If this is expected to be a velocity plot, convert 
        if self._range_velocity[0] != None:
           range_wave[0] = doppler_shift_wave(self._line_wave, self._range_velocity[0])
        if self._range_velocity[1] != None:
           range_wave[1] = doppler_shift_wave(self._line_wave, self._range_velocity[1])

        self._bins_vel  = np.linspace(doppler_shift_vel(range_wave[0], self._line_wave), 
                                      doppler_shift_vel(range_wave[1], self._line_wave),
                                      self._dimensions[0]+1, endpoint=True, dtype=np.float64)
        self._bins_wave  = np.linspace(range_wave[0], range_wave[1],
                                       self._dimensions[0]+1, endpoint=True, dtype=np.float64)                     
        self._bins_delay = np.linspace(range_delay[0], range_delay[1],
                                       self._dimensions[0]+1, endpoint=True, dtype=np.float64)

        self._flux, junk, junk = np.histogram2d(data[:,1], data[:,0], weights=data[:,2], 
                       bins=[self._bins_delay, self._bins_wave], normed=True)                    
        self._count, junk, junk = np.histogram2d(data[:,1], data[:,0], 
                       bins=[self._bins_delay, self._bins_wave], normed=True)
        self._flux_min = np.amin(self._flux[np.nonzero(self._flux)])
        self._flux[self._flux == 0] = 1e-35
        self._flux_spec = np.sum(self._flux, 0)
        self._flux_resp = np.sum(self._flux, 1)

        # Convert delay from seconds to days, speed from m/s to km/s
        self._bins_delay = np.true_divide(self._bins_delay, float(seconds_per_day))
        self._bins_vel = np.true_divide(self._bins_vel, 1000.0)

        self._bins_wave_midp = np.zeros(shape=self._dimensions)
        self._bins_delay_midp = np.zeros(shape=self._dimensions)
        self._bins_vel_midp = np.zeros(shape=self._dimensions)
        for i in range(0,self._dimensions[0]):
            self._bins_wave_midp[i] = (self._bins_wave[i] + self._bins_wave[i+1])/ 2
            self._bins_vel_midp[i] = (self._bins_vel[i] + self._bins_vel[i+1])/ 2
        for i in range(0,self._dimensions[0]):
            self._bins_delay_midp[i] = (self._bins_wave[i] + self._bins_delay[i+1])/ 2

        return self
        
    def flux(self, wave, delay):
        return(self._flux[np.searchsorted(self._bins_delay, delay),
                          np.searchsorted(self._bins_wave, wave)])
    def count(self, wave, delay):
        return(self._count[np.searchsorted(self._bins_delay, delay),
                           np.searchsorted(self._bins_wave, wave)])

    def plot(self, dynamic_range=5):
        cb_max = np.log10(np.amax(self._flux))
        cb_min = cb_max-dynamic_range
        
        fig, ((ax_spec, ax_none), (ax_tf, ax_resp)) = plt.subplots(2,2,sharex='col', sharey='row',
            gridspec_kw={'width_ratios':[3,1], 'height_ratios':[1,3]})
        fig.subplots_adjust(hspace=0, wspace=0)

        cb_units = None
        if self._as_velocity:
            oom = np.log10(np.amax(self._bins_vel))
            oom = oom - oom%3
            bins_x = self._bins_vel/(10**oom)
            bins_x_midp = self._bins_vel_midp/(10**oom)
            ax_tf.set_xlabel('Velocity ($10^{:.0f}$ km s$^{}$)'.format(oom, '{-1}').encode('string-escape'))
            cb_units = r"Log $\psi$($v, \tau$)"

        else:
            bins_x = self._bins_wave
            bins_x_midp = self._bins_wave_midp
            ax_tf.set_xlabel(r'Wavelength $\lambda$ ($\AA$)')
            cb_units = r"Log $\psi$($\lambda, \tau$)"

        tf = ax_tf.pcolor(bins_x, self._bins_delay, np.log10(self._flux),
                         vmin=cb_min, cmap="afmhot_r")
        ax_tf.set_ylabel(r'Delay $\tau$ (days)')
        ax_tf.set_ylim(top=self._bins_delay[-1])
        ax_tf.set_aspect('auto')

        spec = ax_spec.plot(bins_x_midp, self._flux_spec)
        resp = ax_resp.plot(self._flux_resp, self._bins_delay_midp)

        plt.setp(ax_spec.get_xticklabels(), visible=False)
        plt.setp(ax_spec.get_yticklabels(), visible=False)
        plt.setp(ax_resp.get_xticklabels(), visible=False)
        plt.setp(ax_resp.get_yticklabels(), visible=False)
        ax_none.axis('off')
        cbar = plt.colorbar(tf, orientation="vertical")
        cbar.set_label(cb_units)
        plt.savefig("{}.eps".format(s_file))

#s_root = "ngc5548_1e00_obs_2"
s_user = "root"
s_password = "password"

s_file = "test_large"

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

Session = sqlalchemy.orm.sessionmaker(bind=db_engine)
dbc = Session()

start = time.clock()

try:
    dbc.query(Photon.Weight).first()
    # If so, we go with what we've found.
    print("Found existing filled photon database '{}'".format(s_file))
except sqlalchemy.exc.SQLAlchemyError as e:
    print(e)
    # If not, we populate from the delay dump file. This bit is legacy!
    print("No existing filled photon database, reading from file '{}.delay_dump'".format(s_file))
    Base.metadata.create_all(db_engine)

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
    dbc.commit()


print("Input took {:.2f} seconds".format(time.clock()-start))

query = dbc.query(Photon.Wavelength, Photon.Delay, Photon.Weight)
query = query.filter(and_(Photon.Spectrum == 2))
query = query.filter(and_(Photon.Wavelength > 1475, Photon.Wavelength < 1625))
query = query.filter(and_(Photon.Delay < (3 * seconds_per_day)))
TransferFunction(query, [100,100]).line(1550).as_velocity().run().plot()
#query = query.filter(and_(or_(Photon.ResonantScatters > 0, Photon.Origin == 3, Photon.Origin_matom == True)))
#query = query.filter(and_(Photon.Resonance.in_([380, 381, 382, 416])))

