# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import sqlite3
import numpy as np
import astropy as ap
from astropy.table import Table
import time

class Query(object):
    def __init__(self, dbc, spectrum):
        if spectrum < 0:
            print "Error: {} is not a valid observer".format(spectrum)
            self.query = "ERROR"
        else:
            self.spectrum = spectrum
            self.dbc = dbc
            self.as_column_names = as_get_column_names(dbc)
            self.s_query = "SELECT Wavelength, Delay, Weight FROM Photons WHERE Spectrum == {}".format(self.spectrum)

    def clear(self):
            self.query = "SELECT Wavelength, Delay, Weight FROM Photons WHERE Spectrum == {}".format(self.spectrum)
            return self

    def equal(self, s_column, i_val):
        if s_column not in self.as_column_names:
            print "Error: Column '{}' does not exist in file.".format(s_column) 
        else:   
            self.s_query = self.s_query + " AND {} == {}".format(s_column, i_val)
            return self
    def not_equal(self, s_column, i_val):
        if s_column not in self.as_column_names:
            print "Error: Column '{}' does not exist in file.".format(s_column) 
        else:   
            self.s_query = self.s_query + " AND {} != {}".format(s_column, i_val)
            return self

    def in_range(self, s_column, r_val_min, r_val_max):
        if s_column not in self.as_column_names:
            print "Error: Column '{}' does not exist in file.".format(s_column) 
        else:  
            self.s_query = self.s_query + " AND {} BETWEEN {} AND {}".format(s_column, r_val_min, r_val_max)
            return self
    def not_in_range(self, s_column, r_val_min, r_val_max):
        if s_column not in self.as_column_names:
            print "Error: Column '{}' does not exist in file.".format(s_column) 
        else:  
            self.s_query = self.s_query + " AND {} NOT BETWEEN {} AND {}".format(s_column, r_val_min, r_val_max)
            return self

    def in_list(self, s_column, ar_val):
        if s_column not in self.as_column_names:
            print "Error: Column '{}' does not exist in file.".format(s_column) 
        else:  
            self.s_query = self.s_query + " AND {} IN ({})".format(s_column,", ".join(map(str, ar_val)))
            return self
    def not_in_list(self, s_column, ar_val):
        if s_column not in self.as_column_names:
            print "Error: Column '{}' does not exist in file.".format(s_column) 
        else:  
            self.s_query = self.s_query + " AND {} NOT IN ({})".format(s_column,", ".join(map(str, ar_val)))
            return self

    def run(self):
        return dbc.execute(self.s_query).fetchall()

def as_get_column_names(dbc):
    # PRAGMA returns entries like 'key|name|format|flag1|flag2|flag3'
    return [ entry[1] for entry in  dbc.execute("PRAGMA table_info(Photons)").fetchall()]
    
s_root = "ngc5548_1e00_obs_2"
s_file_db = s_root+".db"
s_file_dd = s_root+".delay_dump"

q_table_exists = "SELECT name FROM sqlite_master WHERE type='table' AND name='Photons'"
q_table_add = """CREATE TABLE Photons(Key INT, Wavelength DOUBLE, Weight DOUBLE, X DOUBLE, Y DOUBLE, Z DOUBLE, 
              ContinuumScatters INT, ResonantScatters INT, Delay DOUBLE, Spectrum INT, Origin INT, Resonance INT)"""
q_photon_add = "INSERT INTO Photons VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"


as_names = ["Frequency", "Wavelength", "Weight", "X", "Y", "Z", "Scatters", "Resonant Scatters", "Delay", "Extracted", "Spectrum", "Origin", "Resonance"]
t_photons = Table.read(s_file_dd, format="ascii", names=as_names)
t_photons.remove_columns(["Extracted", "Frequency"])
t_photons["Scatters"] = t_photons["Scatters"] - t_photons["Resonant Scatters"]
t_photons.rename_column('Scatters', "Continuum Scatters")

### TRY OPENING THE DATABASE ###
db = None
try:
    db = sqlite3.connect(s_file_db)
except sqlite3.Error, e:
    print "Error {}: ".format(e.args[0])
    sys.exit(1)

### DOES IT ALREADY EXIST? ###
dbc = db.cursor()
if dbc.execute(q_table_exists).fetchone():
    # If so, we go with what we've found.
    print "Found existing filled photon database '{}'".format(s_file_db)
else:
    # If not, we populate from the delay dump file. This bit is legacy!
    print "No existing filled photon database, reading from file '{}'".format(s_file_dd)
    dbc = db.cursor()
    dbc.execute("DROP TABLE IF EXISTS Photons")
    dbc.execute(q_table_add)

    i_key = 0
    for line in t_photons:
        dbc.execute(q_photon_add, line.insert(i_key,0))
        key=key+1

    db.commit()

start = time.clock()
Query(dbc, 1).equal("Resonance",416).in_range("Wavelength", 1475, 1625).run()
elapsed = time.clock() - start
print "Seconds ",elapsed,"Minutes",elapsed/60,"Hours",elapsed/360
db.close()



# class Photon(object):
#     def __init__(self, vals):
#         print("INIT PHOTON: ",vals)
#         self.freq = vals[1]       
#         self.position = (vals[3], vals[4], vals[5])
#         self.time = vals[6]
#         self.weight = vals[2]
#         self.origin = vals[8]
#         self.observer = vals[7]
#         self.last_scattered_line = vals[9]
#         self.scatters_resonant = vals[4]
#         self.scatters_continuum = vals[5] - vals[4]

# class Input(object):
#     def InputPhotons(filename):    
#         photon_list = list()
#         file = open(filename, 'r')
#         for line in file:
#             print(line)
#             if(line.startswith('#')): 
#                 continue
#             vals = line.split()
#             print("READ PHOTON: ",vals)
#             photon_list.append(Photon(np.asarray(vals,float)))  
#         return(photon_list)

#     Photons = staticmethod(InputPhotons)


# class TransferFunction(object):    
#     def __init__(self,dimensions, range_freq, range_time, observer, 
#                  last_scattered_lines=None, photon_origins=None, 
#                  scatters_resonant=None, scatters_continuum=None):
#         self._bins_freq = np.linspace(range_freq[0], range_freq[1],
#                                       dimensions[0]+1, endpoint=True, dtype=np.float64)
#         self._bins_time = np.linspace(range_time[0], range_time[1],
#                                  dimensions[1]+1, endpoint=True, dtype=np.float64)
#         self._bins_flux = np.zeroes(shape=dimensions, dtype=np.float64)
#         self._bins_count = np.zeroes(shape=dimensions, type=np.intc)
#         self._last_scattered_lines = last_scattered_lines
#         self._scatters_resonant = scatters_resonant
#         self._scatters_continuum = scatters_continuum
#         self._photon_origins = photon_origins
#         self._observer = observer        
        
#     def addPhoton(self, freq, time, weight, origin, observer,
#                   last_scattered_line, scatters_resonant, scatters_continuum):
#         if self._observer != observer:
#             pass
#         elif not self._bins_time[0] < time < self._bins_time[-1]: 
#             pass
#         elif not self._bins_freq[0] < freq < self._bins_freq[-1]:
#             pass
#         elif self._photon_origins != None and \
#             not origin in self._photon_origins:
#             pass
#         elif self._scatters_resonant != None and \
#             not scatters_resonant in self._scatters_resonant:
#             pass
#         elif self._scatters_continuum != None and \
#             not scatters_continuum in self._scatters_continuum:
#             pass
#         elif self._last_scattered_lines != None and \
#             not last_scattered_line in self._last_scattered_lines:
#             pass
#         else:
#             i_freq = np.searchsorted(self._bins_freq, freq)
#             i_time = np.searchsorted(self._bins_time, time)
#             self._bins_flux[i_freq, i_time] = self._bins_flux[i_freq, i_time] + weight
#             self._bins_count[i_freq, i_time] = self._bins_count[i_freq, i_time] + 1
                
#     def flux(self, freq, time):
#         return(self._bins_flux[np.searchsorted(self._bins_freq, freq),
#                                np.searchsorted(self._bins_time, time)])
#     def count(self, freq, time):
#         return(self._bins_count[np.searchsorted(self._bins_freq, freq),
#                                 np.searchsorted(self._bins_time, time)])


# Input.Photons("test.delay_dump")            