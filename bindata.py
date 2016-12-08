# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#import sqlite3 as sql
import numpy as np
import astropy as ap
import time
import sys
import matplotlib.pyplot as plt
import MySQLdb as sql

class Query(object):
    def __init__(self, dbc, spectrum):
        if spectrum < 0:
            print("Error: {} is not a valid observer".format(spectrum))
            self.query = "ERROR"
        else:
            self.dbc = dbc
            # PRAGMA returns entries like 'key|name|format|flag1|flag2|flag3'
            self.column_names = [ entry[1] for entry in  dbc.execute("PRAGMA table_info(Photons)").fetchall()]
            self.query = "SELECT Wavelength, Delay, Weight FROM Photons WHERE Spectrum == {}".format(spectrum)
            self._equal = dict()
            self._not_equal = dict()
            self._in_range = dict()
            self._not_in_range = dict()
            self._in_list = dict()
            self._not_in_list = dict()            
            self._equal["Spectrum"] = spectrum

    def clear(self):
            self.query = "SELECT Wavelength, Delay, Weight FROM Photons WHERE Spectrum == {}".format(self.spectrum)
            spectrum = self._equal["Spectrum"]            
            self._equal.clear()
            self._not_equal.clear()
            self._in_range.clear()
            self._not_in_range.clear()
            self._in_list.clear()
            self._not_in_list.clear()
            self._equal["Spectrum"] = spectrum            
            return self

    def equal(self, column, value):
        if column not in self.column_names:
            print("Error: Column '{}' does not exist in file.".format(column))
        else:   
            self._equal[column] = value
            self.query = self.query + " AND {} == {}".format(column, value)
            return self
    def not_equal(self, column, value):
        if column not in self.column_names:
            print("Error: Column '{}' does not exist in file.".format(column)) 
        else:   
            self._not_equal[column] = value
            self.query = self.query + " AND {} != {}".format(column, value)
            return self

    def in_range(self, column, value_min, value_max):
        if column not in self.column_names:
            print("Error: Column '{}' does not exist in file.".format(column) )
        else:  
            self._in_range[column] = [value_min, value_max]
            self.query = self.query + " AND {} BETWEEN {} AND {}".format(column, value_min, value_max)
            return self
    def not_in_range(self, column, value_min, value_max):
        if column not in self.acolumn_names:
            print("Error: Column '{}' does not exist in file.".format(column))
        else:  
            self._not_in_range[column] = [value_min, value_max]
            self.query = self.query + " AND {} NOT BETWEEN {} AND {}".format(column, value_min, value_max)
            return self

    def in_list(self, column, values):
        if column not in self.column_names:
            print("Error: Column '{}' does not exist in file.".format(column))
        else:  
            self._in_list[column] = values
            self.query = self.query + " AND {} IN ({})".format(column,", ".join(map(str, values)))
            return self
    def not_in_list(self, column, values):
        if column not in self.column_names:
            print("Error: Column '{}' does not exist in file.".format(column))
        else:  
            self._not_in_list[column] = values
            self.query = self.query + " AND {} NOT IN ({})".format(column,", ".join(map(str, values)))
            return self

    def run(self):
        return dbc.execute(self.query).fetchall()


class TransferFunction(Query):    
    def __init__(self,dbc, dimensions, spectrum):
        super(TransferFunction,self).__init__(dbc, spectrum)
        self._dimensions = dimensions
        self._flux = np.zeros(shape=dimensions, dtype=np.float64)
        self._count = np.zeros(shape=dimensions, dtype=np.intc)
        
    def run(self):
        data = np.asarray(super(TransferFunction,self).run())
        # Data returned as Wavelength, Delay, Weight
        if "Wavelength" in self._in_range:
            range_wave = self._in_range["Wavelength"]
        else:
            print("Error: No wavelength range defined, autodetecting...")
            range_wave = [np.amin(data[:,0]), np.amax(data[:,0])]                
        
        if "Delay" in self._in_range:
            range_delay = self._in_range["Delay"]
        else:
            range_delay = [0,np.percentile(data[:,1],99)]
            print("Error: No delay range defined, autodetecting to {:.0f}:{:.2f}".format(*range_delay))
                        
        self._bins_wave  = np.linspace(range_wave[0], range_wave[1],
                                       self._dimensions[0]+1, endpoint=True, dtype=np.float64)
        self._bins_delay = np.linspace(range_delay[0], range_delay[1],
                                       self._dimensions[1]+1, endpoint=True, dtype=np.float64)

        self._flux, junk, junk = np.histogram2d(data[:,1], data[:,0], weights=data[:,2], 
                       bins=[self._bins_delay, self._bins_wave], normed=True)                    
        self._count, junk, junk = np.histogram2d(data[:,0], data[:,1], 
                       bins=[self._bins_delay, self._bins_wave], normed=True)
        self._flux_min = np.amin(self._flux[np.nonzero(self._flux)])
        self._flux[self._flux == 0] = 1e-35
        return self
        
    def flux(self, wave, delay):
        return(self._flux[np.searchsorted(self._bins_delay, delay),
                          np.searchsorted(self._bins_wave, wave)])
    def count(self, wave, delay):
        return(self._count[np.searchsorted(self._bins_delay, delay),
                           np.searchsorted(self._bins_wave, wave)])

    def plot(self, log_range=5):
        cb_max = np.log10(np.amax(self._flux))
        cb_min = np.log10(np.amax(self._flux))-log_range
        print(cb_min, cb_max)
        
        fig = plt.figure()
        axis = fig.add_subplot(111)
        tf = axis.pcolor(self._bins_wave, self._bins_delay, np.log10(self._flux),
                         vmin=cb_min, cmap="afmhot_r")
        axis.set_ylim(top=self._bins_delay[-1])
        axis.set_aspect('auto')
        cbar = plt.colorbar(tf, orientation="vertical")
        fig.show()

s_root = "ngc5548_1e00_obs_2"
s_file_db = s_root+".db"
s_file_dd = s_root+".delay_dump"

q_table_exists = "SELECT name FROM sqlite_master WHERE type='table' AND name='Photons'"
q_table_add = """CREATE TABLE Photons(Key INT, Wavelength DOUBLE, Weight DOUBLE, X DOUBLE, Y DOUBLE, Z DOUBLE, 
              ContinuumScatters INT, ResonantScatters INT, Delay DOUBLE, Spectrum INT, Origin INT, Resonance INT)"""
q_photon_add = "INSERT INTO Photons VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"


as_names = ["Frequency", "Wavelength", "Weight", "X", "Y", "Z", "Scatters", "Resonant Scatters", "Delay", "Extracted", "Spectrum", "Origin", "Resonance"]

### TRY OPENING THE DATABASE ###
db = None
try:
    db = sql.connect(s_file_db)
except sql.Error as e:
    print("Error {}: ".format(e.args[0]))
    sys.exit(1)

### DOES IT ALREADY EXIST? ###
dbc = db.cursor()


start = time.clock()
if dbc.execute(q_table_exists).fetchone():
    # If so, we go with what we've found.
    print("Found existing filled photon database '{}'".format(s_file_db))
else:
    # If not, we populate from the delay dump file. This bit is legacy!
    print("No existing filled photon database, reading from file '{}'".format(s_file_dd))
    dbc = db.cursor()
    dbc.execute(q_table_add)

    key = 0
    delay_dump = open(s_file_dd, 'r')
    for key, line in enumerate(delay_dump):
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
        values.insert(0, key)
        dbc.execute(q_photon_add, values)
    db.commit()

print("Input took {:.2f} seconds".format(time.clock()-start))

start = time.clock()
Query(dbc, 0).equal("Resonance",416).in_range("Wavelength", 1475, 1625).run()
print("Query took {:.2f} seconds".format(time.clock()-start))

start = time.clock()
tf_a = TransferFunction(dbc, [50,50], 0).in_range("Wavelength", 1450, 1650).run()
print("TF took {:.2f} seconds".format(time.clock()-start))
start = time.clock()
tf_b = TransferFunction(dbc, [50,50], 0).in_range("Wavelength", 1450, 1650).in_list("Resonance",[381,416]).run()
print("TF took {:.2f} seconds".format(time.clock()-start))
tf_a.plot()
tf_b.plot()
db.close()