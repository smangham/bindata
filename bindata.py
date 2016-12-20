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
import mysql.connector as sql

class Query(object):
    def __init__(self, dbc, table, answers, limit=-1):
        self._table = table
        self._answers = answers
        self._dbc = dbc
        self._limit = limit
        # PRAGMA returns entries like 'key|name|format|flag1|flag2|flag3'
        self._column_names = [ entry[1] for entry in  dbc.execute("PRAGMA table_info(%s)".format(table)).fetchall()]
        self._equal = dict()
        self._not_equal = dict()
        self._in_range = dict()
        self._not_in_range = dict()
        self._in_list = dict()
        self._not_in_list = dict()            

    def clear(self):
        self._equal.clear()
        self._not_equal.clear()
        self._in_range.clear()
        self._not_in_range.clear()
        self._in_list.clear()
        self._not_in_list.clear()
        return self

    def equal(self, column, value):
        if column not in self.column_names:
            print("Error: Column '{}' does not exist in file.".format(column))
        else:   
            self._equal[column] = value
            return self
    def not_equal(self, column, value):
        if column not in self.column_names:
            print("Error: Column '{}' does not exist in file.".format(column)) 
        else:   
            self._not_equal[column] = value
            return self

    def in_range(self, column, value_min, value_max):
        if column not in self.column_names:
            print("Error: Column '{}' does not exist in file.".format(column) )
        else:  
            self._in_range[column] = [value_min, value_max]
            return self
    def not_in_range(self, column, value_min, value_max):
        if column not in self.acolumn_names:
            print("Error: Column '{}' does not exist in file.".format(column))
        else:  
            self._not_in_range[column] = [value_min, value_max]
            return self

    def in_list(self, column, values):
        if column not in self.column_names:
            print("Error: Column '{}' does not exist in file.".format(column))
        else:  
            self._in_list[column] = values
            return self
    def not_in_list(self, column, values):
        if column not in self.column_names:
            print("Error: Column '{}' does not exist in file.".format(column))
        else:  
            self._not_in_list[column] = values
            self.query = self.query + " AND {} NOT IN ({})".format(column,", ".join(map(str, values)))
            return self

    def run(self):
        query = "SELECT {} FROM {} ".format(", ".join(self._answers), self._table)
        for key, value in self._equal.items():
            query += " AND {} == {}".format(column, value)
        for key, value in self._not_equal.items():
            query += " AND {} != {}".format(column, value)
        for key, value in self._between.items():
            query += " AND {} BETWEEN {} AND {}".format(column, value_min, value_max)
        for key, value in self._not_between.items():
            query += " AND {} NOT BETWEEN {} AND {}".format(column, value_min, value_max)
        for key, value in self._in_list.items():       
            query += " AND {} IN ({})".format(column,", ".join(map(str, values)))
        for key, value in self._in_list.items():       
            query += " AND {} NOT IN ({})".format(column,", ".join(map(str, values)))
        if limit > 0:
            query += " LIMIT {}".format(limit)
        query.replace('AND', 'WITH')
        print(query)
        return _dbc.execute(self.query).fetchall()

class TransferFunction(Query):    
    def __init__(self,dbc, dimensions, spectrum):
        super(TransferFunction,self).__init__(dbc, "Photons")
        self._dimensions = dimensions
        self._flux = np.zeros(shape=dimensions, dtype=np.float64)
        self._count = np.zeros(shape=dimensions, dtype=np.intc)
        self._equal["Spectrum"] = spectrum


    def clear(self):
        super(TransferFunction,self).clear()
        self._equal["Spectrum"] = spectrum

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

    def plot(self, dynamic_range=5):
        cb_max = np.log10(np.amax(self._flux))
        cb_min = cb_max-dynamic_range
        print(cb_min, cb_max)
        
        fig = plt.figure()
        axis = fig.add_subplot(111)
        tf = axis.pcolor(self._bins_wave, self._bins_delay, np.log10(self._flux),
                         vmin=cb_min, cmap="afmhot_r")
        axis.set_ylim(top=self._bins_delay[-1])
        axis.set_aspect('auto')
        cbar = plt.colorbar(tf, orientation="vertical")
        fig.show()

#s_root = "ngc5548_1e00_obs_2"
s_user = "root"
s_password = "password"

s_root = "test"
s_file_db = s_root
s_file_dd = s_root+".delay_dump"

q_table_exists = "SELECT name FROM sqlite_master WHERE type='table' AND name='Photons'"
q_table_add = "CREATE TABLE Photons(ID INT AUTO_INCREMENT, Wavelength DOUBLE, Weight DOUBLE,"\
            " X DOUBLE, Y DOUBLE, Z DOUBLE, ContinuumScatters INT, ResonantScatters INT,"\
            " Delay DOUBLE, Spectrum INT, Origin INT, Resonance INT, PRIMARY KEY (id))"
q_photon_add = "INSERT INTO Photons VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"


as_names = ["Frequency", "Wavelength", "Weight", "X", "Y", "Z", "Scatters", "Resonant Scatters", "Delay", "Extracted", "Spectrum", "Origin", "Resonance"]


### TRY OPENING THE DATABASE ###
print ("Opening database '{}'...".format(s_file_db))

db = None
try:
    db = sql.connect(host="localhost", user=s_user, passwd=s_password, db=s_file_db)
except sql.Error as e:
    print("Error {}: {}".format(e.args[0], e.args[1]))
    sys.exit(1)

### DOES IT ALREADY EXIST? ###
print ("Searching for table 'Photons'...")

dbc = db.cursor()

start = time.clock()
result = dbc.execute("SHOW TABLES LIKE 'Photons'")
print(result)

if dbc.execute("SHOW TABLES LIKE 'Photons'").fetchall():
    # If so, we go with what we've found.
    print("Found existing filled photon database '{}'".format(s_file_db))
else:
    # If not, we populate from the delay dump file. This bit is legacy!
    print("No existing filled photon database, reading from file '{}'".format(s_file_dd))
    dbc.execute(q_table_add)

    key = 0
    delay_dump = open(s_file_dd, 'r')
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