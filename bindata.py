# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np


class Photon(object):
    def __init__(self, vals):
        print("INIT PHOTON: ",vals)
        self.freq = vals[1]       
        self.position = (vals[3], vals[4], vals[5])
        self.time = vals[6]
        self.weight = vals[2]
        self.origin = vals[8]
        self.observer = vals[7]
        self.last_scattered_line = vals[9]
        self.scatters_resonant = vals[4]
        self.scatters_continuum = vals[5] - vals[4]

class Input(object):
    def InputPhotons(filename):    
        photon_list = list()
        file = open(filename, 'r')
        for line in file:
            print(line)
            if(line.startswith('#')): 
                continue
            vals = line.split()
            print("READ PHOTON: ",vals)
            photon_list.append(Photon(np.asarray(vals,float)))  
        return(photon_list)

    Photons = staticmethod(InputPhotons)


class TransferFunction(object):    
    def __init__(self,dimensions, range_freq, range_time, observer, 
                 last_scattered_lines=None, photon_origins=None, 
                 scatters_resonant=None, scatters_continuum=None):
        self._bins_freq = np.linspace(range_freq[0], range_freq[1],
                                      dimensions[0]+1, endpoint=True, dtype=np.float64)
        self._bins_time = np.linspace(range_time[0], range_time[1],
                                 dimensions[1]+1, endpoint=True, dtype=np.float64)
        self._bins_flux = np.zeroes(shape=dimensions, dtype=np.float64)
        self._bins_count = np.zeroes(shape=dimensions, type=np.intc)
        self._last_scattered_lines = last_scattered_lines
        self._scatters_resonant = scatters_resonant
        self._scatters_continuum = scatters_continuum
        self._photon_origins = photon_origins
        self._observer = observer        
        
    def addPhoton(self, freq, time, weight, origin, observer,
                  last_scattered_line, scatters_resonant, scatters_continuum):
        if self._observer != observer:
            pass
        elif not self._bins_time[0] < time < self._bins_time[-1]: 
            pass
        elif not self._bins_freq[0] < freq < self._bins_freq[-1]:
            pass
        elif self._photon_origins != None and \
            not origin in self._photon_origins:
            pass
        elif self._scatters_resonant != None and \
            not scatters_resonant in self._scatters_resonant:
            pass
        elif self._scatters_continuum != None and \
            not scatters_continuum in self._scatters_continuum:
            pass
        elif self._last_scattered_lines != None and \
            not last_scattered_line in self._last_scattered_lines:
            pass
        else:
            i_freq = np.searchsorted(self._bins_freq, freq)
            i_time = np.searchsorted(self._bins_time, time)
            self._bins_flux[i_freq, i_time] = self._bins_flux[i_freq, i_time] + weight
            self._bins_count[i_freq, i_time] = self._bins_count[i_freq, i_time] + 1
                
    def flux(self, freq, time):
        return(self._bins_flux[np.searchsorted(self._bins_freq, freq),
                               np.searchsorted(self._bins_time, time)])
    def count(self, freq, time):
        return(self._bins_count[np.searchsorted(self._bins_freq, freq),
                                np.searchsorted(self._bins_time, time)])


Input.Photons("test.delay_dump")            