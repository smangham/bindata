# -*- coding: utf-8 -*-
import tfpy
import astropy as ap
import astropy.constants as apc
from astropy import units as u
from astropy.units import cds as ucds

import numpy as np
import numpy.random as npr

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as ani

import sys

def interpolation_across_range(x,y,x_int):
    """
    Simple linear interpolation function

    Args:
        x (astropy.Quantity):    X values
        y (astropy.Quantity):    Y values
        x_int (float):           X to find Y for
    Returns:
        float:  Linear interpolation of Y for x_int
    """

    if x_int >= x.quantity[-1]:
        return y.quantity[-1]
    elif x_int <= x.quantity[0] == 0:
        return y.quantity[0]
    else:
        x_max_index = np.searchsorted(x, x_int)
        x_min_index = x_max_index - 1
        x_int_fraction = (x_int.value - x[x_min_index]) / (x[x_max_index] - x[x_min_index])
        y_diff = y.quantity[x_max_index] - y.quantity[x_min_index]
        return y.quantity[x_min_index] + (y_diff * x_int_fraction)

def import_spectra_times(spectra_times_file, time_units=None):
    """
    Imports a ASCII list of spectrum times

    Args:
        spectra_times_file (string):    Name of the file with times
        time_units (astropy.Unit):      Unit the times are in (e.g. u.s, u.d)

    Returns:
        astropy.QTable:                 Single-column table of time in given units
    """
    spectra_times = ap.table.QTable.read(spectra_times_file, format='ascii', names=['time'])
    spectra_times['time'].unit = time_units
    return spectra_times

def import_lightcurve(lightcurve_file, time_units=None, value_units=None,
                    bolometric_correction=None, error_ratio=None, delta_continuum_range=None,
                    target_bolometric_luminosity=None):
    """
    Inputs a two- or three-column ascii file of times, continuum values and errors.

    Args:
        lightcurve_file (string):   Name of the file to read
        time_units (astropy.Unit):  Unit the times are in (e.g. u.s, u.d)
        value_units (astropy.Unit): Unit the values are in (e.g. )
        bolometric_correction(astropy.Quantity):
                                    Conversion factor from e.g. monochromatic flux to bolometric
        target_bolometric_luminosity (astropy.Quantity):
                                    Target mean bolometric luminosity to rescale to.
        error_ratio (float):        F/F_error ratio to create errors at

    Returns:
        astropy.Table:              Two/three column
    """
    assert bolometric_correction is None or target_bolometric_luminosity is None,\
        "Rescale either by correction factor or by giving target luminosity, not both!"
    assert bolometric_correction is not None or target_bolometric_luminosity is not None,\
        "The lightcurve must be converted to bolometric luminosity! Either provide a correction"+\
        "factor, or provide a mean bolometric luminosity to rescale to."

    if error_ratio is None:
        # If we've not been given a S/N ratio, the input file should have one
        try:
            lightcurve = ap.table.Table.read(lightcurve_file, format='ascii', names=['time', 'value', 'error'])
            lightcurve['error'].unit = value_units
        except:
            print("The input file does not have errors! Provide S/N ratio via 'error_ratio' argument")
            sys.exit(1)
    else:
        # Otherwise, construct a fake error from the S/N ratio
        lightcurve = ap.table.Table.read(lightcurve_file, format='ascii', names=['time', 'value'])
        lightcurve['error'] = lightcurve['value']/error_ratio

    lightcurve['time'].unit = time_units
    lightcurve['value'].unit = value_units

    if delta_continuum_range is not None:
        # If we're correcting the continuum range
        lc_min = np.amin(lightcurve['value'])
        lc_max = np.amax(lightcurve['value'])
        lc_mean = np.mean(lightcurve['value'])
        lc_dc_range = (lc_max - lc_min) / (lc_mean * 2)
        print("Rescaling ΔC. Current range: {:.2g}, ({:.3g}:{:.3g}:{:.3g} {})".format(
            lc_dc_range, lc_min, lc_mean, lc_max, lightcurve['value'].unit))
        lightcurve['value'] -= lc_mean
        lightcurve['value'] *= delta_continuum_range / lc_dc_range
        lightcurve['error'] *= delta_continuum_range / lc_dc_range
        lightcurve['value'] += lc_mean
        lc_dc_range = (np.amax(lightcurve['value']) - np.amin(lightcurve['value'])) / (np.mean(lightcurve['value']) * 2)
        print("Rescaled ΔC. New: {:.2g}, ({:.3g}:{:.3g}:{:.3g} {})".format(lc_dc_range,
            np.amin(lightcurve['value']), np.mean(lightcurve['value']), np.amax(lightcurve['value']), lightcurve['value'].unit))

    if bolometric_correction:
        # If we're correcting e.g. from monochromatic to bolometric
        lightcurve['value'] *= bolometric_correction.value
        lightcurve['error'] *= bolometric_correction.value
        lightcurve['value'].unit *= bolometric_correction.unit
        lightcurve['error'].unit *= bolometric_correction.unit
    elif target_bolometric_luminosity:
        # If we're rescaling this to a given bolometric luminosity
        rescale_factor = target_bolometric_luminosity.value / np.mean(lightcurve['value'])
        lightcurve['value'] *= rescale_factor
        lightcurve['error'] *= rescale_factor
        lightcurve['value'].unit = target_bolometric_luminosity.unit
        lightcurve['error'].unit = target_bolometric_luminosity.unit

    # Calculate the bounds of the lightcurve for use later
    lightcurve.meta['min'] = np.amin(lightcurve['value']) * lightcurve['value'].unit
    lightcurve.meta['mean'] = np.mean(lightcurve['value']) * lightcurve['value'].unit
    lightcurve.meta['max'] = np.amax(lightcurve['value']) * lightcurve['value'].unit

    return lightcurve

def import_spectrum(spectrum_file, bins, values, frequency=True, limits=None,
                    wave_units=None, value_units=None, error_ratio=None, errors=None,
                    subtract_continuum_with_mask=None):
    """
    Imports a spectrum, and converts to target units


    Returns:
        astropy.Table:  Table of input file. Key columns are 'wave', 'value' and 'error'

    """

    # Import a spectrum with bins
    spectrum = ap.table.Table.read(spectrum_file, format='ascii')
    assert type(wave_units) is u.Quantity or type(wave_units) is u.Unit, \
        "Please provide the units that wavelength is to be taken or produced in!"

    # Rename value column and assign units
    spectrum.rename_column(values, 'value')
    spectrum['value'].unit = value_units
    if errors:
        # If there are errors, set the value
        spectrum.rename_column(errors, 'error')
    else:
        # If there's no ratio, there's zero
        spectrum['error'] = 0
    spectrum['error'].unit = value_units

    # Sort out which way around the array goes
    if frequency:
        # If we're using frequency, not wavelength
        spectrum.rename_column(bins, 'freq')
        if spectrum['freq'][0] < spectrum['freq'][-1]:
            # We want to go from high frequency to low frequency
            spectrum.reverse()
        if limits:
            # If we're removing data outside of certain limits
            to_remove = []
            for i in range(0, len(spectrum)):
                # For each line
                if spectrum['freq'][i] > limits[1].value or spectrum['freq'][i] < limits[0].value:
                    # If the frequency is outside of our range, remove it
                    to_remove.append(i)
            spectrum.remove_rows(to_remove)

    else:
        spectrum.rename_column(bins, 'wave')
        # We want to go from low wavelength to high wavelength
        if spectrum['wave'][0] > spectrum['wave'][-1]:
            spectrum.reverse()
        if limits:
            # If we're removing data outside of certain limits
            to_remove = []
            for i in range(0, len(spectrum)):
                # For each line
                if spectrum['wave'][i] > limits[1].value or spectrum['wave'][i] < limits[0].value:
                    # If the wavelength is outside of our range, remove it
                    to_remove.append(i)
            spectrum.remove_rows(to_remove)

    # If we're taking in frequency and converting to wavelength
    if frequency is True:
        # Rename to internal names and calculate minima and maxima (a la voronoi?)
        bin_min_data = np.zeros(len(spectrum))
        bin_max_data = np.zeros(len(spectrum))
        for i in range(0, len(spectrum['freq'])-1):
            bin_max_data[i] = (spectrum['freq'][i+1] + spectrum['freq'][i]) / 2
        for i in range(1, len(spectrum['freq'])):
            bin_min_data[i] = (spectrum['freq'][i-1] + spectrum['freq'][i]) / 2

        # Assume end bins either side are symmetrical about the midpoint
        # This is VERY NOT TRUE for log spectra; TODO add log mode that assumes even spacing in logspace
        bin_min_data[0] = spectrum['freq'][0] - (spectrum['freq'][1]-spectrum['freq'][0])
        bin_max_data[-1] = spectrum['freq'][-1] + (spectrum['freq'][-1]-spectrum['freq'][-2])

        # Assign bin bound arrays
        spectrum["freq_min"] = bin_min_data
        spectrum["freq_max"] = bin_max_data

        # Add units to everything
        # Calculate wavelength bins from the frequencies we've been given
        spectrum['freq'].unit = 1/u.s
        spectrum['freq_min'].unit = 1/u.s
        spectrum['freq_max'].unit = 1/u.s
        spectrum['wave'] = (apc.c / spectrum['freq'].quantity).to(wave_units)
        spectrum['wave_max'] = (apc.c / spectrum['freq_min'].quantity).to(wave_units)
        spectrum['wave_min'] = (apc.c / spectrum['freq_max'].quantity).to(wave_units)

    else:
        bin_min_data = np.zeros(len(spectrum))
        bin_max_data = np.zeros(len(spectrum))
        for i in range(0, len(spectrum)-1):
            bin_max_data[i] = (spectrum['wave'][i+1] + spectrum['wave'][i]) / 2
        bin_max_data[-1] = spectrum['wave'][-1] + (spectrum['wave'][-1]-spectrum['wave'][-2])

        for i in range(1, len(spectrum)):
            bin_min_data[i] = (spectrum['wave'][i-1] + spectrum['wave'][i]) / 2
        bin_min_data[0] = spectrum['wave'][0] - (spectrum['wave'][1]-spectrum['wave'][0])
        spectrum["wave_min"] = bin_min_data
        spectrum["wave_max"] = bin_max_data


        spectrum['wave'].unit = wave_units
        spectrum['wave_min'].unit = wave_units
        spectrum['wave_max'].unit = wave_units
        spectrum['freq'] = (apc.c / spectrum['wave'].quantity).to(1/u.s)
        spectrum['freq_max'] = (apc.c / spectrum['wave_min'].quantity).to(1/u.s)
        spectrum['freq_min'] = (apc.c / spectrum['wave_max'].quantity).to(1/u.s)

    if subtract_continuum_with_mask is not None:
        bins = np.array(spectrum['wave'])
        values = np.array(spectrum['value'])
        masked_bins = np.ma.masked_inside(bins, subtract_continuum_with_mask[0].value,
                                                subtract_continuum_with_mask[1].value)
        masked_values = np.ma.array(values, mask=np.ma.getmaskarray(masked_bins), copy=True)
        continuum_fit = np.poly1d(np.ma.polyfit(masked_bins, masked_values, 1))
        spectrum['value'] -= np.around(continuum_fit(spectrum['wave']),2)

        fig, ax = plt.subplots(1, 1)
        ax2 = ax.twinx()
        l_unmod = ax.plot(spectrum['wave'], values, label="Original", c='k')
        l_masked = ax.plot(spectrum['wave'], masked_values, label="Masked original", c='g')
        l_fit = ax.plot(spectrum['wave'], continuum_fit(spectrum['wave']), label="Fit to mask", c='b')
        l_fit_step = ax.plot(spectrum['wave'], np.around(continuum_fit(spectrum['wave']), 2), label="Fit (stepped)", c='b')
        l_mod = ax2.plot(spectrum['wave'], spectrum['value'], label="Subtracted", c='r')
        ax.set_xlabel("Wavelength (Å)")
        ax.set_ylabel("Flux (non-subtracted)")
        ax2.set_ylabel("Flux (subtracted)")

        lns = l_unmod+l_masked+l_fit+l_fit_step+l_mod
        labs = [l.get_label() for l in lns]
        ax.legend(lns, labs)
        plt.show()
        #sys.exit(1)

    if error_ratio:
        # If we're creating errors based on a passed error ratio
        spectrum['error'] = spectrum['value']/error_ratio
        for i in range(0,len(spectrum)):
            if spectrum['error'][i] < 0:
                spectrum['error'][i] = np.abs(spectrum['error'][i])

    return spectrum

def output_CARAMEL(lightcurve, spectra, spectra_times, spectrum, suffix):
    """
    Given a lightcurve, series of spectra and time outputs to CARAMEL format

    Args:
        lightcurve (Table):     Continuum values and times, in seconds and real units
        spectra (Table):        Table of wavelengths and spectra
        spectra_times (Table):  Table of spectrum times
        spectrum (Table):       Table of the start spectra (used for errors)
        suffix (String):        Suffix appended to filename

    Returns:
        None

    Outputs:
        caramel_lightcurve.txt
        caramel_spectra.txt
        caramel_spectra_times.txt
    """

    # Lightcurve file format:
    # Time (s)
    # Value (rescaled from 1-100)
    # Error (rescaled as above)

    time = np.array(lightcurve['time'].quantity.to(u.s))
    # Rescale value and error from 1-100
    value = np.array((lightcurve['value'] - np.amin(lightcurve['value']))
                    / (np.amax(lightcurve['value']) - np.amin(lightcurve['value'])))
    value = (value * 9) + 1
    error = np.array(lightcurve['error']
                    / (np.amax(lightcurve['value']) - np.amin(lightcurve['value'])))
    error = (error * 9)

    np.savetxt('caramel_lightcurve_{}.txt'.format(suffix), np.column_stack((time, value, error)))

    # Spectrum times file format:
    # Time (s)
    # Dummy value (0)
    # Dummy value (0)
    times = np.array(spectra_times['time'].to(u.s))
    dummy = [0] * len(times)

    np.savetxt('caramel_spectra_times_{}.txt'.format(suffix), np.column_stack((times, dummy, dummy)))

    # Spectra file format:
    # First row is # INT, where INT is number of wavelength bins
    # Second row is central pixel wavelength in angstroms
    # Third row is rescaled flux from 1-100 for spectrum 1
    # Fourth row is rescaled flux error for spectrum 1
    to_save = [spectra['wave']]
    spec_min = 999e99
    spec_max = -999e99
    for column in spectra.colnames[1:]:
        if np.amax(spectra[column]) > spec_max:
            spec_max = np.amax(spectra[column])
        if np.amin(spectra[column]) < spec_min:
            spec_min = np.amin(spectra[column])

    for column in spectra.colnames[2:]:
        value = (np.array(spectra[column]) - spec_min) / (spec_max - spec_min)
        value = (value * 9) + 1
        error = np.array(spectrum['error']) / (spec_max - spec_min)
        error = (error * 9)
        to_save.append(value)
        to_save.append(error)

    np.savetxt('caramel_spectra_{}.txt'.format(suffix), to_save, header='{}'.format(len(spectra)))
    return

def output_MEMECHO(lightcurve, spectra, spectra_times, spectrum, suffix):
    """
    Given a lightcurve, series of spectra and time outputs to MEMECHO format

    Args:
        lightcurve (Table):     Continuum values and times, in seconds and real units
        spectra (Table):        Table of wavelengths and spectra
        spectra_times (Table):  Table of spectrum times
        spectrum (Table):       Table of the start spectra (used for errors)
        suffix (String):        Suffix appended to file name

    Returns:
        None

    Outputs:
        prepspec_*.txt
        prepspec_times.txt
        memecho_lightcurve.txt
    """

    names = []
    for iter, column in enumerate(spectra.colnames[2:]):
        np.savetxt('prepspec_{}_{}.txt'.format(suffix, iter), np.column_stack((spectra['wave'], spectra[column], spectrum['error'])))
        names.append('prepspec_{}_{}.txt'.format(suffix, iter))

    with open('prepspec_times_{}.txt'.format(suffix), 'w') as file:
        for name, time in zip(names, spectra_times['time']):
            file.write("{} {}".format(name, time))

    with open('memecho_lightcurve_{}.txt'.format(suffix), 'w') as file:
        for time, value, error in zip(lightcurve['time'], lightcurve['value'], lightcurve['error']):
            file.write("{} {} {}".format(time, value, error))
    return

# ========== SETTINGS ==========
# ---- LIGHTCURVE SETTINGS -----
lightcurve_file = "light_1158.dat"
lightcurve_time_units = ucds.MJD
lightcurve_value_units = 1e-15 * u.erg / (u.s * u.angstrom * u.cm * u.cm)
#lightcurve_bolometric_correction = 9.0 * 5100.0 * u.angstrom
# We want a curve of the observed bolometric luminosity, so we need to rescale to 100 Pc
lightcurve_target_lum =  1.043e44 * u.erg / (u.s * np.pi * np.power(u.parsec.to(u.cm)*100, 2) * u.cm * u.cm)
# ----- SPECTRUM SETTINGS ------
spectrum_file = "qso_100.spec"
spectrum_bins_name = "Lambda"
spectrum_value_name = "A30P0.50"
# python spectra are per cm2 at 100 pc -> we need to divide ΔC by this value too
spectrum_value_units = u.erg / (u.s * u.angstrom * (np.pi * np.power(u.parsec.to(u.cm)*100, 2) * u.cm * u.cm))
#spectrum_value_to_lightcurve_value = (np.pi * np.power(u.parsec.to(u.cm) * 100, 2)) * u.cm * u.cm / 1000
#spectrum_value_units = 1e-14 * u.erg / (u.s * u.angstrom * u.cm * u.cm)
spectrum_wave_units = u.angstrom
spectrum_flux_error_ratio = 50
spectrum_wave_range = [4861 - 250, 4861 + 250] * u.angstrom
spectrum_continuum_subtract_range = [4861 - 175, 4861 + 175] * u.angstrom
# ------ SPECTRA TIMES ---------
spectra_file = "spectra_qso.dat"
spectra_times_file = "spectra_times.dat"
spectra_times_units = ucds.MJD
spectra_fudge_factor = 1
# ---- ANIMATION SETTINGS ------
animation_file = "timeseries_qso.mp4"
is_reversed = False
# --------- RESCALING ----------
delay_max = 30 * u.d
delta_continuum_range = 0.1
# ----------- OUTPUT -----------
suffix = "qso"

# -------- TF SETTINGS ---------
delay_bins = 100
wave = 4861
# line = 31
# scale_mid = 1/20 # SEY
# scale_min = 1/40
# scale_max = 1/40
line = 49
scale_mid = 1/40 # QSO
scale_min = 1/40
scale_max = 1/40
lim_test = 99999999999

# ===============
# Program begins!
# ===============

print("=== tssproduce started! ===")

# Import all data files
print("Importing lightcurve file '{}'...".format(lightcurve_file))
lightcurve = import_lightcurve(lightcurve_file,
                            time_units=lightcurve_time_units,
                            value_units=lightcurve_value_units,
                            target_bolometric_luminosity=lightcurve_target_lum,
                            delta_continuum_range=0.1)

print("Importing spectrum file '{}'...".format(spectrum_file))
spectrum = import_spectrum(spectrum_file, spectrum_bins_name, spectrum_value_name,
                            frequency=False,
                            wave_units=spectrum_wave_units,
                            value_units=spectrum_value_units,
                            error_ratio=spectrum_flux_error_ratio,
                            limits=spectrum_wave_range,
                            subtract_continuum_with_mask=spectrum_continuum_subtract_range)

print("Importing spectra timing file '{}'...".format(spectra_times_file))
spectra_times = import_spectra_times(spectra_times_file,
                            time_units=spectra_times_units)

# Take spectrum and produce full list of wavelength bins for it to feed to the TF
spectrum_bounds = list(spectrum["wave_min"])
spectrum_bounds.append(spectrum["wave_max"][-1])
spectrum_bounds = np.array(spectrum_bounds)

# Produce a TF
print("Generating Ψ...")
db_100 = tfpy.open_database("/Users/amsys/paper_tss/qso_100", "root", "password")
db_090 = tfpy.open_database("/Users/amsys/paper_tss/qso_090", "root", "password")
db_110 = tfpy.open_database("/Users/amsys/paper_tss/qso_110", "root", "password")
tf_mid = tfpy.TransferFunction(db_100, "qso100", continuum=1.043e46,
                                wave_bins=(len(spectrum_bounds)-1), delay_bins=delay_bins)

tf_mid.line(line, wave).wavelength_bins(spectrum_bounds).delay_dynamic_range(2).run(
            scaling_factor=scale_mid, limit=lim_test,verbose=True).plot()
tf_min = tfpy.TransferFunction(db_090, "qso090", continuum=1.043e46*0.9, template=tf_mid).run(
            scaling_factor=scale_min, limit=lim_test).plot()
tf_max = tfpy.TransferFunction(db_110, "qso110", continuum=1.043e46*1.1, template=tf_mid).run(
            scaling_factor=scale_max, limit=lim_test).plot()


tf_mid.response_map_by_tf(tf_min, tf_max).plot(response_map=True, name='resp')
print("Total Ψ = {}".format(np.sum(tf_mid._response)))
print("Total L = {}".format(np.sum(spectrum['value'])))


# ========
# Convolve
# ========

# Set up the spectra array and set to baseline of start spectrum
spectra = ap.table.Table([spectrum['wave'], spectrum['value']])
# For each point in the lightcurve file
for time in spectra_times['time']:
    # Add a base spectrum to the output spectra file
    spectra["value {}".format(time)] = spectra['value']
    spectra["value {}".format(time)].meta['time'] = time #* lightcurve['time'].unit
    #spectra["value {}".format(time)] = 0.0

# We need to evaluate at every bin-width's step
delay_bins = (tf_mid.delay_bins() * u.s).to(lightcurve['time'].unit)
if delay_max:
    # If we need to rescale this to a given maximum
    delay_bins *= (delay_max / delay_bins[-1])
bin_width = (delay_bins[1] - delay_bins[0]).value
times = np.arange(lightcurve['time'][0], lightcurve['time'][-1] + bin_width, bin_width) * lightcurve['time'].unit

# We need continuum in terms of delta from the mean
delta_continuum = np.zeros(len(times)) * lightcurve['value'].unit
for step, time in enumerate(times):
    # calculate the delta continuum from the 'current value and the starting value, hence the pulse contribution to later timesteps
    delta_continuum[step] = interpolation_across_range(lightcurve['time'], lightcurve['value'], time) - lightcurve.meta['mean']

# -------------
# Begin process
# -------------
print("Beginning {} time steps to generate {} spectra...".format(len(times), len(spectra.columns)-2))
# For each timestep, we send out a 'pulse' of continuum
for step, time in enumerate(times):
    # For each time bin this pulse is smeared out over
    print("Step {}: ΔC = {}".format(step, delta_continuum[step]))
    for i in range(0, len(delay_bins)-1):
        # Figure out what the bounds are for the region this pulse affects
        time_min = time + delay_bins[i]
        time_max = time + delay_bins[i+1]
        for column in spectra.colnames[2:]:
            # For each spectrum, if it occurs at a timestamp within this bin
            if time_min < spectra[column].meta['time'] < time_max:
                    # Add this pulse's contribution to it
                    #print("Current spectra:",np.sum(spectra[column])*spectra[column].unit)
                    #print("Continuum, response:",delta_continuum[step], np.sum(tf_mid.response(delay_index=i)))
                    #print("Total:",np.sum((delta_continuum[step] * tf_mid.response(delay_index=i)/spectrum_wave_units).to(spectrum_value_units)))

                    spectra[column] += ((delta_continuum[step] * tf_mid.response(delay_index=i)
                                     * spectra_fudge_factor / spectrum_wave_units)).value

# Now we have the final spectra, add the experimental errors
for column in spectra.colnames[2:]:
    # For each spectrum
    if 'value' in column:
        for i in range(0, len(spectrum_bounds)-1):
            # For each wavelength bin in each spectrum, add a random error
            spectra[column][i] += npr.normal(scale=spectrum['error'][i])


# Write the spectra out to an easily plotted file
ap.io.ascii.write(spectra, spectra_file, overwrite=True)
ap.io.ascii.write(lightcurve, "lightcurve_out.dat", overwrite=True)

output_CARAMEL(lightcurve, spectra, spectra_times, spectrum, suffix)
output_MEMECHO(lightcurve, spectra, spectra_times, spectrum, suffix)

# Now plot the animation
print("Generating animation...")
figure, (axis, ax_dc) = plt.subplots(2, 1, gridspec_kw={'height_ratios':[3,1]})
figure.subplots_adjust(hspace=.3, wspace=0)

# Find the maxima and minima over all the spectra to scale the plot appropriately
ymin = +np.inf
ymax = -np.inf
for column in spectra.colnames[1:]:
    if np.amax(spectra[column]) > ymax:
        ymax = np.amax(spectra[column])
    if np.amin(spectra[column]) < ymin:
        ymin = np.amin(spectra[column])

# For the upper (spectra) axis, label and set limits
axis.set_xlabel("λ ({})".format(spectrum_wave_units))
axis.set_ylabel("L ({})".format(spectrum_value_units))
axis.set_ylim([ymin,ymax])
axis.set_xlim([spectra['wave'][0], spectra['wave'][-1]])

# For the lower (continuum) axis, label and set limits
y_dc_min = np.amin(delta_continuum)
y_dc_max = np.amax(delta_continuum)
ax_dc.set_xlabel("T ({})".format(lightcurve_time_units))
ax_dc.set_ylabel("ΔC ({})".format(lightcurve_value_units))
ax_dc.set_xlim([lightcurve['time'][0], lightcurve['time'][-1]])
ax_dc.set_ylim([y_dc_min.value, y_dc_max.value])
# Put a line for the mean continuum on time-continuum
line_mean = ax_dc.plot((times[0].value, times[-1].value), (0.0, 0.0), '-', c='k', zorder=-len(times))
# Set up a second axis to show % difference
ax_dc_percent = ax_dc.twinx()
ax_dc_percent.set_ylim([y_dc_min.value * 100.0 / lightcurve.meta['mean'].value,
                        y_dc_max.value * 100.0 / lightcurve.meta['mean'].value])
ax_dc_percent.set_ylabel("ΔC (%)")
ax_dc_percent.get_yaxis().get_major_formatter().set_useOffset(False)

# Plot the 'original' spectrum in black
line_start = axis.plot(spectra['wave'], spectra['value'], '-', c='k', zorder=0)
line_colours = [ plt.cm.jet(x) for x in np.linspace(0.0, 1.0, len(spectra_times)) ]

# Artists is the list of objects redrawn each animation (so all plots)
artists = [line_start]
# We evaluate the spectra one-by-one for timestamps, starting with 1st (i.e. spectra column 2)
next_column = 2

def zorder(step, steps, reverse=False):
    # Straightforward- do we want new points over old, or old over new?
    # Latter is useful for showing how a step change propagates
    if reverse:
        return -step
    else:
        return step-steps

def update_figure(step):
    global next_column
    #print("Step {}".format(step))
    time_colour = 'grey'
    if spectra_times['time'][0] <= times[step] <= spectra_times['time'][-1]:
        time_colour = plt.cm.jet( (times[step]-spectra_times['time'][0]) /
                                (spectra_times['time'][-1] - spectra_times['time'][0]) )

    if times[step] > spectra.columns[next_column].meta['time']:
        line_new = axis.plot(spectra['wave'], spectra.columns[next_column], '-', markersize=0.5,
                             c=time_colour, zorder=zorder(step, len(times), is_reversed))
        times_for_line = np.array([times[step].value, times[step].value])
        values_for_line = np.array([y_dc_min.value, y_dc_max.value])
        line_vertical = ax_dc.plot(times_for_line, values_for_line, ':', linewidth=1, c=time_colour, zorder=0)
        artists.append(line_vertical)
        artists.append(line_new)
        next_column += 1
    point_new = ax_dc.plot(times[step], delta_continuum[step], 'o', c=time_colour, zorder=zorder(step, len(times), is_reversed))
    artists.append(point_new)
    return artists,

# Generate animation and save to file
animation = ani.FuncAnimation(figure, update_figure, frames=len(times))
animation.save(animation_file, fps=24)
