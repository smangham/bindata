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

    value_orig = lightcurve['value']

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

    fig, ax1 = plt.subplots(1)
    ax2 = ax1.twinx()
    ax1.set_title("Continuum Rescaling")
    ax1.set_xlabel("Time (MJD)")
    ax1.set_ylabel("Flux (original)")
    ax2.set_ylabel("Flux (rescaled)")
    l_orig = ax1.plot(lightcurve["time"], value_orig, '-', c='r', label='Original')
    l_resc = ax2.plot(lightcurve["time"], lightcurve["value"], '--', c='b', label='Rescaled')
    lns = l_orig+l_resc
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs)
    plt.show()

    return lightcurve


def import_spectrum(spectrum_file, bins, values, frequency=True, limits=None,
                    wave_units=None, value_units=None,
                    errors=None, ,error_absolute=None, error_relative=None,
                    subtract_continuum_with_mask=None, rebin_to=None):
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
    if error:
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

        spectrum.remove_rows(slice(np.searchsorted(spectrum["wave"], subtract_continuum_with_mask[1]) , len(spectrum)+1))
        spectrum.remove_rows(slice(0, np.searchsorted(spectrum["wave"], subtract_continuum_with_mask[0])))

        fig, ax = plt.subplots(1, 1)
        ax2 = ax.twinx()
        ax.set_title("Continuum Subtraction")
        l_unmod = ax.plot(bins, values, label="Original", c='k')
        l_masked = ax.plot(bins, masked_values, label="Masked original", c='g')
        l_fit = ax.plot(spectrum['wave'], continuum_fit(spectrum['wave']), label="Fit to mask", c='b')
        # No longer necessary now we output to many DP
        # l_fit_step = ax.plot(spectrum['wave'], np.around(continuum_fit(spectrum['wave']), 2), label="Fit (stepped)", c='b')
        l_mod = ax2.plot(spectrum['wave'], spectrum['value'], label="Subtracted", c='r')
        ax.set_xlabel("Wavelength (Å)")
        ax.set_ylabel("Flux (non-subtracted)")
        ax2.set_ylabel("Flux (subtracted)")

        lns = l_unmod+l_masked+l_fit+l_mod
        #lns = l_unmod+l_masked+l_fit+l_fit_step+l_mod
        labs = [l.get_label() for l in lns]
        ax.legend(lns, labs)
        plt.show()

        #sys.exit(1)

    if rebin_to:
        # If we're rebinning to X bins
        wave_bounds = np.linspace(spectrum['wave'][0], spectrum['wave'][-1], rebin_to+1)
        wave_midpoints = np.zeros(rebin_to)
        values = np.zeros(rebin_to)
        errors = np.zeros(rebin_to)
        values_min = values - errors
        values_max = values + errors

        for i in range(0, rebin_to):
            wave_midpoints[i] = (wave_bounds[i] + wave_bounds[i+1]) / 2
            full_bin_width = wave_bounds[i+1] - wave_bounds[i]
            for j in range(0, len(spectrum)):
                if spectrum["wave_min"][j] > wave_bounds[i+1] or spectrum["wave_max"][j] < wave_bounds[i]:
                    continue
                elif spectrum["wave_min"][j] > wave_bounds[i] and spectrum["wave_max"][j] < wave_bounds[i+1]:
                    bin_width = spectrum["wave_max"][j] - spectrum["wave_min"][j]
                elif spectrum["wave_min"][j] < wave_bounds[i]:
                    bin_width = spectrum["wave_max"][j] - wave_bounds[i]
                elif spectrum["wave_max"][j] > wave_bounds[i+1]:
                    bin_width = wave_bounds[i+1] - spectrum["wave_min"][j]

                values[i] += spectrum["value"][j] * bin_width / full_bin_width
                errors[i] += spectrum["error"][j] * bin_width / full_bin_width

            if values[i] < 0:
                values[i] = 0

        freq_bounds = (apc.c / (wave_bounds * wave_units)).to(1/u.s).value
        freq_midpoints = (apc.c / (wave_midpoints * wave_units)).to(1/u.s).value
        freq_min = freq_bounds[1:]
        freq_max = freq_bounds[:-1]

        fig, ax = plt.subplots(1)
        ax.set_title("Rebinning from {} to {} bins".format(len(spectrum), rebin_to))
        ax.set_xlabel("Wavelength (Å)")
        ax.set_ylabel("Flux")
        ax.plot(spectrum["wave"], spectrum["value"], '-', c='r', zorder=1, label="Original")
        ax.errorbar(wave_midpoints, values, errors, c='b', label="Rebinned")
        ax.legend()
        plt.show()

        fig, ax = plt.subplots(1)
        ax.set_title("Rebinning checking consistency of wave/freq bins")
        ax.plot(spectrum["freq"], spectrum["wave"], label="Original")
        ax.plot(freq_midpoints, wave_midpoints, 'x', label="Rebinned")
        ax.legend()
        plt.show()

        fig, ax = plt.subplots(1)
        ax.set_title("Rebinning checking width of freq bins")
        ax.errorbar(freq_midpoints, freq_midpoints, fmt='x',
            xerr=[freq_midpoints-freq_min, freq_max-freq_midpoints],
            yerr=[freq_midpoints-freq_min, freq_max-freq_midpoints],
            label="Rebinned")
        ax.errorbar(spectrum["freq"], spectrum["freq"],
            yerr=[spectrum["freq"]-spectrum["freq_min"], spectrum["freq_max"]-spectrum["freq"]],
            label="Original")
        ax.set_xlabel("Min")
        ax.set_ylabel("Max")
        ax.legend()
        plt.show()

        # Replace the existing spectrum with rebinned values and errors
        spectrum.remove_rows(slice(rebin_to, len(spectrum)+1))
        spectrum["value"] = values
        spectrum["error"] = errors
        spectrum["value"].unit = spectrum_value_units
        spectrum["error"].unit = spectrum_value_units

        spectrum["wave"] = wave_midpoints
        spectrum["wave_min"] = wave_bounds[:-1]
        spectrum["wave_max"] = wave_bounds[1:]
        spectrum["wave"].unit = spectrum_wave_units
        spectrum["wave_min"].unit = spectrum_wave_units
        spectrum["wave_max"].unit = spectrum_wave_units

        spectrum["freq"] = freq_midpoints
        spectrum["freq_min"] = freq_min
        spectrum["freq_max"] = freq_max
        spectrum["freq"].unit = 1/u.s
        spectrum["freq_min"].unit = 1/u.s
        spectrum["freq_max"].unit = 1/u.s


    if error_absolute:
        # If we're adding an absolute error, apply
        spectrum['error'] = error_absolute

    if error_relative:
        # If we're adding a relative error scaled to the peak
        index = np.argmax(spectrum['value'])
        spectrum['error'] = spectrum['value'][index] * error_relative

    return spectrum

def output_CARAMEL(lightcurve, spectra, spectra_times, spectrum, suffix):
    """
    Given a lightcurve, series of spectra and time outputs to CARAMEL format

    Args:
        lightcurve (Table):     Continuum values and times, in seconds and real units
        spectra (Table):        Table of wavelengths and spectra. Continuum-subtracted.
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
        spectra (Table):        Table of wavelengths and spectra. Not continuum-subtracted.
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
        np.savetxt('prepspec_{}_{:03d}.txt'.format(suffix, iter), np.column_stack((spectra['wave'], spectra[column], spectrum['error'])))
        names.append('prepspec_{}_{:03d}.txt'.format(suffix, iter))

    with open('prepspec_times_{}.txt'.format(suffix), 'w') as file:
        for name, time in zip(names, spectra_times['time']):
            file.write("{} {}\n".format(name, time.value))

    with open('memecho_lightcurve_{}.txt'.format(suffix), 'w') as file:
        for time, value, error in zip(lightcurve['time'], lightcurve['value'], lightcurve['error']):
            file.write("{} {} {}\n".format(time, value, error))
    return

# ========== SETTINGS ==========
# -------- TF SETTINGS ---------
tf_lim_test = 999999999
tf_delay_bins = 100
tf_wave = 6562.8
# ---- Seyfert Settings ----
# spectrum_file = "sey_100.spec"
# line = 28
# databases = { 'min':{'path':"/home/swm1n12/bindata/sey_090", 'scale':60, 'continuum': 1.043e44*0.9},
#                   'mid':{'path':"/home/swm1n12/bindata/sey_100", 'scale':60, 'continuum': 1.043e44},
#                   'max':{'path':"/home/swm1n12/bindata/sey_110", 'scale':60, 'continuum': 1.043e44*1.1}  }

# ---- QSO Settings ----
spectrum_file = "qso_100.spec"
tf_line = 44
databases = { 'min':{'path':"/home/swm1n12/bindata/qso_090", 'scale':50,  'continuum': 1.043e46*0.9},
              'mid':{'path':"/home/swm1n12/bindata/qso_100", 'scale':100, 'continuum': 1.043e46},
              'max':{'path':"/home/swm1n12/bindata/qso_110", 'scale':50,  'continuum': 1.043e46*1.1}  }

# ---- LIGHTCURVE SETTINGS -----
lightcurve_file = "light_1158.dat"
lightcurve_time_units = ucds.MJD
lightcurve_value_units = 1e-15 * u.erg / (u.s * u.angstrom * u.cm * u.cm)
#l ightcurve_bolometric_correction = 9.0 * 5100.0 * u.angstrom
# We want a curve of the observed bolometric luminosity, so we need to rescale to 100 Pc
lightcurve_target_lum =  tf_continuum * u.erg / (u.s * np.pi * np.power(u.parsec.to(u.cm)*100, 2) * u.cm * u.cm)

# ----- SPECTRUM SETTINGS ------
spectrum_bins_name = "Lambda"
spectrum_value_name = "A40P0.50"
# python spectra are per cm2 at 100 pc -> we need to divide ΔC by this value too
spectrum_value_units = u.angstrom * u.erg / (u.s * u.angstrom * (np.pi * np.power(u.parsec.to(u.cm)*100, 2) * u.cm * u.cm))
#spectrum_value_to_lightcurve_value = (np.pi * np.power(u.parsec.to(u.cm) * 100, 2)) * u.cm * u.cm / 1000
#spectrum_value_units = 1e-14 * u.erg / (u.s * u.angstrom * u.cm * u.cm)
spectrum_wave_units = u.angstrom

# ----- SPECTRUM (FULL) SETTINGS -----
# Used to create the non-continuum-subtracted spectra for MEMECHO
# As much continuum as possible! Can't go further due to the limits of Python.
spectrum_full_wave_range = [tf_wave - 6200, tf_wave + 7000] * u.angstrom
spectrum_full_rebin_to = 100
spectrum_full_continuum_error = 0.01

# ----- SPECTRUM (SUBTRACTED) SETTINGS -----
# Used to create the continuum-subtracted spectra for CARAMEL
# As little continuum as possible! We want just the line to make it faster.
spectrum_line_wave_range = [tf_wave - 250, tf_wave + 250] * u.angstrom
spectrum_line_subtract_range = [tf_wave - 190, tf_wave + 190] * u.angstrom
spectrum_line_rebin_to = 50
spectrum_line_error = 0.02

# ------ SPECTRA TIMES ---------"
spectra_times_file = "spectra_times.dat"
spectra_times_units = ucds.MJD
spectra_fudge_factor = 1

# ---- VISUALIZATION & OUTPUT SETTINGS ------
trailed_spectrogram_file = "out_trailedspectrogram_qso"
animation_file = "out_timeseries_qso.mp4"
is_reversed = False
debug_spectra_file = "out_spectra_qso"
debug_lightcurve_file = "out_lightcurve"

# --------- RESCALING ----------
delay_max = 30 * u.d
delta_continuum_range = 0.1


# ===============
# Program begins!
# ===============

print("=== tssproduce started! ===")

# Import all data files
print("Importing spectrum file '{}'...".format(spectrum_file))
spectrum_line = import_spectrum(spectrum_file, spectrum_bins_name, spectrum_value_name,
                            frequency=False,
                            wave_units=spectrum_wave_units,
                            value_units=spectrum_value_units,
                            error_relative=spectrum_line_error,
                            limits=spectrum_line_wave,
                            subtract_continuum_with_mask=spectrum_line_subtract_range,
                            rebin_to=spectrum_line_rebin_to)
spectrum_full = import_spectrum(spectrum_file, spectrum_bins_name, spectrum_value_name,
                            frequency=False,
                            wave_units=spectrum_wave_units,
                            value_units=spectrum_value_units,
                            error_absolute=spectrum_line['error'][0],
                            limits=spectrum_full_wave_range,
                            rebin_to=spectrum_full_rebin_to)

print("Importing lightcurve file '{}'...".format(lightcurve_file))
lightcurve = import_lightcurve(lightcurve_file,
                            time_units=lightcurve_time_units,
                            value_units=lightcurve_value_units,
                            target_bolometric_luminosity=lightcurve_target_lum,
                            delta_continuum_range=0.1)

print("Importing spectra timing file '{}'...".format(spectra_times_file))
spectra_times = import_spectra_times(spectra_times_file,
                            time_units=spectra_times_units)


# Produce a TF
print("Generating Ψ...")
def generate_tf(databases, spectrum, delay_bins, line, wave, limit=999999999):
    db_mid = tfpy.open_database(databases['mid']['path'], "root", "password")
    db_min = tfpy.open_database(databases['min']['path'], "root", "password")
    db_max = tfpy.open_database(databases['max']['path'], "root", "password")
    # Take spectrum and produce full list of wavelength bins for it to feed to the TF
    bounds = np.array(list(spectrum["wave_min"]).append(spectrum["wave_max"][-1]))

    tf_mid = tfpy.TransferFunction( db_mid, "qso100", continuum=databases['mid']['continuum'],
                wave_bins=(len(bounds)-1), delay_bins=delay_bins)
    tf_mid.line(line, wave).wavelength_bins(bounds).delay_dynamic_range(2).run(
                scaling_factor=databases['mid']['scale'], limit=limit,verbose=True).plot()
    tf_min = tfpy.TransferFunction(db_min, "qso090", continuum=databases['min']['continuum'], template=tf_mid).run(
                scaling_factor=databases['min']['scale'], limit=lim_test).plot()
    tf_max = tfpy.TransferFunction(db_max, "qso110", continuum=databases['max']['continuum'], template=tf_mid).run(
                scaling_factor=databases['max']['scale'], limit=lim_test).plot()
    tf_mid.response_map_by_tf(tf_min, tf_max).plot(response_map=True, name='resp')
    return tf_mid

# Repeat for full and subtracted spectra
tf_full = generate_tf(databases, spectrum_full, tf_delay_bins, tf_line, tf_wave, tf_lim_test)
tf_line = generate_tf(databases, spectrum_line, tf_delay_bins, tf_line, tf_wave, tf_lim_test)

# ========
# Convolve
# ========

# Set up the spectra array and set to baseline of start spectrum
spectra_full = ap.table.Table([spectrum_full['wave'], spectrum_full['value']])
spectra_line = ap.table.Table([spectrum_line['wave'], spectrum_line['value']])
# For each point in the lightcurve file

for time in spectra_times['time']:
    # Add a base spectrum to the output spectra file
    spectra_full["value {}".format(time)] = spectra['value']
    spectra_full["value {}".format(time)].meta['time'] = time
    spectra_line["value {}".format(time)] = spectra['value']
    spectra_line["value {}".format(time)].meta['time'] = time

# We need to evaluate at every bin-width's step
delay_bins = (tf_mid.delay_bins() * u.s).to(lightcurve['time'].unit)
delay_midpoints = np.zeros(len(delay_bins))
for i in range(0, len(delay_bins)-1):
    delay_midpoints = (delay_bins[i] + delay_bins[i+1]) / 2

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
def add_step_response(spectrum, spectra, stepname, delta_continuum, response):
    for j in range(0, len(spectra)):
        # Matching the format of spectra.c line:
        # x *= (freq * freq * 1e-8) / (dfreq * dd * C);
        dfreq = spectrum["freq_max"].quantity[j] - spectrum["freq_min"].quantity[j]
        invwave = (spectrum["freq"].quantity[j] / apc.c).to(1/spectrum_wave_units)
            spectra[column][j] += (delta_continuum * response[j] *
            invwave * spectrum["freq"][j] / dfreq).value

print("Beginning {} time steps to generate {} spectra...".format(len(times), len(spectra.columns)-2))
# For each timestep, we send out a 'pulse' of continuum
for step, time in enumerate(times):
    # For each time bin this pulse is smeared out over
    print("Step {}: {:.2f}%".format(step, 100*step/len(times)))
    for i in range(0, len(delay_bins)-1):
        # Figure out what the bounds are for the region this pulse affects
        time_range = [time + delay_bins[i], time + delay_bins[i+1]]
        for column in spectra_line.colnames[2:]:
            # For each spectrum, if it occurs at a timestamp within this bin
            if time_range[0] < spectra_line[column].meta['time'] < time_range[1]:
                    # Add this pulse's contribution to it
                    add_step_response(  spectrum_line, spectra_line, column, delta_continuum[step],
                                        tf_line.response(delay_index=i)  )
                    add_step_response(  spectrum_full, spectra_full, column, delta_continuum[step],
                                        tf_full.response(delay_index=i) )

# Now we have the final spectra, add the experimental errors
for column in spectra.colnames[2:]:
    # For each spectrum
    if 'value' in column:
        for i in range(0, len(spectrum_bounds)-1):
            # For each wavelength bin in each spectrum, add a random error
            spectra_full[column][i] += npr.normal(scale=spectrum_full['error'][i])
            spectra_line[column][i] += npr.normal(scale=spectrum_line['error'][i])


# Write the spectra out to an easily plotted file
ap.io.ascii.write(spectra_full, debug_spectra_file+'_full.dat', overwrite=True)
ap.io.ascii.write(spectra_line, debug_spectra_file+'_line.dat', overwrite=True)
ap.io.ascii.write(lightcurve, debug_lightcurve_file+'.dat', overwrite=True)

output_CARAMEL(lightcurve, spectra_line, spectra_times, spectrum_line, suffix)
output_MEMECHO(lightcurve, spectra_full, spectra_times, spectrum_full, suffix)


# ==============================================================================
# OUTPUT VISUALIZATIONS
# ==============================================================================
# Generate trailed spectrogram
# ------------------------------------------------------------------------------
print("Generating trailed spectrogram...")

def generate_trailed_spectrogram(spectra, spectra_times, filename):
    fig, (ax_ts, ax_c) = plt.subplots(1, 2, gridspec_kw={'width_ratios':[3,1]})
    pcolor_data = np.zeros([len(spectra_times), len(spectra)])
    for i, column in enumerate(spectra.colnames[2:]):
        pcolor_data[i,:] = spectra[column]

    time_bounds = np.zeros(len(spectra_times)+1)
    for i in range(0, len(spectra_times)-1):
        time_bounds[i+1] = (spectra_times['time'][i] + spectra_times['time'][i+1]).value / 2
    time_bounds[0] = (spectra_times['time'][0] - (spectra_times['time'][1] - spectra_times['time'][0]) / 2).value
    time_bounds[-1] = (spectra_times['time'][-1] + (spectra_times['time'][-1] - spectra_times['time'][-2]) / 2).value

    pcol = ax.pcolor(spectrum_bounds, time_bounds, pcolor_data)
    ax_ts.set_xlabel("λ ({})".format(spectrum_wave_units))
    ax_ts.set_ylabel("L ({})".format(spectra_times_units))
    ax_c.invert_xaxis()
    ax_c.plot(lightcurve['value'], lightcurve['time'], '-', c='m')
    ax_c.set_xlabel("Continuum luminosity ()")

    fig.colorbar(pcol).set_label("Flux ({})".format(spectrum_value_units))
    plt.savefig("{}.eps".format(trailed_spectrogram_file), bbox_inches='tight')

# ------------------------------------------------------------------------------
# Generate animation
# ------------------------------------------------------------------------------
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
        # If this timestep is in the range of spectra, assign it a colour
        time_colour = plt.cm.jet( (times[step]-spectra_times['time'][0]) /
                                (spectra_times['time'][-1] - spectra_times['time'][0]) )

    if times[step] > spectra.columns[next_column].meta['time']:
        # If we've gone past the current column,
        line_new = axis.plot(spectra['wave'], spectra.columns[next_column], '-', markersize=0.5,
                             c=time_colour, zorder=zorder(step, len(times), is_reversed))
        times_for_line = np.array([times[step].value, times[step].value])
        values_for_line = np.array([y_dc_min.value, y_dc_max.value])
        line_vertical = ax_dc.plot(times_for_line, values_for_line, ':', linewidth=1, c=time_colour, zorder=0)
        artists.append(line_vertical)
        artists.append(line_new)
        next_column += 1

    # Plot the continuum brightness for this step in the appropriate colour
    point_new = ax_dc.plot(times[step], delta_continuum[step], 'o', c=time_colour, zorder=zorder(step, len(times), is_reversed))
    artists.append(point_new)
    return artists,

# Generate animation and save to file
animation = ani.FuncAnimation(figure, update_figure, frames=len(times))
animation.save(animation_file, fps=24)
