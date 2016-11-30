#bindata_keplerian.py
import csv
import numpy as np
import scipy as sp
import astropy as ap
from astropy import units as u
from astropy.coordinates import Angle

def calculate_delay(r_angle, r_phase, r_rad, timescale):
	"""Delay relative to continuum observed at angle for emission at radius"""
	# Calculate delay compared to continuum photons
	# Draw plane at r_rad_min out. Find x projection of disk position.
	# Calculate distance travelled to that plane from the current disk position
	# Delay relative to continuum is thus (distance from centre to plane)
	# + distance from centre to point
	vr_disk 	= np.array([r_rad.value*np.cos(r_phase), 0.0]) * u.m
	vr_normal 	= np.array([np.cos(r_angle), np.sin(r_angle)])
	vr_plane 	= r_rad * vr_normal
	return (np.dot((vr_plane - vr_disk), vr_normal) * u.m / ap.constants.c).to(timescale)

def keplerian_velocity(r_mass, r_rad):
	"""Calculates Keplerian velocity at radius"""
	return np.sqrt(ap.constants.G * r_mass / r_rad)

def path_to_delay(r_path):
	"""Converts path to time delay"""
	return r_path / ap.constants.C

def doppler_shift(r_wave, r_vel):
	"""Doppler shifts the passed wavelength"""
	return r_wave * (1 + (r_vel / ap.constants.c))

# INPUTS
#r_line_inp = 1550.7
r_line_inp = 6562.8
r_rad_min_inp = 50
r_rad_max_inp = 100
#r_mass_bh_inp = 1e9
r_mass_bh_inp = 1e7
r_angle_inp = 40
i_steps 	= 100

# PROCESSING
print "=== STARTING ==="

print "--- PARSE INPUT ---"
r_line  	= r_line_inp * u.angstrom
r_angle 	= Angle(r_angle_inp, u.deg)
r_mass_bh 	= r_mass_bh_inp * ap.constants.M_sun
r_rad_grav 	= (6 * ap.constants.G * r_mass_bh / np.power(ap.constants.c, 2))
r_rad_min 	= r_rad_min_inp * r_rad_grav
r_rad_max 	= r_rad_max_inp * r_rad_grav

print "--- SETUP ARRAYS ---"
ar_wave  	= np.zeros(i_steps) * u.angstrom
ar_delay 	= np.zeros(i_steps) * u.s
ar_phase 	= np.linspace(0, np.pi*2, i_steps)
ar_rad		= np.linspace(r_rad_min,r_rad_max*20.0,i_steps)

print "--- SETUP OUTPUT ---"

csvfile = open('bindata_keplerian.csv', 'w')
csvout = csv.writer(csvfile, delimiter=' ', quotechar='"', quoting=csv.QUOTE_MINIMAL)
as_row_titles = ["Wavelength", "Delay"]
csvout.writerow(as_row_titles)


# ITERATE OVER INNER EDGE
r_vel_max = keplerian_velocity( r_mass_bh, r_rad_min )
print "Wave: %.2e, Vel: %.2e" % (doppler_shift(r_line, r_vel_max).value, r_vel_max.value)
print "Wave: %.2e, Vel: %.2e" % (doppler_shift(r_line, r_vel_max*np.sin(r_angle)).value, r_vel_max.value*np.sin(r_angle))
for r_phase, r_wave, r_delay in np.nditer([ar_phase, ar_wave, ar_delay], op_flags=['readwrite']):
	r_vel 			= r_vel_max *  np.sin(r_phase) * np.sin(r_angle)
	r_wave[...]		= doppler_shift( r_line, r_vel )
	r_delay[...]	= calculate_delay( r_angle, r_phase, r_rad_min, u.day)

for r_wave, r_delay in np.nditer([ar_wave, ar_delay]):
	csvout.writerow((r_wave, r_delay))
csvout.writerow([])

# # ITERATE OVER OUTER EDGE
# r_vel_max = keplerian_velocity( r_mass_bh, r_rad_max )
# for r_phase, r_wave, r_delay in np.nditer([ar_phase, ar_wave, ar_delay], op_flags=['readwrite']):
# 	r_vel 			= r_vel_max * np.sin(r_phase) * np.sin(r_angle)
# 	r_wave[...] 	= doppler_shift( r_line, r_vel )
# 	r_delay[...]	= calculate_delay( r_angle, r_phase, r_rad_max, u.day)

# for r_wave, r_delay in np.nditer([ar_wave, ar_delay]):
# 	csvout.writerow((r_wave, r_delay))
# csvout.writerow([])


# ITERATE OVER BLUE BOUND
for r_rad, r_wave, r_delay in np.nditer([ar_rad, ar_wave, ar_delay], op_flags=['readwrite']):
	r_rad 			= r_rad * u.m
	r_vel 			= keplerian_velocity( r_mass_bh, r_rad ) * np.sin(r_angle)
	r_wave[...]		= doppler_shift( r_line, r_vel )
	r_delay[...]	= calculate_delay( r_angle, np.pi/2, r_rad, u.day )

for r_wave, r_delay in np.nditer([ar_wave, ar_delay]):
	csvout.writerow((r_wave, r_delay))
csvout.writerow([])

# ITERATE OVER RED BOUND
for r_rad, r_wave, r_delay in np.nditer([ar_rad, ar_wave, ar_delay], op_flags=['readwrite']):
 	r_rad 			= r_rad * u.m
 	r_vel 			= -keplerian_velocity( r_mass_bh, r_rad ) * np.sin(r_angle)
 	r_wave[...]		= doppler_shift( r_line, r_vel )
 	r_delay[...]	= calculate_delay( r_angle, np.pi/2, r_rad, u.day )

for r_wave, r_delay in np.nditer([ar_wave, ar_delay]):
 	csvout.writerow((r_wave, r_delay))
csvout.writerow([])
