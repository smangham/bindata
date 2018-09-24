import numpy as np
import tfpy
import sqlalchemy
import matplotlib.pyplot as plt
import datetime
from matplotlib.colors import LogNorm

n = 10000
limit = False

qso_100 = tfpy.open_database('qso_100', "root", "password")
Session = sqlalchemy.orm.sessionmaker(bind=qso_100)
session_100 = Session()

print('Query at', datetime.datetime.now())
query = session_100.query(tfpy.Photon.X, tfpy.Photon.Y, tfpy.Photon.Z, tfpy.Photon.Weight)
query.filter(tfpy.Photon.Resonance == 44)

data = None
if limit:
    data = np.asarray(query.limit(n).all())
else:
    data = np.asarray(query.all())

angle = np.radians(40)

bins = 200
factor = 0.193789563

bound_xy = 2e17 * factor
bound_z = 5e16 * factor

scale_xy = 1e16
scale_z = 1e15
scale_xy_str = '10^{16}'
scale_z_str = '10^{15}'

bounds_xy = [-bound_xy/scale_xy, bound_xy/scale_xy]
bounds_z = [-bound_z/scale_z, bound_z/scale_z]
bounds_z = [0, bound_z/scale_z]

x = data[:, 0] * factor / scale_xy
y = data[:, 1] * factor / scale_xy
z = data[:, 2] * factor / scale_z
w = data[:, 3]

min_x = np.amin(x)
min_y = np.amin(y)
min_z = np.amin(z)

max_x = np.amax(x)
max_y = np.amax(y)
max_z = np.amax(z)

bins_x = np.linspace(min_x, max_x, bins, endpoint=True)
bins_y = np.linspace(min_y, max_y, bins, endpoint=True)
bins_z = np.linspace(min_z, max_z, bins, endpoint=True)
bins_z = np.linspace(0, max_z, bins, endpoint=True)

fig, ((ax_xy, ax_dummy), (ax_xz, ax_yz)) = plt.subplots(2, 2, sharey='row', sharex='col')
ax_dummy.set_axis_off()

fig.subplots_adjust(hspace=0, wspace=0)
ax_xz.set_xlabel(r'X (${}$ cm)'.format(scale_xy_str))
ax_yz.set_xlabel(r'Y (${}$ cm)'.format(scale_xy_str))
ax_xz.set_ylabel(r'Z (${}$ cm)'.format(scale_z_str))

ax_xy.set_xlabel(r'X (${}$ cm)'.format(scale_xy_str))
ax_xy.set_ylabel(r'Y (${}$ cm)'.format(scale_xy_str))

(vals_yz, dummy, dummy) = np.histogram2d(y, z, bins=[bins_y, bins_z], weights=w)
(vals_xz, dummy, dummy) = np.histogram2d(x, z, bins=[bins_x, bins_z], weights=w)
(vals_xy, dummy, dummy) = np.histogram2d(x, y, bins=[bins_x, bins_y], weights=w)
cmax = np.max([np.max(vals_yz), np.max(vals_xz), np.max(vals_xy)])

im_xz = ax_xz.pcolor(bins_x, bins_z, vals_xz.T/cmax, vmin=0, vmax=1)
ax_xz.plot([0, np.sin(angle)*100], [0, np.cos(angle)*100], color='red', linewidth=2, alpha=0.6)
im_yz = ax_yz.pcolor(bins_y, bins_z, vals_yz.T/cmax, vmin=0, vmax=1)
ax_yz.plot([0, 0], [0, 100], color='red', linewidth=2, alpha=0.6)
im_xy = ax_xy.pcolor(bins_x, bins_y, vals_xy.T/cmax, vmin=0, vmax=1)
ax_xy.plot([0, 100], [0, 0], color='red', linewidth=2, alpha=0.6)

ax_xz.set_xlim(bounds_xy)
ax_xz.set_ylim(bounds_z)

ax_yz.set_xlim(bounds_xy)
ax_yz.set_ylim(bounds_z)

ax_xy.set_xlim(bounds_xy)
ax_xy.set_ylim(bounds_xy)

ax_cbar = fig.add_axes([0.55, 0.55, 0.05, 0.3])
cbar = fig.colorbar(im_xz, cax=ax_cbar)
cbar.set_label(r'$L/L_{max}$')

fig.savefig('lum-xyz.eps')

# =============

qso_110 = tfpy.open_database('qso_110', "root", "password")
Session = sqlalchemy.orm.sessionmaker(bind=qso_110)
session_110 = Session()

print('Query at', datetime.datetime.now())
query = session_110.query(tfpy.Photon.X, tfpy.Photon.Y, tfpy.Photon.Z, tfpy.Photon.Weight)
query.filter(tfpy.Photon.Resonance == 44)

data_110 = None
if limit:
    data_110 = np.asarray(query.limit(n).all())
else:
    data_110 = np.asarray(query.all())

x_110 = data_110[:, 0] * factor / scale_xy
y_110 = data_110[:, 1] * factor / scale_xy
z_110 = data_110[:, 2] * factor / scale_z
w_110 = data_110[:, 3]

qso_090 = tfpy.open_database('qso_090', "root", "password")
Session = sqlalchemy.orm.sessionmaker(bind=qso_090)
session_090 = Session()

print('Query at', datetime.datetime.now())
query = session_090.query(tfpy.Photon.X, tfpy.Photon.Y, tfpy.Photon.Z, tfpy.Photon.Weight)
query.filter(tfpy.Photon.Resonance == 44)

if limit:
    data_090 = np.asarray(query.limit(n).all())
else:
    data_090 = np.asarray(query.all())

x_090 = data_090[:, 0] * factor / scale_xy
y_090 = data_090[:, 1] * factor / scale_xy
z_090 = data_090[:, 2] * factor / scale_z
w_090 = -data_090[:, 3]

x = np.hstack((x_090, x_110))
y = np.hstack((y_090, y_110))
z = np.hstack((z_090, z_110))
w = np.hstack((w_090, w_110))

fig, ((ax_xy, ax_dummy), (ax_xz, ax_yz)) = plt.subplots(2, 2, sharey='row', sharex='col')
ax_dummy.set_axis_off()

fig.subplots_adjust(hspace=0, wspace=0)
ax_xz.set_xlabel(r'X (${}$ cm)'.format(scale_xy_str))
ax_yz.set_xlabel(r'Y (${}$ cm)'.format(scale_xy_str))
ax_xz.set_ylabel(r'Z (${}$ cm)'.format(scale_z_str))

ax_xy.set_xlabel(r'X (${}$ cm)'.format(scale_xy_str))
ax_xy.set_ylabel(r'Y (${}$ cm)'.format(scale_xy_str))

for tk in ax_yz.get_xticklabels():
        tk.set_visible(True)

(vals_yz, dummy, dummy) = np.histogram2d(y, z, bins=[bins_y, bins_z], weights=w)
(vals_xz, dummy, dummy) = np.histogram2d(x, z, bins=[bins_x, bins_z], weights=w)
(vals_xy, dummy, dummy) = np.histogram2d(x, y, bins=[bins_x, bins_y], weights=w)
cmax = np.max([np.max(np.abs(vals_yz)), np.max(np.abs(vals_xy)), np.max(np.abs(vals_xz))])

im_xz = ax_xz.pcolor(bins_x, bins_z, vals_xz.T/cmax, vmin=-1, vmax=1, cmap='RdBu_r')
ax_xz.plot([0, np.sin(angle)*100], [0, np.cos(angle)*100], color='black', linewidth=2, alpha=0.6)
im_yz = ax_yz.pcolor(bins_y, bins_z, vals_yz.T/cmax, vmin=-1, vmax=1, cmap='RdBu_r')
ax_yz.plot([0, 0], [0, 100], color='black', linewidth=2, alpha=0.6)
im_xy = ax_xy.pcolor(bins_x, bins_y, vals_xy.T/cmax, vmin=-1, vmax=1, cmap='RdBu_r')
ax_xy.plot([0, 100], [0, 0], color='black', linewidth=2, alpha=0.6)

ax_xz.set_xlim(bounds_xy)
ax_xz.set_ylim(bounds_z)

ax_yz.set_xlim(bounds_xy)
ax_yz.set_ylim(bounds_z)

ax_xy.set_xlim(bounds_xy)
ax_xy.set_ylim(bounds_xy)

ax_cbar = fig.add_axes([0.55, 0.55, 0.05, 0.3])
cbar = fig.colorbar(im_xz, cax=ax_cbar)
cbar.set_label(r'$\Delta L/\Delta L_{max}$')
fig.savefig('lum-xyz_resp.eps')

# =============

r = np.sqrt(x*x + y*y)
bins_r = np.linspace(0, max_x, bins, endpoint=True)


fig_rz, ax_rz = plt.subplots()
ax_rz.set_xlabel(r'R (${}$ cm)'.format(scale_xy_str))
ax_rz.set_ylabel(r'Z (${}$ cm)'.format(scale_z_str))

(vals_rz, dummy, dummy) = np.histogram2d(r, z, bins=[bins_r, bins_z], weights=w)
cmax = np.max([np.max(np.abs(vals_rz)), np.max(np.abs(vals_rz))])

im_rz = ax_rz.pcolor(bins_r, bins_z, vals_rz.T/cmax, vmin=-1, vmax=1, cmap='RdBu_r')
ax_rz.plot([0, np.sin(angle)*100], [0, np.cos(angle)*100], color='black', linewidth=2, alpha=0.6)

ax_rz.set_xlabel(r'Radius (${}$ cm)'.format(scale_xy_str))
ax_rz.set_ylabel(r'Z (${}$ cm)'.format(scale_z_str))

ax_rz.set_xlim([0, bounds_xy[-1]])
ax_rz.set_ylim(bounds_z)

fig_rz.colorbar(im_xz).set_label(r'$\Delta L/\Delta L_{max}$')
fig_rz.savefig('lum-rz_resp.png')
