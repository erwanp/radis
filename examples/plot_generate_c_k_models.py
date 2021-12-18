# -*- coding: utf-8 -*-
"""
==============================
Generate a grid of C-K spectra
==============================

Generate a grid of spectra for various temperature ``Tgas``, pressure ``pressure_mbar``,
to be used in atmosphere calculations

Then uses :py:mod:`exo_k` to derive the cross-section tables, then the
k-tables.

"""

import matplotlib.pyplot as plt
import numpy as np

from radis import SpectrumFactory

sf = SpectrumFactory(
    wavenum_min=2900,
    wavenum_max=3200,
    molecule="OH",
    truncation=5,  # cm-1
    verbose=0,  # more for more details
    wstep="auto",
    # wstep=0.003,
)
sf.fetch_databank("hitemp")

# This line is very important :
# it initiate a SpecDatabase to automatically retrieve spectra if they
# exist already, or compute them if they don't

sf.init_database(
    "Highres_spectra_database",
    add_info=[
        "Tgas",
        "pressure_mbar",
    ],  #  info that appear in file names. All conditions are stored in metadata anyway
    # autoretrieve='force'  # make sure we do not recompute. Useful for debugging.
)

db = sf.SpecDatabase
if db:
    print(db)  # print database if not empty already initialized

#%%
# Generate a database of Spectra:

p_grid = [1e-3, 1e-2, 1e-1, 1]  # bar
T_grid = [1000, 1250, 1600]  # K

for Tgas in T_grid:
    for p in p_grid:
        sf.eq_spectrum(Tgas=Tgas, pressure=p)  # automatically stores in database


# %%
# Plot the spectra precomputed in database
db.plot_cond("Tgas", "pressure_mbar")

# Generate an (array) grid of filenames :
file_grid = db.create_fname_grid(["pressure_mbar", "Tgas"])
print(file_grid)

#%% Note : if we wanted to interpolate all spectra on the same range, we could do :

# # Resample all on spectrum of minimum wstep
# s_wstep_min = db.get(wstep=float(db.see("wstep").min()))[0]
# db.map(lambda s: s.resample(s_wstep_min))
# # Export to a new database:
# db.compress_to(db.path+'_interp', if_exists_then='replace')


# %% Now send it to exo-k to generate a C-K table

# Inspired from http://perso.astrophy.u-bordeaux.fr/~jleconte/exo_k-doc/examples-exo_k.html#Creating-k-coefficients-for-a-new-species-not-in-ExoMol-from-high-resolution-spectra-from-the-petitRADTRANS-database


import exo_k as xk

Hires_spectra = xk.hires_to_xtable(
    path=db.path,
    filename_grid=file_grid,
    logpgrid=[np.log10(float(p)) for p in p_grid],
    tgrid=T_grid,
    mol="OH",
    grid_p_unit="bar",
    p_unit="bar",
    binary=True,
    mass_amu=17.0,
)

## We create a custom g-grid with 8 gauss legendre points between 0 and 0.95
##   and 8 points between 0.95 and 1 (as usual with petitRADTRANS data)
weights, ggrid, gedges = xk.split_gauss_legendre(order=16, g_split=0.95)

wnedges = Hires_spectra.wnedges  # reduce (??)

ktab = xk.Ktable(
    xtable=Hires_spectra, wnedges=wnedges, weights=weights, ggrid=ggrid, p_unit="bar"
)  # , remove_zeros=True)


fig, ax = plt.subplots(1, 1, sharey=False, figsize=(8, 4))
# Hires_spectra.plot_spectrum(ax, p=1.e3, t=1300., xscale='log', yscale='log', label='p=1 mbar')
Hires_spectra.plot_spectrum(ax, p=1, t=1300.0, yscale="log", label="p=1 bar", lw=2)
ktab.plot_spectrum(ax, p=1, t=1300.0, g=1, label="ktab, g=1, R=1000, p=1 bar")
ax.legend(loc="upper right")


# %% Compare

plt.figure(figsize=(16, 6))
db.get_unique(pressure_mbar=1000, Tgas=1250).plot(
    "xsection",
    wunit="cm-1",
    yscale="log",
    lw=2,
    label="RADIS LBL",
    Iunit=Hires_spectra.kdata_unit,
    nfig="same",
)
Hires_spectra.plot_spectrum(
    plt.gca(), p=1, t=1250.0, label="Exo-K Hi-res", x_axis="cm-1"
)
assert Hires_spectra.p_unit == "bar"
plt.legend()


# %% Compare perf

# (disconnect SpecDatabase first)
sf.autoupdatedatabase = False
sf.autoretrievedatabase = False


I_ck = ktab.spectrum_to_plot(t=1300, p=1, g=1)
s = sf.eq_spectrum(Tgas=1300, pressure=1)

"""
%timeit ktab.spectrum_to_plot(t=1300, p=1, g=1)
57.9 ms ± 2.27 ms per loop (mean ± std. dev. of 7 runs, 10 loops each)

%timeit sf.eq_spectrum(Tgas=1300, pressure=1)
23.2 ms ± 931 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)
"""

plt.figure(figsize=(16, 6))
s.plot(
    "xsection",
    lw=2,
    label="RADIS LBL",
    yscale="log",
    Iunit=ktab.kdata_unit,
    nfig="same",
)
plt.plot(ktab.wns, I_ck, lw=1.5, label="Exo-K C-K")
plt.legend()
