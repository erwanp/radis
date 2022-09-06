# -*- coding: utf-8 -*-
"""
Created on Sun Aug 22 13:58:04 2021

@author: erwan
"""

# TODO: add labellines

import astropy.units as u
import matplotlib.pyplot as plt

from radis import Spectrum, calc_spectrum, get_residual
from radis.test.utils import getValidationCase

#%% CO 1900 - 2400 cm-1 ; 70 bar
title='CO-CO 1900-2400 cm-1  297K 70bar'

s_exp = Spectrum.from_txt(getValidationCase("Brodbeck1994_CO_1900-2400cm1_70bar_297K_transmittance.txt"),
                          'transmittance',
                          'cm-1', unit='', delimiter=',')
s_exp.plot(lw=2, color='k')

res = {}
for truncation in [10, 30, 50, 75, 100, 200, 300, 500]:

    s = calc_spectrum(1900 / u.cm, 2400 / u.cm,
                      molecule='CO',
                      isotope='1,2,3',
                      pressure=70 * u.bar,
                      Tgas=297,           # K
                      mole_fraction=1,
                      path_length=7.1 * u.cm,
                      databank='hitemp',  # or use 'hitemp',
                      truncation=truncation,
                      )
    s.apply_slit(0.2, "cm-1")  # cf [Brodbeck1994]
    s.name = f' truncation: {s.conditions["truncation"]}cm-1'

    res[truncation] = get_residual(s_exp, s, 'transmittance', ignore_nan=True)

    # plot_diff(s_exp, s)
    s.plot("transmittance", nfig='same')

plt.figure()
plt.title(title)
plt.plot(res.keys(), res.values(), 'ok')
plt.xlabel('Truncation')
plt.ylabel('Residual')

#%%
title='[CO-CO 1900-2400 cm-1 691K 24.4bar'

s_exp = Spectrum.from_txt(getValidationCase("Brodbeck1994_CO_1900-2400cm1_26.4bar_691K_transmittance.txt"),
                          'transmittance',
                          'cm-1', unit='', delimiter=',')
s_exp.plot(lw=2, color='k')

res = {}
for truncation in [10, 30, 50, 75, 100, 200, 300, 500]:

    s = calc_spectrum(1900 / u.cm, 2400 / u.cm,
                      molecule='CO',
                      isotope='1,2,3',
                      pressure=26.4 * u.bar,
                      Tgas=691,           # K
                      mole_fraction=1,
                      path_length=7.1 * u.cm,
                      databank='hitemp',  # or use 'hitemp',
                      truncation=truncation,
                      )
    s.apply_slit(0.2, "cm-1")  # cf [Brodbeck1994]
    s.name = f' truncation: {s.conditions["truncation"]}cm-1'

    res[truncation] = get_residual(s_exp, s, 'transmittance')

    # plot_diff(s_exp, s)
    s.plot("transmittance", nfig='same')


plt.figure()
plt.title(title)
plt.plot(res.keys(), res.values(), 'ok')
plt.xlabel('Truncation')
plt.ylabel('Residual')
