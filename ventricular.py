#!/usr/bin/env python3
#
# Relative contributions of the major ionic currents in human ventricular
# models.
#
import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import myokit
import myokit.lib.plots as mp

import shared


# Update matplotlib styles
matplotlib.rcParams['axes.spines.right'] = False
matplotlib.rcParams['axes.spines.top'] = False
matplotlib.rcParams['mathtext.default'] = 'regular'

# Current colors
cmap = matplotlib.cm.get_cmap('tab20')
current_colours = {
    'I_Kr': 0,
    'I_Ks': 1,
    'I_to': 2,
    'I_Kb': 3,
    # I_f
    # I_Kur
    'I_K1': 4,
    'I_NaK': 5,
    'I_Na': 16,
    'I_NaL': 17,
    'I_CaL': 10,
    'I_NaCa': 11,
    'I_Na,B': 12,
    'I_Ca,B': 13,
    'I_ClCa': 6,
    'I_Cl,B': 7,
    'I_Ca,P': 14,
    'I_K,ACh': 18,
    'I_K,ATP': 19,
}

# Human atrial models
model_names = {
    #'priebe': 'priebe-1998.mmt',
    #'iyer': 'iyer-2004.mmt',
    #'grandi': 'grandi-2010.mmt',
    #'tnnp': 'tentusscher-2004.mmt',
    #'tp': 'tentusscher-2006.mmt',
    #'ohara': 'ohara-2011.mmt',
    #'cipa': 'ohara-cipa-v1-2017.mmt',
    'tomek': 'tomek-2020-chloride-epi.mmt',
}

fancy_names = {
    'priebe': 'Priebe & Beuckelmann, 1998',
    'iyer': 'Iyer et al., 2004 (epi)',
    'grandi': 'Grandi et al., 2010 (epi)',
    'tnnp': 'Ten Tusscher et al., 2004 (epi)',
    'tp': 'Ten Tusscher & Panfilov 2006 (epi)',
    'ohara': 'O\'Hara et al., 2011 (epi)',
    'cipa': 'O\'Hara et al., 2017 CiPA (epi)',
    'tomek': 'Tomek et al., 2020 (epi)',
}


def current_variables(model, colours=False):
    """ Returns an ordered list of transmembrane current variable names. """
    name = model.name().lower()
    if 'priebe' in name:
        currents = {
            'I_NaCa': 'inaca.i_NaCa',
            'I_to': 'ito.i_to',
            'I_Ks': 'iks.i_Ks',
            'I_Kr': 'ikr.i_Kr',
            'I_K1': 'ik1.i_K1',
            'I_NaK': 'inak.i_NaK',
            'I_CaL': 'ica.i_Ca',
            'I_Ca,B': 'icab.i_b_Ca',
            'I_Na,B': 'inab.i_b_Na',
            'I_Na': 'ina.i_Na',
        }
    elif 'iyer' in name:
        currents = {
            'I_NaCa': 'inaca.inaca',
            'I_to': 'ito.Ito1',
            'I_Ks': 'iks.iks',
            'I_Kr': 'ikr.ikr',
            'I_Ca,P': 'ipca.ipca',
            'I_K1': 'ik1.ik1',
            'I_NaK': 'inak.inak',
            'I_CaL': 'ical.ICa_total',
            'I_Ca,B': 'icab.icab',
            'I_Na,B': 'inab.inab',
            'I_Na': 'ina.ina',
        }
    elif 'grandi' in name:
        currents = {
            'I_Cl,B': 'iclb.IClb',
            'I_ClCa': 'iclca.iclca',
            'I_to': 'ito.ito',
            'I_Kb': 'ikp.I_kp',
            'I_Ks': 'iks.I_ks',
            'I_Kr': 'ikr.I_kr',
            'I_Ca,P': 'ipca.I_pca',
            'I_K1': 'ik1.I_k1',
            'I_NaK': 'inak.I_nak',
            'I_CaL': 'ical.I_CaL',
            'I_NaCa': 'incx.I_ncx',
            'I_Ca,B': 'icabk.I_cabk',
            'I_Na,B': 'inab.I_nabk',
            'I_Na': 'ina.I_Na',
        }
    elif 'tusscher_2004' in name:
        currents = {
            'I_NaCa': 'inaca.INaCa',
            'I_to': 'ito.Ito',
            'I_Kb': 'ipk.IpK',
            'I_Ks': 'iks.IKs',
            'I_Kr': 'ikr.IKr',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_Ca,P': 'ipca.IpCa',
            'I_CaL': 'ical.ICaL',
            'I_Ca,B': 'icab.ICab',
            'I_Na,B': 'inab.INab',
            'I_Na': 'ina.INa',
        }
    elif 'tusscher-2006' in name:
        currents = {
            'I_to': 'ito.Ito',
            'I_Kb': 'ipk.IpK',
            'I_Ks': 'iks.IKs',
            'I_Kr': 'ikr.IKr',
            'I_K1': 'ik1.IK1',
            'I_NaCa': 'inaca.INaCa',
            'I_NaK': 'inak.INaK',
            'I_Ca,P': 'ipca.IpCa',
            'I_CaL': 'ical.ICaL',
            'I_Ca,B': 'icab.ICab',
            'I_Na,B': 'inab.INab',
            'I_Na': 'ina.INa',
        }
    elif 'ohara-2011' in name:
        currents = {
            'I_to': 'ito.Ito',
            'I_Kb': 'ikb.IKb',
            'I_Ks': 'iks.IKs',
            'I_Kr': 'ikr.IKr',
            'I_Ca,P': 'ipca.IpCa',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_CaL': 'ical.ICaL_total',
            'I_NaL': 'inal.INaL',
            'I_NaCa': 'inaca.INaCa_total',
            'I_Ca,B': 'icab.ICab',
            'I_Na,B': 'inab.INab',
            'I_Na': 'ina.INa',
        }
    elif 'cipa' in name:
        currents = {
            'I_to': 'ito.Ito',
            'I_Kb': 'ikb.IKb',
            'I_Ks': 'iks.IKs',
            'I_Kr': 'ikr.IKr',
            'I_Ca,P': 'ipca.IpCa',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_CaL': 'ical.ICaL_total',
            'I_NaL': 'inal.INaL',
            'I_NaCa': 'inaca.INaCa_total',
            'I_Ca,B': 'icab.ICab',
            'I_Na,B': 'inab.INab',
            'I_Na': 'ina.INa',
        }
    elif 'torord' in name:
        currents = {
            'I_Cl,B': 'ICl.IClb',
            'I_ClCa': 'ICl.IClCa',
            'I_to': 'Ito.Ito',
            'I_Kb': 'IKb.IKb',
            'I_Ks': 'IKs.IKs',
            'I_Kr': 'IKr.IKr',
            'I_K,ATP': 'I_katp.I_katp',
            'I_Ca,P': 'IpCa.IpCa',
            'I_K1': 'IK1.IK1',
            'I_NaK': 'INaK.INaK',
            'I_CaL': 'ICaL.ICaL',
            'I_NaL': 'INaL.INaL',
            'I_NaCa': 'INaCa.INaCa',
            'I_Ca,B': 'ICab.ICab',
            'I_Na,B': 'INab.INab',
            'I_Na': 'INa.INa',
        }

    else:
        currents = shared.guess_currents(model)
        print('\n'.join(currents))
        print(len(currents))
        raise NotImplementedError('Unknown model: ' + model.name())

    if colours:
        colours = [cmap(current_colours[x]) for x in currents.keys()]
        currents = list(currents.values())
        return currents, colours
    return list(currents.values())


# Create protocol
cl = 1000
protocol = myokit.pacing.blocktrain(cl, duration=0.5, offset=50)

# Load and prepare models
models = {}
for name, fname in model_names.items():
    pre_pace = True
    if 'priebe' in name or 'iyer' in name:
        pre_pace = False
    model = myokit.load_model(os.path.join('models', 'ventricular', fname))
    shared.prepare_model(model, protocol, current_variables(model), pre_pace)
    models[name] = model

# Maximum time to show in plots
tmax = 800

def text(ax, x, y, t, c='w'):
    ax.text(x, y, t, color=c, transform=ax.transAxes, fontweight='bold',
            horizontalalignment='right', verticalalignment='center')

# Create figure
fig = plt.figure(figsize=(9, 9))
fig.subplots_adjust(0.075, 0.05, 0.98, 0.97, hspace=0.35, wspace=0.2)
grid = GridSpec(3, 3)

#
# Top row: Various models
#
# Priebe & Beuckelmann 1998
code = 'priebe'
model = models[code]
currents, colours = current_variables(model, True)
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
d = s.run(tmax)
ax = fig.add_subplot(grid[0, 0])
ax.set_title(fancy_names[code])
ax.set_xlabel('Time (s)')
ax.set_ylabel('Relative contribution')
ax.set_xlim(0, tmax)
ax.set_ylim(-1.02, 1.02)
mp.cumulative_current(d, currents, ax, colors=colours, normalise=True)

# Iyer et al. 2004
code = 'iyer'
model = models[code]
currents, colours = current_variables(model, True)
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
d = s.run(tmax)
ax = fig.add_subplot(grid[0, 1])
ax.set_title(fancy_names[code])
ax.set_xlabel('Time (s)')
ax.set_yticklabels([])
ax.set_xlim(0, tmax)
ax.set_ylim(-1.02, 1.02)
mp.cumulative_current(d, currents, ax, colors=colours, normalise=True)

# Grandi 2010
code = 'grandi'
model = models[code]
currents, colours = current_variables(model, True)
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
d = s.run(tmax)
ax = fig.add_subplot(grid[0, 2])
ax.set_title(fancy_names[code])
ax.set_xlabel('Time (s)')
ax.set_yticklabels([])
ax.set_xlim(0, tmax)
ax.set_ylim(-1.02, 1.02)
mp.cumulative_current(d, currents, ax, colors=colours, normalise=True)

#
# Middle row: Ten Tusscher & Panfilov models
#
# TNNP 2004
code = 'tnnp'
model = models[code]
currents, colours = current_variables(model, True)
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
d = s.run(tmax)
ax = fig.add_subplot(grid[1, 0])
ax.set_title(fancy_names[code])
ax.set_xlabel('Time (s)')
ax.set_ylabel('Relative contribution')
ax.set_xlim(0, tmax)
ax.set_ylim(-1.02, 1.02)
mp.cumulative_current(d, currents, ax, colors=colours, normalise=True)


# TP 2006
code = 'tp'
model = models[code]
currents, colours = current_variables(model, True)
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
d = s.run(tmax)
ax = fig.add_subplot(grid[1, 1])
ax.set_title(fancy_names[code])
ax.set_xlabel('Time (s)')
ax.set_yticklabels([])
ax.set_xlim(0, tmax)
ax.set_ylim(-1.02, 1.02)
mp.cumulative_current(d, currents, ax, colors=colours, normalise=True)

#
# Bottom row: O'Hara models
#
# O'Hara et al. 2011
code = 'ohara'
model = models[code]
currents, colours = current_variables(model, True)
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
d = s.run(tmax)
ax = fig.add_subplot(grid[2, 0])
ax.set_title(fancy_names[code])
ax.set_xlabel('Time (s)')
ax.set_ylabel('Relative contribution')
ax.set_xlim(0, tmax)
ax.set_ylim(-1.02, 1.02)
mp.cumulative_current(d, currents, ax, colors=colours, normalise=True)

# O'Hara et al. 2017 CiPA update
code = 'cipa'
model = models[code]
currents, colours = current_variables(model, True)
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
d = s.run(tmax)
ax = fig.add_subplot(grid[2, 1])
ax.set_title(fancy_names[code])
ax.set_xlabel('Time (s)')
ax.set_yticklabels([])
ax.set_xlim(0, tmax)
ax.set_ylim(-1.02, 1.02)
mp.cumulative_current(d, currents, ax, colors=colours, normalise=True)

# Tomek et al. 2020
code = 'tomek'
model = models[code]
currents, colours = current_variables(model, True)
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
d = s.run(tmax)
ax = fig.add_subplot(grid[2, 2])
ax.set_title(fancy_names[code])
ax.set_xlabel('Time (s)')
ax.set_yticklabels([])
ax.set_xlim(0, tmax)
ax.set_ylim(-1.02, 1.02)
mp.cumulative_current(d, currents, ax, colors=colours, normalise=True)


#
# Legend
#
ax = fig.add_subplot(grid[1, 2])
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.set_frame_on(False)
lines = []
for current, i in current_colours.items():
    lines.append(matplotlib.lines.Line2D([0], [0], color=cmap(i), lw=5))
labels = [x.replace('_', '') for x in current_colours]
ax.legend(lines, labels, loc=(0.05, 0.05), ncol=2)
#ax.legend(lines, labels, loc=(0.05, -0.7), ncol=1)


# Show / store
plt.savefig('ventricular.png')
plt.savefig('ventricular.pdf')
print('Done')
