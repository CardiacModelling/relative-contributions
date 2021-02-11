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
current_names = [
    'I_ClCa',    # 0
    'I_NaCa',    # 1
    'I_to',      # 2
    'I_Kb',      # 3
    'I_Kr',      # 4
    'I_Ks',      # 5
    'I_Cl,B',    # 6
    'I_K,ATP',   # 7
    'I_K1',      # 8
    'I_NaK',     # 9
    'I_CaL',     # 10
    'I_NaL',     # 11
    'I_Na,B',    # 12
    'I_Ca,B',    # 13
    'I_Na',      # 14
    'I_Ca,P',    # 15
]

# Human atrial models
model_names = {
    'priebe': 'priebe-1998.mmt',
    #'iyer': 'iyer-2004.mmt',
    'grandi': 'grandi-2010.mmt',
    #'tnnp1': 'tentusscher-2004.mmt',
    #'tnnp2': 'tentusscher-2006.mmt',
    #'ohara': 'ohara-2011.mmt',
    #'cipa': 'ohara-cipa-v1-2017.mmt',
    'tomek': 'tomek-2020-chloride-epi.mmt',
}

fancy_names = {
    'priebe': 'Priebe & Beuckelmann, 1998',
    'iyer': 'Iyer et al., 2004 (epi)',
    'grandi': 'Grandi et al., 2010 (epi)',
    'tnnp1': 'Ten Tusscher et al., 2004 (epi)',
    'tnnp2': 'Ten Tusscher & Panfilov 2006 (epi)',
    'ohara': 'O\'Hara et al., 2011 (epi)',
    'cipa': 'O\'Hara et al., 2017 CiPA (epi)',
    'tomek': 'Tomek et al., 2020 (epi)',
}


def current_variables(model, colours=False):
    """ Returns an ordered list of transmembrane current variable names. """
    name = model.name().lower()
    if 'priebe' in name:
        currents = {
            'inaca.i_NaCa': 1,
            'ito.i_to': 2,
            'ikr.i_Kr': 4,
            'iks.i_Ks': 5,
            'ik1.i_K1': 8,
            'inak.i_NaK': 9,
            'ica.i_Ca': 10,
            'inab.i_b_Na': 12,
            'icab.i_b_Ca': 13,
            'ina.i_Na': 14,
        }
    elif 'grandi' in name:
        currents = {
            'iclca.iclca': 0,
            'incx.I_ncx': 1,
            'ito.ito': 2,
            'ikp.I_kp': 3,
            'ikr.I_kr': 4,
            'iks.I_ks': 5,
            'iclb.IClb': 6,
            'ik1.I_k1': 8,
            'inak.I_nak': 9,
            'ical.I_CaL': 10,
            'inab.I_nabk': 12,
            'icabk.I_cabk': 13,
            'ina.I_Na': 14,
            'ipca.I_pca': 15,
        }
    elif 'torord' in name:
        currents = {
            'ICl.IClCa': 0,
            'INaCa.INaCa': 1,
            'Ito.Ito': 2,
            'IKb.IKb': 3,
            'IKr.IKr': 4,
            'IKs.IKs': 5,
            'ICl.IClb': 6,
            'I_katp.I_katp': 7,
            'IK1.IK1': 8,
            'INaK.INaK': 9,
            'ICaL.ICaL': 10,
            'INaL.INaL': 11,
            'INab.INab': 12,
            'ICab.ICab': 13,
            'INa.INa': 14,
            # TODO
            'IpCa.IpCa': 15,


        }

    else:
        currents = shared.guess_currents(model)
        print('\n'.join(currents))
        print(len(currents))
        raise NotImplementedError('Unknown model: ' + model.name())

    if colours:
        colours = [cmap(x) for x in currents.values()]
        currents = list(currents.keys())
        return currents, colours
    return list(currents.keys())


# Create protocol
cl = 1000
protocol = myokit.pacing.blocktrain(cl, duration=0.5, offset=50)

# Load and prepare models
models = {}
for name, fname in model_names.items():
    pre_pace = True
    if 'priebe' in name:
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
# Priebe 1998
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
'''
text(ax, 0.950, 0.870, 'INaK')
text(ax, 0.950, 0.640, 'IK1')
text(ax, 0.950, 0.420, 'INaCa')
text(ax, 0.950, 0.230, 'ICa,B')
text(ax, 0.950, 0.060, 'INa,B')
text(ax, 0.210, 0.530, 'IKur')
text(ax, 0.220, 0.450, 'ICaL')
text(ax, 0.300, 0.880, 'Ito')
ax.arrow(0.20, 0.88, -0.11, -0.08, transform=ax.transAxes, zorder=999)
text(ax, 0.400, 0.06, 'INa')
ax.arrow(0.28, 0.06, -0.09, -0.03, transform=ax.transAxes, zorder=999)
'''
'''
# Maleckar 2009
code = 'maleckar'
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
text(ax, 0.950, 0.870, 'INaK')
text(ax, 0.950, 0.620, 'IK1')
text(ax, 0.950, 0.445, 'INaCa')
text(ax, 0.950, 0.320, 'ICa,B')
text(ax, 0.950, 0.125, 'INa,B')
text(ax, 0.210, 0.530, 'IKur')
text(ax, 0.220, 0.450, 'ICaL')
text(ax, 0.300, 0.900, 'Ito')
ax.arrow(0.20, 0.90, -0.11, -0.08, transform=ax.transAxes, zorder=999)
text(ax, 0.380, 0.06, 'INa')
ax.arrow(0.26, 0.06, -0.09, -0.03, transform=ax.transAxes, zorder=999)
'''

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
'''
text(ax, 0.950, 0.880, 'INaK')
text(ax, 0.950, 0.640, 'IK1')
text(ax, 0.950, 0.440, 'INaCa')
text(ax, 0.950, 0.260, 'ICa,B')
text(ax, 0.950, 0.080, 'INa,B')
text(ax, 0.210, 0.530, 'IKur')
text(ax, 0.215, 0.450, 'ICaL')
text(ax, 0.300, 0.880, 'Ito')
ax.arrow(0.20, 0.88, -0.11, -0.08, transform=ax.transAxes, zorder=999)
text(ax, 0.360, 0.053, 'INa')
ax.arrow(0.24, 0.053, -0.09, -0.03, transform=ax.transAxes, zorder=999)
'''
'''
#
# Middle row: Courtemanche models
#
# Courtemanche 1998
code = 'courtemanche'
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
text(ax, 0.950, 0.950, 'ICa,P')
text(ax, 0.950, 0.838, 'INaK')
text(ax, 0.950, 0.630, 'IK1')
text(ax, 0.950, 0.440, 'INaCa')
text(ax, 0.950, 0.236, 'ICa,B')
text(ax, 0.950, 0.055, 'INa,B')
text(ax, 0.210, 0.530, 'IKur')
text(ax, 0.210, 0.450, 'ICaL')
text(ax, 0.270, 0.840, 'Ito')
ax.arrow(0.172, 0.840, -0.088, -0.120, transform=ax.transAxes, zorder=999)
text(ax, 0.440, 0.660, 'IKs')
ax.arrow(0.332, 0.660, -0.093, -0.060, transform=ax.transAxes, zorder=999)
text(ax, 0.420, 0.580, 'IKr')
ax.arrow(0.320, 0.580, -0.086, -0.030, transform=ax.transAxes, zorder=999)

# Ni 2017
code = 'ni'
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
text(ax, 0.950, 0.950, 'ICa,P')
text(ax, 0.950, 0.838, 'INaK')
text(ax, 0.950, 0.635, 'IK1')
text(ax, 0.950, 0.420, 'INaCa')
text(ax, 0.950, 0.230, 'ICa,B')
text(ax, 0.950, 0.055, 'INa,B')
text(ax, 0.215, 0.550, 'IKur')
text(ax, 0.225, 0.440, 'ICaL')
text(ax, 0.280, 0.850, 'Ito')
ax.arrow(0.182, 0.850, -0.088, -0.120, transform=ax.transAxes, zorder=999)
text(ax, 0.500, 0.720, 'IKs')
ax.arrow(0.392, 0.720, -0.093, -0.060, transform=ax.transAxes, zorder=999)
text(ax, 0.480, 0.640, 'IKr')
ax.arrow(0.380, 0.640, -0.086, -0.030, transform=ax.transAxes, zorder=999)

#
# Bottom row: O'Hara models
#
# Grandi 2011
code = 'grandi'
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
text(ax, 0.950, 0.865, 'INaK')
text(ax, 0.950, 0.625, 'IK1')
text(ax, 0.950, 0.434, 'ICl,B')
text(ax, 0.950, 0.316, 'INaCa')
text(ax, 0.950, 0.190, 'ICa,B')
text(ax, 0.950, 0.065, 'INa,B')
text(ax, 0.207, 0.542, 'IKur')
text(ax, 0.230, 0.440, 'ICaL')
text(ax, 0.300, 0.900, 'Ito')
ax.arrow(0.202, 0.900, -0.098, -0.120, transform=ax.transAxes, zorder=999)

# Voigt-Heijman 2013
code = 'voigt'
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
text(ax, 0.950, 0.865, 'INaK')
text(ax, 0.950, 0.625, 'IK1')
text(ax, 0.950, 0.423, 'ICl,B')
text(ax, 0.950, 0.306, 'INaCa')
text(ax, 0.950, 0.189, 'ICa,B')
text(ax, 0.950, 0.061, 'INa,B')
text(ax, 0.207, 0.542, 'IKur')
text(ax, 0.230, 0.440, 'ICaL')
text(ax, 0.300, 0.900, 'Ito')
ax.arrow(0.202, 0.900, -0.098, -0.120, transform=ax.transAxes, zorder=999)
text(ax, 0.510, 0.640, 'IKr')
ax.arrow(0.410, 0.640, -0.086, -0.030, transform=ax.transAxes, zorder=999)
'''

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
'''
text(ax, 0.950, 0.865, 'INaK')
text(ax, 0.950, 0.625, 'IK1')
text(ax, 0.950, 0.423, 'ICl,B')
text(ax, 0.950, 0.306, 'INaCa')
text(ax, 0.950, 0.189, 'ICa,B')
text(ax, 0.950, 0.061, 'INa,B')
text(ax, 0.207, 0.542, 'IKur')
text(ax, 0.230, 0.440, 'ICaL')
text(ax, 0.300, 0.900, 'Ito')
ax.arrow(0.202, 0.900, -0.098, -0.120, transform=ax.transAxes, zorder=999)
text(ax, 0.510, 0.640, 'IKr')
ax.arrow(0.410, 0.640, -0.086, -0.030, transform=ax.transAxes, zorder=999)
'''

#
# Legend
#
ax = fig.add_subplot(grid[1, 2])
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.set_frame_on(False)
lines = []
for i, current in enumerate(current_names):
    lines.append(matplotlib.lines.Line2D([0], [0], color=cmap(i), lw=5))
labels = [x.replace('_', '') for x in current_names]
ax.legend(lines, labels, loc=(0.05, 0.05), ncol=2)
#ax.legend(lines, labels, loc=(0.05, -0.7), ncol=1)


# Show / store
plt.savefig('ventricular.png')
plt.savefig('ventricular.pdf')
print('Done')
