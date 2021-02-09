#!/usr/bin/env python3
#
# Relative contributions of the major ionic currents in human atrial models.
# all models.
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
    'I_Kur',    # 0
    'I_to',     # 1
    'I_Cl,B',   # 2
    'I_CaL',    # 3
    'I_NaCa',   # 4
    'I_Kr',     # 5
    'I_Ks',     # 6
    'I_K1',     # 7
    'I_NaK',    # 8
    'I_Ca,P',   # 9
    'I_Ca,B',   # 10
    'I_Na,B',   # 11
    'I_KACh',   # 12
    'I_ClCa',   # 13
    'I_Kp',     # 14
    'I_f',      # 15
    'I_Na',     # 16
    'I_NaL',    # 17
]

# Human atrial models
model_names = {
    'courtemanche': 'courtemanche-1998.mmt',
    'grandi': 'grandi-2011.mmt',
    'koivumaki': 'koivumaki-2011.mmt',
    'maleckar': 'maleckar-2008.mmt',
    'ni': 'ni-2017.mmt',
    'nygren': 'nygren-1998.mmt',
    'voigt': 'voigt-heijman-2013.mmt',
}

fancy_names = {
    'courtemanche': 'Courtemanche et al., 1998',
    'grandi': 'Grandi-Pandit-Voigt et al., 2011',
    'koivumaki': 'Koivumaki et al., 2011',
    'maleckar': 'Maleckar et al., 2008',
    'ni': 'Ni et al., 2017',
    'nygren': 'Nygren et al., 1998',
    'voigt': 'Voigt-Heijman et al., 2013',
}


def current_variables(model, colours=False):
    """ Returns an ordered list of transmembrane current variable names. """
    name = model.name().lower()
    if 'courtemanche-1998' in name:
        currents = {
            'ikur.i_Kur': 0,
            'ito.i_to': 1,
            'ical.i_Ca_L': 3,
            'inaca.i_NaCa': 4,
            'ikr.i_Kr': 5,
            'iks.i_Ks': 6,
            'ik1.i_K1': 7,
            'inak.i_NaK': 8,
            'ipca.i_PCa': 9,
            'ib.i_B_Ca': 10,
            'ib.i_B_Na': 11,
            'ina.i_Na': 16,
        }
    elif 'grandi-2011' in name:
        currents = {
            'ikur.IKur': 0,
            'ito.Ito': 1,
            'iclb.IClB': 2,
            'ical.ICaL': 3,
            'inaca.INaCa': 4,
            'ikr.IKr': 5,
            'iks.IKs': 6,
            'ik1.IK1': 7,
            'inak.INaK': 8,
            'ipca.IpCa': 9,
            'icab.ICaB': 10,
            'inab.INaB': 11,
            'iclca.IClCa': 13,
            'ikp.IKp': 14,
            'ina.INa': 16,
            'inal.INaL': 17,
        }
    elif 'koivumaki' in name:
        currents = {
            'ikur.IKur': 0,
            'it.It': 1,
            'ical.ICaL': 3,
            'inaca.INaCa': 4,
            'ikr.IKr': 5,
            'iks.IKs': 6,
            'ik1.IK1': 7,
            'inak.INaK': 8,
            'icap.ICaP': 9,
            'icab.ICab': 10,
            'inab.INab': 11,
            'if.If': 15,
            'ina.INa': 16,
        }
    elif 'maleckar-' in name:
        currents = {
            'ikur.i_Kur': 0,
            'it.i_t': 1,
            'ical.i_Ca_L': 3,
            'inaca.i_NaCa': 4,
            'ikr.i_Kr': 5,
            'iks.i_Ks': 6,
            'ik1.i_K1': 7,
            'inak.i_NaK': 8,
            'icap.i_CaP': 9,
            'ib.i_B_Ca': 10,
            'ib.i_B_Na': 11,
            'ikach.i_KACh': 12,
            'ina.i_Na': 16,
        }
    elif 'ni-' in name:
        currents = {
            'ikur.IKur': 0,
            'ito.Ito': 1,
            'ical.ICaL': 3,
            'inaca.INaCa': 4,
            'ikr.IKr': 5,
            'iks.IKs': 6,
            'ik1.IK1': 7,
            'inak.INaK': 8,
            'icap.ICap': 9,
            'ibca.IbCa': 10,
            'ibna.IbNa': 11,
            'ina.INa': 16,
        }
    elif 'nygren' in name:
        currents = {
            'isus.i_sus': 0,
            'it.i_t': 1,
            'ical.iCaL': 3,
            'inaca.i_NaCa': 4,
            'ikr.i_Kr': 5,
            'iks.i_Ks': 6,
            'ik1.i_K1': 7,
            'inak.i_NaK': 8,
            'icap.i_CaP': 9,
            'ib.i_B_Ca': 10,
            'ib.i_B_Na': 11,
            'ina.i_Na': 16,
        }
    elif 'voigt' in name:
        currents = {
            'ikur.IKur': 0,
            'ito.Ito': 1,
            'iclb.IClB': 2,
            'ical.ICaL': 3,
            'inaca.INaCa': 4,
            'ikr.IKr': 5,
            'iks.IKs': 6,
            'ik1.IK1': 7,
            'inak.INaK': 8,
            'ipca.IpCa': 9,
            'icab.ICaB': 10,
            'inab.INaB': 11,
            'ikach.IKACh': 12,
            'iclca.IClCa': 13,
            'ikp.IKp': 14,
            'ina.INa': 16,
            'inal.INaL': 17,
        }
    else:
        # Get all currents
        def rec(parent, currents=set()):
            for var in parent.refs_to():
                if 'tot' in var.name():
                    rec(var, currents)
                else:
                    currents.add(var)
            return currents

        currents = rec(model.labelx('cellular_current'))
        currents = [x.qname() for x in currents]
        currents.sort()
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
    if 'koiv' in name:
        pre_pace = False
    model = myokit.load_model(os.path.join('models', fname))
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
# Top row: Nygren models
#
# Nygren 1998
code = 'nygren'
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

# Koivumaki 2011
code = 'koivumaki'
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
# Bottom row: Grandi models
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
#ax.legend(lines, labels, loc=(0.05, 0.05), ncol=2)
ax.legend(lines, labels, loc=(0.05, -0.7), ncol=1)


# Show / store
plt.savefig('atrial.png')
plt.savefig('atrial.pdf')
print('Done')
