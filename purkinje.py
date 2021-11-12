#!/usr/bin/env python3
#
# Relative contributions of the major ionic currents in human atrial models.
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
current_colours = dict(shared.current_colours)
#del(current_colors['I_f'])
#del(current_colors['I_Kur'])
#del(current_colors['I_CaT'])

# Human atrial models
model_names = {
    'sampson': 'sampson-2010.mmt',
    'stewart': 'stewart-2009.mmt',
    'trovato': 'trovato-2020.mmt',
}

fancy_names = {
    'sampson': 'Sampson-Iyer et al., 2010',
    'stewart': 'Stewart et al., 2009',
    'trovato': 'Trovato et al., 2020',
}


def current_variables(model, colours=False):
    """ Returns an ordered list of transmembrane current variable names. """
    name = model.name().lower()
    if 'sampson' in name:
        currents = {
            'I_Kr': 'ikr.IKr',
            'I_Ks': 'iks.IKs',
            'I_to': 'ito.Ito_total',
            #'I_Kb': '',
            'I_f': 'ihcn.IHCN',
            #'I_Kur': '',
            'I_K1': 'ik1.IK1',
            'I_NaK': 'inak.INaK',
            'I_Na': 'nav15.INa',
            'I_NaL': 'nav11.INa1',
            'I_CaL': 'ical.ICaL_total',
            'I_CaT': 'icat.I',
            'I_NaCa': 'inaca.INaCa',
            #'I_Na,B': '',
            #'I_Ca,B': '',
            #'I_ClCa': '',
            #'I_Cl,B': '',
            'I_Ca,P': 'ipca.IpCa',
            #'I_K,ACh': '',
            #'I_K,ATP': '',
        }
    elif 'stewart' in name:
        currents = {
            'I_Kr': 'ikr.i_Kr',
            'I_Ks': 'iks.i_Ks',
            'I_to': 'ito.i_to_total',
            'I_Kb': 'ipk.i_p_K',
            'I_f': 'if.i_f_total',
            #'I_Kur': '',
            'I_K1': 'ik1.i_K1',
            'I_NaK': 'inak.i_NaK',
            'I_Na': 'ina.i_Na',
            #'I_NaL': '',
            'I_CaL': 'ical.i_CaL',
            #'I_CaT': '',
            'I_NaCa': 'inaca.i_NaCa',
            'I_Na,B': 'ibna.i_b_Na',
            'I_Ca,B': 'ibca.i_b_Ca',
            #'I_ClCa': '',
            #'I_Cl,B': '',
            'I_Ca,P': 'ipca.i_p_Ca',
            #'I_K,ACh': '',
            #'I_K,ATP': '',
        }
    elif 'trovato' in name:
        currents = {
            'I_Kr': 'IKr.IKr',
            'I_Ks': 'IKs.IKs',
            'I_to': 'Ito.Ito_total',
            #'I_Kb': '',
            'I_f': 'If.If',
            #'I_Kur': '',
            'I_K1': 'IK1.IK1',
            'I_NaK': 'INaK.INaK',
            'I_Na': 'INa.INa',
            'I_NaL': 'INaL.INaL',
            'I_CaL': 'ICaL.ICaL_total',
            'I_CaT': 'ICaT.ICaT',
            'I_NaCa': 'INaCa_i.INaCa_total',
            'I_Na,B': 'INab.INab',
            'I_Ca,B': 'ICab.ICab',
            #'I_ClCa': '',
            #'I_Cl,B': '',
            'I_Ca,P': 'IpCa.IpCa',
            #'I_K,ACh': '',
            #'I_K,ATP': '',
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
    model = myokit.load_model(os.path.join('models', 'purkinje', fname))
    if 'stewart' in name:
        c = model.get('ito')
        v = c.add_variable('i_to_total')
        v.set_unit(c.get('i_to').unit())
        v.set_rhs('ito.i_to + isus.i_sus')
        c = model.get('if')
        v = c.add_variable('i_f_total')
        v.set_unit(c.get('i_f_Na').unit())
        v.set_rhs('i_f_Na + i_f_K')

    elif 'sampson' in name:
        c = model.get('ito')
        v = c.add_variable('Ito_total')
        v.set_unit(c.get('Ito1').unit())
        v.set_rhs('ito.Ito1 + isus.Isus')
        c = model.get('ical')
        v = c.add_variable('ICaL_total')
        v.set_unit(c.get('ICa').unit())
        v.set_rhs('ICa + ICaK')

    elif 'trovato' in name:
        c = model.get('Ito')
        v = c.add_variable('Ito_total')
        v.set_unit(c.get('Ito').unit())
        v.set_rhs('Ito.Ito + Isus.Isus')
        c = model.get('ICaL')
        v = c.add_variable('ICaL_total')
        v.set_unit(c.get('ICaL').unit())
        v.set_rhs('ICaL + ICaK + ICaNa')
        c = model.get('INaCa_i')
        v = c.add_variable('INaCa_total')
        v.set_unit(c.get('INaCa_i').unit())
        v.set_rhs('INaCa_i.INaCa_i + INaCa_ss.INaCa_ss')

    pre_pace = True
    if 'stewart' in name:
        pre_pace = False

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

# Stewart 2009
code = 'stewart'
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

# Sampson 2010
code = 'sampson'
model = models[code]
currents, colours = current_variables(model, True)
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
d = s.run(tmax)
ax = fig.add_subplot(grid[1, 0])
ax.set_title(fancy_names[code])
ax.set_xlabel('Time (s)')
ax.set_yticklabels([])
ax.set_xlim(0, tmax)
ax.set_ylim(-1.02, 1.02)
mp.cumulative_current(d, currents, ax, colors=colours, normalise=True)

# Trovato 2020
code = 'trovato'
model = models[code]
currents, colours = current_variables(model, True)
s = myokit.Simulation(model, protocol)
s.set_tolerance(1e-8, 1e-8)
d = s.run(tmax)
ax = fig.add_subplot(grid[2, 0])
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
labels = [shared.current_names[x] for x in current_colours]
#ax.legend(lines, labels, loc=(0.05, 0.05), ncol=2)
ax.legend(lines, labels, loc=(0.05, -0.7), ncol=1)

# Show / store
plt.savefig('purkinje.png')
plt.savefig('purkinje.pdf')
print('Done')
