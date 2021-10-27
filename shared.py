#!/usr/bin/env python3
#
# Shared code for model current "relative contribution" graphs.
#
import myokit
import numpy as np


current_colours = {
    'I_Kr': 0,
    'I_Ks': 1,
    'I_to': 2,
    'I_Kb': 3,
    'I_f': 8,
    'I_Kur': 9,
    'I_K1': 4,
    'I_NaK': 5,
    'I_Na': 16,
    'I_NaL': 17,
    'I_CaL': 10,
    'I_CaT': 11,
    'I_NaCa': 12,
    'I_Na,B': 14,
    'I_Ca,B': 15,
    'I_ClCa': 6,
    'I_Cl,B': 7,
    'I_Ca,P': 13,
    'I_K,ACh': 18,
    'I_K,ATP': 19,
}

current_names = {
    'I_Kr': 'IKr',
    'I_Ks': 'IKs',
    'I_to': 'Ito (+Isus)',
    'I_Kb': 'IKb (IbK)',
    'I_f': 'If',
    'I_Kur': 'IKur',
    'I_K1': 'IK1',
    'I_NaK': 'INaK',
    'I_Na': 'INa',
    'I_NaL': 'INaL',
    'I_CaL': 'ICaL',
    'I_CaT': 'ICaT',
    'I_NaCa': 'INaCa (INCX)',
    'I_Na,B': 'INa,B',
    'I_Ca,B': 'ICa,B',
    'I_ClCa': 'IClCa',
    'I_Cl,B': 'ICl,B',
    'I_Ca,P': 'ICa,P (IpCa)',
    'I_K,ACh': 'IK,ACh',
    'I_K,ATP': 'IK,ATP',
}


def prepare_model(model, protocol, currents, pre_pace=True):
    """
    Prepares a model by setting the desired units, adding a voltage-clamp
    switch, and pre-pacing.

    The following variables / labels guaranteed to exist after this:

    - A time variable in ms
    - ``membrane_potential`` in mV
    - ``membrane_capacitance`` (units unchanged)
    - All variables in ``currents``, in A/F

    Pre-pacing can be disabled by setting ``pre_pace=False``.
    """
    # Get model variables
    t = model.timex()
    v = model.labelx('membrane_potential')

    # Check variables have units
    for var in (v, t):
        if var.unit() is None:
            raise ValueError(
                'No unit set for ' + str(var) + ' in ' + str(model))

    # Optionally get capacitance
    helpers = []
    C = model.label('membrane_capacitance')
    if C is not None:
        if C.unit() is None:
            raise ValueError(
                'No unit set for ' + str(C) + ' in ' + str(model))
        helpers.append(C.rhs())

    # Convert variable units
    i_unit = 'A/F'
    v.convert_unit('mV')
    for qname in currents:
        var = model.get(qname)
        if var.unit() is None:
            raise ValueError(
                'No unit set for ' + str(var) + ' in ' + str(model))
        var.convert_unit(i_unit, helpers=helpers)
    t.convert_unit('ms')

    # Pre-pace
    if pre_pace and not 'koiv' in model.name():
        print('Pre-pacing: ' + model.name())
        model.set_state(limit_cycle(model, protocol))
        print(model.format_state(model.state()))
    else:
        print('NOT Pre-pacing: ' + model.name())


def demote(var):
    """
    Changes a state variable to a non-state variable.

    1. Replaces any references to its derivative with an inlined equation
    2. Sets the variable's RHS to its state value
    3. Demotes the variable.
    """
    if not var.is_state():
        return

    # Replace references to dot(x) by inlined equation
    subst = {var.lhs(): var.rhs()}
    for ref in list(var.refs_by()):
        ref.set_rhs(ref.rhs().clone(subst=subst))

    # Demote
    var.set_rhs(var.state_value())
    var.demote()


def limit_cycle(model, protocol, cl=None, rel_tol=1e-5, max_beats=20000,
                max_period=10, path=None):
    """
    Pre-paces a model to periodic orbit ("steady state").

    Arguments
    ``model``
    ``protocol``
    ``cl``
    ``rel_tol``
    ``max_beats``
    ``max_period``
    ``path``
    """

    # Create simulation
    s = myokit.Simulation(model, protocol)
    s.set_tolerance(1e-9, 1e-9)
    if cl is None:
        cl = protocol.characteristic_time()

    # Load steady-state from file, if given
    try:
        loaded = myokit.load_state(path)
        s.set_state(loaded)
    except Exception:
        loaded = None
        pass
    else:
        print('Loaded state from ' + str(path))

    # Get scale of each state
    states = list(model.states())
    d = s.run(cl, log=myokit.LOG_STATE)
    x = np.array([d[var] for var in states])
    scale = np.max(x, axis=1) - np.min(x, axis=1)
    scale[scale==0] = 1

    # Check if already at steady-state
    d = s.run(2 * cl, log_interval=cl, log=myokit.LOG_STATE)
    x = np.array([d[var] for var in states]).T
    dx = np.abs(x[0] - x[1]) / scale
    if np.max(dx) < rel_tol:
        return s.state() if loaded is None else loaded

    beats = 0
    period = 0
    duration = max_period * cl
    while beats < max_beats:

        # Run and capture a number of beats
        d = s.run(duration, log_interval=cl, log=myokit.LOG_STATE)
        beats += max_period
        x = np.array([d[var] for var in states]).T

        # Check to see if any of the beats are a repeat of the first beat
        period = 0
        for i in range(1, max_period):
            dx = np.abs(x[0] - x[i]) / scale
            if np.max(dx) < rel_tol:
                period = i
                break

        # Repeat found! Check if there's a second repeat
        if period > 0:
            # Simulate more beats if needed
            if 2 * period > max_period:
                d = s.run(duration, log_interval=cl, log=myokit.LOG_STATE)
                beats += max_period
                x = np.concatenate((x, np.array([d[var] for var in states]).T))

            dx = np.abs(x[period] - x[2 * period]) / scale
            if np.max(dx) < rel_tol:
                print('Terminating after ' + str(beats) + ' beats')
                break
            else:
                period = 0

    # Save state to file
    if path is not None:
        print('Saving final state to ' + str(path))
        myokit.save_state(path, s.state())

    if period > 1:
        print('WARNING: Detected alternans with period ' + str(period) + '.')
    elif period == 0:
        print('WARNING: Terminating after maximum number of beats.')
        dx = np.abs(x[0] - x[i]) / scale
        print('Final dx: ' + str(np.max(dx)))

    return s.state()


def guess_currents(model):
    """ Guess all transmembrane currents in a given ``model``. """
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
    return currents

