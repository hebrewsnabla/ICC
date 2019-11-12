"""
version 1.1
author: Felix Plasser
usage:Package for drawing energy schemes using the python matplotlib library.
license: The code may be freely used and copied.
    About distribution or modifications please contact the author:
    felix.plasser@gmx.at
"""

import pylab

class energy_scheme:
    def __init__(self):
        self.energy_levels = {}
    
    def add_energy_level(self, name, *en_lev_args, **extra_plot_args):
        self.energy_levels[name] = energy_level(name, *en_lev_args, **extra_plot_args)

    def plot_energy_levels(self, **extra_plot_args):
        """
        Plot all the energy levels.
        """
        pylab.ylabel('energy (eV)')
        # load default arguments and assign extra arguments
        used_plot_args = {'linestyle':'-', 'color':'black', 'linewidth':2}
        for arg, val in extra_plot_args.iteritems():
            used_plot_args[arg] = val
        try:
            label = used_plot_args.pop('label')        # extract the label
        except KeyError:
            label = ''
            
        for name, en_lev in self.energy_levels.iteritems():
            en_lev.plot(**used_plot_args)
        pylab.plot([0],[0], label=label, **used_plot_args)
        

    def plot_connection(self, en_lev1_name, en_lev2_name, **extra_plot_args):
        """
        Plot the connection line between two energy levels.
        en_lev1 is to the left of en_lev2
        """
        en_lev1 = self.energy_levels[en_lev1_name]
        en_lev2 = self.energy_levels[en_lev2_name]

        # load default arguments and assign extra arguments
        used_plot_args = {'linestyle':'--', 'color':'grey'}
        for arg, val in extra_plot_args.iteritems():
            used_plot_args[arg] = val
            
        pylab.plot([en_lev1.pos + .25, en_lev2.pos -.25], [en_lev1.energy, en_lev2.energy], **used_plot_args)

    def level_text(self, level, text, hdisp=0., vdisp=0.):
        """
        Write <text> to level <level> with horizontal/vertical displacement <hdisp>/<vdisp>.
        """
        pylab.text(self.energy_levels[level].pos+hdisp, self.energy_levels[level].energy+vdisp, text, horizontalalignment='center')

class energy_level:
    def __init__(self, name, energy, pos, text=None, **extra_plot_args):
        """
        <name> internal recognition
        <energy> in unit chosen
        <pos> horizontal position
        <text> text to be written in a standard way (see also energy_scheme.level_txt)
        """
        self.name = name
        self.energy = energy
        self.pos = pos
        self.text = text
        self.extra_plot_args = extra_plot_args

    def plot(self, **extra_plot_args):
        # load default arguments and assign extra arguments
        used_plot_args = {'linestyle':'-', 'color':'black', 'linewidth':2}
        for arg, val in extra_plot_args.iteritems(): # arguments passed to plot
            used_plot_args[arg] = val
        for arg, val in self.extra_plot_args.iteritems(): # arguments of the energy level
            used_plot_args[arg] = val
            
        pylab.plot([self.pos - .25, self.pos + .25], [self.energy, self.energy], **used_plot_args)
        if not self.text == None:
            pylab.text(self.pos, self.energy+.01, self.text, horizontalalignment='center')
