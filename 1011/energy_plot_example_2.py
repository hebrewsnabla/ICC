"""
version 1.0
author: Felix Plasser
usage: Example for the usage of the energy_plot.py package.
"""

import pylab
import energy_plot
        

# Define classes that contain the connections between energy schemes.
class my_energy_scheme(energy_plot.energy_scheme):
    """
    Definition of connections between the structures on the energy scheme.
    """
    def plot_connections(self):
        self.plot_connection('Struc1', 'TS1')
        self.plot_connection('TS1', 'Struc2')
        self.plot_connection('Struc2', 'TS2')
        self.plot_connection('TS2', 'Struc3')

# Define the specific plot commands with all the parameters.
def plot_calc1():
    """
    Parameteres for plotting
    """
    pylab.title('calc1')
    en_sch = my_energy_scheme()
    
    en_sch.add_energy_level('Struc1', 3.2, -2) # name, energy, position
    en_sch.add_energy_level('TS1', 3.6, -1, '3.6') # optional label in the plot
    en_sch.add_energy_level('Struc2', 3, 0)
    en_sch.add_energy_level('FC', 3.9, 0, linestyle = '--', label='FC1') # optional matplotlib commands
    en_sch.add_energy_level('TS2', 3.2, 1)
    en_sch.add_energy_level('Struc3', 2.7, 2)

    en_sch.plot_connections()
    
    en_sch.plot_energy_levels(linestyle = '-', color='green', label='calc1')
    #pylab.axis('off')
    
def plot_calc2():
    pylab.title('calc2')
    en_sch = my_energy_scheme()
    
    en_sch.add_energy_level('Struc1', 3.1, -2) # name, energy, position
    en_sch.add_energy_level('TS1', 3.5, -1)
    en_sch.add_energy_level('Struc2', 3.2, 0)
    en_sch.add_energy_level('FC', 3.8, 0, linestyle = '--', label='FC2')
    en_sch.add_energy_level('TS2', 3.4, 1)
    en_sch.add_energy_level('Struc3', 2.6, 2)

    en_sch.plot_connections()
    
    en_sch.plot_energy_levels(linestyle = '-', color='red', label='calc2')

    en_sch.level_text('TS1', '3.5', vdisp=-.06) # add a label at a non-standard position


if __name__ == '__main__':
    # plot the energy scheme
    plot_calc1()
    plot_calc2()
    
    
    pylab.xticks([-2,0,2],['Struc1', 'Struc2','Struc3']) # write the names of the structures on the x-axis.
    pylab.title('calc1 and calc2')
    pylab.legend()
    pylab.axis([-2.5, 2.5, 2.5, 4.1]) # (x_low, x_high, y_low, y_high)
        
    

# show everything
    pylab.show()
