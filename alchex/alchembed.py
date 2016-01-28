'''
Alchembed Methods
'''
from alchex.gromacs_interface import GromacsEditableTOPFile



class AlchembedAccumulator(object):
    def __init__(self, 
        input_gro,
        input_top,
        simulation_container,
        working_directory="/alchembed-accumulator"):
    self.simulations = simulation_container
    simulation_container.cd(working_directory)
    

