'''
Alchembed Methods
'''
from alchex.gromacs_interface import GromacsEditableTOPFile, SimulationContainer, GromacsWrapper, GromacsMDPFile
from alchex.config import default_configuration

config = default_configuration()

def alchembed(simulation_container, input_gro, input_top, parameters={}):
    pass



class AlchembedAccumulator(object):
    def __init__(self, 
        input_gro,
        input_top,
        simulation_container,
        working_directory="/alchembed-accumulator"):
        self.simulations = simulation_container
        simulation_container.cd(working_directory)

s = SimulationContainer("alchembed-test", GromacsWrapper("/sbcb/packages/opt/Linux_x86_64/gromacs/5.1/bin/gmx_sse"))

em_mdp = GromacsMDPFile()
em_mdp.attrs = config.grompp_parameters["alchembed"]
em_mdp.to_file(s.resolve_path("alchembed.mdp"))

s.gromacs.grompp(kwargs={"-f":"alchembed.mdp", "-p":"cox1-cg.top"})