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

s = SimulationContainer("alchembed-test", GromacsWrapper("gmx"))

em_mdp = GromacsMDPFile()
em_mdp.attrs = config.grompp_parameters["alchembed"]
em_mdp.to_file(s.resolve_path("alchembed.mdp"))


s.makedirs("em")

s.gromacs.grompp(
	kwargs={
	"-f":"em-cg.mdp", 
	"-p":"cox1-cg.top",
	"-c":"cox1-cg.pdb",
	"-o":"em/em.tpr",
	"-maxwarn":"1"})

s.cd("em")

print s.gromacs.mdrun(kwargs={"-deffnm":"em"})[1]

s.cd("../")

s.makedirs("alchembed")

s.gromacs.grompp(
	kwargs={
	"-f":"alchembed.mdp", 
	"-p":"cox1-cg.top",
	"-c":"em/em.gro",
	"-o":"alchembed/alchembed.tpr",
	"-maxwarn":"1"})

s.cd("alchembed")

print("--------------")

print s.gromacs.mdrun(kwargs={"-deffnm":"alchembed"})[1]
