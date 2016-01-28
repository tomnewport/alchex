'''
Alchembed Methods
'''
from alchex.gromacs_interface import GromacsEditableTOPFile, SimulationContainer, GromacsWrapper, GromacsMDPFile
from alchex.config import default_configuration
from copy import deepcopy

config = default_configuration()

def alchembed(simulation_container, input_gro, input_top, mdp_file, name="alchembed", parameters={}):
    mdp_file = deepcopy(mdp_file)
    for k, v in parameters.items():
        mdp_file.attrs[k] = v
    alch_mdp.to_file(s.resolve_path(name + ".mdp"))
    simulation_container.gromacs.grompp(
        kwargs={
        "-f": name+".mdp", 
        "-p": input_top,
        "-c": input_gro,
        "-o": name+".tpr",
        "-maxwarn":"1"})
    return simulation_container.gromacs.mdrun(kwargs={"-deffnm":name})


class AlchembedAccumulator(object):
    def __init__(self, 
            alchex_config,
            input_gro,
            input_top,
            simulation_container,
            parameters={},
            working_directory="/alchembed-accumulator"):
        '''
        
        '''

        self.alchex_config = alchex_config
        self.simulations = simulation_container
        simulation_container.makedirs(working_directory)
        simulation_container.cd(working_directory)
        self.working_directory = working_directory

        self.input_gro = input_gro
        self.input_top = input_top

        self.simulations.copy(input_gro, "start.gro")

        alch_mdp = GromacsMDPFile()
        alch_mdp.attrs = deepcopy(self.alched_config.grompp_parameters["alchembed"])

        for k, v in parameters.items():
            alch_mdp.attrs[k] = v

        alch_mdp.to_file(s.resolve_path("alchembed.mdp"))






s = SimulationContainer("alchembed-test", GromacsWrapper("/sbcb/packages/opt/Linux_x86_64/gromacs/5.1/bin/gmx_sse_d"))

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

alch_mdp = GromacsMDPFile()
alch_mdp.attrs = config.grompp_parameters["alchembed"]


print alchembed(s, "em/em.gro", "cox1-cg.top", alch_mdp, name="messing_around", parameters={"couple-moltype":"PROTEIN"})[1]
