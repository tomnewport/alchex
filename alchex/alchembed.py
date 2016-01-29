'''
Alchembed Methods
'''
from alchex.gromacs_interface import GromacsEditableTOPFile, SimulationContainer, GromacsWrapper, GromacsMDPFile, GromacsEditableTOPFile
from alchex.config import default_configuration
from alchex.workarounds import WAEditableGrofile
from copy import deepcopy
from os import path
from collections import OrderedDict

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



s = SimulationContainer("alchembed-combine-test", GromacsWrapper("/sbcb/packages/opt/Linux_x86_64/gromacs/5.1/bin/gmx_sse"))

def alchembed_combine(
    simulation_container,
    alchex_config,
    system_1_structure,
    system_1_topology,
    system_2_structure,
    system_2_topology,
    parameters={},
    working_directory="/alchembed-combine"
    ):
    #simulation_container.delete_folder(working_directory)
    simulation_container.makedirs(working_directory)

    simulation_container.cd(working_directory)

    alch_mdp = GromacsMDPFile()
    alch_mdp.attrs = deepcopy(config.grompp_parameters["alchembed"])
    alch_mdp.attrs["couple-moltype"] = "ALCHEX2"
    alch_mdp.to_file(simulation_container.resolve_path("alch.mdp"))

    prep_mdp = GromacsMDPFile()
    prep_mdp.attrs = deepcopy(config.grompp_parameters["alchembed"])
    del prep_mdp.attrs["couple-moltype"]
    del prep_mdp.attrs["couple-lambda0"]
    del prep_mdp.attrs["couple-lambda1"]
    prep_mdp.to_file(simulation_container.resolve_path("preprocess.mdp"))

    status, message = simulation_container.gromacs.grompp(kwargs={
        "-c" : simulation_container.resolve_path(system_1_structure),
        "-p" : simulation_container.resolve_path(system_1_topology),
        "-f" : "preprocess.mdp",
        "-pp": simulation_container.resolve_path("alch1.top"),
        "-maxwarn" : 1
    })
    if status != 0:
        return False

    status, message = simulation_container.gromacs.grompp(kwargs={
        "-c" : simulation_container.resolve_path(system_2_structure),
        "-p" : simulation_container.resolve_path(system_2_topology),
        "-f" : "preprocess.mdp",
        "-pp": simulation_container.resolve_path("alch2.top"),
        "-maxwarn" : 1
    })
    if status != 0:
        return False

    simulation_container.copy_file(
        system_1_structure, 
        working_directory, 
        rename="alch1.gro"
        )
    simulation_container.copy_file(
        system_2_structure, 
        working_directory, 
        rename="alch2.gro"
        )

    system1_top = GromacsEditableTOPFile()
    system1_top.from_file(
        simulation_container.resolve_path("alch1.top")
        )

    system1_unified_top = GromacsEditableTOPFile()
    system1_unified_top.add_moltype(
        system1_top.unify_moltype("ALCHEX1")
        )
    system1_unified_top.to_file(
        simulation_container.resolve_path("alch1.itp")
        )

    system2_top = GromacsEditableTOPFile()
    system2_top.from_file(
        simulation_container.resolve_path("alch2.top")
        )

    system2_unified_top = GromacsEditableTOPFile()
    system2_unified_top.add_moltype(
        system2_top.unify_moltype("ALCHEX2")
        )
    system2_unified_top.to_file(
        simulation_container.resolve_path("alch2.itp")
        )

    with open(
        simulation_container.resolve_path("alch1.top"),
        "r") as read_file_handle:
        with open(
            simulation_container.resolve_path("alch.itp"),
            "w"
            ) as write_file_handle:
            for line in read_file_handle:
                if "system" in line and "[" in line:
                    break
                else:
                    write_file_handle.write(line)



    alch_top = GromacsEditableTOPFile()
    alch_top.includes.append("alch.itp")
    alch_top.includes.append("alch1.itp")
    alch_top.includes.append("alch2.itp")

    alch_top.tables["system"] = [
        OrderedDict([("system_name", "Alchex Alchembed")])
    ]
    alch_top.tables["molecules"] = [
        OrderedDict([("moltype", "ALCHEX1"),
                     ("count", "1")]),
        OrderedDict([("moltype", "ALCHEX2"),
                     ("count", "1")])
    ]

    alch_top.to_file(
        simulation_container.resolve_path("alch.top")
        )

    system1_structure = WAEditableGrofile()
    system1_structure.from_file(
        simulation_container.resolve_path("alch1.gro")
        )
    system2_structure = WAEditableGrofile()
    system2_structure.from_file(
        simulation_container.resolve_path("alch2.gro")
        )

    system1_structure.combine(system2_structure)
    system1_structure.to_file(
        simulation_container.resolve_path("alch.gro")
        )

    print simulation_container.gromacs.grompp(kwargs={
        "-c" : "alch.gro",
        "-p" : "alch.top",
        "-f" : "alch.mdp",
        "-o" : "alch-run.tpr",
        "-maxwarn" : "1"
    })[1]

    print simulation_container.gromacs.mdrun(kwargs={
        "-deffnm" : "alch-run",
        "-rdd" : 1
        })[1]
    # assuming this has worked, we now need to
    # 1. Renumber residues to match the original
    #    .gro files (assuming they dovetail nicely)
    # 2. Sort by residue id
    # 3. Build a new topology file

print alchembed_combine(
    s,
    config,
    "/inputs/alchex_component.0.em.gro",
    "/inputs/alchex_component.0.top",
    "/inputs/alchex_component.1.em.gro",
    "/inputs/alchex_component.1.top"
    )

'''
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
'''