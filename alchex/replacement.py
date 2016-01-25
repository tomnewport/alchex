from alchex.gromacs_interface import GromacsWrapper, SimulationContainer, GromacsTOPFile, GromacsMDPFile
from alchex.config import default_configuration
from alchex.geometry import PointCloud
from os import path

class ReplacementSpecification(object):
    def __init__(self, 
        selection=None, 
        parsimonious=None, 
        composition=None, 
        from_file=None):
        if from_file is not None:
            raise NotImplementedError()
        else:
            self.selection = selection
            self.parsimonious = parsimonious
            self.composition = composition
    def normalised_composition(self):
        sum_values = 1.0 * sum(self.composition.values())
        return {k:v/sum_values for k, v in self.composition.items()}


class ReplaceableEntity(object):
    def __init__(self, replacement_system, exchange_model, residue_ids, new_residue_ids):
        self.replacement_system = replacement_system
        self.residue_ids = residue_ids
        self.exchange_model = exchange_model




class ReplacementSystem(object):
    def __init__(self, 
        alchex_config=default_configuration(),
        input_structure_filename="conf.gro", 
        input_topology_filename="topology.top", 
        output_structure_filename="alchex.gro",
        output_topology_filename="alchex.top",
        em_parameters_filename=None,
        alchembed_parameters_filename=None,
        include_from=None,
        replacement_specifications=[],
        root_folder="./alchex-replacement"):
        self.input_structure_filename = input_structure_filename
        self.input_topology_filename = input_topology_filename
        self.output_structure_filename = output_structure_filename
        self.output_topology_filename = output_topology_filename
        self.replacement_specifications = replacement_specifications
        self.em_parameters_filename = em_parameters_filename
        self.alchembed_parameters_filename = alchembed_parameters_filename
        self.root_folder = root_folder
        self.alchex_config = alchex_config
        self.include_from = include_from
        if include_from is None:
            self.include_from = path.split(self.input_topology_filename)[0]
        self.setup()
        self.preprocess()
        self.replace()
        self.energy_minimise()
        self.alchembed()
    def setup(self):
        # Add specified files:
        self.simulations = SimulationContainer(self.root_folder, GromacsWrapper(self.alchex_config.gromacs_executable))
        self.simulations.delete_folder()
        self.simulations.add_file(self.input_structure_filename, rename="input.gro")
        self.simulations.add_file(self.input_topology_filename, rename="input.top")
        # Load extra parameters if set:
        if self.em_parameters_filename is not None:
            mdp_reader = GromacsMDPFile()
            mdp_reader.from_file(self.em_parameters_filename)
            self.additional_em_parameters = mdp_reader.attrs
        else:
            self.additional_em_parameters = {}
        if self.alchembed_parameters_filename is not None:
            mdp_reader = GromacsMDPFile()
            mdp_reader.from_file(self.alchembed_parameters_filename)
            self.additional_alchembed_parameters = mdp_reader.attrs
        else:
            self.additional_alchembed_parameters = {}
        # Read topology and get any additional ITP files
        topology = GromacsTOPFile()
        topology.from_file(self.input_topology_filename)
        for line_id, line, included in topology.get_includes():
            self.simulations.add_file(path.join(self.include_from, included))
        self.simulations.copy_folder("/", "original")
    def preprocess(self):
        input_universe = self.simulations.universe("input.gro")
        # Build a set of replacable entities from the input structure file
        for replacement_specification in self.replacement_specifications:
            input_selection    = input_universe.select_atoms(replacement_specification.selection)
            new_composition    = replacement_specification.normalised_composition()
            groupings = None
            baseunits = {}
            # Check to see if any mappings require residues to be grouped
            for from_resname in set([x.atoms[0].resname for x in input_selection.residues]):
                for to_resname in new_composition.keys():
                    exchange_map = self.alchex_config.get_exchange_map(from_resname, to_resname)
                    baseunits[from_resname] = exchange_map.to_count
                    if exchange_map.grouping is not None:
                        if groupings is not None:
                            # Cannot currently do more than one kind of grouping
                            raise NotImplementedError()
                        else:
                            groupings = exchange_map.grouping
                            groupings["from"] = exchange_map.from_count
                            groupings["to"]   = exchange_map.to_count

            # Select atoms matching the source selection
            replacement_source_selection = input_universe.select_atoms(replacement_specification.selection)

            # Pair residues up if required:
            if groupings is not None:
                group_centres_selection = replacement_source_selection.select_atoms("name " + groupings["atom"])
                group_centres = PointCloud(3)
                group_centres.add_points(group_centres_selection.coordinates())
                paired, unpaired, max_distance = group_centres.group_points(2, d_tol=10, d_max=11)
                group_centres.lines = paired
                # paired can now be used to build exchangable entities

            # Now we need to re-normalise the composition to take into account differences
            # caused by merging residues

            # Determine pairings
            # Add to list of replacable entities
    def replace(self):
        # Carry out replacement of each replacable entity
        pass
    def energy_minimise(self):
        # Determine non-clashing groups of replaceable entites
        # Energy minimise non-clashing groups of replacable entities
        pass
    def alchembed(self):
        # Alchemically add in non-clashing groups of replaceable entites
        pass


d = ReplacementSystem(
    input_structure_filename="data/popc_patch/sample.gro",
    input_topology_filename="data/popc_patch/sample.top",
    replacement_specifications=[
                                    ReplacementSpecification(
                                        selection="resname POPC", 
                                        parsimonious=False, 
                                        composition={"DPPC" : 1, "CDL0" : 1})
                                ]
    )