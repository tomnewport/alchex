#!/usr/bin/env python
# -*- coding: utf-8 -*- 

from alchex.gromacs_interface import GromacsWrapper, SimulationContainer, GromacsTOPFile, GromacsMDPFile, GromacsEditableTOPFile
from alchex.config import default_configuration
from alchex.geometry import PointCloud
from alchex.errors import GromacsInterfaceError
from alchex.residue import ResidueStructure, MultiResidueStructure
from alchex.workarounds import WAEditableGrofile
from alchex.lipid_analysis import find_bilayer_leaflets
from os import path
from collections import Counter
import logging
from random import shuffle
import networkx as nx

logging.basicConfig(format=' ⚗ %(levelname)s : %(message)s', level=logging.INFO)

def norm_dict(dictionary):
    sum_values = 1.0 * sum(dictionary.values())
    return {k:v/sum_values for k, v in dictionary.items()}

def calculate_swap(graph, a, b):
    sp = nx.shortest_path(graph, a,b)
    r = 1
    for n1, n2 in zip(sp[:-1], sp[1:]):
        r *= graph.get_edge_data(n1, n2)["multiply"]
    return r

class ReplacementSpecification(object):
    def __init__(self, 
        alchex_config,
        selection=None, 
        parsimonious=None, 
        composition=None, 
        from_file=None):
        self.alchex_config = alchex_config
        if from_file is not None:
            raise NotImplementedError()
        else:
            self.selection = selection
            self.parsimonious = parsimonious
            self.composition = composition
    def normalised_composition(self):
        return norm_dict(self.composition)
    def map_counts(self, input_counts):
        # Accepts a dictionary of input residue counts and
        # returns a dictionary of output residue counts such
        # that creating v new instances of k will use up all
        # inputs and result in a final ratio as close as possible
        # to the self.composition ratio.
        swap_ratios = nx.DiGraph()
        to_counts = { k:0 for k in self.composition.keys() }
        sizes = {}
        for from_resname, from_count in input_counts.items():
            ratios = {}
            if from_resname not in swap_ratios:
                swap_ratios.add_node(from_resname)
            for to_resname in self.composition.keys():
                ex_map = self.alchex_config.get_exchange_map(from_resname, to_resname)
                if to_resname not in sizes:
                    sizes[to_resname] = []
                sizes[to_resname].append(ex_map.from_count/ex_map.to_count)
                if to_resname not in swap_ratios:
                    swap_ratios.add_node(to_resname)
                swap_ratios.add_edge(from_resname, to_resname, multiply=ex_map.to_count/ex_map.from_count)
                swap_ratios.add_edge(to_resname, from_resname, multiply=ex_map.from_count/ex_map.to_count)

                ratios[to_resname] = self.composition[to_resname] * (1.0*ex_map.from_count/ex_map.to_count)
            ratios = norm_dict(ratios)
            for to_resname, to_ratio in ratios.items():
                ex_map = self.alchex_config.get_exchange_map(from_resname, to_resname)
                to_counts[to_resname] = (from_count * to_ratio) * (1.0*ex_map.to_count/ex_map.from_count)
        # Ideal counts obtained, we should modify them to ensure whole molecules:
        real_counts = {k:v for k, v in to_counts.items()}
        ideal_counts = {k:v for k, v in to_counts.items()}
        for resname, ideal_count in sorted(ideal_counts.items(), key=lambda x : -sum(sizes[x[0]])):
            ideal_count = to_counts[resname]
            del to_counts[resname]
            real_count = int(round(ideal_count))
            if real_count == 0:
                real_count += 1
            real_counts[resname] = real_count
            # Need to redistribute a fraction of this molecule amongst others
            if len(to_counts) > 0:
                redistribute = (ideal_count - real_count) / len(to_counts)
                for rdist_resname, count in to_counts.items():
                    to_counts[rdist_resname] = count + redistribute * calculate_swap(swap_ratios, resname, rdist_resname)
        return real_counts

class ReplaceableEntity(object):
    def __init__(self, exchange_model, residue_ids, new_residue_ids, replacement_system=None):
        self.replacement_system = replacement_system
        self.residue_ids = residue_ids
        self.new_residue_ids = new_residue_ids
        self.exchange_model = exchange_model
    def setup(self):
        from_parameters = self.replacement_system.alchex_config.parameters[self.exchange_model.from_resname]
        to_parameters = self.replacement_system.alchex_config.parameters[self.exchange_model.to_resname]
        residues = []
        mda_universe = self.replacement_system.simulations.universe("input.gro")
        for from_resid in self.residue_ids:
            residues.append(ResidueStructure(mda_universe.select_atoms("resid " + str(from_resid)), from_parameters))
        if len(residues) == 1:
            self.from_residue = residues[0]
        else:
            self.from_residue = MultiResidueStructure(residues)
    def replace(self):
        to_residue = self.replacement_system.alchex_config.get_reference_structure(self.exchange_model.to_resname)
        self.replacement = self.exchange_model.run(self.from_residue, to_residue)
        self.replacement.resid = self.new_residue_ids
        self.replacement.moltype = self.exchange_model.to_moltype
        return self.replacement

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
        root_folder="./alchex-replacement"):
        self.input_structure_filename = input_structure_filename
        self.input_topology_filename = input_topology_filename
        self.output_structure_filename = output_structure_filename
        self.output_topology_filename = output_topology_filename
        self.em_parameters_filename = em_parameters_filename
        self.alchembed_parameters_filename = alchembed_parameters_filename
        self.root_folder = root_folder
        self.alchex_config = alchex_config
        self.include_from = include_from
        if include_from is None:
            self.include_from = path.split(self.input_topology_filename)[0]
        self.setup()
        self._custom_groups = {}
    def custom_group(self, group_name):
        universe = self.simulations.universe("input.gro")
        if group_name not in self._custom_groups:
            if group_name in ["leaflet upper", "leaflet lower"]:  
                lower_leaflet_resids, upper_leaflet_resids = find_bilayer_leaflets(universe)
                self._custom_groups["leaflet lower"] = "("+" or ".join(
                    ["resid "+str(resid) for resid in lower_leaflet_resids]
                    )+")"
                self._custom_groups["leaflet upper"] = "("+" or ".join(
                    ["resid "+str(resid) for resid in upper_leaflet_resids]
                    )+")"
        return self._custom_groups[group_name]
    def preprocess_selection(self, string):
        custom_groups = ["leaflet upper", "leaflet lower"]
        for cg in custom_groups:
            if cg in string:
                string = string.replace(cg, self.custom_group(cg))
        return string
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
    def build_replaceable_entities(self, map_counts, replacement_specification):
        # Accepts a list of residue names and counts
        input_universe       = self.simulations.universe("input.gro")
        replaceable_universe = input_universe.select_atoms(self.preprocess_selection(replacement_specification.selection))
        available_entities   = replaceable_universe.residues
        resid_to_name        = {r.id : r.name for r in available_entities}
        from_resnames = Counter([r.name for r in available_entities])


        from_group_sizes = {}
        groupings = {}

        for to_resname in map_counts:
            from_group_sizes[to_resname] = {}
            for from_resname in from_resnames.keys():
                ex_map = self.alchex_config.get_exchange_map(from_resname, to_resname)
                from_group_sizes[to_resname][from_resname] = ex_map.from_count
                if from_resname not in groupings or ex_map.from_count not in groupings[from_resname]:
                    group = []
                    if ex_map.from_count == 1:
                        group = [{r.id,} for r in available_entities if r.name==from_resname]
                    elif ex_map.from_count == 2:
                        group_centres_selection = replaceable_universe.select_atoms("name " + ex_map.grouping["atom"])
                        group_centres = PointCloud(3)
                        group_centres.add_points(group_centres_selection.coordinates())
                        logging.info("Pairing points using atom name " + ex_map.grouping["atom"] + ". This may take some time...")
                        paired, unpaired, max_distance = group_centres.group_points(2, d_tol=10, d_max=11)
                        logging.info("Pairing completed:\n  {len_unpaired} were not paired, maximum distance was {d_max}Å".format(len_unpaired=len(unpaired), d_max=max_distance))
                        # Return to residue names:
                        for p1, p2 in paired:
                            group.append({group_centres_selection[p1].resid, group_centres_selection[p2].resid})
                    if from_resname not in groupings:
                        groupings[from_resname] = {}
                    groupings[from_resname][ex_map.from_count] = group
        
        replaceable_entities = []
        for to_resname, final_count in sorted(map_counts.items(), key=lambda x : -sum(from_group_sizes[x[0]].values())):
            # We want to randomly pick final_count groups to become
            # to_resname.
            available = []
            for from_resname in from_resnames.keys():
                ex_map = self.alchex_config.get_exchange_map(from_resname, to_resname)
                available += groupings[from_resname][ex_map.from_count]
            shuffle(available)
            used = available[:final_count]
            used_resids = set([y for x in used for y in x])

            # Remove used residue ids
            for resname, bycount in groupings.items():
                for c in bycount:
                    groupings[resname][c] = [x for x in groupings[resname][c] if len(used_resids.intersection(x)) == 0]
            for from_group in used:
                new_resid = min(from_group)
                from_resname = resid_to_name[min(from_group)]
                ex_map = self.alchex_config.get_exchange_map(from_resname, to_resname)
                replaceable_entities.append(
                    ReplaceableEntity(
                        ex_map, 
                        from_group, 
                        new_resid)
                    )

        return replaceable_entities
    def _replace(self, replaceable_entities, name="alchex"):
        '''
        Given a set of replaceable_entities, performs replacement in the specified directory
        '''
        logging.info("Preparing to replace molecules...")
        self.simulations.copy_folder("/original", "/" + name)
        self.simulations.cd("/"+name)
        original = WAEditableGrofile()
        original.from_file(self.simulations.resolve_path("/"+name+"/input.gro"))
        original_topology = GromacsEditableTOPFile()
        original_topology.from_file(self.simulations.resolve_path("/"+name+"/input.top"))
        em_mdp = GromacsMDPFile()
        em_mdp.attrs = self.alchex_config.grompp_parameters["em"]
        em_mdp.to_file(self.simulations.resolve_path("em.mdp"))
        exitcode, message = self.simulations.gromacs.grompp(kwargs={
                    "-f" : "em.mdp", 
                    "-c" : "input.gro", 
                    "-p" : "input.top",
                    "-o" : "preprocess",
                    "-pp": "input-pp.top",
                    "-maxwarn":"1"})
        if exitcode != 0:
            raise GromacsInterfaceError(message)
        full_topology = GromacsEditableTOPFile()
        full_topology.from_file(self.simulations.resolve_path("/"+name+"/input-pp.top"))
        original.add_topology(full_topology)
        replaced = WAEditableGrofile()
        logging.info("Performing replacements...")
        todo = len(replaceable_entities)
        done = 0
        for entity in replaceable_entities:
            to_resname = entity.exchange_model.to_resname
            to_moltype = entity.exchange_model.to_moltype
            if to_moltype not in original.moltypes:
                insert_idx = original.moltypes.index(entity.exchange_model.from_moltype)
                original.moltypes.insert(insert_idx, entity.exchange_model.to_moltype)
            entity.replacement_system = self
            entity.setup()
            new_residues = entity.replace()
            original.delete_by_resids(entity.residue_ids)
            original.residues.append(new_residues)
            done += 1
            print(done, todo)
        original_topology.modify_molecules(original.top_molecules())
        original_topology.to_file(self.simulations.resolve_path("replaced.top"))
        original.sort_residues()
        original.renumber_moltypes()
        logging.info("Performing energy minimisation simulation...")
        #print(len(original.declash_moltypes(2)))
        # Energy minimise the new system
        original.to_file(self.simulations.resolve_path("replaced.gro"))   
             
        statuscode, message = self.simulations.gromacs.grompp(kwargs={
                    "-f" : "em.mdp", 
                    "-c" : "replaced.gro", 
                    "-p" : "replaced.top",
                    "-o" : "all-em.tpr",
                    "-maxwarn":"1"})
        if statuscode != 0:
            raise GromacsInterfaceError(message)

        statuscode, message = self.simulations.gromacs.mdrun(kwargs={
            "-deffnm" : "all-em"
            })
        if statuscode != 0:
            raise GromacsInterfaceError(message)

        original.reload(self.simulations.resolve_path("all-em.gro"), same_atom_names=False)
        #print(len(original.declash_moltypes(2)))
        # Look for clashes
        # Use alchembed to resolve clashes
    def replace(self, replaceable_entities, name="alchex"):
        '''
        Old method, see _replace instead.
        Given a set of replaceable_entities, performs replacement in the specified directory
        '''
        self.simulations.copy_folder("/original", "/" + name)
        self.simulations.cd("/"+name)
        replacement_groups = {}
        '''
        for entity in replaceable_entities:
            to_resname = entity.exchange_model.to_resname
            if to_resname not in replacement_groups:
                replacement_groups[to_resname] = []
            entity.replacement_system = self
            entity.setup()
            replacement_groups[to_resname].append(entity.replace())
        '''
        original = WAEditableGrofile()
        original.from_file(self.simulations.resolve_path("/"+name+"/input.gro"))
        original_topology = GromacsTOPFile()
        original_topology.from_file(self.simulations.resolve_path("/"+name+"/input.top"))
        em_mdp = GromacsMDPFile()
        em_mdp.attrs = self.alchex_config.grompp_parameters["em"]

        '''
        for p, v in self.additional_em_parameters.items():
            em_mdp.attrs[p] = v
        for to_resname, residues in replacement_groups.items():
            self.simulations.cd("/"+name)
            self.simulations.copy_folder("/original", to_resname)
            self.simulations.cd(to_resname)
            em_mdp.to_file(self.simulations.resolve_path("em.mdp"))
            structure = WAEditableGrofile()
            structure.residues = residues
            structure.sysname = "Alchex replacement component " + to_resname
            structure.box_vector = original.box_vector
            structure.to_file("debug-"+to_resname+".gro")
            structures = structure.declash(2)
            for idx, dc_structure in enumerate(structures):
                self.simulations.cd("/" + name + "/" + to_resname)
                base_filename = "alchex_component."+str(idx)
                em_dir = base_filename+".em"
                dc_structure.to_file(self.simulations.resolve_path(base_filename+".gro"))
                original_topology.modify_molecules(dc_structure.top_molecules())
                original_topology.to_file(self.simulations.resolve_path(base_filename+".top"))
                self.simulations.makedirs(em_dir)
                print self.simulations.gromacs.grompp(kwargs={
                    "-f":"em.mdp", 
                    "-c": base_filename+".gro", 
                    "-p": base_filename+".top",
                    "-o": em_dir+"/"+em_dir+".tpr",
                    "-maxwarn":"1"})
                self.simulations.cd(em_dir)
                grompp_result, grompp_message = self.simulations.gromacs.mdrun(kwargs={"-deffnm":em_dir})
                if "Steepest Descents converged to Fmax" in grompp_message:
                    logging.info("Successfully energy minimised " + to_resname + "/" +  base_filename + ".")
        '''
    def auto_replace(self, *groups):
        '''Shorthand to create replaceable entities and then run replacement'''
        replaceable_entities = []
        for group in groups:
            selection = group["selection"]
            composition = group["composition"]

            replacement_specification = ReplacementSpecification(self.alchex_config,
                selection=selection, 
                composition=composition, 
                parsimonious=False)
            input_universe       = self.simulations.universe("input.gro")

            replaceable_universe = input_universe.select_atoms(self.preprocess_selection(selection))
            # Determine how many of each residue are going to be created:
            map_counts = replacement_specification.map_counts(
                Counter(
                [x.atoms[0].resname for x in replaceable_universe.residues]
                )
                )
            replaceable_entities += self.build_replaceable_entities(map_counts, replacement_specification)
        self._replace(replaceable_entities)
    def energy_minimise(self):
        # Determine non-clashing groups of replaceable entites
        # Energy minimise non-clashing groups of replacable entities
        pass
    def alchembed(self):
        # Alchemically add in non-clashing groups of replaceable entites
        pass

'''
d = ReplacementSystem(
    input_structure_filename="data/testpatch/sample.gro",
    input_topology_filename="data/testpatch/sample.top")

d.auto_replace(selection="resname POPC", composition={"DPPC" : 50, "CDL0" : 24})


d = ReplacementSystem(
    input_structure_filename="pip-test/original/centred.gro",
    input_topology_filename="pip-test/original/input.top",
    root_folder="pip-test-run"
    )
print d.preprocess_selection("hello world")
print d.preprocess_selection("leaflet upper")
print d.preprocess_selection("leaflet lower and resid 10")
'''