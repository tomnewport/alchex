# coding: utf-8

from itertools import product
import networkx as nx
from cgswap.geometry import PointCloud, TransformationMatrix, plot_3d
import matplotlib.pyplot as plt
from workarounds import WAEditableResidue
from cgswap.residue import MultiResidueStructure


class ExchangeMap(object):
    def __init__(self):
        self.from_resname = ""
        self.from_count = 0
        self.to_resname = ""
        self.to_count = 0
        self.actions = []
        self.methods_implemented = [
            "direct_overlay", 
            "fragment_align", 
            "molecule_align", 
            "simple_bridge", 
            "scale_match_vector"]
    def used_actions(self):
        return set([x["method"] for x in self.actions])
    def scorecard(self):
        sc = ""
        ua = self.used_actions()
        for m in self.methods_implemented:
            if m in ua:
                sc += "☑"
            else:
                sc += "☐"
        return sc
    def _new_map_fragment_alignment(self, remaining, mapped, mapping):
        # Given a list of the remaining <to> and mapped <to>
        to_from_map = {v:k for k, v in mapping.items()}
        to_g = self.to_itp.graph()
        remaining_g = to_g.subgraph(remaining)
        fragment_subgraphs = nx.connected_component_subgraphs(remaining_g)
        fragments = []
        for fragment_subgraph in fragment_subgraphs:
            fragment = {"method":"fragment_align"}
            fragment["nodes"] = fragment_subgraph.nodes()
            connected = list()
            for node in fragment["nodes"]:
                c_nodes = to_g.neighbors(node)
                for c_node in c_nodes:
                    if c_node in mapped:
                        connected.append(c_node)
            while len(connected) < 3:
                cc_nodes = set()
                for c_node in connected:
                    cc_nodes.update(set(to_g.neighbors(c_node)))
                cc_nodes = list(set(cc_nodes).intersection(mapped) - set(connected))
                connected += cc_nodes[:3 - len(connected)]
            fragment["reference"] = [[to_from_map[to_node], to_node] for to_node in connected]
            fragment["centroid_weighting"] = [1,0,0]
            fragments.append(fragment)
        return fragments
    def _new_map_distort_cg_lipid_tails(self, overlay_map):
        # First find lipid tails
        to_tails = self.to_itp.find_lipid_tails()
        from_tails = self.from_itp.find_lipid_tails()
        if self.from_count > 1:
            from_tails = [
                        [    str(i+1) +"."+atom_id
                         for atom_id in tail_atoms
                        ] 
                    for i in range(self.from_count) 
                for tail_atoms in from_tails
            ]
        # Determine which corresponds to which:
        distortions = []
        for from_tail in from_tails:
            for from_atom in from_tail:
                if from_atom in overlay_map:
                    to_atom = overlay_map[from_atom]
                    for to_tail in to_tails:
                        if to_atom in to_tail:
                            distortions.append(
                                {"from"  :from_tail, 
                                 "to"    :to_tail,
                                 "method":"scale_match_vector",
                                 "pins" : {0:0, -1:-1}
                                }
                            )
                            break
                    break
        return distortions
      
    def _new_martini_lipid(self, draw=False):
        direct_overlay = self.from_itp.mcs(self.to_itp)[0]
        # These atoms can be directly replaced
        aligned = set(direct_overlay.values())
        remaining = set(self.to_itp.atoms.keys()) - aligned
        fragments = self._new_map_fragment_alignment(remaining, aligned, direct_overlay)
        self.actions = [{"method":"direct_overlay", "map" : direct_overlay}] + fragments
        self.actions += self._new_map_distort_cg_lipid_tails(overlay_map=direct_overlay)
    def _new_martini_lipid_to_card(self, draw=False):
        card  = self.to_itp
        lipid = self.from_itp
        # There should be two clear analogues of lipid within card
        candidates = lipid.mcs(card, threshold=0.7)
        card_lipid_1 = candidates[0]
        candidates = [x for x in candidates if 
                      len(set(x.values()).intersection(set(card_lipid_1.values()))
                         ) == 0]
        card_lipid_2 = candidates[0]
        lipid_card_1 = {v:k for k, v in card_lipid_1.items()}
        lipid_card_2 = {v:k for k, v in card_lipid_2.items()}
        if draw:
            print("Find two lipid motifs in cardiolipin:")
            card.draw(orange=card_lipid_1.values(), cyan=card_lipid_2.values())
        card_g = card.graph()
        most_central = [x[0] for x in sorted(nx.eccentricity(card_g).items(), key=lambda x : x[1])]
        if draw:
            print("Find central nodes of cardiolipin:")
            card.draw(
                red=most_central[0], 
                orange=most_central[1:3], 
                yellow=most_central[3:5],
                green=most_central[5:9],
                blue=most_central[5:-8]
            )
        lipid_1_bridge = [x for x in most_central if x in card_lipid_1.values()][0]
        lipid_2_bridge = [x for x in most_central if x in card_lipid_2.values()][0]
        l1_bridge_id = str(lipid_card_1[lipid_1_bridge])
        l1_bridge_name = self.from_itp.atoms[l1_bridge_id].attrs["atom"]
        l2_bridge_id = str(lipid_card_2[lipid_2_bridge])
        l2_bridge_name = self.from_itp.atoms[l2_bridge_id].attrs["atom"]
        bridge = nx.shortest_path(card_g, lipid_1_bridge, lipid_2_bridge)
        if draw:
            card.draw(orange=[lipid_1_bridge], cyan=[lipid_2_bridge], yellow=bridge[1:-1])
        combined_map = {k:v for k,v in [("2."+x[0],x[1]) for x in card_lipid_2.items()] + [("1."+x[0],x[1]) for x in card_lipid_1.items()]}
        self.actions = []
        self.actions.append({
                "method"  : "direct_overlay",
                "map"    : combined_map
            })
        self.actions.append({
                "method"   : "simple_bridge",
                "from"  : combined_map["1." + l1_bridge_id],
                "to"    : combined_map["2." + l2_bridge_id],
                "from_selector" : "resname " + self.from_itp.resname + " and name " + l1_bridge_name, 
                "to_selector"   : "resname " + self.from_itp.resname + " and name " + l2_bridge_name, 
                "atoms" : bridge
            })
        remaining = [x for x in card.atoms if x not in combined_map.values() and x not in bridge]
        fragments = self._new_map_fragment_alignment(remaining, combined_map.values() + bridge, combined_map)
        self.actions += self._new_map_fragment_alignment(remaining, combined_map.values() + bridge, combined_map)
        self.actions += self._new_map_distort_cg_lipid_tails(overlay_map = combined_map)
        #print(remaining)
        #print(self)
    def _new_martini_static_planar_alignment(self, clusters, draw=False):
        id_clusters = []
        for from_cluster, to_cluster, weight in clusters:
            id_clusters.append([
                    self.from_itp.atoms_to_ids(from_cluster), 
                    self.to_itp.atoms_to_ids(to_cluster), 
                    weight])
        self.actions = [
            {"method":"molecule_align",
             "clusters" : id_clusters}
        ]
    def __repr__(self):
        return " 〘 ⚗ Change {from_count}x {from_resname} to {to_count}x {to_resname} {sc} 〙 ".format(sc=self.scorecard(), **self.__dict__)
    def new(self, from_itp=None, to_itp=None, method="martini.lipid", draw=False, **kwargs):
        self.from_itp = from_itp
        self.to_itp = to_itp
        self.from_resname = self.from_itp.atoms.values()[0].attrs["residue"]
        self.from_count = 1
        self.to_resname = self.to_itp.atoms.values()[0].attrs["residue"]
        self.to_count = 1
        if method == "martini.lipid":
            return self._new_martini_lipid(draw=draw)
        if method == "martini.lipid_to_card":
            self.from_count = 2
            return self._new_martini_lipid_to_card(draw=draw)
        if method == "martini.static_planar_alignment":
            return self._new_martini_static_planar_alignment(draw=draw, **kwargs)
    def _run_molecule_align(self, action, from_residue, to_residue, new_residue):
        from_pointcloud = PointCloud(3)
        to_pointcloud = PointCloud(3)
        new_residue.import_mdanalysis_atoms(to_residue)
        for from_cluster, to_cluster, weight in action["clusters"]:
            from_pointcloud.add_points([from_residue.position(from_cluster).mean(axis=0)])
            to_pointcloud.add_points([to_residue.position(to_cluster).mean(axis=0)])
        transformation, rmse, aligned = from_pointcloud.paired_3d_align(to_pointcloud, inv=False)
        new_residue.transform(transformation)
    def run(self, from_residue, to_residue):
        # Takes two MDAnalysis residues and replaces one with the other
        new_residue = WAEditableResidue(resname=to_residue.resname, resid=to_residue.resid)
        for action in self.actions:
            if action["method"] == "molecule_align":
                self._run_molecule_align(action, from_residue, to_residue, new_residue)
            elif action["method"] == "direct_overlay":
                new_residue.overlay(
                    from_residue=from_residue, 
                    to_residue=to_residue, 
                    mapping = action["map"])
            elif action["method"] == "fragment_align":
                new_residue.align_fragment(
                    from_residue = from_residue,
                    to_residue = to_residue,
                    params = action
                    )
            elif action["method"] == "scale_match_vector":
                pass
                #new_residue.scale_match_vector(from_residue, action)
            elif action["method"] == "simple_bridge":
                new_residue.build_bridge(from_residue, to_residue, action)
            else:
                print(action["method"])
        #plot_3d(from_residue.point_cloud(), to_residue.point_cloud(), new_residue.coordinates)
        return new_residue