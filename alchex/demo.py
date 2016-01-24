import sys
sys.path.append("../")

from alchex.alchex_config import AlchexConfig
from alchex.box import SimulationBox

config = AlchexConfig()

config.load_itp_file("../data/DLPG.itp", "DLPG")
config.load_itp_file("../data/DVPE.itp", "DVPE")
config.load_itp_file("../data/CDL.itp", "CDL0")
config.load_itp_file("../data/chol.itp", "CHOL")
config.load_itp_file("../data/dppc.itp", "DPPC")
config.load_itp_file("../data/popc_patch/popc_ac.itp", "POPC")

config.build_exchange_map(
    from_resname="DPPC", 
    to_resname = "CHOL", 
    exchange_model="martini.static_planar_alignment", 
    draw=False, 
    clusters=[
        [
            ["PO4", "NC3"],
            ["ROH"],
            1
        ],
        [
            ["C1B", "C2A"],
            ["R2", "R4"],
            1
        ],
        [
            ["C4B", "C4A"],
            ["C2"],
            1
        ]
    ])

config.build_exchange_map(
    from_resname="DPPC",
    to_resname="DLPG",
    exchange_model="martini.lipid")

config.build_exchange_map(
    from_resname="POPC",
    to_resname="CDL0",
    exchange_model="martini.lipid_to_card")

config.build_exchange_map(
    from_resname="POPC",
    to_resname="DLPG",
    exchange_model="martini.lipid")

config.add_reference_structure("DLPG","../data/DLPG-em.gro")

config.add_composition( "all_dlpg", DLPG=1)

popc_patch = SimulationBox(config)

popc_patch.load_universe("../data/popc_patch/sample.gro")
popc_patch.add_replacement("resname POPC", "all_dlpg")
popc_patch.perform_replacement("../output.gro")
