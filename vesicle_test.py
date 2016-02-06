from alchex.replacement import ReplacementSystem


d = ReplacementSystem(
    input_structure_filename="vesicle/pope_vesicle.gro",
    input_topology_filename="vesicle/pope_vesicle.top",
    root_folder="vesicle-test-run"
    )

d.auto_replace({
	"selection": "resname POPE",
	"composition" :{
		"DPPC" : 1,
		"DLPG" : 1
		}})