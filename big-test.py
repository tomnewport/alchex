from alchex.replacement import ReplacementSystem

d = ReplacementSystem(
    input_structure_filename="big/equil_memb.gro",
    input_topology_filename="big/equil_memb.top",
    root_folder="big-test-run"
    )

'''
d.auto_replace(
		{
			"selection"  :"leaflet upper",
		 	"composition": {
		 	    "DLPG" : 30
		 	    }
		},
		{	
			"selection"  :"leaflet lower",
		 	"composition": {
		 	    "PI3P" : 1, 
		 	    "CDL0" : 1,
		 	    "CHOL" : 1,
		 	    "DLPG" : 1,
		 	    "DPPC" : 1
		 	    }}
    )
'''

d.auto_replace({
	"selection": "leaflet upper",
	"composition" :{
		"DPPC" : 1,
		"CDL0" : 1
		}},
	{
	"selection": "leaflet lower",
	"composition" :{
		"DPPC" : 1,
		"CDL0" : 1
		}}
	)




'''
d = ReplacementSystem(
    input_structure_filename="vesicle/pope_vesicle.gro",
    input_topology_filename="vesicle/pope_vesicle.top",
    root_folder="vesicle-test-run"
    )

d.auto_replace({
	"selection": "resname POPE",
	"composition" :{
		"DPPC" : 1,
		"DLPG" : 1,
		"CDL0" : 1
		}})
'''