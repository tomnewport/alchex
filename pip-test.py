from alchex.replacement import ReplacementSystem

d = ReplacementSystem(
    input_structure_filename="pip-test/original/centred.gro",
    input_topology_filename="pip-test/original/input.top",
    root_folder="pip-test-run"
    )

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
