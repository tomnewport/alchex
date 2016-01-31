from alchex.replacement import ReplacementSystem

d = ReplacementSystem(
    input_structure_filename="pip-test/original/input.gro",
    input_topology_filename="pip-test/original/input.top",
    root_folder="pip-test-run"
    )

d.auto_replace(
    selection="resname DPPC", 
    composition={"PI3P" : 100}
    )
