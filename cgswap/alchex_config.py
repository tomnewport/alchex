from cgswap.gromacs import GromacsITPFile
from cgswap.exchange_map import ExchangeMap

class AlchexConfig(object):
    def __init__(self, folder="alchex_configuration"):
        self.folder = folder
        self.parameters = {}
        self.exchange_maps = {}
        self.compositions = {}
    def load_itp_file(self, filename, resname):
        itp        = GromacsITPFile(filename)
        parameters = itp.read_residue(resname)
        self.parameters[resname] = parameters
    def build_exchange_map(self, from_resname, to_resname, exchange_model, draw=False, **kwargs):
        from_parameters = self.parameters[from_resname]
        to_parameters = self.parameters[to_resname]
        newmap = ExchangeMap()
        newmap.new(
            from_itp=from_parameters, 
            to_itp=to_parameters, 
            method=exchange_model, 
            draw=False, 
            **kwargs)
        if from_parameters not in self.exchange_maps:
            self.exchange_maps[from_parameters] = {}
        self.exchange_maps[from_parameters][to_parameters] = newmap
        print(newmap)
    def add_composition(self, name, **resname_fractions):
        self.compositions[name] = resname_fractions
