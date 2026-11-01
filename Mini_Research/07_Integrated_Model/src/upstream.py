"""Read-only, uniquely named imports; bytecode writes disabled by entry scripts."""
from pathlib import Path
import importlib.util
import importlib
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
BASE = '493dbdc213f3789b7195d93760bd5c0c867cb98c'

def load(stage, module):
    name = 'mr7_readonly_' + stage[:2]
    if name not in sys.modules:
        directory = ROOT.parent / stage / 'src'
        spec = importlib.util.spec_from_file_location(name, directory/'__init__.py',
                   submodule_search_locations=[str(directory)])
        package = importlib.util.module_from_spec(spec)
        sys.modules[name] = package
        spec.loader.exec_module(package)
    return importlib.import_module(name+'.'+module)

thermo = load('04_Temperature_Barrier', 'thermodynamics')
mr4 = load('04_Temperature_Barrier', 'model')
optics = load('05_Plasmonic_Heating', 'plasmonics')
Parameters = thermo.Parameters
T0 = thermo.T0
RU = thermo.RU
