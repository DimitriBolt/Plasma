"""Read accepted sibling modules under unique package names; never rewrite them."""
from pathlib import Path
import importlib.util
import importlib
import sys

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
BASE = '0204ce7a8429406b7ec235b02cec8581274fc8c6'

def load(stage, module):
    name = 'mr6_readonly_' + stage[:2]
    if name not in sys.modules:
        directory = ROOT.parent / stage / 'src'
        spec = importlib.util.spec_from_file_location(name, directory/'__init__.py',
                                                     submodule_search_locations=[str(directory)])
        package = importlib.util.module_from_spec(spec)
        sys.modules[name] = package
        spec.loader.exec_module(package)
    return importlib.import_module(name+'.'+module)

optics = load('05_Plasmonic_Heating', 'plasmonics')
thermo = load('04_Temperature_Barrier', 'thermodynamics')
saved = load('03_In_Bubble_Thermodynamics', 'saved_mechanics')
Parameters = thermo.Parameters
