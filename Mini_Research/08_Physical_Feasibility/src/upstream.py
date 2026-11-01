"""Read-only imports of accepted code under distinct package names."""
from pathlib import Path
import importlib.util
import importlib
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]

def load(stage, module):
    name = 'mr8_readonly_' + stage[:2]
    if name not in sys.modules:
        directory = ROOT.parent / stage / 'src'
        spec = importlib.util.spec_from_file_location(name, directory/'__init__.py',
                    submodule_search_locations=[str(directory)])
        package = importlib.util.module_from_spec(spec)
        sys.modules[name] = package
        spec.loader.exec_module(package)
    return importlib.import_module(name+'.'+module)

mechanics = load('02_Compressible_Bubble_Dynamics', 'model')
thermo = load('04_Temperature_Barrier', 'thermodynamics')
