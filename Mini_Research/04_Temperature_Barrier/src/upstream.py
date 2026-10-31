"""Read-only imports of accepted MR3, under a distinct package name."""
from pathlib import Path
import importlib.util
import importlib
import sys

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
MR3 = ROOT.parent / '03_In_Bubble_Thermodynamics'
BASE = '1a12fac295dc904569ab2a57ce4f68a8ea75f0c4'
name = 'mr4_readonly_mr3'
spec = importlib.util.spec_from_file_location(name, MR3/'src/__init__.py',
                                             submodule_search_locations=[str(MR3/'src')])
module = importlib.util.module_from_spec(spec)
sys.modules[name] = module
spec.loader.exec_module(module)
def old(module_name):
    return importlib.import_module(name+'.'+module_name)
Parameters = old('bubble_thermo').Parameters
RU = old('bubble_thermo').RU
SavedMechanics = old('saved_mechanics').SavedMechanics
