"""Independent high-amplitude variational sign check on archived Radau mechanics."""
from .stable_pairs import archived
from .diagnostics import extrema,evaluate
from .experiments import table

def main():
    b,check=archived('extra_high_1.15_Radau_base');e=extrema(b)
    s=b.model.sensitivity(b,e['tc_s']-5e-6)
    r,z=evaluate(s,'high_Radau_1.15_tangent',do_quad=False)
    table('high_Radau_sensitivity',[r]);table('high_Radau_interpolation',[check])
    print('Independent 1.15 bar Radau derivative',r['dTmax_dF_K_m2_J'],flush=True)

if __name__=='__main__':main()
