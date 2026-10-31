"""Make every scientific illustration from saved CSV/JSON, without reintegration."""
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
COLORS = ["#777777", "#0072B2", "#009E73", "#D55E00", "#CC79A7", "#111111"]


def rows(name):
    with (ROOT / "results" / name).open() as f:
        return list(csv.DictReader(f))


def write_report_tables(meta, summary, convergence):
    conv = [row for row in convergence if float(row["Pa_Pa"]) == 1.3e5]
    conv_one = [row for row in convergence if float(row["Pa_Pa"]) == 1.0e5]
    """Keep report numbers mechanically tied to the saved data."""
    out = ROOT / "results"
    def sci(x):
        if float(x) == 0:
            return "0"
        mantissa, exponent = f"{float(x):.2e}".split("e")
        return mantissa+r"\\times 10^{"+str(int(exponent))+"}"
    lines = [r"\begin{tabular}{rrrrrrr}",r"\toprule",
             r"$P_a$ (bar) & $R_{\max}$ (\textmu m) & $R_{\min}$ (\textmu m) & $U_+$ (m/s) & $|U_-|$ (m/s) & $M_{w,\max}$ & $N_+/N_-$ \\",r"\midrule"]
    for row in summary:
        values=[f"{float(row['Pa_bar']):.1f}",f"{float(row['Rmax_m'])*1e6:.4f}",
                f"{float(row['Rmin_m'])*1e6:.5f}",f"{float(row['Umax_positive_m_s']):.4g}",
                f"{float(row['Umax_inward_m_s']):.4g}",f"{float(row['max_wall_Mach']):.3g}",
                row['local_maxima_count']+"/"+row['local_minima_count']]
        if values[0]=='0.0':
            values[3:6]=['0','0','0']
        lines.append(" & ".join(values)+r" \\")
    lines.extend([r"\bottomrule",r"\end{tabular}"])
    (out/'report_amplitude_table.tex').write_text('\n'.join(lines)+'\n')
    for selection, name in ((conv, 'report_convergence_table.tex'),
                            (conv_one, 'report_convergence_1bar_table.tex')):
        lines=[r"\begin{tabular}{lrrrrr}",r"\toprule",
               r"Method / rtol & $R_{\max}$ (\textmu m) & $R_{\min}$ (\textmu m) & $t_c$ (\textmu s) & $|U|_{\max}$ (m/s) & Max. rel. diff. \\",r"\midrule"]
        for row in selection:
            err=max(float(row['relative_error_'+key]) for key in ('Rmax_m','Rmin_m','collapse_time_s','max_abs_U_m_s'))
            lines.append(row['method']+' / $'+sci(row['rtol']).replace('\\\\','\\')+'$ & '+
                         f"{float(row['Rmax_m'])*1e6:.6f} & {float(row['Rmin_m'])*1e6:.7f} & "+
                         f"{float(row['collapse_time_s'])*1e6:.7f} & {float(row['max_abs_U_m_s']):.3f} & $"+
                         sci(err).replace('\\\\','\\')+r"$ \\")
        lines.extend([r"\bottomrule",r"\end{tabular}"])
        (out/name).write_text('\n'.join(lines)+'\n')
    st=summary[-1]
    vals=dict(eqError=sci(meta['equilibrium']['relative_radius_error']),
              eqSpeed=sci(meta['equilibrium']['max_abs_U_m_s']),
              eqResidual=sci(meta['equilibrium']['pressure_residual_Pa']),
              pgInitial=f"{meta['derived']['pg0_Pa']:.3f}",
              naturalFreq=f"{meta['derived']['natural_frequency_Hz']/1000:.2f}",
              strongCollapseTime=f"{float(st['collapse_time_s'])*1e6:.6f}",
              strongCollapseDuration=f"{float(st['collapse_duration_s'])*1e6:.3f}",
              strongExpansionDuration=f"{float(st['preceding_expansion_duration_s'])*1e6:.3f}",
              strongMach=f"{float(st['max_wall_Mach']):.2f}",
              linearLow=sci(meta['linear_checks'][0]['max_radius_error_over_R0']),
              linearHigh=sci(meta['linear_checks'][1]['max_radius_error_over_R0']),
              convergedDiff=sci(max(float(conv[2]['relative_error_'+k]) for k in ('Rmax_m','Rmin_m','collapse_time_s','max_abs_U_m_s'))),
              radauDiff=sci(max(float(conv[4]['relative_error_'+k]) for k in ('Rmax_m','Rmin_m','collapse_time_s','max_abs_U_m_s'))),
              convergedOnebarDiff=sci(max(float(conv_one[2]['relative_error_'+k]) for k in ('Rmax_m','Rmin_m','collapse_time_s','max_abs_U_m_s'))),
              radauOnebarDiff=sci(max(float(conv_one[4]['relative_error_'+k]) for k in ('Rmax_m','Rmin_m','collapse_time_s','max_abs_U_m_s'))),
              minStrongStep=sci(meta['solver_runs'][-1]['min_step_s']))
    (out/'report_values.tex').write_text('\n'.join(
        "\\newcommand{\\"+k+"}{"+v.replace('\\\\','\\')+"}" for k,v in vals.items())+'\n')


def main():
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":9,
                         "axes.titlesize":10, "axes.labelsize":9,
                         "legend.fontsize":8, "lines.linewidth":1.15,
                         "axes.spines.top":False, "axes.spines.right":False,
                         "axes.grid":True, "grid.alpha":.18, "pdf.fonttype":42})
    figdir = ROOT / "figures"
    figdir.mkdir(exist_ok=True)
    meta = json.loads((ROOT / "results/run_metadata.json").read_text())
    p = meta["parameters_SI"]
    T, R0, c = 1/p["f"], p["R0"], p["c"]
    data = np.genfromtxt(ROOT / "results/trajectories.csv", delimiter=",", names=True)
    summary = rows("amplitude_summary.csv")
    amps = np.array([float(row["Pa_bar"]) for row in summary])
    curves = {bar:data[np.isclose(data["Pa_Pa"],bar*1e5)] for bar in amps}
    strong = curves[1.3]
    sr = summary[-1]
    tc = float(sr["collapse_time_s"])
    tm = float(sr["preceding_maximum_time_s"])

    def finish(fig,name):
        fig.savefig(figdir / (name+".pdf"),bbox_inches="tight")
        fig.savefig(figdir / (name+".png"),dpi=300,bbox_inches="tight")
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.8,3.0),layout="constrained")
    ax.set(xlim=(-2.8,3.2),ylim=(-1.35,1.4),aspect="equal")
    ax.set_facecolor("#eef6fa")
    ax.add_patch(Circle((0,0),1,facecolor="white",edgecolor="black",linewidth=1.5))
    ax.annotate("",xy=(.71,.71),xytext=(0,0),arrowprops=dict(arrowstyle="->",color="black"))
    ax.text(.18,.52,r"$R(t)$",rotation=45)
    ax.text(-.61,-.25,"Noncondensable gas\n"+r"$p_g(R)$ + vapor $p_v$",fontsize=9)
    ax.text(1.2,.75,"Water at 20 °C\n"+r"$\rho,\ \mu,\ \sigma$",fontsize=9)
    for x,y,xe,ye in [(-2,0,-1.1,0),(2,0,1.1,0),(0,1.35,0,1.06),(0,-1.35,0,-1.06)]:
        ax.annotate("",xy=(xe,ye),xytext=(x,y),arrowprops=dict(arrowstyle="->",color=COLORS[1]))
    ax.text(-2.7,.75,r"$p_\infty(t)=p_0-P_a\sin(2\pi f t)$",fontsize=9)
    ax.text(-2.7,-1.17,r"Spherical interface; $R_0=6\ \mu$m; $f=30$ kHz",fontsize=9)
    ax.axis("off")
    finish(fig,"01_bubble_schematic")

    fig,axes = plt.subplots(2,1,figsize=(6.8,4.1),sharex=True,layout="constrained")
    for ax,bar in zip(axes,[1.,1.3]):
        d = curves[bar]
        ax.plot(d["t_s"]/T,d["R_m"]/R0,color=COLORS[3] if bar==1 else COLORS[-1])
        ax.set(ylabel=r"$R/R_0$",title=f"$P_a={bar:.1f}$ bar; $R_0=6$ µm; $f=30$ kHz")
    axes[-1].set_xlabel(r"Time $t/T$ (one period = 33.333 µs)")
    finish(fig,"02_radius_vs_time")

    fig,axes=plt.subplots(2,1,figsize=(6.8,4.1),sharex=True,layout="constrained")
    for ax,bar in zip(axes,[1.,1.3]):
        d=curves[bar]
        t,U=d["t_s"]/T,d["U_m_s"]
        ax.plot(t,U,color="black",linewidth=.7)
        ax.fill_between(t,U,0,where=U>=0,color=COLORS[1],alpha=.4,label="Expansion: U > 0")
        ax.fill_between(t,U,0,where=U<0,color=COLORS[3],alpha=.4,label="Collapse: U < 0")
        ax.set(ylabel="U (m/s)", title=f"$P_a={bar:.1f}$ bar; $R_0=6$ µm; $f=30$ kHz")
    axes[0].legend(loc="lower left",ncol=2,fontsize=7)
    axes[1].set_yscale("symlog",linthresh=1.)
    axes[1].set_ylabel("U (m/s; symlog)")
    axes[1].set_xlabel(r"Time $t/T$")
    finish(fig,"03_wall_velocity_vs_time")

    fig,axes=plt.subplots(1,2,figsize=(6.8,3.1),layout="constrained")
    for ax,bar in zip(axes,[1.,1.3]):
        d=curves[bar]
        ax.plot(d["R_m"]/R0,d["U_m_s"],color=COLORS[1] if bar==1 else "black",linewidth=.7)
        ax.plot(1,0,"o",ms=4,color=COLORS[3],label="Initial state")
        ax.set(xlabel=r"$R/R_0$",ylabel=r"$U$ (m/s)",title=f"$P_a={bar:.1f}$ bar")
    axes[0].legend()
    axes[1].set_yscale("symlog",linthresh=1.)
    axes[1].set_ylabel("U (m/s; symlog)")
    fig.suptitle(r"Phase trajectories over five cycles; $R_0=6$ µm, $f=30$ kHz",fontsize=10)
    finish(fig,"04_phase_portrait")

    fig,ax=plt.subplots(figsize=(6.8,3.4),layout="constrained")
    for bar,color in zip(amps,COLORS):
        d=curves[bar]
        ax.plot(d["t_s"]/T,d["R_m"]/R0,color=color,label=f"{bar:.1f} bar",linewidth=.9)
    ax.set(xlabel=r"Time $t/T$",ylabel=r"$R/R_0$",
           title=r"Acoustic-amplitude comparison; $R_0=6$ µm, $f=30$ kHz")
    ax.legend(ncol=3,loc="upper right")
    finish(fig,"05_amplitude_comparison")

    fig,axes=plt.subplots(1,2,figsize=(6.8,3.2),layout="constrained")
    d=strong
    ta=max(0,tm-float(sr["preceding_expansion_duration_s"]))
    selected=(d["t_s"]>=ta)&(d["t_s"]<=tc+.5e-6)
    axes[0].plot(d["t_s"][selected]*1e6,d["R_m"][selected]/R0,color="black")
    axes[0].axvspan(ta*1e6,tm*1e6,color=COLORS[1],alpha=.12,label="Expansion interval")
    axes[0].axvspan(tm*1e6,tc*1e6,color=COLORS[3],alpha=.12,label="Collapse interval")
    axes[0].set(xlabel="Time (µs)",ylabel=r"$R/R_0$",title="Selected expansion and collapse")
    axes[0].legend(fontsize=7)
    selected=abs(d["t_s"]-tc)<.4e-12
    axes[1].plot((d["t_s"][selected]-tc)*1e12,d["R_m"][selected]/R0,color="black")
    twin=axes[1].twinx()
    twin.plot((d["t_s"][selected]-tc)*1e12,d["U_m_s"][selected],color=COLORS[3])
    twin.set_ylabel("U (m/s)",color=COLORS[3])
    twin.grid(False)
    axes[1].set(xlabel=r"Time from minimum $t-t_c$ (ps)",ylabel=r"$R/R_0$",title="Formal final 0.8 ps window")
    fig.suptitle(r"Formal stress test: $P_a=1.3$ bar, $R_0=6$ µm, $f=30$ kHz",fontsize=10)
    finish(fig,"06_collapse_zoom")

    fig,axes=plt.subplots(1,2,figsize=(6.8,3.2),layout="constrained")
    axes[0].plot(amps,[float(x["Rmax_over_R0"]) for x in summary],"o-",color=COLORS[1],label=r"$R_{max}/R_0$")
    axes[0].plot(amps,[float(x["Rmin_over_R0"]) for x in summary],"s-",color=COLORS[3],label=r"$R_{min}/R_0$")
    axes[0].set(ylabel="Radius ratio")
    axes[1].semilogy(amps[1:],[float(x["Umax_positive_m_s"]) for x in summary[1:]],"o-",color=COLORS[1],label="Maximum outward U")
    axes[1].semilogy(amps[1:],[float(x["Umax_inward_m_s"]) for x in summary[1:]],"s-",color=COLORS[3],label="Maximum inward |U|")
    axes[1].set_ylabel("Wall speed (m/s; log scale)")
    for ax in axes:
        ax.set_xlabel(r"Acoustic amplitude $P_a$ (bar)")
        ax.axvspan(p["p0"]/1e5,1.34,color="gray",alpha=.12)
        ax.legend(fontsize=7)
    fig.suptitle(r"Five-cycle extrema; shading: $p_\infty<0$ during part of a cycle",fontsize=10)
    finish(fig,"07_amplitude_sweep_metrics")

    conv=rows("convergence.csv")
    write_report_tables(meta,summary,conv)
    fig,axes=plt.subplots(1,2,figsize=(6.8,3.1),layout="constrained")
    x=np.arange(4)
    for ax,bar in zip(axes,[1.,1.3]):
        selection=[row for row in conv if float(row["Pa_Pa"]) == bar*1e5]
        for key,label,color in [("Rmax_m",r"$R_{max}$",COLORS[1]),("Rmin_m",r"$R_{min}$",COLORS[2]),
                                ("collapse_time_s",r"$t_c$",COLORS[3]),("max_abs_U_m_s",r"$|U|_{max}$",COLORS[4])]:
            vals=[float(selection[j]["relative_error_"+key]) for j in (0,1,2,4)]
            ax.semilogy(x,np.maximum(vals,1e-16),"o-",color=color,label=label,markersize=3)
        ax.set_xticks(x,["DOP853\n1e-7","DOP853\n1e-9","DOP853\n1e-10","Radau\n1e-10"],fontsize=6)
        ax.set_title(f"$P_a={bar:.1f}$ bar")
        ax.set_ylabel("Relative difference from reference",fontsize=8)
    axes[0].legend(ncol=2,fontsize=7)
    fig.suptitle("Convergence; reference: DOP853 rtol 1e-11, step ≤ T/2000",fontsize=10)
    finish(fig,"08_numerical_convergence")

    fig,axes=plt.subplots(1,2,figsize=(6.8,3.1),layout="constrained")
    for bar,color in zip(amps[1:],COLORS[1:]):
        d=curves[bar]
        axes[0].plot(d["t_s"]/T,d["wall_Mach"],color=color,label=f"{bar:.1f} bar",linewidth=.8)
    axes[0].axhline(.1,color=COLORS[3],ls="--",label="Caution: Mw = 0.1")
    axes[0].axhline(1.,color="gray",ls=":",label="Mw = 1")
    axes[0].set(xlabel=r"Time $t/T$",ylabel=r"Wall Mach number $|U|/c$",yscale="log",ylim=(1e-6,1e3))
    axes[0].legend(fontsize=6,ncol=2,loc="upper left")
    phase=np.linspace(0,1,500)
    for bar,color in zip([1.,1.3],[COLORS[4],"black"]):
        axes[1].plot(phase,(p["p0"]-bar*1e5*np.sin(2*np.pi*phase))/1e5,color=color,label=f"{bar:.1f} bar")
    axes[1].axhline(0,color=COLORS[3],ls="--")
    axes[1].set(xlabel=r"Acoustic phase $t/T$",ylabel=r"Absolute far-field pressure (bar)")
    axes[1].legend()
    fig.suptitle(r"Validity diagnostics; $c=1482$ m/s, $R_0=6$ µm, $f=30$ kHz",fontsize=10)
    finish(fig,"09_validity_diagnostics")
    print("Generated nine figures, each in vector PDF and 300 dpi PNG.")


if __name__ == "__main__":
    main()
