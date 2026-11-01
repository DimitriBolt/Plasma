"""Finite Au sphere + spherical bubble, axisymmetric linear heat equation in SI.

This module solves a FROZEN geometry, not a moving-boundary flow problem.
Temperatures are laser-tagged perturbations with the gas perturbation held zero.
Positive q means liquid -> gas, n_g points from gas into liquid. The export is
a history-dependent small-increment closure: -MR4 loss + delta q, not two full
boundary fluxes. Geometry/optical/thermal admission must be checked by the caller.
No clipping of temperatures, energies, radius, or flux. Contact is rejected.
"""
from dataclasses import dataclass
import json
import numpy as np
from scipy.sparse import coo_matrix, diags, bmat, csc_matrix
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.special import erf
import triangle
from .upstream import optics

MODEL_ID = 'axisymmetric-P1-lumped-mass-frozen-spheres-v1'
SOURCE_CONTRACT_ID = 'frozen-reference-baseline-plus-laser-v2'

@dataclass(frozen=True)
class Geometry:
    radius_m: float = 6e-6
    particle_radius_m: float = 20e-9
    gap_m: float = 20e-9
    outer_padding_m: float = 1e-6

    def __post_init__(self):
        if not all(np.isfinite(x) and x > 0 for x in
                   (self.radius_m, self.particle_radius_m, self.gap_m, self.outer_padding_m)):
            raise ValueError('Positive SI geometry required; d<=0 is contact/intersection')

def particle_trajectory(radius_m, particle_radius_m, initial_gap_m, mode='passive'):
    """r^3-R^3=const follows from 3 r^2 rdot = 3 R^2 Rdot.

    First radius defines the initial placement. 'capture' is an imposed gap,
    requiring an unspecified retaining force, never a passive trajectory.
    Returns center radius, SURFACE gap and explicit contact status at every time.
    """
    R = np.asarray(radius_m, dtype=float)
    if R.ndim != 1 or len(R) == 0 or np.any(R <= 0) or particle_radius_m <= 0 or initial_gap_m <= 0:
        raise ValueError('Positive radius sequence, a and initial gap required')
    initial = R[0]+particle_radius_m+initial_gap_m
    if mode == 'passive':
        r = np.cbrt(R**3+initial**3-R[0]**3)
    elif mode == 'fixed':
        r = np.full_like(R, initial)
    elif mode == 'capture':
        r = R+particle_radius_m+initial_gap_m
    else:
        raise ValueError('Unknown kinematic scenario')
    gap = r-R-particle_radius_m
    return r, gap, np.where(gap <= 0, 'contact_or_intersection', 'separated')

def gaussian_source(t_s, energy_J, center_s, sigma_s, interval_s):
    """Finite-interval normalized diagnostic, not a diffusion kernel. W in SI."""
    a,b = interval_s
    if not (b > a and sigma_s > 0 and energy_J >= 0):
        raise ValueError('Invalid Gaussian normalization')
    z = .5*(erf((b-center_s)/(np.sqrt(2)*sigma_s))-erf((a-center_s)/(np.sqrt(2)*sigma_s)))
    if z <= 1e-14:
        raise ValueError('Gaussian interval has unresolved probability mass')
    t = np.asarray(t_s)
    return np.where((t >= a)&(t <= b), energy_J/(np.sqrt(2*np.pi)*sigma_s*z)
                    *np.exp(-.5*((t-center_s)/sigma_s)**2), 0.)

def point_sphere_fraction(t_s, R_m, center_distance_m, alpha_m2_s):
    """Exact absorbing sphere impulse benchmark in infinite liquid; no finite Au.

    The hitting fraction R/r erfc[(r-R)/sqrt(4 alpha t)] follows from the radial
    heat equation after w=r*T. It is an analytic control, not the main model.
    """
    from scipy.special import erfc
    if not (center_distance_m > R_m > 0 and alpha_m2_s > 0):
        raise ValueError('Source must be outside sphere')
    t = np.asarray(t_s, dtype=float)
    out = np.zeros_like(t)
    positive = t > 0
    out[positive] = R_m/center_distance_m*erfc((center_distance_m-R_m)/np.sqrt(4*alpha_m2_s*t[positive]))
    return out

def point_robin_fraction(t_s,R_m,center_distance_m,alpha_m2_s,k_W_mK,h_W_m2K):
    """Infinite-water point impulse -> finite-R Robin sphere, analytic Laplace inversion.

    For w=r*T, w_r=(1/R+h/k)w at R. The integrated hitting transform is
    (R/r)*(h/k)*exp[-(r-R)*sqrt(s/alpha)]/[s*(sqrt(s/alpha)+1/R+h/k)].
    The erfcx expression avoids cancellation/overflow of exp(beta^2 alpha t).
    """
    from scipy.special import erfc,erfcx
    if not (center_distance_m>R_m>0 and alpha_m2_s>0 and k_W_mK>0 and h_W_m2K>=0):
        raise ValueError('Invalid Robin impulse geometry/properties')
    t=np.asarray(t_s,dtype=float);out=np.zeros_like(t);mask=t>0
    beta=1/R_m+h_W_m2K/k_W_mK
    x=(center_distance_m-R_m)/np.sqrt(4*alpha_m2_s*t[mask])
    out[mask]=R_m/center_distance_m*(h_W_m2K/k_W_mK)/beta*(erfc(x)-np.exp(-x*x)*erfcx(x+beta*np.sqrt(alpha_m2_s*t[mask])))
    return out

def _mesh(geometry, refinement=1., isolated=False, planar=False):
    """Constrained Delaunay mesh of exact sphere boundaries approximated by chords.

    Coordinates normalized by a for robust triangulation. Axis is x=0; sphere
    apex at z=0; full bubble center at z=-R, NP center H=a+d. All azimuths are
    included through the 2*pi*x weight. Outer domain is a finite cylinder.
    """
    g = geometry; a=g.particle_radius_m; R=g.radius_m/a; H=0. if isolated else 1+g.gap_m/a
    pad=g.outer_padding_m/a
    width = pad if (isolated or planar) else R+pad
    bottom = -pad if (isolated or planar) else -2*R-pad
    if planar: bottom=0.
    top = H+1+pad
    vertices=[]; markers=[]
    def add(point, preceding_marker):
        vertices.append(point); markers.append(preceding_marker)
    add([0,bottom],3);add([width,bottom],1 if planar else 2)
    add([width,top],2);add([0,top],2);add([0,H+1],3)
    theta=np.linspace(0,np.pi,int(48*refinement)+1)
    # Near-contact refinement resolves the finite sphere's lower cap.
    scale=min(.4, np.sqrt(g.gap_m/a))
    if not isolated:
        theta=np.unique(np.r_[theta,np.pi-np.geomspace(scale/30,scale, max(6,int(12*refinement)))])
    for th in theta[1:]: add([float(np.sin(th)),float(H+np.cos(th))],4)
    if not isolated and not planar:
        add([0,0],3)
        small=min(.15,1/(R*refinement))
        th=np.unique(np.r_[0.,np.geomspace(small/6,.3,max(12,int(30*refinement))),
                           np.linspace(.3,np.pi,max(12,int(30*refinement)))])
        for angle in th[1:]: add([float(R*np.sin(angle)),float(-R+R*np.cos(angle))],1)
    nboundary=len(vertices)
    segments=[[i,(i+1)%nboundary] for i in range(nboundary)]
    segment_markers=[markers[(i+1)%nboundary] for i in range(nboundary)]
    # Logarithmic rings provide nested local resolution without a uniform vast grid.
    nr=max(12,int(28*refinement)); nt=max(24,int(48*refinement))
    distances=np.geomspace(1.02, max(width,top-H,H-bottom)*1.2, nr)
    for s in distances:
        for th in np.linspace(0,np.pi,nt+1):
            x,z=s*np.sin(th),H+s*np.cos(th)
            if not (x>1e-9 and x<width and bottom<z<top): continue
            if not isolated and not planar and x*x+(z+R)**2 <= R*R*(1+1e-12): continue
            vertices.append([float(x),float(z)])
    result=triangle.triangulate(dict(vertices=np.asarray(vertices),segments=np.asarray(segments),
                                    segment_markers=np.asarray(segment_markers)[:,None]),'pq20Q')
    return (np.asarray(result['vertices'])*a,np.asarray(result['triangles']),
            np.asarray(result['segments']),np.asarray(result['segment_markers']).ravel())

class FrozenHeatModel:
    """Energy-conservative axisymmetric P1 FEM, diagonal exact row-sum capacity.

    Robin Au/water exchange uses boundary quadrature and the SAME conductance
    in both ledgers. Bubble Robin exchange uses fixed zero gas perturbation.
    No advection or moving-domain work is present; geometry is fixed throughout.
    """
    def __init__(self, geometry=Geometry(), h_W_m2K=4400., props=None, refinement=1.,
                 isolated=False, planar=False, outer_bc='bath', interface_W_m2K=None,
                 mesh=None):
        self.geometry=geometry; self.props=optics.ThermalProperties() if props is None else props
        if not np.isfinite(h_W_m2K) or h_W_m2K < 0 or refinement <= 0 or outer_bc not in ('bath','insulated'):
            raise ValueError('Nonnegative h, positive refinement, bath/insulated required')
        self.h=h_W_m2K; self.isolated=isolated; self.planar=planar
        self.G=self.props.interface_W_m2K if interface_W_m2K is None else interface_W_m2K
        if not np.isfinite(self.G) or self.G<0: raise ValueError('Nonnegative Au-water G required')
        if mesh is None:
            xy,tri,edges,markers=_mesh(geometry,refinement,isolated,planar)
        else:
            # An archived mesh is an explicit numerical input. No remesher or
            # cross-machine assertion about Triangle's node choices is involved.
            xy,tri,edges,markers=[np.asarray(x).copy() for x in mesh]
            if (xy.ndim!=2 or xy.shape[1]!=2 or not np.isfinite(xy).all()
                    or np.any(xy[:,0]<0) or tri.ndim!=2 or tri.shape[1]!=3
                    or edges.ndim!=2 or edges.shape[1]!=2 or markers.shape!=(len(edges),)
                    or tri.dtype.kind not in 'iu' or edges.dtype.kind not in 'iu'
                    or np.any(tri<0) or np.any(tri>=len(xy))
                    or np.any(edges<0) or np.any(edges>=len(xy))
                    or not np.isin(markers,[1,2,3,4]).all()):
                raise ValueError('Invalid archived meridional mesh')
        self.xy,self.tri,self.edges,self.markers=xy,tri,edges,markers
        n=len(xy); points=xy[tri]; r=points[:,:,0]
        u=points[:,1]-points[:,0];v=points[:,2]-points[:,0]
        twice=u[:,0]*v[:,1]-u[:,1]*v[:,0]
        if np.any(twice==0): raise ValueError('Degenerate archived triangle')
        area=np.abs(twice)/2
        gradients=np.stack([np.roll(points[:,:,1],-1,axis=1)-np.roll(points[:,:,1],-2,axis=1),
                            np.roll(points[:,:,0],-2,axis=1)-np.roll(points[:,:,0],-1,axis=1)],axis=-1)/twice[:,None,None]
        local=self.props.k_water_W_mK*2*np.pi*area[:,None,None]*r.mean(axis=1)[:,None,None]*np.einsum('tik,tjk->tij',gradients,gradients)
        rows=np.broadcast_to(tri[:,:,None],local.shape).ravel(); cols=np.broadcast_to(tri[:,None,:],local.shape).ravel()
        K=coo_matrix((local.ravel(),(rows,cols)),shape=(n+1,n+1)).tocsc()
        capacity=np.zeros(n+1)
        for j in range(3):
            np.add.at(capacity,tri[:,j],self.props.rho_water_kg_m3*self.props.cp_water_J_kgK*2*np.pi*area*(r.sum(axis=1)+r[:,j])/12)
        capacity[-1]=optics.thermal_scales(geometry.particle_radius_m,self.props)['capacity_J_K']
        self.full_capacity=capacity
        surface={tag:np.zeros(n+1) for tag in (1,2,4)}
        for edge,tag in zip(edges,markers):
            if tag not in surface: continue
            length=np.linalg.norm(xy[edge[1]]-xy[edge[0]])
            rr=xy[edge,0]
            surface[tag][edge]+=2*np.pi*length*(rr.sum()+rr)/6
        self.au_area=surface[4].sum(); self.bubble_area=surface[1].sum()
        au=surface[4]*self.G
        # Off-diagonal links connect water nodes to the one gold node.
        ai=np.flatnonzero(au)
        K=K+diags(au+surface[1]*self.h,format='csc')
        K=K+coo_matrix((np.r_[-au[ai],-au[ai],au.sum()],
                       (np.r_[ai,np.full(len(ai),n),n],np.r_[np.full(len(ai),n),ai,n])),shape=K.shape).tocsc()
        bath_nodes=np.unique(edges[markers==2]) if outer_bc=='bath' else np.array([],dtype=int)
        self.active=np.setdiff1d(np.arange(n+1),bath_nodes)
        self.capacity=capacity[self.active]
        self.K=K[self.active][:,self.active].tocsc()
        self.bubble_g=(surface[1]*self.h)[self.active]
        # Sum of diffusion rows to eliminated fixed-temperature boundary nodes.
        self.bath_g=np.asarray(K[self.active][:,bath_nodes].sum(axis=1)).ravel()*-1 if len(bath_nodes) else np.zeros(len(self.active))
        self.au_g=au[self.active]
        self.gold_index=int(np.flatnonzero(self.active==n)[0])
        self.A=-diags(1/self.capacity)@self.K
        self.forcing=np.zeros(len(self.active));self.forcing[self.gold_index]=1/self.capacity[self.gold_index]
        self.outer_bc=outer_bc
        self.wall_weights=surface[1]

    @classmethod
    def from_saved_field(cls,path,geometry=None,props=None):
        """Reassemble the heat operator on the archived mesh, in SI.

        Field archives already contain nodes, triangles, boundary edges/markers
        and h. Geometry is read from metadata or explicitly supplied for the
        original reference archive, which predates geometry metadata. Material
        properties must match the source run; defaults are the MR6 reference.
        """
        with np.load(path,allow_pickle=False) as z:
            meta=json.loads(str(z['metadata_json'])) if 'metadata_json' in z else {}
            recorded=Geometry(**meta['geometry']) if 'geometry' in meta else None
            if geometry is None: geometry=recorded
            if geometry is None: raise ValueError('Reference field requires explicit Geometry')
            if recorded is not None and geometry!=recorded:
                raise ValueError('Geometry differs from archived field metadata')
            mesh=tuple(z[key].copy() for key in ('xy_m','triangles','boundary_edges','boundary_markers'))
            h=float(z['h_W_m2K'])
        return cls(geometry,h,props,mesh=mesh)

    def pulse(self, absorbed_energy_J, duration_s, time_s, rtol=2e-7, method='BDF', max_step_s=np.inf, save_field=False):
        """Start a top-hat at t=0. Return temperatures [K], signed flux [W], J.

        Integrated gas and bath fluxes are independent ODE states. The energy
        residual never defines a missing channel. Time domain has no extrapolation.
        """
        t=np.asarray(time_s,dtype=float)
        if not (np.all(np.isfinite(t)) and len(t)>1 and t[0]>=0 and np.all(np.diff(t)>0)
                and absorbed_energy_J>=0 and duration_s>0): raise ValueError('Invalid pulse/time domain')
        n=len(self.capacity); unit=1e-15
        # Ledger states are scaled J/unit; integration time is scaled to ns.
        ledger=np.vstack([self.bubble_g,self.bath_g])/unit
        operator=bmat([[self.A,None],[csc_matrix(ledger),csc_matrix((2,2))]],format='csc')*1e-9
        laser=np.r_[self.forcing,np.zeros(2)]*absorbed_energy_J/duration_s*1e-9
        initial=np.zeros(n+2); solutions=[];start=0.
        for end,on in ((min(duration_s,t[-1]),True),(t[-1],False)):
            if end<=start: continue
            sol=solve_ivp(lambda u,y:operator@y+(laser if on else 0.),(start/1e-9,end/1e-9),initial,
                          method=method,jac=operator,rtol=rtol,atol=rtol*1e-4,dense_output=True,
                          max_step=max_step_s/1e-9)
            if not sol.success: raise RuntimeError(sol.message)
            solutions.append((start,end,sol));initial=sol.y[:,-1];start=end
        states=np.empty((n+2,len(t)))
        for start,end,sol in solutions:
            mask=(t>=start)&(t<=end);states[:,mask]=sol.sol(t[mask]/1e-9)
        temp=states[:n]; egold=self.capacity[self.gold_index]*temp[self.gold_index]
        ein=absorbed_energy_J*np.minimum(t/duration_s,1.)
        ewater=self.capacity@temp-egold
        qgas=self.bubble_g@temp; qbath=self.bath_g@temp
        full=np.zeros((len(self.full_capacity),len(t)));full[self.active]=temp
        au_weights=self.au_g.copy();au_weights[self.gold_index]=0
        water_nearest=np.argmin(np.linalg.norm(self.xy-np.array([0.,self.geometry.gap_m]),axis=1))
        out=dict(time_s=t,gold_rise_K=temp[self.gold_index],nearest_water_rise_K=full[water_nearest].copy(),
                 mean_Au_surface_water_rise_K=(au_weights@temp/au_weights.sum() if au_weights.sum()>0 else np.zeros_like(t)),
                 q_wall_laser_increment_W=qgas,q_bath_W=qbath,absorbed_so_far_J=ein,
                 gold_energy_J=egold,water_energy_J=ewater,gas_energy_J=states[n]*unit,
                 bath_energy_J=states[n+1]*unit,
                 balance_residual_J=ein-egold-ewater-states[n]*unit-states[n+1]*unit,
                 minimum_water_rise_K=np.min(full[:-1],axis=0))
        if save_field:
            out.update(water_rise_K=full[:-1],xy_m=self.xy,triangles=self.tri,
                       boundary_edges=self.edges,boundary_markers=self.markers,
                       wall_node_area_m2=self.wall_weights[:-1],h_W_m2K=np.array(self.h))
        return out

class ThermalHistory:
    """Stored reference-state diagnostic and frozen laser increment, SI.

    The baseline belongs to the archived reference gas-temperature trace. Its
    sum with the laser increment is NOT the full physical flux for another gas
    state. No q_wall_total_W is exported. An ODE with its own current-state MR4
    loss must add ONLY laser_increment_W(), exactly once. The liquid response
    remains frozen with zero gas-temperature perturbation; no MR7 model is implied.
    """
    def __init__(self, time_s, laser_W, baseline_W, metadata=None,reference_temperature_K=None,
                 reference_radius_m=None):
        self.t=np.asarray(time_s,dtype=float)
        if len(self.t)<2 or np.any(np.diff(self.t)<=0): raise ValueError('Strictly ordered times required')
        self.laser=np.asarray(laser_W,dtype=float);self.baseline=np.asarray(baseline_W,dtype=float)
        if self.laser.shape!=self.t.shape or self.baseline.shape!=self.t.shape or not np.isfinite(np.r_[self.t,self.laser,self.baseline]).all(): raise ValueError('Finite matching histories required')
        self.metadata={} if metadata is None else dict(metadata)
        self.metadata['source_contract_id']=SOURCE_CONTRACT_ID
        self.metadata['baseline_scope']='stored_reference_state_only'
        self.metadata['gas_ode_addition']='q_wall_laser_increment_W_only'
        self.reference_temperature=None if reference_temperature_K is None else np.asarray(reference_temperature_K,dtype=float)
        if self.reference_temperature is not None and (self.reference_temperature.shape!=self.t.shape
                or not np.isfinite(self.reference_temperature).all() or np.any(self.reference_temperature<=0)):
            raise ValueError('Finite positive reference gas temperatures must match time grid')
        self.temperature_function=(None if self.reference_temperature is None else
                                   PchipInterpolator(self.t,self.reference_temperature,extrapolate=False))
        self.reference_radius=None if reference_radius_m is None else np.asarray(reference_radius_m,dtype=float)
        if self.reference_radius is not None and (self.reference_radius.shape!=self.t.shape
                or not np.isfinite(self.reference_radius).all() or np.any(self.reference_radius<=0)):
            raise ValueError('Finite positive reference radii must match time grid')
        self.radius_function=(None if self.reference_radius is None else
                              PchipInterpolator(self.t,self.reference_radius,extrapolate=False))
        self.functions=[PchipInterpolator(self.t,q,extrapolate=False) for q in (self.laser,self.baseline)]

    def thermal_source(self,t_s,bubble_state=None):
        if not np.isfinite(t_s) or t_s<self.t[0] or t_s>self.t[-1]: raise ValueError('Outside serialized thermal history; extrapolation prohibited')
        if bubble_state is not None and 'frozen_radius_m' in self.metadata:
            if not np.isclose(bubble_state['R_m'],self.metadata['frozen_radius_m'],rtol=self.metadata.get('radius_tolerance',.01),atol=0):
                raise ValueError('Geometry differs from precomputed response; solve stateful liquid equation')
        laser,base=[float(f(t_s)) for f in self.functions];reference_sum=laser+base
        return dict(q_wall_reference_sum_W=reference_sum,q_wall_reference_baseline_W=base,
                    q_wall_laser_increment_W=laser,q_wall_reference_positive_W=max(reference_sum,0.),
                    q_wall_reference_negative_W=min(reference_sum,0.),
                    reference_gas_temperature_K=(None if self.temperature_function is None else float(self.temperature_function(t_s))),
                    reference_bubble_radius_m=(None if self.radius_function is None else float(self.radius_function(t_s))),
                    current_state_total_flux_available=False,source_contract_id=SOURCE_CONTRACT_ID,
                    model_id=MODEL_ID,status=self.metadata.get('status','conditional_frozen_geometry'))

    def laser_increment_W(self,t_s,bubble_state=None):
        """The only additive term for a gas ODE already evaluating its heat loss.

        Does not re-evaluate the loss at a supplied T_g_K or include gas feedback
        into the liquid. Geometry/time restrictions remain those of this history.
        """
        return self.thermal_source(t_s,bubble_state)['q_wall_laser_increment_W']

    def save(self,path):
        arrays=dict(time_s=self.t,q_wall_laser_increment_W=self.laser,
                    q_wall_reference_baseline_W=self.baseline,q_wall_reference_sum_W=self.laser+self.baseline,
                    metadata_json=np.array(json.dumps(self.metadata,sort_keys=True)))
        if self.reference_temperature is not None:arrays['reference_gas_temperature_K']=self.reference_temperature
        if self.reference_radius is not None:arrays['reference_bubble_radius_m']=self.reference_radius
        np.savez_compressed(path,**arrays)

    @classmethod
    def load(cls,path):
        with np.load(path,allow_pickle=False) as z:
            key='q_wall_reference_baseline_W' if 'q_wall_reference_baseline_W' in z else 'q_wall_baseline_W'
            meta=json.loads(str(z['metadata_json']))
            if key=='q_wall_baseline_W':meta['legacy_reference_fields_loaded']=True
            T=z['reference_gas_temperature_K'] if 'reference_gas_temperature_K' in z else None
            radius=z['reference_bubble_radius_m'] if 'reference_bubble_radius_m' in z else None
            return cls(z['time_s'],z['q_wall_laser_increment_W'],z[key],meta,T,radius)

def thermal_source(t_s,history,bubble_state=None):
    """Stateless facade over an explicitly supplied, history-bearing object."""
    return history.thermal_source(t_s,bubble_state)

def integrated_parts(time_s,flux_W):
    t=np.asarray(time_s);q=np.asarray(flux_W)
    return dict(signed_J=float(np.trapezoid(q,t)),positive_J=float(np.trapezoid(np.maximum(q,0),t)),
                negative_J=float(np.trapezoid(np.minimum(q,0),t)))
