###############################################################################
# Some functions for looking at sampels of cooked planets
###############################################################################

import sys, os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import TOFPlanet
import ahelpers as ah

def load_planets(fname):
    with open(fname, 'rb') as f:
        planets = pickle.load(f)
        print(f"Found {len(planets)} planets in {fname}.")
    return planets

def _get_canvas(newfig):
    if newfig or not plt.get_fignums():
        plt.figure(figsize=(8,6))
    else:
        plt.figure(plt.gcf().number)

def hist_mass(planets, newfig=False, bins='auto', density=True, **kwargs):
    # Prepare the data
    planets = load_planets(planets) if type(planets) is str else planets
    M = np.array([p.M for p in planets])/1e24

    # Prepare the canvas
    _get_canvas(newfig)

    # Plot the histogram
    plt.hist(M, bins=bins, density=density, **kwargs)

    # Style, annotate, and show
    plt.xlabel(r'Integrated mass $M$ [$10^{24}$ kg]')
    # plt.xlim(M.mean()-2*M.std(),M.mean()+2*M.std())
    plt.yticks([])
    plt.show(block=False)

def hist_params(planets,rol=False,dims='all',trans=None,**kwargs):
    # Prepare the data
    planets = np.loadtxt(planets) if type(planets) is str else planets
    try:
        C = np.array([p._params for p in planets])
    except AttributeError:
        print("Planet object did not store sample parameter.")
        return
    if trans is not None:
        C = np.array([trans(row) for row in C])
    # Inspect histograms
    if dims == 'all':
        dims = tuple(range(C.shape[-1]))
    if type(dims) is int:
        dims = [dims]
    for dim in dims:
        x = C[:,dim]
        if rol:
            x = x[np.abs(x - x.mean()) < 3*x.std()]
        plt.figure()
        plt.hist(x, density=True, **kwargs)
        plt.ylabel('x{}'.format(dim))
        plt.annotate('mean={:g}'.format(
            np.mean(x)), (0.01,0.9), xycoords='axes fraction');
        plt.show(block=False)
        pass

def hist_moi(planets,rol=False,newfig=False,density=True,show=True,**kwargs):
    # Prepare the data
    planets = load_planets(planets) if type(planets) is str else planets
    ice = np.array([p.NMoI for p in planets])
    mu, sig = ice.mean(), ice.std()
    if rol:
        ice = ice[np.abs(ice - mu) < 3*sig]

    # Prepare the canvas
    _get_canvas(newfig)

    # Plot the histogram
    plt.hist(ice,density=density,**kwargs)
    plt.vlines(mu, *plt.ylim(), ls='--', color='r', label=f"{mu:.4f}")

    # Style, annotate, and show
    plt.xlabel(r'Normalized moment of inertia, $I/Ma_0^2$')
    plt.yticks([])
    plt.legend()
    if show:
        plt.show(block=False)

def hist_rho_c(planets, newfig=False, bins='auto', density=True, **kwargs):
    # Prepare the data
    planets = load_planets(planets) if type(planets) is str else planets
    rcs = np.array([p.rhoi[-1] for p in planets])

    # Prepare the canvas
    _get_canvas(newfig)

    # Plot the histogram
    plt.hist(rcs, bins=bins, density=density, **kwargs)

    # Style, annotate, and show
    plt.xlabel(r'$\rho_c$ [1000 kg/m$^3$]')
    plt.yticks([])
    plt.show(block=False)

def hist_rho0(planets, newfig=False, bins='auto', density=True, **kwargs):
    # Prepare the data
    planets = load_planets(planets) if type(planets) is str else planets
    r0s = np.array([p.rhoi[0] for p in planets])

    # Prepare the canvas
    _get_canvas(newfig)

    # Plot the histogram
    plt.hist(r0s, bins=bins, density=density, **kwargs)

    # Style, annotate, and show
    plt.xlabel(r'$\rho_0$ [1000 kg/m$^3$]')
    plt.yticks([])
    plt.show(block=False)

def ensemble_of_profs(planets, newfig=False, nlines=20, alfa=0.4, **kwargs):
    # Prepare the data
    planets = load_planets(planets) if type(planets) is str else planets
    profs = np.array([p.rhoi for p in planets]).T
    rcs = profs[-1,:]
    ind = np.argsort(rcs)
    profs = profs[:,ind]

    # Prepare the canvas
    _get_canvas(newfig)

    # Plot the lines
    x = planets[0].si/planets[0].s0
    x = np.append(x, 0)
    skip = int(np.ceil(len(planets)/nlines))
    for k in range(0,len(planets),skip):
        y = profs[:,k]
        y = np.append(y,y[-1])
        plt.plot(x,y,alpha=alfa,**kwargs)

    # Style, annotate, and show
    plt.xlabel(r'Level surface radius, $s/R_m$')
    plt.ylabel(r'$\rho$ [1000 kg/m$^3$]')
    plt.show(block=False)

def ensemble_of_barotropes(planets, newfig=False, nlines=20, alfa=0.4, **kwargs):
    # Prepare the data
    planets = load_planets(planets) if type(planets) is str else planets
    pees = 1e-11*np.array([p.Pi for p in planets]).T
    rhos = 1e-3*np.array([p.rhoi for p in planets]).T

    rcs = rhos[-1,:]
    ind = np.argsort(rcs)
    pees = pees[:,ind]
    rhos = rhos[:,ind]

    # Prepare the canvas
    _get_canvas(newfig)

    # Plot the lines
    skip = int(np.ceil(len(planets)/nlines))
    for k in range(0,len(planets),skip):
        x = pees[:,k]
        y = rhos[:,k]
        plt.plot(x,y,alpha=alfa,**kwargs)

    # Style, annotate, and show
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(left=1e-5)
    plt.ylim(1e-3, 100)
    plt.xlabel(r'$p$ [Mbar]')
    plt.ylabel(r'$\rho$ [1000 kg/m$^3$]')
    plt.show(block=False)

def density_envelope(planets, newfig=False, prctile=2, show=True, **kwargs):
    # Prepare the data
    planets = load_planets(planets) if type(planets) is str else planets
    profs = np.array([p.rhoi for p in planets]).T
    prcs_lo = prctile
    prcs_hi = 100 - prcs_lo
    x = planets[0].si/planets[0].s0
    ylo = np.percentile(profs, prcs_lo, axis=1)/1000
    yhi = np.percentile(profs, prcs_hi, axis=1)/1000

    # Prepare the canvas
    _get_canvas(newfig)

    # Plot the shaded regions
    plt.fill_between(x, ylo, yhi, **kwargs)

    # Style, annotate, and show
    plt.xlabel(r'Level surface radius, $s/R_m$')
    plt.ylabel(r'$\rho$ [1000 kg/m$^3$]')
    plt.ylim((0,plt.ylim()[1]))
    if show:
        plt.show(block=False)

def barotrope_envelope(planets, newfig=False, prctile=2, show=True, **kwargs):
    from scipy.interpolate import interp1d
    # Prepare the data
    planets = load_planets(planets) if type(planets) is str else planets
    pees = 1e-11*np.array([p.Pi for p in planets]).T
    rhos = 1e-3*np.array([p.rhoi for p in planets]).T
    x = np.logspace(-6, np.ceil(np.log10(pees.max())), 1024)
    y = np.nan*np.ones((1024, len(planets)))
    for k in range(len(planets)):
        rhoofp = interp1d(pees[:,k], rhos[:,k], kind='cubic')
        ind = x < pees[-1,k]
        y[ind,k] = rhoofp(x[ind])

    prcs_lo = prctile
    prcs_hi = 100 - prcs_lo
    ylo = np.nanpercentile(y, prcs_lo, axis=1)
    yhi = np.nanpercentile(y, prcs_hi, axis=1)

    # Prepare the canvas
    _get_canvas(newfig)

    # Plot the shaded regions
    ind = ~np.isnan(ylo)
    plt.fill_between(x[ind], ylo[ind], yhi[ind], **kwargs)

    # Style, annotate, and show
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(left=1e-5)
    plt.ylim(1e-3, 100)
    plt.vlines(pees[-1,:].min(), *plt.ylim(), linestyle=':', color='r')
    plt.xlabel(r'$p$ [Mbar]')
    plt.ylabel(r'$\rho$ [1000 kg/m$^3$]')
    if show:
        plt.show(block=False)

def plot_profile(s, rho, newfig=False, **kwargs):
    # Prepare the data
    x = np.append(s, 0)/s[0]
    y = np.append(rho, rho[-1])/1000

    # Prepare the canvas
    _get_canvas(newfig)

    # Plot the lines
    plt.plot(x, y, **kwargs)

    # Style, annotate, and show
    plt.xlabel(r'Level surface radius, $s/R_m$')
    plt.ylabel(r'$\rho$ [1000 kg/m$^3$]')
    plt.show(block=False)

def hist_loss(planets,obs,Jmax,newfig=False,bins='auto',density=True,**kwargs):
    # Prepare the data
    planets = load_planets(planets) if type(planets) is str else planets
    L = ah.lossify_planets(planets,obs,Jmax)

    # Prepare the canvas
    _get_canvas(newfig)

    # Plot the histogram
    plt.hist(L, bins=bins, density=density, **kwargs)

    # Style, annotate, and show
    plt.xlabel('Loss function value')
    plt.yticks([])
    plt.show(block=False)

def hist_mass_err(planets,obs,newfig=False,bins='auto',density=True,**kwargs):
    # Prepare the data
    planets = load_planets(planets) if type(planets) is str else planets
    M_err = (np.array([p.M for p in planets]) - obs.M)/obs.dM

    # Prepare the canvas
    _get_canvas(newfig)

    # Plot the histogram
    plt.hist(M_err, bins=bins, density=density, **kwargs)

    # Style, annotate, and show
    plt.xlabel(r'$(M-M^\mathrm{obs})/\sigma{M}^\mathrm{obs}$')
    # plt.xlim(M.mean()-2*M.std(),M.mean()+2*M.std())
    plt.yticks([])
    plt.show(block=False)

def hist_J(planets, n, newfig=False, bins='auto', density=True, **kwargs):
    # Prepare the data
    planets = load_planets(planets) if type(planets) is str else planets
    J = 1e6*np.array([p.Js[n//2] for p in planets])

    # Prepare the canvas
    _get_canvas(newfig)

    # Plot the histogram
    plt.hist(J, bins=bins, density=density, **kwargs)

    # Style, annotate, and show
    plt.xlabel(f'$J_{n}$X$10^6$')
    plt.yticks([])
    plt.show(block=False)

def hist_J_err(planets,n,obs,newfig=False,bins='auto',density=True,**kwargs):
    # Prepare the data
    planets = load_planets(planets) if type(planets) is str else planets
    J = np.array([p.Js[n//2] for p in planets])
    J_err = 1e6*(J - obs.Js[n//2])

    # Prepare the canvas
    _get_canvas(newfig)

    # Plot the histogram
    plt.hist(J_err, bins=bins, density=density, **kwargs)

    # Style, annotate, and show
    plt.xlabel(f'$(J_{n}-J_{n}^*)$X$10^6$')
    # plt.xlim(M.mean()-2*M.std(),M.mean()+2*M.std())
    plt.yticks([])
    plt.show(block=False)

def hist_rotation(planets,obs=None,newfig=False,bins='auto',density=True,**kwargs):
    # Prepare the data
    planets = load_planets(planets) if type(planets) is str else planets
    s0 = np.array([p.s0 for p in planets])
    mrot = np.array([p.mrot for p in planets])
    GM = np.array([p.GM for p in planets])
    w_rot = np.sqrt(mrot*GM/s0**3)
    P_rot = 2*np.pi/w_rot
    if obs is None:
        P_fid = P_rot.mean()
    else:
        P_fid = obs.P

    # Prepare the canvas
    _get_canvas(newfig)

    # Plot the histogram
    plt.hist(P_rot - P_fid, bins=bins, density=density, **kwargs)

    # Style, annotate, and show
    fid_h = int(P_fid/3600)
    fid_m = int((P_fid%3600)/60)
    fid_s = np.round((P_fid%3600)%60,1)
    plt.xlabel(f'rotation period  - {fid_h}:{fid_m}:{fid_s}')
    plt.yticks([])
    plt.show(block=False)

def hist_attr(planets, attr,
    obs=None,newfig=False,bins='auto',density=True,**kwargs):
    # Prepare the data
    planets = load_planets(planets) if type(planets) is str else planets
    vals = np.array([getattr(p,attr,np.nan) for p in planets])
    if obs is None:
        v_mu = vals.mean()
        v_sig = vals.std()
    else:
        v_mu = getattr(obs,attr)
        v_sig = getattr(obs, 'd'+attr)

    # Prepare the canvas
    _get_canvas(newfig)

    # Plot the histogram
    plt.hist(vals, bins=bins, density=density, **kwargs)
    plt.vlines(v_mu, *plt.ylim(), ls='--', lw=2, color='r')
    plt.vlines([v_mu-v_sig,v_mu+v_sig], *plt.ylim(), ls=':', lw=1, color='k')

    # Style, annotate, and show
    plt.xlabel(attr)
    plt.yticks([])
    plt.show(block=False)

def k2_v_I(planets,newfig=False,rol=False,show=True,**kwargs):
    # Prepare the data
    planets = load_planets(planets) if type(planets) is str else planets
    I = np.array([p.NMoI for p in planets])
    K = np.array([p.k2 for p in planets])

    # Prepare the canvas
    _get_canvas(newfig)

    # Plot
    plt.plot(I, K, '+', **kwargs)

    # Style, annotate, and show
    plt.xlabel(r'Normalized MoI, $C/MR^2$')
    plt.ylabel(r'Tidal $k_2$')
    if show:
        plt.show(block=False)

def I_v_k2(planets,newfig=False,rol=False,show=True,**kwargs):
    # Prepare the data
    planets = load_planets(planets) if type(planets) is str else planets
    I = np.array([p.NMoI for p in planets])
    K = np.array([p.k2 for p in planets])

    # Prepare the canvas
    _get_canvas(newfig)

    # Plot
    plt.plot(K, I, '+', **kwargs)

    # Style, annotate, and show
    plt.xlabel(r'Tidal $k_2$')
    plt.ylabel(r'Normalized MoI, $C/MR^2$')
    if show:
        plt.show(block=False)

def dI_v_dk2(planetss,newfig=False,rol=False,show=True,**kwargs):
    # Prepare the data
    dI = []
    dK = []
    for planets in planetss:
        I = np.array([p.NMoI for p in planets])
        K = np.array([p.k2 for p in planets])
        dI.append(I.std())
        dK.append(K.std())

    # Prepare the canvas
    _get_canvas(newfig)

    # Plot
    plt.plot(dK, dI, '+', **kwargs)

    # Style, annotate, and show
    plt.xscale('log')
    # plt.xticks(dK, [f"{x:0.3e}" for x in dK])
    plt.xticks([],minor=True)
    plt.xlabel(r'$\sigma\left(k_2\right)$')
    plt.ylabel(r'$\sigma\left(C/MR^2\right)$')
    if show:
        plt.show(block=False)

    if True:
        pass
