################################################################################
# Some helper functions
################################################################################

import sys, os
import numpy as np
import TOFPlanet
import tof4, tof7
from scipy.stats import norm, chi2
cout = sys.stdout.write

def cook_ppbs_planet(x, obs, mdl, transform, **kwargs):
    """Create a planet object from sample-space parameters."""
    opts = {**{'fix_rot':True,
               'fix_mass':True,
               'no_spin':False,
               'with_k2':False,
               'toforder':4,
               'toflevels':4096,
               'xlevels':128,
               'savess':False},
            **kwargs}

    if opts['fix_rot']:
        Prot = obs.P
        x = x
    else:
        Prot = x[-1]*obs.dP/2 + obs.P
        x = x[:-1]
    y = transform(x)
    tp = mdl(opts['toflevels'], y, obs, opts['toforder'], opts['xlevels'])
    if opts['no_spin']:
        tp.period = np.inf
    tp.relax_to_barotrope(fixmass=opts['fix_mass'])

    tofun = tof4.tof4 if opts['toforder'] == 4 else tof7.tof7
    _, out = tofun(tp.si,tp.rhoi,tp.mrot,xlevels=opts['xlevels'],calc_moi=True)
    tp.NMoI = out.NMoI

    if opts['with_k2']:
        tp.k2 = lovek2(tp.si, tp.rhoi)

    if not opts['savess']:
        tp.ss = None
        tp.SS = None

    return tp

def cook_ppwd_planet(x, obs, mdl, transform, **kwargs):
    """Create a planet object from sample-space parameters."""
    opts = {**{'fix_rot':True,
               'fix_mass':True,
               'preserve_period':True,
               'no_spin':False,
               'with_k2':False,
               'toforder':4,
               'toflevels':4096,
               'xlevels':256,
               'savess':False},
            **kwargs}

    if opts['fix_rot']:
        Prot = obs.P
        x = x
    else:
        Prot = x[0]*obs.dP/2 + obs.P
        x = x[1:]
    y = transform(x, obs)
    svec, dvec = mdl(opts['toflevels'], y, obs.rho0)
    tp = TOFPlanet.TOFPlanet(obs)
    tp.si = svec*obs.s0
    tp.rhoi = dvec
    tp.period = Prot
    tp.opts['toforder'] = opts['toforder']
    tp.opts['xlevels'] = opts['xlevels']
    if opts['no_spin']:
        tp.period = np.inf
    if opts['preserve_period']:
        tp.relax_to_rotation(opts['fix_mass'])
    else:
        tp.mrot = (2*np.pi/tp.period)**2*tp.s0**3/tp.GM
    tp.relax_to_HE(fixmass=opts['fix_mass'],moi=True,pressure=True)

    if opts['with_k2']:
        tp.k2 = lovek2(tp.si, tp.rhoi)
    if not opts['savess']:
        tp.ss = None
        tp.SS = None
    return tp

def _ecdf(a):
    x, counts = np.unique(a, return_counts=True)
    return x, np.cumsum(counts)/len(a)

def lossify_planets(planets, obs, Jmax, rho0=True, mass=False, period=False):
    import losses
    J_ord = [2*n for n in range(1,int(Jmax/2)+1)]
    L = []
    for p in planets:
        el = losses.euclid_Jnm(p.Js, obs, J_ord)
        if rho0:
            el = np.sqrt(el**2 + losses.rho0((p.si,p.rhoi),obs)**2)
        if mass:
            el = np.sqrt(el**2 + losses.mass((p.si,p.rhoi),obs)**2)
        if period:
            mu, sig = obs.P, obs.dP/2
            el = np.sqrt(el**2 + ((p.period - mu)/sig)**2)
        L.append(el)
    return np.array(L)

def load_planets(fname):
    import pickle
    import TOFPlanet
    with open(fname, 'rb') as f:
        planets = pickle.load(f)
        print(f"Found {len(planets)} planets in {fname}.")
    return planets

def load_profiles(fname):
    planets = load_planets(fname)
    profs = np.array([p.rhoi for p in planets]).T
    return profs

def density_prctiles(profs, prcs):
    """High and low density values of symmetric percentile in sample."""
    prcs_lo = prcs
    prcs_hi = 100 - prcs_lo
    YLO = np.percentile(profs, prcs_lo, axis=1)
    YHI = np.percentile(profs, prcs_hi, axis=1)
    return (YLO, YHI)

def _sig2mahal(sig,p):
    """Return Mahalanobis distance in p dims equivalent to Z-score sig."""
    return np.sqrt(chi2.ppf(norm.cdf(sig) - norm.cdf(-sig), p))

def _mahal2sig(d, p):
    """Return z-score equivalent of Mahalanobis distance d in p dims."""
    return -norm.ppf((1 - chi2.cdf(d**2,p))/2)

def winners(C,L,df,prior,rotflag,sig=2.5):
    """Find walkers with endstate likelihood above sigma-equivalent."""

    L = L[:,-1]
    X = C[:,-1,:] if rotflag else C[:,-1,1:]
    P = [prior(x) for x in X]
    thresh = -0.5*_sig2mahal(sig, df)**2
    wins = np.squeeze(np.argwhere((L - P) > thresh))
    return wins

def filter_walker_acceptance(fname,t=1.0):
    """Return logical array where walker ar > than mu - t*sig."""
    C,_ = load_chain(fname)
    A = acceptance_ratio(C)
    return A > (A.mean() - t*A.std())

def _find_outwalkers(chain,thresh=3.0):
    """Try to detect outlier walkers in ensemble sampler chain."""
    C = chain
    t = thresh
    outwalkers = np.zeros(C.shape[0]) > 0

    # Mark with shame the lowest performers
    ars = acceptance_ratio(C)
    outwalkers = outwalkers | ((ars - ars.mean())**2 > t*ars.var())

    # Then check for odd parameter values
    means = [np.mean(C[:,:,dim], axis=1) for dim in range(C.shape[-1])]
    mus = [m.mean() for m in means]
    sigs = [m.var() for m in means]
    for dim in range(C.shape[-1]):
        outwalkers = outwalkers | ((means[dim] - mus[dim])**2 > t*sigs[dim])

    return outwalkers

def mass_integral(svec, dvec):
    """Return approximate mass integral."""
    from scipy.integrate import trapz
    return -4*np.pi*trapz(dvec*svec**2, x=svec)

def lovek2(zvec, dvec):
    """ Tidal Love number k2 from density profile."""
    if zvec[1] - zvec[0] < 0:
        zvec = zvec[::-1]
        dvec = dvec[::-1]
    zvec = zvec/zvec[-1]
    m = dvec[0]*zvec[0]**3; # starting mass
    rhom = m/zvec[0]**3; # starting mean density
    eta = 0
    for k in range(zvec.size - 1):
        s1 = (6 - 6*(dvec[k]/rhom)*(eta + 1) + eta - eta**2)/zvec[k];
        zhalf = zvec[k] + 0.5*(zvec[k+1] - zvec[k]);
        dhalf = dvec[k] + 0.5*(dvec[k+1] - dvec[k]);
        mhalf = m + dhalf*(zhalf**3 - zvec[k]**3);
        rhalf = mhalf/zhalf**3;
        ehalf = eta + s1*(zhalf - zvec[k]);
        s2 = (6 - 6*(dhalf/rhalf)*(ehalf + 1) + ehalf - ehalf**2)/zhalf;
        eta = eta + s2*(zvec[k+1] - zvec[k]);
        m = mhalf + dvec[k+1]*(zvec[k+1]**3 - zhalf**3);
        rhom = m/zvec[k+1]**3;
    k2 = (3 - eta)/(2 + eta)
    return k2

def load_chain(fname):
    """Load a chain from dump file and return as np arrays C and L."""

    C = np.load(fname)['C']
    try:
        L = np.load(fname)['L']
    except:
        L = np.zeros(C.shape[:-1])
    return (C,L)

def _load_from_backend(fname):
    """Load chain and logprob from HDF backend and return in standard form."""

    import emcee.backends
    be = emcee.backends.HDFBackend(fname)
    C = np.swapaxes(be.get_chain(), 0, 1)
    L = np.swapaxes(be.get_log_prob(), 0, 1)
    return (C,L)

def plot_traces(fname,ws='all',dims='all',burn=0,**kwargs):
    """Visually inspect chain traces."""
    
    # Load the sample(s)
    C,_ = load_chain(fname)
    if (type(ws) is str) and (ws == 'all'):
        ws = tuple(range(C.shape[0]))
    if type(ws) is int:
        ws = (ws,)
    C = (C[ws,burn:,:])
    print('C.shape =', C.shape)
    
    # Visual inspection
    if dims == 'all':
        dims = range(C.shape[-1])
    plot_compare_chains(C, dims=dims, **kwargs)

def plot_correlations(fname,ws=0,dims='all',burn=0,skip=1,showcorn=False,ml=100):
    """Visually inspect chain correlations."""
    
    # Load the sample(s)
    C,_ = load_chain(fname)
    C = (C[ws,burn::skip,:])
    print('C.shape =', C.shape)
    assert type(ws) is int
    
    # Inspect autocorrelations
    import matplotlib.pyplot as plt
    if dims == 'all':
        dims = tuple(range(C.shape[-1]))
    for dim in dims:
        plt.figure()
        x = C[:,dim]
        plt.acorr(x - np.mean(x), usevlines=True, maxlags=ml)
        plt.ylabel("dim{}".format(dim))
        plt.show(block=False)

    # Inspect corners
    if showcorn:
        import corner
        corner.corner(C);
        plt.show(block=False)

def plot_runningmeans(fname,ws=(0,),dims='all',burn=0,skip=1,**kwargs):
    """Visually inspect the running means."""
    
    # Load the sample(s)
    C,_ = load_chain(fname)
    if (type(ws) is str) and (ws == 'all'):
        ws = tuple(range(C.shape[0]))
    C = (C[ws,burn::skip,:])
    print('C.shape =', C.shape)
    
    # Inspect running means
    import matplotlib.pyplot as plt
    runmean = lambda x: np.cumsum(x)/(np.arange(x.shape[0]) + 1)
    if dims == 'all':
        dims = tuple(range(C.shape[-1]))
    for dim in dims:
        plt.figure()
        for w in range(C.shape[0]):
            x = C[w,:,dim]
            plt.plot(range(x.shape[0]), runmean(x), ',', **kwargs)
        plt.ylabel('x{}'.format(dim))
        if len(ws) == 1:
            plt.annotate('mean={:g}'.format(
                np.mean(x)), (0.01,0.9), xycoords='axes fraction');
        plt.show(block=False)
        pass

def gelmanrubin(fname,ws='all',burn=0,skip=1):
    """Calculate Gelman-Rubin score from output file."""

    # Load the sample(s)
    C,_ = load_chain(fname)
    if (type(ws) is str) and (ws == 'all'):
        ws = tuple(range(C.shape[0]))
    C = C[ws,burn::skip,:]

    # The Gelman-Rubin score
    return _grubin(C)

def _grubin(C):
    """Gelman Rubin score of chain."""

    M = C.shape[0]
    hats = np.mean(C, axis=1)
    hat = np.mean(hats, axis=0)
    B = 1/(M-1)*np.sum((hats - hat)**2, axis=0)
    W = np.mean(np.var(C, axis=1), axis=0)
    return 1 + (M + 1)/M*B/W

def hist_compare_chains(C1,C2,ws=...,dims='all',burn=0,**kwargs):
    """Compare dimensional histograms."""

    if dims == 'all':
        dims = tuple(range(C1.shape[-1]))
    if type(dims) is int:
        dims = (dims,)

    # Inspect histograms for each dim
    import matplotlib.pyplot as plt
    for dim in dims:
        plt.figure()
        x1 = C1[ws,burn:,dim].flat
        x2 = C2[ws,burn:,dim].flat
        plt.hist(x1, histtype='step', density=True, **kwargs)
        plt.hist(x2, histtype='step', density=True, **kwargs)
        plt.ylabel('x{}'.format(dim))
        plt.show(block=False)

def compare_walker_histograms(fname,ws='all',dims='all',burn=0,skip=1,**kwargs):
    """Compare single-walker histrograms."""

    # Load the sample(s)
    C,_ = load_chain(fname)
    if (type(ws) is str) and (ws == 'all'):
        ws = tuple(range(C.shape[0]))
    if dims == 'all':
        dims = tuple(range(C.shape[-1]))
    C = C[ws,burn::skip,:] # note no vstack, all ws all dims for now
    print('C.shape = ', C.shape)

    # Inspect histograms for each dim
    import matplotlib.pyplot as plt
    for dim in dims:
        plt.figure()
        for w in range(C.shape[0]):
            x = C[w,:,dim]
            plt.hist(x, 'auto', histtype='step', **kwargs)
        plt.ylabel('x{}'.format(dim))
        plt.show(block=False)
        pass

def compare_subhists(fname,subs,ws='all',dims='all',burn=0,skip=1,**kwargs):
    """Inspect histograms taken on random subsamples."""
    
    # Load the sample(s)
    C,_ = load_chain(fname)
    if (type(ws) is str) and (ws == 'all'):
        ws = tuple(range(C.shape[0]))
    C = np.vstack(C[ws,burn::skip,:])
    print('C.shape =', C.shape)
    
    # Inspect histograms
    import matplotlib.pyplot as plt
    if dims == 'all':
        dims = tuple(range(C.shape[-1]))
    for dim in dims:
        X = C[:,dim]
        plt.figure()
        for k in range(subs[0]):
            x = X[np.random.permutation(X.size)[:subs[1]]]
            plt.hist(x, 'auto', histtype='step', **kwargs)
        plt.ylabel('x{}'.format(dim))
        plt.show(block=False)
        pass

def plot_histograms(fname,ws='all',dims='all',burn=0,skip=1,trans=None,**kwargs):
    """Visually inspect individual histograms."""
    
    # Load the sample(s)
    C,_ = load_chain(fname)
    if (type(ws) is str) and (ws == 'all'):
        ws = tuple(range(C.shape[0]))
    C = np.vstack(C[ws,burn::skip,:])
    if trans is not None:
        C = np.array([trans(x) for x in C])
    print('C.shape =', C.shape)
    
    # Inspect histograms
    import matplotlib.pyplot as plt
    if dims == 'all':
        dims = tuple(range(C.shape[-1]))
    if type(dims) is int:
        dims = [dims]
    for dim in dims:
        x = C[:,dim]
        plt.figure()
        plt.hist(x, 'auto', **kwargs)
        plt.ylabel('x{}'.format(dim))
        plt.annotate('mean={:g}'.format(
            np.mean(x)), (0.01,0.9), xycoords='axes fraction');
        plt.show(block=False)
        pass

class ProgressBar:
    def __init__(self, iterations):
        self.iterations = iterations
        self.prog_bar = '[]'
        self.fill_char = '*'
        self.width = 50
        self.__update_amount(0)

    def animate(self, iter):
        print('\r', self, end='', file=sys.stderr)
        sys.stdout.flush()
        self.update_iteration(iter + 1)

    def update_iteration(self, elapsed_iter):
        self.__update_amount((elapsed_iter / float(self.iterations)) * 100.0)
        self.prog_bar += '  %d of %s complete' % (elapsed_iter, self.iterations)

    def __update_amount(self, new_amount):
        percent_done = int(round((new_amount / 100.0) * 100.0))
        all_full = self.width - 2
        num_hashes = int(round((percent_done / 100.0) * all_full))
        self.prog_bar = '[' + self.fill_char * num_hashes + ' ' * (all_full - num_hashes) + ']'
        pct_place = (len(self.prog_bar) // 2) - len(str(percent_done))
        pct_string = '%d%%' % percent_done
        self.prog_bar = self.prog_bar[0:pct_place] + \
            (pct_string + self.prog_bar[pct_place + len(pct_string):])

    def __str__(self):
        return str(self.prog_bar)

def _running_acceptance(chain):
    nwalk = chain.shape[0]
    nsamp = chain.shape[1]
    A = np.zeros((nwalk,nsamp))
    pb = ProgressBar(nwalk)
    for j in range(nwalk):
        pb.animate(j+1)
        nu = 1
        for k in range(1,nsamp):
            if not np.all(chain[j,k,:] == chain[j,k-1,:]):
                nu += 1
            A[j,k] = nu/(k+1)
            pass
        pass
    return A

def acceptance_ratio(C):
    """Walkers' acceptance."""
    a = np.zeros(C.shape[0])
    n = C.shape[1]
    for k in range(len(a)):
        a[k] = np.unique(C[k,:,:], axis=0).shape[0]/n
    return a

def plot_compare_chains(chains,dims=None,hardcopy=False,walkers=None, *args, **kwargs):
    import matplotlib.pyplot as plt
    if type(chains) is not tuple:
        chains = (chains,)
    if dims is None:
        dims = range(chains[0].shape[2])
    if type(dims) is int:
        dims = [dims]
    if dims is Ellipsis:
        dims = [dims]
    if walkers is None:
        walkers = range(chains[0].shape[0])
    if type(walkers) is int:
        walkers = [walkers]

    styles = ['b,','r,','g,']
    for dim in dims:
        plt.figure()
        try:
            if len(chains) == 1:
                plt.plot(chains[0][walkers,:,dim].T, ',', alpha=0.2, *args, **kwargs)
            else:
                for k in range(len(chains)):
                    plt.plot(chains[k][walkers,:,dim].T, styles[k%len(styles)], alpha=0.2, *args, **kwargs)
            plt.xlabel('step')
            plt.ylabel('dim {}'.format(dim))
            plt.show(block=False)
        except:
            plt.close()
            raise
        if hardcopy:
            plt.savefig('dim{}.png'.format(dim))
            plt.close()
            pass
        pass
    return

def reseed(fname):
    import os
    chain = np.load(fname)['C']
    seeds = chain[:,-1,:]
    outname = os.path.join(*os.path.split(fname)[0:-1],'reseeds.txt')
    np.savetxt(outname,seeds)
    print("seeds saved in {}.".format(outname))

def like_traces(fname,ws=(0,),burn=0,skip=1,**kwargs):
    """Visually inspect log-likelihood traces."""
    
    # Load the sample(s)
    _,L = load_chain(fname)
    if (type(ws) is str) and (ws == 'all'):
        ws = tuple(range(L.shape[0]))
    L = L[ws,burn::skip]
    print('L.shape =', L.shape)
    
    # Visual inspection
    import matplotlib.pyplot as plt
    plot_compare_chains(L, dims=..., **kwargs)
    plt.ylabel('log-likelihood')
    plt.annotate('max log-like = {:g}'.format(L.max()),
        (0.01,0.9), xycoords='axes fraction')

def _samplify_chain(fname,ws='all',burn=0,skip=1):
    """Cull a EnsembleSampler chain into a usable sample."""

    # Load the sample(s)
    C,_ = load_chain(fname)
    if (type(ws) is str) and (ws == 'all'):
        ws = tuple(range(C.shape[0]))
    C = C[ws,burn::skip,:]
    Z = C.reshape((-1,C.shape[-1]))

    print("Z.shape = {}".format(Z.shape))
    return Z

def samplify_chain(fname,ws='all',burn=0,skip=1):
    """Call _samplify_chain and save to consistently named output."""

    Z = _samplify_chain(fname, ws=ws, burn=burn, skip=skip)
    chain_dir, cf = os.path.split(fname)
    out_dir = os.path.join(chain_dir, 'samplified')
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    else:
        cout("WARNING: directory {} already exists, ".format(out_dir))
        cout("files may be overwritten.\n")
    sf = os.path.join(out_dir, cf)
    print("saving to", sf)
    np.savez_compressed(sf,Z=Z)
    af = os.path.join(out_dir, os.path.basename(cf)[:-3]+'txt')
    print("saving to", af)
    np.savetxt(af, Z)

def _fixpoly(x,rho0):
    return np.hstack((x,0,rho0-np.polyval(np.hstack((x,0,0)),1)))

def _append(fname,Z):
    Z0 = np.loadtxt(fname)
    Z = np.vstack((Z0,Z))
    np.savetxt(fname,Z)

def _hydrostatic_pressure(zvec, dvec, obs, ss, SS, flipum=False):
    if flipum:
        zvec = np.flipud(zvec)
        dvec = np.flipud(dvec)

    # let's get the potential
    from tof4 import Upu
    N = len(ss[0]) - 1
    s0 = ss[0][N]; s2 = ss[1][N]; s4 = ss[2][N]; s6 = ss[3][N]; s8 = ss[4][N]
    aos = 1 + s0 - (1/2)*s2 + (3/8)*s4 - (5/16)*s6 + (35/128)*s8
    U = Upu(ss, SS, obs.m)
    U = np.flipud(U)  # WARNING: U was bottom up b/c of ss and SS
    U = U*obs.G*obs.M*aos/obs.a0*zvec**2

    # need density in real units now
    svec = zvec*obs.a0/aos
    drho = np.hstack((dvec[0],np.diff(dvec)))
    m = (4*np.pi/3)*sum(drho*svec**3)
    rvec = dvec*obs.M/m

    # integrate hydrostatic equilibrium top down
    rho = 0.5*(rvec[:-1] + rvec[1:])
    P = np.zeros(dvec.shape)
    P[0] = obs.P0
    P[1:] = P[0] + np.cumsum(-rho*np.diff(U))

    return P

if __name__ == "__main__":
    print("alo world")
