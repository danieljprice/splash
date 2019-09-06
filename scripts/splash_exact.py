import os.path
from pathlib import Path

from ctypes import *

####        Add to your script
#
# import sys
# sys.path.append('/path/to/splash/scripts') (directory, which contains this script)
# import splash_exact
# import numpy
#
####        After those lines, all splash exact solutions are available for use
#
# x = numpy.linspace(-0.25, 0.25, num=1000)
# x, y, ierr = splash_exact.shock(x)
# if (ierr != 0):
#     plt.plot(x,y)
#     plt.show()
#
####

def prerun():
    libexactpath = Path(__file__).parent.parent.joinpath('build/libexact.so')
    if (not libexactpath.exists()):
        print("Could not find `libexact.so")
        print("Please compile with `make libexact` first.")
        exit(1)
    else:
        lib = cdll.LoadLibrary(str(libexactpath))
        lib.check_argcv_f()
        return lib

def checkfmt(str):
    return str.lower().replace(' ', '_').replace('-', '_').replace('/', '_')

def shock(
    x,
    plot         = 'density',
    time         = 0.1,
    gamma        = 5./3.,
    xshock       = 0.,
    rho_L        = 1.,
    rho_R        = 0.125,
    p_L          = 1.,
    p_R          = 0.1,
    v_L          = 0,
    v_R          = 0.,
    rdust_to_gas = 0.):

    pin = checkfmt(plot)
    if   (pin == 'density'):
        iplot = 1
    elif (pin == 'pressure'):
        iplot = 2
    elif (pin == 'velocity'):
        iplot = 3
    elif (pin == 'uthermal'):
        iplot = 4
    elif (pin == 'deltav'):
        iplot = 5
    elif (pin == 'dustfrac'):
        iplot = 6
    else:
        print("Splash Exact: Unrecognised plot type.")
        print("Splash Exact: plot = density || pressure || velocity || uthermal || deltav || dustfrac")
        exit(1)

    exactdll = prerun()
    ierr = 0
    c_x = (c_float*len(x))()
    c_y = (c_float*len(x))()
    c_x[:] = x[:]
    exactdll.shock_(
        byref(c_int(iplot)),
        byref(c_int(len(x))),
        byref(c_float(time)),
        byref(c_float(gamma)),
        byref(c_float(xshock)),
        byref(c_float(rho_L)),
        byref(c_float(rho_R)),
        byref(c_float(p_L)),
        byref(c_float(p_R)),
        byref(c_float(v_L)),
        byref(c_float(v_R)),
        byref(c_float(rdust_to_gas)),
        byref(c_x),
        byref(c_y),
        byref(c_int(ierr))
    )
    return c_x[:], c_y[:],ierr

def shock_sr(
    x,
    plot    = 'density',
    time    = 0.1,
    gamma   = 5./3.,
    rho_L   = 1.,
    rho_R   = 0.125,
    p_L     = 1.0,
    p_R     = 0.1,
    v_L     = 0.,
    v_R     = 0.):

    pin = checkfmt(plot)
    if   (pin == 'density'):
        iplot = 1
    elif (pin == 'pressure'):
        iplot = 2
    elif (pin == 'velocity'):
        iplot = 3
    elif (pin == 'uthermal'):
        iplot = 4
    elif (pin == 'density*'):
        iplot = 5
    else:
        print("Splash Exact: Unrecognised plot type.")
        print("Splash Exact: plot = density || pressure || velocity || uthermal || density*")
        exit(1)

    exactdll = prerun()
    ierr = 0
    c_x = (c_float*len(x))()
    c_y = (c_float*len(x))()
    c_x[:] = x[:]
    exactdll.shock_sr_(
        byref(c_int(iplot)),
        byref(c_int(len(x))),
        byref(c_float(time)),
        byref(c_float(gamma)),
        byref(c_float(rho_L)),
        byref(c_float(rho_R)),
        byref(c_float(p_L)),
        byref(c_float(p_R)),
        byref(c_float(v_L)),
        byref(c_float(v_R)),
        byref(c_x),
        byref(c_y),
        byref(c_int(ierr))
    )
    return c_x[:],c_y[:],ierr

def sedov(
    r,
    plot    = 'density',
    time    = 0.1,
    gamma   = 5./3.,
    rhozero = 1.,
    energy  = 1.,
    rmax    = 1.):

    pin = checkfmt(plot)
    if   (pin == 'density'):
        iplot = 1
    elif (pin == 'pressure'):
        iplot = 2
    elif (pin == 'uthermal'):
        iplot = 3
    elif (pin == 'kinetic_energy'):
        iplot = 4
    elif (pin == 'velocity'):
        iplot = 5
    else:
        print("Splash Exact: Unrecognised plot type.")
        print("Splash Exact: plot = density || pressure || uthermal || kinetic_energy ||  velocity ")
        exit(1)

    exactdll = prerun()
    ierr = 0
    c_r = (c_float*len(r))()
    c_y = (c_float*len(r))()
    c_r[:] = r[:]
    exactdll.sedov_(
        byref(c_int(iplot)),
        byref(c_int(len(r))),
        byref(c_float(time)),
        byref(c_float(gamma)),
        byref(c_float(rhozero)),
        byref(c_float(energy)),
        byref(c_float(rmax)),
        byref(c_r),
        byref(c_y),
        byref(c_int(ierr))
    )
    return c_r[:],c_y[:],ierr

def polytrope(
    r,
    gamma   = 5./3.,
    polyk   = 1.,
    totmass = 1.):

    exactdll = prerun()
    ierr = 0
    c_r = (c_float*len(r))()
    c_y = (c_float*len(r))()
    c_r[:] = r[:]
    npartout = len(r)
    nout = c_int(0)
    exactdll.polytrope_(
        byref(c_int(len(r))),
        byref(c_float(gamma)),
        byref(c_float(polyk)),
        byref(c_float(totmass)),
        byref(c_r),
        byref(c_y),
        byref(c_nout),
        byref(c_int(ierr))
    )
    return c_r[:nout.value],c_y[:nout.value],ierr

def toystar1D(
    x,
    plot    = 'density',
    time    = 1.,
    gamma   = 5./3.,
    H0      = 1.,
    A0      = 1.,
    C0      = 1.,
    sigma   = 0.,
    norder  = -1):

    pin = checkfmt(plot)
    if   (pin == 'density'):
        iplot = 1
    elif (pin == 'pressure'):
        iplot = 2
    elif (pin == 'uthermal'):
        iplot = 3
    elif (pin == 'velocity_x'):
        iplot = 4
    elif (pin == 'mag_field_y'):
        iplot = 5
    elif (pin == 'ac_plane'):
        iplot = 7
    else:
        print("Splash Exact: Unrecognised plot type.")
        print("Splash Exact: plot = density || pressure || uthermal || velocity_x ||  mag_field_y || ac_plane ")
        exit(1)

    exactdll = prerun()
    ierr = 0
    c_x = (c_float*len(x))()
    c_y = (c_float*len(x))()
    c_x[:] = x[:]
    exactdll.toystar1d_(
        byref(c_int(iplot)),
        byref(c_int(len(x))),
        byref(c_float(time)),
        byref(c_float(gamma)),
        byref(c_float(H0)),
        byref(c_float(A0)),
        byref(c_float(C0)),
        byref(c_float(sigma)),
        byref(c_int(norder)),
        byref(c_x),
        byref(c_y),
        byref(c_int(ierr))
    )
    return c_x[:],c_y[:],ierr

def toystar2D(
    x,
    plot    = 'density',
    time    = 1.,
    gamma   = 5./3.,
    polyk   = 1.,
    totmass = 1.,
    A0      = 1.,
    H0      = 1.,
    C0      = 1.,
    jorder  = -1,
    morder  = 0,
    V11     = 0,
    V22     = 0,
    V12     = 0,
    V21     = 0):

    pin = checkfmt(plot)
    if   (pin == 'density'):
        iplot = 1
    elif (pin == 'pressure'):
        iplot = 2
    elif (pin == 'uthermal'):
        iplot = 3
    elif (pin == 'velocity_x'):
        iplot = 4
    elif (pin == 'velocity_y'):
        iplot = 5
    elif (pin == 'x_vs_y'):
        iplot = 0
    else:
        print("Splash Exact: Unrecognised plot type.")
        print("Splash Exact: plot = density || pressure || uthermal || velocity_x ||  velocity_y || x_vs_y ")
        exit(1)

    exactdll = prerun()
    ierr = 0
    c_x = (c_float*len(x))()
    c_y = (c_float*len(x))()
    c_x[:] = x[:]

    exactdll.toystar2d_(
        byref(c_int(iplot)),
        byref(c_int(len(x))),
        byref(c_float(time)),
        byref(c_float(gamma)),
        byref(c_float(polyk)),
        byref(c_float(totmass)),
        byref(c_float(A0)),
        byref(c_float(H0)),
        byref(c_float(C0)),
        byref(c_int(jorder)),
        byref(c_int(morder)),
        byref(c_float(V11)),
        byref(c_float(V22)),
        byref(c_float(V12)),
        byref(c_float(V21)),
        byref(c_x),
        byref(c_y),
        byref(c_int(ierr))
    )
    return c_x[:],c_y[:],ierr

def gresho(
    x,
    plot = 'velocity_phi'):

    pin = checkfmt(plot)
    if   (pin == 'velocity_phi'):
        iplot = 1
    elif (pin == 'pressure'):
        iplot = 2
    else:
        print("Splash Exact: Unrecognised plot type.")
        print("Splash Exact: plot = velocity_phi || pressure ")
        exit(1)

    exactdll = prerun()
    ierr = 0
    c_x = (c_float*len(x))()
    c_y = (c_float*len(x))()
    c_x[:] = x[:]

    exactdll.gresho_(
        byref(c_int(iplot)),
        byref(c_int(len(x))),
        byref(c_x),
        byref(c_y),
        byref(c_int(ierr))
    )
    return c_x[:],c_y[:],ierr

def mhdshock(
    x,
    plot     = 'density',
    solution = '7 jump',
    time     = 0.2,
    gamma    = 5./3.,
    xmin     = -1,
    xmax     = 1,
    xshock   = 0):

    pin = checkfmt(plot)
    if   (pin == 'density'):
        iplot = 1
    elif (pin == 'pressure'):
        iplot = 2
    elif (pin == 'velocity_x'):
        iplot = 3
    elif (pin == 'velocity_y'):
        iplot = 4
    elif (pin == 'velocity_z'):
        iplot = 5
    elif (pin == 'mag_field_y'):
        iplot = 6
    elif (pin == 'mag_field_z'):
        iplot = 7
    elif (pin == 'uthermal'):
        iplot = 8
    elif (pin == 'Bxzero'):
        iplot = 9
    else:
        print("Splash Exact: Unrecognised plot type.")
        print("Splash Exact: plot = density || pressure || velocity_x ||  velocity_y || velocity_z || mag_field_y || mag_field_z || uthermal || Bxzero ")
        exit(1)

    ishk = 0
    pin = checkfmt(solution)
    if   (pin == 'brio_wu'):
        ishk = 1
    elif (pin == 'fast_slow'):
        ishk = 2
    elif (pin == '7_jump'):
        ishk = 3
    elif (pin == 'isothermal'):
        ishk = 4
    elif (pin == 'rarefaction'):
        ishk = 5
    elif (pin == 'mach_25'):
        ishk = 6
    elif (pin == 'toth'):
        ishk = 7
    else:
        print("Splash Exact: Unrecognised solution type.")
        print("Splash Exact: solution = Brio/Wu || fast/slow || 7_jump || isothermal ||  rarefaction || Mach_25 || Toth ")
        exit(1)

    exactdll = prerun()
    ierr = 0
    c_x = (c_float*len(x))()
    c_y = (c_float*len(x))()
    c_x[:] = x[:]
    c_nout = c_int(0)
    exactdll.mhdshock_(
        byref(c_int(iplot)),
        byref(c_int(len(x))),
        byref(c_int(ishk)),
        byref(c_float(time)),
        byref(c_float(gamma)),
        byref(c_float(xmin)),
        byref(c_float(xmax)),
        byref(c_float(xshock)),
        byref(c_x),
        byref(c_y),
        byref(c_nout),
        byref(c_int(ierr))
    )
    return c_x[:c_nout.value],c_y[:c_nout.value],ierr

def rhoh(
    x,
    plot        = 'density',
    ndim        = 3,
    hfact       = 1.2,
    pmassval    = 1.0):

    pin = checkfmt(plot)
    if   (pin == 'h'):
        iplot = 1
    elif (pin == 'density'):
        iplot = 2
    else:
        print("Splash Exact: Unrecognised plot type.")
        print("Splash Exact: plot = h || density ")
        exit(1)

    exactdll = prerun()
    ierr = 0
    c_x = (c_float*len(x))()
    c_y = (c_float*len(x))()
    c_x[:] = x[:]

    exactdll.rhoh_(
        byref(c_int(iplot)),
        byref(c_int(len(x))),
        byref(c_int(ndim)),
        byref(c_float(hfact)),
        byref(c_float(pmassval)),
        byref(c_x),
        byref(c_y),
        byref(c_int(ierr))
    )
    return c_x[:],c_y[:],ierr

def densityprofiles(
    x,
    plot    = 'density',
    profile = 'Plummer',
    Msphere = [1.0,0.0],
    rsoft   = [1.0,0.1]):

    pin = checkfmt(plot)
    if   (pin == 'density'):
        iplot = 1
    elif (pin == 'potential'):
        iplot = 2
    elif (pin == 'force'):
        iplot = 3
    else:
        print("Splash Exact: Unrecognised plot type.")
        print("Splash Exact: plot = density || potential || force ")
        exit(1)

    pin = checkfmt(profile)
    if   (pin == 'plummer'):
        iprofile = 1
    elif (pin == 'hernquist'):
        iprofile = 2
    else:
        print("Splash Exact: Unrecognised plot type.")
        print("Splash Exact: profile = Plummer || Hernquist ")
        exit(1)

    exactdll = prerun()
    ierr = 0
    c_x = (c_float*len(x))()
    c_y = (c_float*len(x))()
    c_x[:] = x[:]

    exactdll.densityprofiles_(
        byref(c_int(iplot)),
        byref(c_int(len(x))),
        byref(c_int(iprofile)),
        byref(c_float(Msphere[0])),
        byref(c_float(Msphere[1])),
        byref(c_float(rsoft[0])),
        byref(c_float(rsoft[1])),
        byref(c_x),
        byref(c_y),
        byref(c_int(ierr))
    )
    return c_x[:],c_y[:],ierr

def torus(
    x,
    plot        = 'density',
    torus       = 'default',
    Mstar       = 1.0,
    Rtorus      = 1.0,
    polyk       = 0.0764,
    distortion  = 1.1,
    gamma       = 5./3.):

    pin = checkfmt(plot)
    if   (pin == 'density'):
        iplot = 1
    elif (pin == 'pressure'):
        iplot = 2
    elif (pin == 'uthermal'):
        iplot = 3
    elif (pin == 'Btheta'):
        iplot = 4
    elif (pin == 'Jphi_current'):
        iplot = 5
    else:
        print("Splash Exact: Unrecognised plot type.")
        print("Splash Exact: plot = density || pressure || uthermal || Btheta || Jphi_current")
        exit(1)

    itorus = 0
    pin = checkfmt(torus)
    if   (pin == 'default'):
        itorus = 1
    elif (pin == 'tokamak'):
        itorus = 2
    else:
        print("Splash Exact: Unrecognised torus type.")
        print("Splash Exact: torus = Default || Tokamak ")
        exit(1)

    exactdll = prerun()
    ierr = 0
    c_x = (c_float*len(x))()
    c_y = (c_float*len(x))()
    c_x[:] = x[:]

    exactdll.torus_(
        byref(c_int(iplot)),
        byref(c_int(len(x))),
        byref(c_int(itorus)),
        byref(c_float(Mstar)),
        byref(c_float(Rtorus)),
        byref(c_float(polyk)),
        byref(c_float(distortion)),
        byref(c_float(gamma)),
        byref(c_x),
        byref(c_y),
        byref(c_int(ierr))
    )
    return c_x[:],c_y[:],ierr

def ringspread(
    x,
    plot   = 'density',
    time   = 1.0,
    Mdisk  = 1.0,
    Rdisk  = 1.0,
    viscnu = 1.e-3):

    pin = checkfmt(plot)
    if   (pin == 'density'):
        iplot = 1
    else:
        print("Splash Exact: Unrecognised plot type.")
        print("Splash Exact: plot = density")
        exit(1)

    exactdll = prerun()
    ierr = 0
    c_x = (c_float*len(x))()
    c_y = (c_float*len(x))()
    c_x[:] = x[:]

    exactdll.ringspread_(
        byref(c_int(iplot)),
        byref(c_int(len(x))),
        byref(c_float(time)),
        byref(c_float(Mdisk)),
        byref(c_float(Rdisk)),
        byref(c_float(viscnu)),
        byref(c_x),
        byref(c_y),
        byref(c_int(ierr))
    )
    return c_x[:],c_y[:],ierr

def dustywave(
    x,
    plot       = 'gas_density',
    time       = 1.0,
    ampl       = 1.0,
    cs         = 1.0,
    Kdrag      = 1.0,
    lambdacoef = 1.0,
    x0         = 1.0,
    rhog0      = 1.0,
    rhod0      = 1.0):

    pin = checkfmt(plot)
    if   (pin == 'gas_velocity'):
        iplot = 1
    elif (pin == 'dust_velocity'):
        iplot = 2
    elif (pin == 'gas_density'):
        iplot = 3
    elif (pin == 'dust_density'):
        iplot = 4
    else:
        print("Splash Exact: Unrecognised plot type `" + plot + "`")
        print("Splash Exact: plot = gas_velocity || dust_velocity || gas_density || dust_density")
        exit(1)

    exactdll = prerun()
    ierr = 0
    c_x = (c_float*len(x))()
    c_y = (c_float*len(x))()
    c_x[:] = x[:]

    exactdll.dustywave_(
        byref(c_int(iplot)),
        byref(c_int(len(x))),
        byref(c_float(time)),
        byref(c_float(ampl)),
        byref(c_float(cs)),
        byref(c_float(Kdrag)),
        byref(c_float(lambdacoef)),
        byref(c_float(x0)),
        byref(c_float(rhog0)),
        byref(c_float(rhod0)),
        byref(c_x),
        byref(c_y),
        byref(c_int(ierr))
    )
    return c_x[:],c_y[:],ierr

def rochelobe(
    x,
    primatypos = [0.,0.],
    secondarypos = [1.,0.],
    primarymass =  1.,
    secondarymass = 1.):

    exactdll = prerun()
    ierr = 0
    c_x = (c_float*len(x))()
    c_y = (c_float*len(x))()
    c_x[:] = x[:]

    exactdll.rochelobe_(
        byref(c_int(len(x))),
        byref(c_float(primatypos[0])),
        byref(c_float(primatypos[1])),
        byref(c_float(secondarypos[0])),
        byref(c_float(secondarypos[0])),
        byref(c_float(primarymass)),
        byref(c_float(secondarymass)),
        byref(c_x),
        byref(c_y),
        byref(c_int(ierr))
    )
    return c_x[:],c_y[:],ierr

def cshock(
    x,
    plot='density',
    time=0.2,
    gamma=5./3.,
    machs = 50.,
    macha = 5.,
    xmin=-0.25,
    xmax= 0.25,):

    pin = checkfmt(plot)
    if   (pin == 'density'):
        iplot = 1
    elif (pin == 'mag_field_y'):
        iplot = 2
    elif (pin == 'velocity_x'):
        iplot = 3
    elif (pin == 'velocity_y'):
        iplot = 4
    elif (pin == 'mag_field_x'):
        iplot = 5
    else:
        print("Splash Exact: Unrecognised plot type `" + plot + "`")
        print("Splash Exact: plot = density || mag_field_x || mag_field_y || velocity_x || velocity_y")
        exit(1)

    exactdll = prerun()
    ierr = 0
    c_x = (c_float*len(x))()
    c_y = (c_float*len(x))()
    c_x[:] = x[:]

    exactdll.cshock_(
        byref(c_int(iplot)),
        byref(c_int(len(x))),
        byref(c_float(time)),
        byref(c_float(gamma)),
        byref(c_float(machs)),
        byref(c_float(macha)),
        byref(c_float(xmin)),
        byref(c_float(xmax)),
        byref(c_x),
        byref(c_y),
        byref(c_int(ierr))
    )
    return c_x[:],c_y[:],ierr

def check_spiral_params(i1,i2,j1,j2,nparams,nsolutions):
    if (i1 > i2):
        print("Splash Exact: Wrong Spiral Params. i1 > i2. Should define blocks [i1,i2,j1,j2,val].")
        exit(1)
    if (j1 > j2):
        print("Splash Exact: Wrong Spiral Params. j1 > j2. Should define blocks [i1,i2,j1,j2,val].")
        exit(1)
    if (i2 >= nparams):
        print("Splash Exact: Wrong Spiral Params. i1 > "+str(nparams)+". Should define blocks [i1,i2,j1,j2,val].")
        exit(1)
    if (j2 >= nsolutions):
        print("Splash Exact: Wrong Spiral Params. i1 > "+str(nsolutions)+". Should define blocks [i1,i2,j1,j2,val].")
        exit(1)

def planetdisc(
    x,
    plot='phi/r plane',
    spiral='Ogilvie/Rafikov',
    time=0.2,
    HonR = 0.05,
    rplanet = 1.,
    q_index = 0.25,
    narms = 1,
    spiral_params=[[1,1,0,6,360]]):
    # same as filling block i1 to i2, j1 to j2 with x
    # spiral_params = 0.
    # spiral_params(2,:) = 360.
    pin = checkfmt(plot)
    if   (pin == 'phi_r_plane'):
        iplot = 1
    elif (pin == 'x_y_plane'):
        iplot = 2
    else:
        print("Splash Exact: Unrecognised plot type `" + plot + "`")
        print("Splash Exact: plot = phi-r plane || x-y plane")
        exit(1)

    ispiral = 0
    pin = checkfmt(spiral)
    if   (pin == 'ogilvie_rafikov'):
        ispiral = 1
    elif (pin == 'spiral_arm_fitting'):
        ispiral = 2
    else:
        print("Splash Exact: Unrecognised plot type `" + plot + "`")
        print("Splash Exact: plot = Ogilvie/Rafikov || spiral_arm_fitting")
        exit(1)

    exactdll = prerun()
    ierr = 0
    c_x = (c_float*len(x))()
    c_y = (c_float*len(x))()
    c_x[:] = x[:]
    nparams = 7
    nsolutions = 10
    c_params = (c_float*nparams*nsolutions)()
    for e in spiral_params:
        i1 = e[0]
        i2 = e[1]
        j1 = e[2]
        j2 = e[3]
        val = e[4]
        check_spiral_params(i1,i2,j1,j2,nparams,nsolutions)
        for i in range(i1,i2+1):
            for j in range(j1,j2+1):
                c_params[i][j] = val

    exactdll.planetdisc_(
        byref(c_int(iplot)),
        byref(c_int(len(x))),
        byref(c_int(ispiral)),
        byref(c_float(time)),
        byref(c_float(HonR)),
        byref(c_float(rplanet)),
        byref(c_float(q_index)),
        byref(c_int(narms)),
        byref(c_params),
        byref(c_x),
        byref(c_y),
        byref(c_int(ierr))
    )
    return c_x[:],c_y[:],ierr

def bondi(
    x,
    plot = 'density',
    time = 10,
    gamma = 5./3.,
    const1 = 8.,
    const2 = 1.,
    Mstar = 1,
    relativistic  = True,
    geodesic_flow = False,
    is_wind       = True):

    pin = checkfmt(plot)
    if   (pin == 'velocity_x'):
        iplot = 1
    elif (pin == 'uthermal'):
        iplot = 2
    elif (pin == 'density'):
        iplot = 2
    else:
        print("Splash Exact: Unrecognised plot type `" + plot + "`")
        print("Splash Exact: plot = velocity_x || uthermal || density")
        exit(1)

    exactdll = prerun()
    ierr = 0
    c_x = (c_float*len(x))()
    c_y = (c_float*len(x))()
    c_x[:] = x[:]

    exactdll.bondi_(
        byref(c_int(iplot)),
        byref(c_int(len(x))),
        byref(c_float(time)),
        byref(c_float(gamma)),
        byref(c_float(const1)),
        byref(c_float(const2)),
        byref(c_float(Mstar)),
        byref(c_bool(relativistic)),
        byref(c_bool(geodesic_flow)),
        byref(c_bool(is_wind)),
        byref(c_x),
        byref(c_y),
        byref(c_int(ierr))
    )
    return c_x[:],c_y[:],ierr
