import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def xcorl(star,temp,range1,*args,**kwargs):
    #12-Jun-92 JAV	Added minchi parameter and logic.
    #17-Jun-92 JAV	Added "Max. allowable range" error message.
    #24-Aug-92 JAV	Supressed output of blank line when print keyword absent.
    #3-Jan-97 GB    Added "fine" (# pixs finely resolved) and "mult" options
    #  these give finer resolution to the peak, and use products instead of diffs.
    #8-Jan-97 GB	Added "fshft" to force one of a double peak
    #23-Oct-01 GB   Added /full keyword to simplify the call
    #28-Feb-13 CAT   Ported to Python
    #16-Jun-16 AYK   Added to hammer code

    # Set the defaults
    pr = 0
    fine = 0
    mult=0
    fshft=0
    full=0
    ff = 0

    # Read the arguments
    for arg in args:
        if arg.lower() == 'fine': fine = 1
        if arg.lower() == 'full': full = 1

    # Read the keywords
    for key in kwargs:
        if key.lower() == 'mult':
            mult = kwargs[key]
        if key.lower() == 'fshft':
            fshft = kwargs[key]

    if "plot" in kwargs:
        plot = kwargs["plot"]
    else:
        plot = False

    ln = len(temp)
    ls = len(star)
    length = np.min([ln, ls])
    if range1 > (length-1)/2:
        print( 'Maximum allowable "range" for this case is' + str((length-1)/2))
    newln = length - 2*range1  # Leave "RANGE" on ends for overhang.
    newend = range1 + newln - 1

    #Normalize obs and template spectra by their means.
    te = temp/(np.mean(temp))
    st = star/(np.mean(star))

    x = np.arange(-range1, range1+1)
    chi = np.zeros(len(x))

    if full == 1:
        pr=1

    #Shift observed spectrum across template spectrum and fill chi array.
    for j in range(-range1, range1+1):    # Goose step, baby
        if mult == 1:
            dif = te[range1:newend+1] * st[range1+j:newend+j+1]
            chi[j+range1] = np.sum(abs(dif))
        else:
            dif = te[range1:newend+1] - st[range1+j:newend+j+1]
            chi[j+range1] = np.sum(dif*dif)
    xcr = chi


    length = len(x) * 100
    xl = np.arange(length)
    xl = xl/100. - range1
    xp = xl[0:length-99]
    function2 = interp1d(x, chi, kind='cubic')
    cp = function2(xp)
    if plot:
        fig, ax = plt.subplots()
        ax.set_xlabel("shft")
        ax.set_ylabel("Cross Correlation")
        ax.plot(x, chi, color="green")
        ax.plot(xp, cp, color="mediumseagreen")
    if mult != 0:
        minchi = np.max(cp)
        mm = np.where(cp == minchi)
    else:
        minchi = np.min(cp)
        mm = np.where(cp == minchi)
    shft = xp[mm[0]]
    if pr != 0:
        print( 'XCORL: The shift is: %10.2f'%(shft))
    if abs(shft) > range1:
        ff=1
        return
    if fshft != 0:
        shft = fshft

    if fine != 0:
        nf = fine*20+1
        rf = fine*10.
        nc = -1
        fchi = np.zeros(nf)
        xl = np.zeros(nf)
        for j in range(int(-rf), int(rf+1)):
            xl[nc+1] = shft + j/10.
            nst = self.shfour(st, -xl[nc+1])
            nc += 1
            if mult == 1:
                dif = nst[range1:newend+1] * te[range1:newend+1]
                fchi[nc] = np.sum(abs(dif))
            else:
                dif = nst[range1:newend+1] - te[range1:newend+1]
                fchi[nc] = np.sum(np.real(dif*dif))
        xp = np.arange( (nf-1) * 100 + 1) / 1000. + shft - fine
        function3 = interp1d(xl, fchi, kind='cubic')
        cp = function3(xp)
        if mult == 1:
            minchi = np.max(cp)
            mm = np.where(cp == minchi)
        else:
            minchi = np.min(cp)
            mm = np.where(cp == minchi)
        fshft = xp[mm]

        if pr != 0:
            print( 'XCORL: The final shift is: %10.2f'%(fshft))
    else:
        fshft=shft
    shft=fshft
    return shft
