import numpy as np
from fitstools import manage_dtype, display, mask_fits, row_avg
import matplotlib.pyplot as plt
plt.ion()

@manage_dtype()
def make_throughput_map(fmask, mflat):
    #Generate a blank throughput map
    throughput_map = np.ones_like(fmask)

    #Obtain a list of the used fibers from the fiber mask.
    fnums = [int(n) for n in set(fmask.flatten()) if n != 0]
    
    #Correct for the profile of each flat fiber
    #for fnum in fnums:
    #    flat_spec = row_avg(mask_fits(mflat, fmask, maskval=fnum))
    #    medianf = np.median(flat_spec)
    #    for i,counts in enumerate(flat_spec):
    #        throughput_map[i][np.where(fmask[i]==fnum)] *= flat_spec[i]/medianf


    #Correct for fiber-to-fiber throughput.
    profile_corrected_mflat = mflat/throughput_map
    flat_fiber_avg_vals = []
    for fnum in fnums:
        flat_spec = row_avg(mask_fits(profile_corrected_mflat, fmask, maskval=fnum))
        avg_val = np.mean(flat_spec)
        flat_fiber_avg_vals.append(avg_val)
    medianf = np.median(flat_fiber_avg_vals)
    for f,fnum in zip(flat_fiber_avg_vals,fnums):
        throughput_map[np.where(fmask==fnum)] *= f/medianf

    display(throughput_map)

    plot=True
    if plot:
        throughput_corrected = mflat/throughput_map
        fig, ax = plt.subplots()
        for fnum in fnums:
            flat_spec = row_avg(mask_fits(throughput_corrected, fmask, maskval=fnum))
            ax.plot(flat_spec)
    return throughput_map
