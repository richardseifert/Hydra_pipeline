from ensure_path import ensure_path
import numpy as np
from fitstools import combine
from astropy.io import fits
from flat import find_fibers, make_fiber_profile_map, make_throughput_map
from thar import wvlsol, wvlsolver
from extract import fibers, optimal_extraction, robust_mean_extraction
from spectra import spectrum, median_spectra
from sys import stdout
import os
from os.path import exists
import shutil
import time
from calib import calibrate
from make_recipe import load_recipes
from skyflat import convolve_gaussian, flatten_spec, get_rv

np.warnings.filterwarnings('ignore')

class processor(object):
    cp_fnames = {'master_bias':'master_bias.fits',
                     'master_flat':'master_flat.fits',
                     'fiber_mask':'fiber_mask.fits',
                     'profile_map':'fiber_profile_map.fits',
                     'throughput_map':'throughput_map.fits',
                     'wavelength_solution':'wvlsol.fits',
                     'master_sky_spec':'master_sky_spec.fits'}
    def __init__(self, dname, recipes=None, output_log=None, plotter=None):
		self.dname=dname
		if recipes==None:
			recipes = 'recipes/'+self.dname+'.recipe'
		if os.path.exists(recipes):
			self.recipes = load_recipes(recipes)
		else:
			raise OSError("Recipe file not found: "+recipes)
		self.output_log=output_log
		self.init_dirs()
                self.output_log.set_log_path(self.outdata+'/output.log')
                self.plotter = plotter
    def init_dirs(self):
                self.indata = 'indata/'+self.dname
                if not os.path.exists(self.indata):
                    raise OSError("Indata folder "+self.indata+" not found.")
                self.calib = 'calib/'+self.dname
                ensure_path(self.calib+'/')
                self.calib_dirs = {}
                self.outdata = 'outdata/'+self.dname
                ensure_path(self.outdata+'/')
                self.outdata_dirs = {}
                self.plot_dirs = {}
                for pnum in list({r.pnum for r in self.recipes if r.rtype!='zero'}):
                    self.calib_dirs[pnum] = self.calib+'/P'+str(pnum)
                    ensure_path(self.calib_dirs[pnum])
                    self.outdata_dirs[pnum] = self.outdata+'/P'+str(pnum)
                    ensure_path(self.outdata_dirs[pnum])
                    self.plot_dirs[pnum] = self.calib_dirs[pnum]+'/plots'
                    ensure_path(self.plot_dirs[pnum]+'/')

    def load_calib_products(self, pnum):
        for product in self.cp_fnames.keys():
    		fpath = None
    		if os.path.exists(self.calib_dirs[pnum]+'/'+self.cp_fnames[product]):
    			fpath = self.calib_dirs[pnum]+'/'+self.cp_fnames[product]
    		elif os.path.exists(self.calib+'/'+self.cp_fnames[product]):
    			fpath = self.calib+'/'+self.cp_fnames[product]
    		if fpath != None:
    			f = fits.open(fpath)
    			if product != 'master_sky_spec':
    				d = f[0].data
    			else:
    				d = spectrum(f[1].data, f[0].data, f[2].data)
    			f.close()
    			setattr(self, product, d)

    def output(self, message=None, progress=None, **kwargs):
		if self.output_log != None:
                    if message!=None:
                        self.output_log.edit_message(message, **kwargs)
                    if progress!=None:
                        self.output_log.edit_progress(progress)
    def get_recipes(self, pnum=None, rtype=None):
		return [r for r in self.recipes if (pnum == None or pnum == r.pnum) and (rtype == None or rtype == r.rtype)]
    def make_master_bias(self):
		self.output('Generating master bias frame.')
		filenames = [fname for r in self.get_recipes(rtype='zero') for fname in r.filenames]
		if len(filenames) > 0:
                    #Load files and create master bias
                    biases = [fits.open(self.indata+'/'+fname) for fname in filenames]
                    self.master_bias = combine(*biases, method='mean')
                    self.master_bias.writeto(self.calib+'/master_bias.fits', clobber=True)
                    for f in biases:
                            f.close()
    def get_master_bias(self):
		try:
			return self.master_bias
		except AttributeError:
			mb_path = self.calib+'/master_bias.fits'
			if os.path.exists(mb_path):
				mb = fits.open(mb_path)
				self.master_bias = mb[0].data
				mb.close()
                        elif len(self.get_recipes(rtype='zero')) == 0:
                            #Load the super master bias, since no bias recipes were given.
                            
                            #Determine the readtime by looking at a sample target frame.
                            sample_fits = fits.open( self.indata+'/'+(fname for r in self.get_recipes(rtype='flat') for fname in r.filenames).next() )

                            readtime = sample_fits[0].header['READTIME']
                            sample_fits.close()

                            #Load the appropriate master bias.
                            if readtime == 34:
				mb = fits.open('calib/master_calib/master_bias34.fits')
				self.master_bias = mb[0].data
				mb.close()
                            elif readtime == 45:
				mb = fits.open('calib/master_calib/master_bias45.fits')
				self.master_bias = mb[0].data
				mb.close()
                            elif readtime == 67:
				mb = fits.open('calib/master_calib/master_bias67.fits')
				self.master_bias = mb[0].data
				mb.close()
                            else:
                                raise OSError("Unable to find master bias with a readtime of "+str(readtime)+".")
                                
			else:
				self.make_master_bias()
			return self.master_bias
    
class process_flat(processor):
        def __init__(self, dname, recipes=None, output_log=None, plotter=None):
            processor.__init__(self, dname, recipes=recipes, output_log=output_log, plotter=plotter)
	
        def dispatch(self):
            master_bias = self.get_master_bias()
            for i,r in enumerate(self.get_recipes(rtype='flat')):
                self.output(message='Generating master flat frame.', progress='Processing '+self.dname+' flats: '+str(i+1)+'/'+str(len(self.get_recipes(rtype='flat')))+' |')

                #Generate master flat frame, save to calib directory.
                flats = [fits.open(self.indata+'/'+filename) for filename in r.filenames]
                master_flat = combine(*flats)
                for flat in flats:
                    flat.close()
                mf_path = self.calib_dirs[r.pnum]+'/'+self.cp_fnames['master_flat']
                master_flat.writeto(mf_path, clobber=True)
    		
                #Find fibers, make fiber mask, save to calib directory.
                self.output('Locating fibers.')
                fiber_mask = find_fibers(master_flat, r.fibers)
                fm_path = self.calib_dirs[r.pnum]+'/'+self.cp_fnames['fiber_mask']
                fits.writeto(fm_path, fiber_mask, clobber=True)
                self.output('Fiber mask saved at '+fm_path)
                
                #Generate a fiber profile map, save to calib directory.
                self.output('Generating fiber profile map.')
                profile_map = make_fiber_profile_map(master_flat, fiber_mask)
                fp_path = self.calib_dirs[r.pnum]+'/'+self.cp_fnames['profile_map']
                fits.writeto(fp_path, profile_map, clobber=True)
                self.output('Fiber profile map saved at '+fp_path)
    		
                #Calibrate master flat.
                self.output('Bias correcting master flat frame.')
                master_flat = calibrate(master_flat, master_bias, fiber_mask, lacosmic=False)
    		
                #Generate a fiber thoughput map, save to calib directory.
                self.output('Generating throughput map.')
                throughput_map = make_throughput_map(fiber_mask, master_flat)
                tm_path = self.calib_dirs[r.pnum]+'/'+self.cp_fnames['throughput_map']
                fits.writeto(tm_path, throughput_map, clobber=True)
                self.output('Throughput map saved at '+tm_path)


class process_thar(processor):
        def __init__(self, dname, recipes=None, output_log=None, fast=False, plotter=None):
            processor.__init__(self, dname, recipes=recipes, output_log=output_log, plotter=plotter)
            self.fast=fast
	
        def dispatch(self):
                master_bias = self.get_master_bias()
                for n,r in enumerate(self.get_recipes(rtype='comp')):
                    self.output(message='', progress='Processing '+self.dname+' comps: '+str(n+1)+'/'+str(len(self.get_recipes(rtype='flat')))+' |',add_to_log=False)

                    #Load previously generated calibration products.
                    self.load_calib_products(r.pnum)

                    #Set plotter to save plots in the correct plot directory.
                    self.plotter.set_rootpath(self.plot_dirs[r.pnum])

                    #Generate wavelength solutions for each comp frame, save to calib directory.
                    wvlsol_maps = []
                    summary = {}
                    for i,fname in enumerate(r.filenames):
                        self.output('Finding wavelength solution from '+fname)
                        comp = fits.open(self.indata+'/'+fname)
                        comp = calibrate(comp, master_bias, self.fiber_mask, lacosmic=False)
                        comp[0].data = comp[0].data / self.throughput_map
                        w = wvlsolver(comp, self.fiber_mask, r.fibers, self.profile_map, fast=self.fast, output=self.output_log, plotter=self.plotter)
                        w.set_path(self.calib_dirs[r.pnum]+'/'+fname.split('.fits')[0]+'_summary.dat')
                        w.solve()

                        #Add info to wavelength solution summary
                        npeaks_summary = w.get_fiber_npeaks()
                        for fnum in npeaks_summary.keys():
                            if not fnum in summary:
                                summary[fnum] = []
                            summary[fnum].append(npeaks_summary[fnum])

                        #Save wavelength solution
                        wvlsol_map = w.get_wvlsol_map()
                        wvlsol_maps.append(wvlsol_map)
                        ws_path = self.calib_dirs[r.pnum]+'/wvlsol_'+fname.split('.')[0]+'.fits'
                        fits.writeto(ws_path, wvlsol_map, clobber=True)
                        self.output('Wavelength solution derived from '+fname+' saved at '+ws_path)

                    #Save summary of wavelength solutions
                    f = open(self.calib_dirs[r.pnum]+'/wvlsol_summary.dat', 'w')
                    for fnum in sorted(summary.keys(), key=lambda n: int(n)):
                        f.write(','.join([str(n) for n in ([fnum]+summary[fnum])])+'\n')
                    f.close()

                    #Generate master wavelength solution as average of the individual wavelength solutions,
                    # save to calib directory
                    master_wvlsol = np.mean(wvlsol_maps, axis=0)
                    mws_path = self.calib_dirs[r.pnum]+'/'+self.cp_fnames['wavelength_solution']
                    fits.writeto(mws_path, master_wvlsol, clobber=True)
                    self.output('Master Wavelength solution saved at '+mws_path)
		
class process_skyflat(processor):
    def __init__(self, dname, recipes=None, output_log=None, plotter=None):
        processor.__init__(self, dname, recipes=recipes, output_log=output_log, plotter=plotter)

    def dispatch(self):
        master_bias = self.get_master_bias()
        for n,r in enumerate(self.get_recipes(rtype='skyflat')):
            self.output(message='', progress='Processing '+self.dname+' skyflats: '+str(n+1)+'/'+str(len(self.get_recipes(rtype='skyflat')))+' |')

            #Set plotter to save plots in the correct plot directory
            self.plotter.set_rootpath(self.plot_dirs[r.pnum])

            #If an uncorrected wavelength solution exists, start using it.
            ws_path = self.calib_dirs[r.pnum]+'/'+self.cp_fnames['wavelength_solution']
            if os.path.exists(ws_path[:-5]+"_uncorr.fits"):
                shutil.copyfile(ws_path[:-5]+"_uncorr.fits", ws_path)

            #Load previously generated calibration products
            self.load_calib_products(r.pnum)

            #Extract skyflat spectra from each frame and group by fiber number.
            skyflat_frames = []
            for fname in r.filenames:
                self.output('Loading in '+fname)
                skyflat_frame = calibrate(fits.open(self.indata+'/'+fname), master_bias, lacosmic=False)
                skyflat_frame[0].data = skyflat_frame[0].data / self.throughput_map
                skyflat_frames.append(skyflat_frame)
            self.output('Extracting skyflat spectra.')
            skyflat_fibers = robust_mean_extraction(skyflat_frames, self.fiber_mask, self.profile_map, self.wavelength_solution, r.fibers)
            skyflat_fibers.scale_spectra()

            #Determine fiber RV offests and apply offsets to wavelength solution.
            # Load template solar spectrum.
            solar_w, solar_f = np.loadtxt("calib/master_calib/bass2000_6000_7000.txt", unpack=True, dtype=np.float32)
            solar_f /= 10000.0
            solar_f /= 0.5*(np.nanmax(solar_f)+np.nanmedian(solar_f))
            solar_f_smooth = convolve_gaussian(solar_w, solar_f, sig=.16) #IMPORT convolve_gaussian

            # Make dictionary to store rv offsets
            rv_dict = {}
            for fnum in skyflat_fibers.get_fiber_numbers():
                #Get the skyflat spectrum for this fiber.
                sp = skyflat_fibers.get_spectrum(fnum)
                wav = sp.get_wavelength()
                flux = sp.get_flux()
                wav_flat, flux_flat = flatten_spec(wav, flux, twav=solar_w, tflux=solar_f_smooth)

                #Cross correlate with solar spectrum to get rv offset.
                rv = get_rv(wav_flat, flux_flat, solar_w, solar_f_smooth, logspace=True)
                if abs(rv) > 5. or fnum==50:
                    wav_flat, flux_flat = flatten_spec(wav, flux, twav=solar_w, tflux=solar_f_smooth, plot=True)
                rv_dict[fnum]=-rv #km/s that is applied to the observed skyflat.

            #Go through all pointings and apply the RV offsets.
            for rec in self.get_recipes(rtype="comp"):
                self.apply_offsets(rec, rv_dict)

            #Save the RV offsets.
            rv_path = self.calib_dirs[r.pnum]+'/rv_offsets.txt'
            np.savetxt(rv_path, rv_dict.items(), fmt=["%d","%.8f"])
            
            #Save the skyflat spectra.
            skyflat_fibers.save(self.calib_dirs[r.pnum]+'/skyflat_spec.fits')
            for sp in skyflat_fibers.get_spectra():
                self.plotter.clear()
                sp.plot(p=self.plotter)

    def apply_offsets(self, r, rv_dict):
        #If an uncorrected wavelength solution exists, start using it.
        ws_path = self.calib_dirs[r.pnum]+'/'+self.cp_fnames['wavelength_solution']
        if os.path.exists(ws_path[:-5]+"_uncorr.fits"):
            shutil.copyfile(ws_path[:-5]+"_uncorr.fits", ws_path)
        self.load_calib_products(r.pnum)

        #Apply RV shifts
        wvlsol_corr = self.wavelength_solution.copy()
        for fnum, rv in rv_dict.items():
            wvlsol_corr[self.fiber_mask == fnum] *= (1+rv/3.e5)

        #Save new wavelength solution.
        ws_path = self.calib_dirs[r.pnum]+'/'+self.cp_fnames['wavelength_solution'] 
        if not os.path.exists(ws_path[:-5]+"_uncorr.fits"):
            shutil.copyfile(ws_path, ws_path[:-5]+"_uncorr.fits")
        fits.writeto(ws_path, wvlsol_corr, clobber=True)



class process_sky(processor):
        def __init__(self, dname, recipes=None, output_log=None, plotter=None):
            processor.__init__(self, dname, recipes=recipes, output_log=output_log, plotter=plotter)
	
        def dispatch(self):
            master_bias = self.get_master_bias()

            for i,r in enumerate(self.get_recipes(rtype='sky')):
                self.output(message='', progress='Processing '+self.dname+' flats: '+str(i+1)+'/'+str(len(self.get_recipes(rtype='flat')))+' |')

                #Set plotter to save plots in the correct plot directory
                self.plotter.set_rootpath(self.plot_dirs[r.pnum])

                #Check that there are actually sky fibers in need of reducing
                if len(r.filenames) == 0:
                    self.output('No sky fibers to reduce. Moving on.')
                    continue

                #Load previously generated calibration products
                self.load_calib_products(r.pnum)

                #Make master sky frame, save to calib directory
                self.output('Loading sky frames.')

                #Load in and calibrate frames.
                sky_frames = []
                sky_headers = []
                for fname in r.filenames:
                    sky_frame = fits.open(self.indata+'/'+fname)
                    sky_frame = calibrate(sky_frame, master_bias, self.fiber_mask, lacosmic=False)
                    sky_frame[0].data /= self.throughput_map
                    sky_frames.append(sky_frame)

                #Extract sky spectra, save to outdata directory
                sky_fibers = robust_mean_extraction(sky_frames, self.fiber_mask, self.profile_map, self.wavelength_solution, r.fibers)
                #sky_fibers.scale_spectra()
                for sp in sky_fibers.get_spectra():
                    sp.plot(p=self.plotter, color='gray')

                #Make master sky spectrum by median combining all fibers. Save to calib directory.
                self.output('Producing master sky spectrum')
                master_sky_spec = median_spectra(sky_fibers.get_spectra())
                master_sky_spec.plot(p=self.plotter)
                self.plotter.save('master_sky.pdf')
                mss_path = self.calib_dirs[r.pnum]+'/'+self.cp_fnames['master_sky_spec']
                master_sky_spec.save(mss_path)
                self.output('Master sky spectrum saved at '+mss_path)

class process_target(processor):
        def __init__(self, dname, recipes=None, output_log=None, plotter=None):
            processor.__init__(self, dname, recipes=recipes, output_log=output_log, plotter=plotter)
	
        def dispatch(self):
            master_bias = self.get_master_bias()
            for i,r in enumerate(self.get_recipes(rtype='object')):
                self.output(message='', progress='Processing '+self.dname+' targets: '+str(i+1)+'/'+str(len(self.get_recipes(rtype='object')))+' |')

                #Set plotter to save plots in the outdata directory.
                self.plotter.set_rootpath(self.outdata_dirs[r.pnum])

                #Check that there are actually target fibers in need of reducing.
                if len(r.filenames) == 0:
                    self.output('No target fibers to reduce. Moving on.')
                    continue

                #Load previously generated calibration products.
                self.load_calib_products(r.pnum)

                #Load in target frames.
                self.output('Loading target frames.')
                target_frames = []
                for fname in r.filenames:
                    target_frame = fits.open(self.indata+'/'+fname)
                    target_frame = calibrate(target_frame, master_bias, self.fiber_mask, lacosmic=False)
                    target_frame[0].data /= self.throughput_map
                    target_frames.append(target_frame)

                #Extract the 1D spectrum.
                target_fibers = robust_mean_extraction(target_frames, self.fiber_mask, self.profile_map, self.wavelength_solution, r.fibers)

                #Save output to the outdata directory.
                ts_path = self.outdata_dirs[r.pnum]+'/target_spectra.fits'
                target_fibers.save(ts_path)
                self.output('Target spectrum saved as '+ts_path)

                #Do sky subtraction if master sky spectrum exists (was created in process_sky).
                #Save sky-subtracted spectrum in the outdata directory.
                try:
                    self.output('Subtracting sky spectrum from target spectra.')
                    for fnum, spec in target_fibers.get_spectra(as_dict=True).items():
                        spec = target_fibers.get_spectrum(fnum) - self.master_sky_spec

                        self.plotter.clear()
                        target_fibers.get_spectrum(fnum).plot(p=self.plotter, line_width=2)
                        spec.plot(p=self.plotter, line_color='red', line_width=2)

                        target_fibers.set_spectrum(fnum, spec)
                        self.plotter.save('spec_fib'+str(fnum)+'.html')
                    ts_path = self.outdata_dirs[r.pnum]+'/target_spectra_nosky.fits'
                    target_fibers.save(ts_path)
                    self.output('Target spectrum saved as '+ts_path)
                except AttributeError:
                    self.output('No master sky spectrum found. Skipping sky subtraction.')
                    for fnum, spec in target_fibers.get_spectra(as_dict=True).items():
                        self.plotter.clear()
                        target_fibers.get_spectrum(fnum).plot(p=self.plotter, line_width=2)
                        self.plotter.save('spec_fib'+str(fnum)+'.html')
                
def get_recipes(dname, recipe=None, pnum=None, rtype=None):
    if recipe == None:
        recipe = 'recipes/'+dname+'.recipe'

    #Extract lines from recipe file that pertain to flat fields.
    recipe_lines = [line for line in filter(None, open(recipe).read().split('\n')) if line[0] != '#']
    if rtype != None:
        recipe_lines = [line for line in recipe_lines if len(line.split(',')) > 2 and line.split(',')[1] == rtype]
    if pnum != None:
        pnum = int(pnum)
        recipe_lines = [line for line in recipe_lines if len(line.split(',')) > 1 and int(line.split(',')[0]) == pnum]
    return recipe_lines

def make_master_bias(dname, recipe=None, output=stdout):
    bias_recipes = get_recipes(dname, recipe, rtype='zero')

    calib_dir = 'calib/'+dname
    ensure_path(calib_dir+'/')

    if not os.path.exists(calib_dir+'/'+'master_bias.fits'):
        filenames = []
        for r in bias_recipes:
            info = r.split(',')
            fnames = info[2].split(' ')
            filenames.extend(fnames)

        if len(filenames) > 0:
            #Load files and create master bias
            biases = [fits.open('indata/'+dname+'/'+fname) for fname in filenames]
            master_bias = combine(*biases, method='mean')
            master_bias.writeto(calib_dir+'/master_bias.fits', clobber=True)
            for f in biases:
                f.close()

'''
def process_flat(dname, recipe=None, output=stdout):
    flat_recipes = get_recipes(dname, recipe, rtype='flat')
    num_r = len(flat_recipes)

    for i,r in enumerate(flat_recipes):
        output.edit_progress('Processing '+dname+' flats: '+str(i+1)+'/'+str(num_r)+' |')
        output.edit_message("", add_to_log=False)
        flat_dorecipe(r, dname, recipe, output)
'''
def flat_dorecipe(r, dname, recipe, output=None):
    #Unpack info from recipe line.
    output.edit_message('Reading recipe.')
    info = r.split(',')
    group_num = info[0]
    rtype = info[1]
    filenames = info[2].split(' ')
    use_fibers = [int(f_num) for f_num in info[3].split(' ')]

    #Make directory in calib for this group number
    indata_dir = 'indata/'+dname
    calib_dir = 'calib/'+dname+'/group'+group_num
    ensure_path(calib_dir+'/')

    #Obtain master bias frame
    make_master_bias(dname, recipe)
    mb_path = 'calib/'+dname+'/master_bias.fits'
    if os.path.exists(mb_path):
        mb = fits.open(mb_path)
        master_bias = mb[0].data
        mb.close()
    else:
        master_bias = None

    #Load flat frames, bias correct, and median combine flats, save to calib directory
    output.edit_message('Loading flat frames.')
    flats = [fits.open(indata_dir+'/'+filename) for filename in filenames]
    output.edit_message('Median combining flat frames.')
    master_flat = combine(*flats)
    for flat in flats:
        flat.close()
    mf_path = calib_dir+'/master_flat.fits'
    master_flat.writeto(mf_path, clobber=True)
    output.edit_message('Master flat frame saved at '+mf_path)

    #Find fibers, make fiber mask, save to directory
    output.edit_message('Locating fibers.')
    fiber_mask = find_fibers(master_flat, use_fibers)
    fm_path = calib_dir+'/fiber_mask.fits'
    fits.writeto(fm_path, fiber_mask, clobber=True)
    output.edit_message('Fiber mask saved at '+fm_path)

    #Calibrate master flat
    output.edit_message('Bias correcting master flat frame.')
    master_flat = calibrate(master_flat, master_bias, fiber_mask, lacosmic=False)

    #Generate a fiber thoughput map.
    output.edit_message('Generating throughput map.')
    throughput_map = make_throughput_map(fiber_mask, master_flat)
    tm_path = calib_dir+'/throughput_map.fits'
    fits.writeto(tm_path, throughput_map, clobber=True)
    output.edit_message('Throughput map saved at '+tm_path)
    
    #Generate a fiber profile map.

'''
def process_thar(dname, recipe=None, output=stdout, **kwargs):
    thar_recipes = get_recipes(dname, recipe, rtype='comp')

    num_r = len(thar_recipes)
    for i,r in enumerate(thar_recipes):
        output.edit_progress('Processing '+dname+' thar: '+str(i+1)+'/'+str(num_r)+' |')
        output.edit_message('')
        thar_dorecipe(r, dname, output, **kwargs)
'''

def thar_dorecipe(r, dname, output=None, **kwargs):
    output.edit_message('Reading thar recipe.')
    #Unpack info from recipe line.
    info = r.split(',')
    group_num = info[0]
    rtype = info[1]
    filenames = info[2].split(' ')
    use_fibers = [int(f_num) for f_num in info[3].split(' ')]

    #Make directory in calib for this group number
    indata_dir = 'indata/'+dname
    calib_dir = 'calib/'+dname+'/group'+group_num
    ensure_path(calib_dir+'/')

    #Obtain master bias frame
    mb_path = 'calib/'+dname+'/master_bias.fits'
    if os.path.exists(mb_path):
        mb = fits.open(mb_path)
        master_bias = mb[0].data
        mb.close()
    else:
        master_bias = None

    tm = fits.open(calib_dir+'/throughput_map.fits')
    throughput_map = tm[0].data
    tm.close()

    fm = fits.open(calib_dir+'/fiber_mask.fits')
    fiber_mask = fm[0].data
    fm.close()

    wvlsol_maps = []
    for i,fname in enumerate(filenames):
        output.edit_message('Finding wavelength solution from '+fname)
        comp = fits.open(indata_dir+'/'+fname)
        comp = calibrate(comp, master_bias, fiber_mask)
        comp[0].data = comp[0].data / throughput_map
        wvlsol_map = wvlsol(comp, fiber_mask, use_fibers, **kwargs) #NEED OPTIMAL EXTRACTION IN THERE
        wvlsol_maps.append(wvlsol_map)
        ws_path = calib_dir+'/wvlsol_'+fname.split('.')[0]+'.fits'
        fits.writeto(ws_path, wvlsol_map, clobber=True)
        output.edit_message('Wavelength solution derived from '+fname+' saved at '+ws_path)
    master_wvlsol = np.mean(wvlsol_maps, axis=0)
    mws_path = calib_dir+'/wvlsol.fits'
    fits.writeto(mws_path, master_wvlsol, clobber=True)
    output.edit_message('Master Wavelength solution saved at '+mws_path)

'''
def process_sky(dname, recipe=None, output=stdout):
    #output = output_log(log_path='calib/'+dname+'/output.log')
    sky_recipes = get_recipes(dname, recipe, rtype='sky')

    num_r = len(sky_recipes)
    for i,r in enumerate(sky_recipes):
        output.edit_progress('Processing '+dname+' sky: '+str(i+1)+'/'+str(num_r)+' |')
        output.edit_message('')
        sky_dorecipe(r, dname, output)
'''

def sky_dorecipe(r, dname, output=None):
    #Unpack info from the recipe line.
    output.edit_message('Reading sky recipe.')
    info = r.split(',')
    group_num = info[0]
    rtype = info[1]
    filenames = info[2].split(' ')
    use_fibers = [int(f_num) for f_num in info[3].split(' ')]

    #If there are no fibers to be reduced, stop.
    if len(use_fibers) == 0:
        output.edit_message('No fibers to be reduced. Stopping.')
        return

    #Make directory in calib for this group number
    indata_dir = 'indata/'+dname
    calib_dir = 'calib/'+dname+'/group'+group_num
    ensure_path(calib_dir+'/')

    #Obtain master bias frame, if it exists.
    mb_path = 'calib/'+dname+'/master_bias.fits'
    if os.path.exists(mb_path):
        mb = fits.open(mb_path)
        master_bias = mb[0].data
        mb.close()
    else:
        master_bias = None

    #Load in the fiber mask.
    fm = fits.open(calib_dir+'/fiber_mask.fits')
    fiber_mask = fm[0].data
    fm.close()

    #Load in the master flat frame.
    mf = fits.open(calib_dir+'/master_flat.fits')
    master_flat = mf[0].data
    mf.close()

    #Load in the throughput map.
    tm = fits.open(calib_dir+'/throughput_map.fits')
    throughput_map = tm[0].data
    tm.close()

    #Load in the wavelength solution.
    ws = fits.open(calib_dir+'/wvlsol.fits')
    wvlsol_map = ws[0].data
    ws.close()

    #Make master sky frame
    output.edit_message('Loading sky frames.')
    skys = [fits.open(indata_dir+'/'+filename) for filename in filenames]
    output.edit_message('Median combining sky frames.')
    master_sky = combine(*skys)
    master_sky.writeto('plots/lacosmic/master_sky.fits', clobber=True)
    for sky in skys:
        sky.close()
    output.edit_message('Bias correcting master sky frame.')
    master_sky = calibrate(master_sky, master_bias, fiber_mask)
    output.edit_message('Throughput correcting master sky frame.')
    #master_sky[0].data /= throughput_map
    ms_path = calib_dir+'/master_sky.fits'
    master_sky.writeto(ms_path, clobber=True)
    output.edit_message('Master sky frame saved at '+ms_path)

    output.edit_message('Extracting sky spectra.')
    #sky_spec = extract(fiber_mask, fnum, master_sky, wvlsol_map) # NEED OPTIMAL EXTRACTION IN THERE
    sky_fibers = optimal_extraction(master_sky, fiber_mask, use_fibers, master_flat, wvlsol_map)
    output.edit_message('Producing master sky spectrum')
    ss_path = outdata_dir+'/sky_spectra.fits'
    sky_fibers.save(ss_path)
    output.edit_message('Sky spectra saved at '+ss_path)
    master_sky_spec = median_spectra(sky_fibers.get_spectra())

    mss_path = calib_dir+'/master_sky_spec.dat'
    master_sky_spec.save(mss_path)
    output.edit_message('Master sky spectrum saved at '+mss_path)

'''
def process_target(dname, recipe=None, output=stdout):
    #output = output_log(log_path='calib/'+dname+'/output.log')
    tar_recipes = get_recipes(dname, recipe, rtype='object')

    num_r = len(tar_recipes)
    for i,r in enumerate(tar_recipes):
        output.edit_progress('Processing '+dname+' targets: '+str(i+1)+'/'+str(num_r)+' |')
        output.edit_message('')
        target_dorecipe(r, dname, output)
'''

def target_dorecipe(r, dname, output=None):
    #Unpack info from the recipe line.
    output.edit_message('Reading target recipe.')
    info = r.split(',')
    group_num = info[0]
    rtype = info[1]
    filenames = info[2].split(' ')
    use_fibers = [int(f_num) for f_num in filter(None, info[3].split(' '))]

    #If there are no fibers to be reduced, stop.
    if len(use_fibers) == 0:
        output.edit_message('No fibers to be reduced. Stopping.')
        return

    #Define paths to indata, calib, and outdata
    indata_dir = 'indata/'+dname
    calib_dir = 'calib/'+dname+'/group'+group_num
    ensure_path(calib_dir+'/')
    outdata_dir = 'outdata/'+dname+'/group'+group_num
    ensure_path(outdata_dir+'/')

    #Obtain master bias frame, if it exists.
    mb_path = 'calib/'+dname+'/master_bias.fits'
    if os.path.exists(mb_path):
        mb = fits.open(mb_path)
        master_bias = mb[0].data
        mb.close()
    else:
        master_bias = None

    #Load in the fiber mask
    fm = fits.open(calib_dir+'/fiber_mask.fits')
    fiber_mask = fm[0].data
    fm.close()

    #Load in the wavelength solution.
    ws = fits.open(calib_dir+'/wvlsol.fits')
    wvlsol_map = ws[0].data
    ws.close()

    #Load in the throughput map.
    tm = fits.open(calib_dir+'/throughput_map.fits')
    throughput_map = tm[0].data
    tm.close()

    #Load in the master flat frame.
    mf = fits.open(calib_dir+'/master_flat.fits')
    master_flat = mf[0].data
    mf.close()

    #master_sky_spec.plot()

    #Make master target frame
    output.edit_message('Loading target frames.')
    tars = [fits.open(indata_dir+'/'+filename) for filename in filenames]
    output.edit_message('Median combining target frames.')
    master_tar = combine(*tars)
    for f in tars:
        f.close()
    output.edit_message('Bias correcting master target frame.')
    master_tar = calibrate(master_tar, master_bias)
    output.edit_message('Throughput correcting master target frame.')
    #master_tar[0].data /= throughput_map
    mt_path = calib_dir+'/master_target_frame.fits'
    master_tar.writeto(mt_path, clobber=True)
    output.edit_message('Master target frame saved at '+mt_path)

    header = master_tar[0].header
    output.edit_message('Extracting target spectra.')
    target_fibers = optimal_extraction(master_tar, fiber_mask, use_fibers, master_flat, wvlsol_map)
    ts_path = outdata_dir+'/target_spectra.fits'
    target_fibers.save(ts_path)
    output.edit_message('Target spectrum saved as '+ts_path)
    if os.path.exists(calib_dir+'/master_sky_spec.dat'):
    	#Load in master sky spectrum.
    	mss = fits.open(calib_dir+'/master_sky_spec.dat')
    	master_sky_spec = spectrum(mss[1].data, mss[0].data, mss[2].data)
    	for i in range(len(target_fibers.get_spectra())):
        	spec = target_fibers[i] - master_sky_spec
        	target_fibers[i] = spec
    	ts_path = outdata_dir+'/target_spectra_nosky.fits'
    	target_fibers.save(ts_path)
    	output.edit_message('Target spectrum saved as '+ts_path)
