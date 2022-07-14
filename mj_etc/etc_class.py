import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const
import pandas as pd

import os 

package_directory = os.path.dirname(os.path.abspath(__file__)) + '/'

def mag2flux(mag):
	"""
	Convert magnitude to flux
	
	Input:
		mag: magnitude

	Output:
		flux: flux in erg/s/cm^2/A
	"""
	flux = 10**(-(mag+48.6)/2.5)
	flux = flux * u.erg / u.s / u.cm**2 / u.Hz
	return flux

class ETC():

	def __init__(self,bandpass,detector,telescope,source_mag=np.nan,moon_phase='dark',seeing=2):
		"""
		Class to calculate the signal to noise ratio of a source in a given bandpass for telescopes and detectors at Mt John.

		Input:
			bandpass: name of the bandpass
			detector: name of the detector
			telescope: name of the telescope
			source_mag: magnitude of the source
			moon_phase: phase of the moon (dark, grey, or bright)
			seeing: seeing of the telescope

		"""
		self.detector = self._set_detector(detector,bandpass)
		self.telescope = self._set_telescope(telescope)
		self.source_mag = source_mag
		self.seeing = seeing
		self.aper = self.aperture_size()
		self.aper_pix = self.aperture_size_pix()
		self.sky = Sky(moon_phase, self.aper, self.detector, self.telescope)
		self.source = Source(source_mag, self.detector, self.telescope)
		self.exp_time = np.nan
		# calculated
		self.noise = np.nan
		self.signal = np.nan
		self.snr = np.nan


	def _get_detector_params(self,name):
		"""
		Get the detector parameters from the detectors.csv file.
		
		Input:
			name: name of the detector
		"""
		table = pd.read_csv(package_directory + 'camera_params.csv')
		ind = np.where(name.upper() == table['name'].values)[0]
		if len(ind) == 0:
			m = 'Detector {} does not exist'.format(name)
			raise ValueError(m)
		params = table.iloc[ind]
		return params

	def _get_bandpass(self,bandpass):
		"""
		Set the bandpass object.
		
		Input:
			bandpass: name of the bandpass
		"""
		b = Bandpass(bandpass)
		return b

	def _set_detector(self,name,bandpass):
		"""
		Set the detector object.

		Input:
			name: name of the detector
			bandpass: name of the bandpass
		"""
		p = self._get_detector_params(name)
		b =  self._get_bandpass(bandpass)
		d = Detector(name = p['name'].values[0], dark_current = p['dark_current'].values[0] / u.s, 
					 read_noise = p['read_noise'].values[0], gain = p['gain'].values[0], 
					 pixel_size = p['pixel_size'].values[0] * u.um, bandpass = b)
		return d

	def _get_telescope_params(self,name):
		"""
		Get the telescope parameters from the telescopes.csv file.

		Input:
			name: name of the telescope
		"""
		table = pd.read_csv(package_directory + 'telescope_params.csv')
		try:
			ind = np.where(name == table['name'].values)[0][0]
		except:
			raise ValueError('Telescope {} does not exist'.format(name))
		params = table.iloc[ind]
		return params
	
	def _set_telescope(self,name):
		"""
		Set the telescope object.
		"""
		p = self._get_telescope_params(name)
		t = Telescope(name = p['name'], diameter = p['diameter'], 
					  throughput = p['throughput'], f_num = p['f_num'])
		return t

	def print_status(self):
		"""
		Print the status of the ETC object.
		"""
		print('Detector: {}'.format(self.detector.name))
		print('Filter: {}'.format(self.detector.bandpass.name))
		print('Detector qe: {}'.format(self.detector.qe))
		print('Telescope: {}'.format(self.telescope.name))
		print('Diameter: {}'.format(self.telescope.diameter))
		print('Telescope throughput: {}'.format(self.telescope.throughput))
		print('Source: {}'.format(self.source.mag))
		print('Aperture: {}'.format(self.aper))
		print('Sky: {}'.format(self.sky.sky_brightness))
		print('Sky photons: {}'.format(self.sky.photons))
		print('Source photons: {}'.format(self.source.photons))
		print('\n')


	def calc_sky_photons(self):
		"""
		Calculate the expected number of photons from the sky.
		"""
		self.sky.sky_signal()

	def calc_source_photons(self):
		"""
		Calculate the expected number of photons from the input source.
		"""
		self.source.source_signal()

	def aperture_size_pix(self):
		"""
		Calculate the aperture size of the telescope in the given bandpass, telescope & detector setup.
		"""
		ps = self.telescope.platescale * self.detector.pix_size.to(u.mm).value
		pix2 = self.aper / ps**2
		return pix2

	def aperture_size(self):
		"""
		Calculate the aperture size of the telescope in the given bandpass, telescope & detector setup.

		Output:
			pix2: aperture size in pixels^2
		"""
		pix2 = np.pi*(self.seeing * 1.5)**2
		return pix2

	def calculate_noise(self):
		"""
		Calculate the noise of the source in the given bandpass, telescope & detector setup.
		"""
		self.sky.sky_signal()
		self.source.source_signal()
		self._dc_noise = np.sqrt(self.detector.dark_current * self.exp_time * self.aper_pix
			 					 * self.detector.gain)
		self._sky_noise = np.sqrt(self.sky.photons * self.telescope.throughput 
								  * self.exp_time * self.detector.qe * self.detector.gain)
		self._source_noise = np.sqrt(self.source.photons * self.telescope.throughput 
				 			    	 * self.exp_time * self.detector.qe * self.detector.gain)
		self._read_noise = self.detector.read_noise * np.sqrt(self.aper_pix)
		noise = np.sqrt(self._dc_noise**2+self._sky_noise**2+self._source_noise**2+self._read_noise**2)
		self.noise = noise

	def calculate_signal(self):
		"""
		Calculate the expected signal from the input source.
		"""
		self.source.source_signal()
		source = (self.source.photons * self.telescope.throughput 
				  * self.exp_time * self.detector.qe)
		self.signal = source

	def calculate_SNR(self):
		"""
		Calculate the signal to noise ratio of the source in the given bandpass, telescope & detector setup.
		"""
		self.calculate_noise()
		self.calculate_signal()
		self.snr = self.signal / self.noise


	def time_for_snr(self,snr,mag=None,plot=False):
		"""
		Calculate the time needed to observe a given SNR.

		Inputs:
			snr: SNR to observe
			mag: magnitude of the source
			plot: if True, plot the time vs. SNR
		
		Outputs:
			time: time needed to observe the given SNR

		"""
		if mag is not None:
			self.source = Source(mag, self.detector, self.telescope)
		time = np.arange(.1,1e4,1) * u.s
		# should be updated with a loop
		self.exp_time = time 
		self.calculate_SNR()

		ind = np.argmin(abs(self.snr - snr))
		if ind < len(time)-10:
			if plot:
				t = time < 1.2*time[ind]
				self.snr_plot(time,snr,t,ind)
			#print('Time needed to reach $SNR= {}$'.format(snr)+' is ' + str(time[ind]))
			return time[ind]
		else:
			print('!!! No time found for SNR = {} !!!'.format(snr))
			return None

	def snr_plot(self,time,snr,t,ind):
		"""
		Plot the SNR vs time
		
		Inputs:
			time: time array in seconds
			snr: SNR array
			t: boolean array for time array
			ind: index of the time array where the SNR is closest to the desired SNR
		"""
		plt.figure()
		plt.plot(time[t],self.snr[t])
		plt.axvline(time[ind].value,color='k',linestyle=':')
		plt.text(0.3,0.2,'Time to SNR = {}'.format(time[ind]),transform=plt.gca().transAxes,fontsize=12)
		plt.axhline(snr,ls='--',color='r',label='$SNR={}$'.format(snr))
		plt.ylabel('Signal-to-noise ratio',fontsize=15)
		plt.xlabel('Exposure time [s]',fontsize=15)

	def snr_for_time(self,time,mag=None):
		"""
		Calculate the SNR for a given time.
		
		Input:
			time: time in seconds
			mag: source magnitude
		"""
		if mag is not None:
			self.source = Source(mag, self.detector, self.telescope)
		self.exp_time = time 
		self.calculate_SNR()
		return self.snr




class Telescope():
	def __init__(self,name,diameter,throughput,f_num):
		"""
		Class for telescope parameters.

		Inputs:
			name: string
			diameter: float
			throughput: float
			f_num: float
		"""
		self.name = name
		self.f_num = f_num
		self.diameter = diameter * u.m
		self.throughput = throughput
		self.area = np.pi*(self.diameter/2)**2
		self.platescale = self.calculate_plate_scale()

	def calculate_plate_scale(self):
		"""
		Calculate the plate scale of the telescope.
		"""
		p = 206265 / (self.f_num * self.diameter.to(u.mm).value) # 
		return p


class Sky():
	def __init__(self,moon_phase,aper,detector,telescope):
		"""
		Class for the sky background.

		Inputs:
			moon_phase: Moon phase in 'bright', 'grey', 'dark'
			aper: Aperture for photometry
			detector: Detector object
			telescope: Telescope object
		"""
		self.detector = detector
		self.telescope = telescope
		self.aper = aper
		self.moon_phase = moon_phase
		self.sky_brightness = self._get_sky()
		self.sky_signal()


	def _get_band_brightness(self):
		"""
		Table values taken from: 
		https://www.mso.anu.edu.au/~pfrancis/reference/reference/node4.html
		"""
		tab = pd.read_csv(package_directory + 'sky_brightness.csv')
		t = tab.iloc[self.detector.bandpass._tab_ind]
		return t

	def _get_sky(self):
		"""
		Calculate the sky brightness in the given bandpass.
		"""
		tab = self._get_band_brightness()
		sky = tab[self.moon_phase]
		return sky

	def sky_signal(self):
		"""
		Calculate the sky signal in the given bandpass in photons per second.
		"""
		#ps = self.telescope.platescale * self.detector.pix_size.to(u.mm).value
		# multiply by plate scale
		sky = self.sky_brightness  # this needs to be 
		# convert to mag
		flux = mag2flux(sky) * self.aper
		#convert to SI units 
		flux = flux.to(u.J/u.s/u.Hz/u.m**2)
		# convert to Lambda [J/s/m/m^2/pix]
		flux_lam = flux * const.c / (self.detector.bandpass.wavelength.to(u.m)**2)
		#flux_lam = flux_lam.to(u.J/ u.s/ u.m /u.m**2)
		# multiply by bandwidth [J/s/m^2/pix]
		flux_lam = flux_lam * self.detector.bandpass.bandwidth.to(u.m)
		# multiply by size of telescope [J/s/pix]
		flux_lam = flux_lam * self.telescope.area
		# divide by average energy of a photon 
		sky_photons = flux_lam / ((const.h * const.c)/self.detector.bandpass.wavelength.to(u.m))
		sky_photons = sky_photons.to(1/u.s)
		self.photons = sky_photons



class Source():
	def __init__(self,mag,detector,telescope):
		"""
		Class to calculate the source signal for a given magnitude.

		Inputs:
			mag: magnitude of the source
			detector: detector object
			telescope: telescope object
		"""
		# defined
		self.mag = mag 
		self.detector = detector
		self.telescope = telescope
		# calculated
		self.photons = np.nan

	def source_signal(self):
		"""
		Calculate the source signal in the given bandpass in photons per second.
		"""
		source = mag2flux(self.mag)
		#convert to SI units 
		flux = source.to(u.J/u.s/u.Hz/u.m**2)
		# convert to Lambda [J/s/m/m^2]
		flux_lam = flux * const.c / (self.detector.bandpass.wavelength.to(u.m)**2)
		# multiply by bandwidth [J/s/m^2]
		flux_lam = flux_lam * self.detector.bandpass.bandwidth.to(u.m)
		# multiply by size of pixel [J/s]
		flux_lam = flux_lam * self.telescope.area
		photons = flux_lam / ((const.h * const.c)/self.detector.bandpass.wavelength.to(u.m))
		photons = photons.to(1/u.s)
		self.photons = photons






class Detector():
	def __init__(self,name,dark_current,read_noise,gain,pixel_size,bandpass):
		'''
		Class defining the detector properties.

		Inputs:
			name: name of the detector
			dark_current: dark current in e-/pixel
			read_noise: read noise in e-/pixel
			gain: gain in e-/ADU
			pixel_size: pixel size in m
			bandpass: bandpass object
		'''
		self.name = name
		self.bandpass = bandpass
		self.dark_current = dark_current
		self.read_noise = read_noise
		self.qe = self._get_qe()
		self.gain = gain
		self.pix_size = pixel_size

		self.check_inputs()

	
	def check_inputs(self):
		"""
		Check the inputs to the detector class.
		"""
		if type(self.dark_current) is not type(1/u.s):
			m = 'dark_current must be a float or integer'
			raise ValueError(m)
		#if (type(self.read_noise) is not type(np.int64)) | (type(self.read_noise) is not int):
		#	m = 'read_noise must be a float or integer'
		#	raise ValueError(m)
		if type(self.pix_size) is not type(1*u.um):
			m = 'pixel_size must have an astropy unit'
			raise ValueError(m)
		#if type(self.bandpass) is not type(Bandpass):
		#	m = 'bandpass must be the Bandpass class'
		#	raise ValueError(m)


	def integrated_dark_current(self,time):
		"""
		Calculate the integrated dark current in the given bandpass.

		Inputs:
			time: time in seconds

		Outputs:
			dark_current: dark current in e-
		"""
		dc = self.dark_current * time 
		return dc 

	def _get_qe(self):
		"""
		Calculate the quantum efficiency of the detector.

		Outputs:
			qe: quantum efficiency
		"""
		tab = self.bandpass._table.iloc[self.bandpass._tab_ind]
		qe = tab[self.name + '_qe']
		if qe ==0:
			raise ValueError('!!! Filter is not available for this detector !!!')
		return qe



class Bandpass():
	def __init__(self,name):
		"""
		Class defining the bandpass properties.

		Inputs:
			name: name of the bandpass
		"""
		self.name = name
		self._table = self._load_table()
		self._tab_ind = self._get_index()

		self.bandwidth = self._get_bandwidth()
		self.wavelength = self._get_wavelength()

	def _load_table(self):
		"""
		Load the bandpass table.

		Outputs:
			table: pandas dataframe
		"""
		tab = pd.read_csv(package_directory + 'bandpass_params.csv')
		return tab

	def _get_index(self):
		"""
		Get the index of the bandpass in the table.

		Outputs:
			tab_ind: index of the bandpass in the table
		"""
		try:
			ind = np.where(self.name == self._table['name'].values)[0][0]
		except:
			raise ValueError('No such filter!')
		return ind

	def _get_bandwidth(self):
		"""
		Get the bandwidth of the bandpass.
		
		Outputs:
			bandwidth: bandwidth in nm
		"""
		bw = self._table.iloc[self._tab_ind]['bandwidth'] * u.nm
		return bw 

	def _get_wavelength(self):
		"""
		Get the wavelength of the bandpass.

		Outputs:
			wavelength: wavelength in nm
		"""
		wav = self._table.iloc[self._tab_ind]['wavelength'] * u.nm
		return wav
