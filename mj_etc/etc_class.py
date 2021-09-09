import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const

def mag2flux(mag):
	flux = 10**(-(mag+48.6)/2.5)
	flux = flux * u.erg / u.s / u.cm**2 / u.Hz
	return flux

class ETC():

	def __init__(self,detector,telescope,source_mag=np.nan,moon_phase=0,seeing=2):
		self.detector = detector
		self.telescope = telescope
		self.seeing = seeing
		self.sky = Sky(moon_phase, detector, telescope)
		self.source = Source(source_mag, detector, telescope)
		self.exp_time = np.arange(.1,1000,1) * u.s
		# calculated
		self.pixels = np.nan
		self.noise = np.nan
		self.signal = np.nan
		self.snr = np.nan

	def calc_sky_photons(self):
		self.sky.sky_signal()

	def calc_source_photons(self):
		self.source.source_signal()

	def aperture_size(self):
		pix = ceil(self.seeing / self.detector.platescale)
		pix2 = pix**2
		self.pixels = pix2

	def calculate_noise(self):
		dc = self.detector.dark_current * self.exp_time * self.pixels**2
		sky = np.sqrt(self.sky.sky_photons * self.telescope.throughput * self.exp_time)
		source = np.sqrt(self.source.source_photons * self.telescope.throughput * self.exp_time)
		read = self.detector.read_noise * self.pixels**2

		noise = np.sqrt(dc**2+sky**2+source**2+read**2)
		self.noise = noise

	def calculate_signal(self):
		source = self.source.source_photons * self.telescope.throughput * self.exp_time
		self.signal = source

	def calculate_SNR(self):
		self.calculate_noise()
		self.calculate_signal()
		self.snr = self.signal / self.noise


	def time_for_snr(self,snr,plot=False):
		time = np.arange(.1,1000,1) * u.s
		# should be updated with a loop
		self.time = time 
		self.calculate_SNR()

		ind = np.argmin(abs(self.snr - snr))
		if plot:
			self.snr_plot(snr)
		#print('Time needed to reach $SNR= {}$'.format(snr)+' is ' + str(time[ind]))
		return time[ind]

	def snr_plot(self,snr):
		plt.figure()
		plt.plot(self.time,self.snr)
		plt.axhline(snr,ls='--')
		plt.ylabel('SNR')
		plt.xlabel('Exposure time [s]')

	def snr_for_time(self,time):
		self.exp_time = time 
		self.calculate_SNR()
		return self.snr









class Telescope():
	def __init__(self,name,diameter,throughput):
		self.name = name
		self.diameter = diameter
		self.throughput = throughput
		self.area = np.pi*(diameter/2)**2



sky_brightness_table = 'path'
phases = np.array([])

class Sky():
	def __init__(self,moon_phase,detector,telescope):
		self.detector = detector
		self.moon_phase = moon_phase
		self.sky_brightness = self._get_sky()
		self.photons = np.nan


	def _get_band_brightness():
		tab = pd.read_csv(sky_brightness_table)
		ind = np.where(self.detector.bandpass.name == tab['name'].values)[0]
		if len(ind) == 0:
			m = 'No such filter!'
			raise ValueError(m)
		if len(ind) > 0:
			m = 'Multiple filter definitions!'
			raise ValueError(m)
		return tab.iloc[ind]

	def _get_sky(self):
		tab = self._get_band_brightness()
		ind = np.argmin(abs(self.moon-phases))
		return tab[ind]

	def sky_signal(self):
		sky = mag2flux(mag)
		# multiply by plate scale
		flux = sky * self.detector.platescale**2
		#convert to SI units 
		flux = flux.to(u.J/u.s/u.Hz/u.m**2)
		# convert to Lambda [J/s/m/m^2/pix]
		flux_lam = flux * const.c / (self.detector.bandpass.wavelength.to(u.m)**2)
		# multiply by bandwidth [J/s/m^2/pix]
		flux_lam = flux_lam * self.detector.bandpass.bandwidth.to(u.m)
		# multiply by size of telescope [J/s/pix]
		flux_lam = flux_lam * self.telescope.area
		sky_photons = flux_lam / ((const.h * const.c)/self.detector.bandpass.wavelength.to(u.m))
		self.photons = sky_photons



class Source():
	def __init__(self,mag,detector,telescope):
		# defined
		self.mag = mag 
		self.detector = detector
		self.telescope = telescope
		# calculated
		self.photons = np.nan

	def source_signal(self):
		source = mag2flux(self.mag)
		#convert to SI units 
		flux = flux.to(u.J/u.s/u.Hz/u.m**2)
		# convert to Lambda [J/s/m/m^2]
		flux_lam = flux * const.c / (self.detector.bandpass.wavelength.to(u.m)**2)
		# multiply by bandwidth [J/s/m^2]
		flux_lam = flux_lam * self.detector.bandpass.bandwidth.to(u.m)
		# multiply by size of pixel [J/s]
		flux_lam = flux_lam * self.telescope.area
		photons = flux_lam / ((const.h * const.c)/self.detector.bandpass.wavelength.to(u.m))
		self.photons = photons






class Detector():
	def __init__(self,dark_current,read_noise,gain,qe,platescale,pixel_size,bandpass):
		self.dark_current = dark_current
		self.read_noise = read_noise
		self.qe = qe
		self.gain = gain
		self.pix_size = pixel_size
		self.bandpass = bandpass


		self.check_inputs()

	
	def check_inputs(self):
		if (type(self.dark_current) is not float) | (type(self.dark_current) is not int):
			m = 'dark_current must be a float or integer'
			raise ValueError(m)
		if (type(self.read_noise) is not float) | (type(self.read_noise) is not int):
			m = 'read_noise must be a float or integer'
			raise ValueError(m)
		if type(self.pix_size) is not type(1*u.s):
			m = 'pixel_size must have an astropy unit'
			raise ValueError(m)
		if type(self.bandpass) is not type(Bandpass):
			m = 'bandpass must be the Bandpass class'
			raise ValueError(m)


	def integrated_dark_current(self,time):
		dc = self.dark_current * time 
		return dc 



class Bandpass():
	def __init__(self,name):
		self.name = name
		self.bandwidth = self._get_bandwidth()
		self.wavelength = self._get_wavelength()

	def _get_index():
		bandpass_table = 'path'
		tab = pd.read_csv(bandpass_table)
		ind = np.where(self.name == tab['name'].values)[0]
		if len(ind) == 0:
			m = 'No such filter!'
			raise ValueError(m)
		if len(ind) > 0:
			m = 'Multiple filter definitions!'
			raise ValueError(m)
		return tab.iloc[ind]

	def get_bandwidth(self):
		tab = self._get_entry()
		bw = tab['bandwidth'] * u.nm
		return bw 

	def get_wavelength(self):
		tab = self._get_entry()
		wav = tab['wavelength'] * nm
		return wav

