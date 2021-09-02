import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const

def mag2flux(mag):
	flux = 10**(-(mag+48.6)/2.5)
	flux = flux * u.erg / u.s / u.cm**2 / u.Hz
	return flux

class ETC():

	def __init__(self,mag,detector,telescope,moon_phase=0):
		





class Telescope():
	def __init__(self,diameter,throughput):
		self.diameter = diameter
		self.throughput = throughput
		self.area = np.pi*(diameter/2)**2



sky_brightness_table = 'path'
phases = np.array([])

class Sky():
	def __init__(self,detector,telescope,moon_phase):
		self.detector = detector
		self.moon_phase = moon_phase
		self.sky_brightness = self._get_sky()


	def _get_band_brightness():
		tab = pd.read_csv(sky_brightness_table)
		ind = np.where(self.detector.bandpass.bandpass == tab['name'].values)[0]
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
		# convert to Lambda [J/s/m/m^2]
		flux_lam = flux * const.c / (self.detector.bandpass.wavelength**2)
		# multiply by bandwidth [J/s/m^2]
		flux_lam = flux_lam * self.detector.bandpass.bandwidth
		# multiply by size of pixel [J/s]
		flux_lam = flux_lam * self.telescope.area
		sky_photons = flux_lam / ((const.h * const.c)/self.detector.bandpass.wavelength)
		return sky_photons



class Source():
	def __init__(self,mag,detector,telescope):
		self.mag = mag 
		self.detector = detector
		

	def source_signal(self):
		sky = mag2flux(self.mag)
		# multiply by plate scale
		flux = sky * self.detector.platescale**2
		#convert to SI units 
		flux = flux.to(u.J/u.s/u.Hz/u.m**2)
		# convert to Lambda [J/s/m/m^2]
		flux_lam = flux * const.c / (self.detector.bandpass.wavelength**2)
		# multiply by bandwidth [J/s/m^2]
		flux_lam = flux_lam * self.detector.bandpass.bandwidth
		# multiply by size of pixel [J/s]
		flux_lam = flux_lam * self.telescope.area
		sky_photons = flux_lam / ((const.h * const.c)/self.detector.bandpass.wavelength)
		return sky_photons






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

bandpass_table = 'path'

class Bandpass():
	def __init__(self,bandpass):
		self.bandpass = bandpass
		self.bandwidth = self._get_bandwidth()
		self.wavelength = self._get_wavelength()

	def _get_index():
		tab = pd.read_csv(bandpass_table)
		ind = np.where(self.bandpass == tab['name'].values)[0]
		if len(ind) == 0:
			m = 'No such filter!'
			raise ValueError(m)
		if len(ind) > 0:
			m = 'Multiple filter definitions!'
			raise ValueError(m)
		return tab.iloc[ind]

	def _get_bandwidth(self):
		tab = self._get_entry()
		bw = tab['bandwidth'] * u.Angstrom
		return bw 

	def _get_wavelength(self):
		tab = self._get_entry()
		wav = tab['wavelength'] * u.Angstrom
		return wav
