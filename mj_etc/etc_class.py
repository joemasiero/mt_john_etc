import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const
import pandas as pd

import os 

package_directory = os.path.dirname(os.path.abspath(__file__)) + '/'

def mag2flux(mag):
	flux = 10**(-(mag+48.6)/2.5)
	flux = flux * u.erg / u.s / u.cm**2 / u.Hz
	return flux

class ETC():

	def __init__(self,bandpass,detector,telescope,source_mag=np.nan,moon_phase='dark',seeing=2):
		self.detector = self._set_detector(detector,bandpass)
		self.telescope = self._set_telescope(telescope)
		self.source_mag = source_mag
		self.seeing = seeing
		self.sky = Sky(moon_phase, self.detector, self.telescope)
		self.source = Source(source_mag, self.detector, self.telescope)
		self.exp_time = np.nan
		# calculated
		self.pixels = self.aperture_size()
		self.noise = np.nan
		self.signal = np.nan
		self.snr = np.nan


	def _get_detector_params(self,name):
		table = pd.read_csv(package_directory + 'camera_params.csv')
		ind = np.where(name.upper() == table['name'].values)[0]
		if len(ind) == 0:
			m = 'Detector {} does not exist'.format(name)
			raise ValueError(m)
		params = table.iloc[ind]
		return params

	def _get_bandpass(self,bandpass):
		b = Bandpass(bandpass)
		return b

	def _set_detector(self,name,bandpass):
		p = self._get_detector_params(name)
		b =  self._get_bandpass(bandpass)
		d = Detector(name = p['name'][0], dark_current = p['dark_current'][0] / u.s, 
					 read_noise = p['read_noise'][0], gain = p['gain'][0], 
					 pixel_size = p['pixel_size'][0] * u.um, bandpass = b)
		return d

	def _get_telescope_params(self,name):
		table = pd.read_csv(package_directory + 'telescope_params.csv')
		try:
			ind = np.where(name == table['name'].values)[0][0]
		except:
			raise ValueError('Telescope {} does not exist'.format(name))
		params = table.iloc[ind]
		return params
	
	def _set_telescope(self,name):
		p = self._get_telescope_params(name)
		t = Telescope(name = p['name'], diameter = p['diameter'], 
					  throughput = p['throughput'], f_num = p['f_num'])
		return t


	def calc_sky_photons(self):
		self.sky.sky_signal()

	def calc_source_photons(self):
		self.source.source_signal()

	def aperture_size(self):
		ps = self.telescope.platescale / self.detector.pix_size.to(u.mm).value
		pix = np.ceil(self.seeing / ps)
		pix2 = pix**2
		return pix2

	def calculate_noise(self):
		self.sky.sky_signal()
		self.source.source_signal()
		dc = self.detector.dark_current * self.exp_time * self.pixels**2
		sky = np.sqrt(self.sky.photons * self.telescope.throughput 
					  * self.exp_time * self.detector.qe)
		source = np.sqrt(self.source.photons * self.telescope.throughput 
				     	 * self.exp_time * self.detector.qe)
		read = self.detector.read_noise * self.pixels**2
		#print(dc)
		#print(sky)
		#print(source)
		#print(read)
		noise = np.sqrt(dc**2+sky**2+source**2+read**2)
		self.noise = noise

	def calculate_signal(self):
		self.source.source_signal()
		source = (self.source.photons * self.telescope.throughput 
				  * self.exp_time * self.detector.qe)
		self.signal = source

	def calculate_SNR(self):
		self.calculate_noise()
		self.calculate_signal()
		self.snr = self.signal / self.noise


	def time_for_snr(self,snr,mag=None,plot=False):
		if mag is not None:
			self.source = Source(mag, self.detector, self.telescope)
		time = np.arange(.1,1e5,1) * u.s
		# should be updated with a loop
		self.exp_time = time 
		self.calculate_SNR()

		ind = np.argmin(abs(self.snr - snr))
		if plot:
			self.snr_plot(snr)
		#print('Time needed to reach $SNR= {}$'.format(snr)+' is ' + str(time[ind]))
		return time[ind]

	def snr_plot(self,snr):
		plt.figure()
		plt.plot(self.exp_time,self.snr)
		plt.axhline(snr,ls='--',color='r')
		plt.ylabel('SNR')
		plt.xlabel('Exposure time [s]')

	def snr_for_time(self,time,mag=None):
		if mag is not None:
			self.source = Source(mag, self.detector, self.telescope)
		self.exp_time = time 
		self.calculate_SNR()
		return self.snr




class Telescope():
	def __init__(self,name,diameter,throughput,f_num):
		self.name = name
		self.f_num = f_num
		self.diameter = diameter * u.m
		self.throughput = throughput
		self.area = np.pi*(self.diameter/2)**2
		self.platescale = self.calculate_plate_scale()

	def calculate_plate_scale(self):
		p = 206265 / self.f_num # 
		return p


class Sky():
	def __init__(self,moon_phase,detector,telescope):
		self.detector = detector
		self.telescope = telescope
		self.moon_phase = moon_phase
		self.sky_brightness = self._get_sky()
		self.photons = np.nan


	def _get_band_brightness(self):
		"""
		Table values taken from: 
		https://www.mso.anu.edu.au/~pfrancis/reference/reference/node4.html
		"""
		tab = pd.read_csv(package_directory + 'sky_brightness.csv')
		t = tab.iloc[self.detector.bandpass._tab_ind]
		return t

	def _get_sky(self):
		tab = self._get_band_brightness()
		sky = tab[self.moon_phase]
		return sky

	def sky_signal(self):
		sky = mag2flux(self.sky_brightness)
		# multiply by plate scale
		flux = sky * self.telescope.platescale**2
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
		flux = source.to(u.J/u.s/u.Hz/u.m**2)
		# convert to Lambda [J/s/m/m^2]
		flux_lam = flux * const.c / (self.detector.bandpass.wavelength.to(u.m)**2)
		# multiply by bandwidth [J/s/m^2]
		flux_lam = flux_lam * self.detector.bandpass.bandwidth.to(u.m)
		# multiply by size of pixel [J/s]
		flux_lam = flux_lam * self.telescope.area
		photons = flux_lam / ((const.h * const.c)/self.detector.bandpass.wavelength.to(u.m))
		self.photons = photons






class Detector():
	def __init__(self,name,dark_current,read_noise,gain,pixel_size,bandpass):
		self.name = name
		self.bandpass = bandpass
		self.dark_current = dark_current
		self.read_noise = read_noise
		self.qe = self._get_qe()
		self.gain = gain
		self.pix_size = pixel_size

		#self.check_inputs()

	
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

	def _get_qe(self):
		tab = self.bandpass._table.iloc[self.bandpass._tab_ind]
		qe = tab[self.name + '_qe']
		return qe



class Bandpass():
	def __init__(self,name):
		self.name = name
		self._table = self._load_table()
		self._tab_ind = self._get_index()

		self.bandwidth = self._get_bandwidth()
		self.wavelength = self._get_wavelength()

		

	def _load_table(self):
		tab = pd.read_csv(package_directory + 'bandpass_params.csv')
		return tab

	def _get_index(self):
		try:
			ind = np.where(self.name == self._table['name'].values)[0][0]
		except:
			raise ValueError('No such filter!')
		return ind

	def _get_bandwidth(self):
		
		bw = self._table.iloc[self._tab_ind]['bandwidth'] * u.nm
		return bw 

	def _get_wavelength(self):
		wav = self._table.iloc[self._tab_ind]['wavelength'] * u.nm
		return wav
