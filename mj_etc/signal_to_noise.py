import numpy as np


def mag(d,M):
    return 5*np.log10(d)-5+M

def Sky(mag):
    Flux=10**(-(mag+48.6)/2.5)
    Flux=Flux*10**-7*10**4 # convert to SI [J/s/Hz/m^2/arcsec^2]
    Flux=Flux*3*10**8/(Wavelength)**2 # Converting to wavelength [J/s/m/m^2/arcsec^2]
    Flux=Flux*Bandwidth # multiply by bandwidth [J/s/m^2/arcsec^2]
    Flux=Flux*Area # [J/s/arcsec^2]
    Flux=Flux*QE*Gain #Scalling by instrument [J/s/arcsec^2]
    Flux=Flux*Platescale**2 #[J/s/pix]
    photon=(Flux/((6.626*10**-34*3*10**8)/Wavelength)) #[count/pix/s]
    return photon
    
#H=Sky(Topguess)
#T=Sky(lowguess)
#K=Sky(Keck)
#S=Sky(Sidingspring)

#print(H,K,S)

def SN(mag,skymag,time,Wavelength,Bandwidth,Area,QE,Tau,Gain,Darkcurrent,Apareapix,Readnoise):
    #sky brightness
    Flux=10**(-(skymag+48.6)/2.5)
    Flux=Flux*10**-7*10**4 # convert to SI [J/s/Hz/m^2/arcsec^2]
    Flux=Flux*3*10**8/((Wavelength)**2) # Converting to wavelength [J/s/m/m^2/arcsec^2]
    Flux=Flux*Bandwidth # multiply by bandwidth [J/s/m^2/arcsec^2]
    Flux=Flux*Area # [J/s/arcsec^2]
    Flux=Flux*QE*Tau*Gain #Scalling by instrument [J/s/arcsec^2]
    Flux=Flux*Platescale**2 #[J/s/pix]
    sky=(Flux/((6.626*10**-34*3*10**8)/Wavelength))*time
    
    
    
    #Signal to noise 
    Flux=10**(-(mag+48.6)/2.5)
    Flux=Flux*10**-7*10**4 # convert to SI [J/s/Hz/m^2]
    Flux=Flux*3*10**8/(Wavelength)**2 # Converting to wavelength [J/s/m/m^2
    Flux=Flux*Bandwidth # multiply by bandwidth [J/s/m^2]
    Flux=Flux*Area # [J/s]
    Flux=Flux*QE*Tau*Gain #Scalling by instrument [J/s]
    Counts=(Flux/((6.626*10**-34*3*10**8)/Wavelength)) #[counts/s]
    Counts=Counts*time # [counts/expo/app] SIGNAL
    #errors 
    Dc=np.sqrt(Darkcurrent*time*Apareapix*Gain) # Dark current [count/exposure/app]
    Readnoise=nrd*np.sqrt(Apareapix) # e/app
    Source=np.sqrt(Counts)
    Skyerror=np.sqrt(sky*Apareapix) # [count/expo/app]
    Noise=np.sqrt(Dc**2+Readnoise**2+Source**2+Skyerror**2)
    
    return Counts, Noise


def DistSN(SN,Mag,skymag,exptime):
    #sky brightness
    Flux=10**(-(skymag+48.6)/2.5)
    Flux=Flux*10**-7*10**4 # convert to SI [J/s/Hz/m^2/arcsec^2]
    Flux=Flux*3*10**8/((Wavelength)**2) # Converting to wavelength [J/s/m/m^2/arcsec^2]
    Flux=Flux*Bandwidth # multiply by bandwidth [J/s/m^2/arcsec^2]
    Flux=Flux*Area # [J/s/arcsec^2]
    Flux=Flux*QE*Tau*Gain #Scalling by instrument [J/s/arcsec^2]
    Flux=Flux*Platescale**2 #[J/s/pix]
    sky=(Flux/((6.626*10**-34*3*10**8)/Wavelength))*exptime
    Skyerror=np.sqrt(sky*Apareapix) # [count/expo/app]
    #detector noise
    Dc=np.sqrt(Darkcurrent*exptime*Apareapix*Gain) # Dark current [count/exposure/app]
    Readnoise=nrd*np.sqrt(Apareapix) # e/app
    #calcualte the number of counts required for a provided S\N
    Counts=(-(SN**2)+np.sqrt(SN**4+4*SN**2*(Dc**2+Readnoise**2+Skyerror**2)))/2
    print(SN**4+4*(Dc**2+Readnoise**2+Skyerror**2))
    print(Counts)
    print(np.sqrt(SN**4+4*SN**2*(Dc**2+Readnoise**2+Skyerror**2)))
    Countspert=Counts/exptime
    Flux=Countspert*((6.626*10**-34*3*10**8)/Wavelength) #[J/s]
    Flux=Flux/(QE*Tau*Gain) # dividing by instrument efficiency to recover incident value
    Flux=Flux/Area #divide by area to fet flux units [J/s/m^2]
    Flux=Flux/Bandwidth #[J/s/m^2/m]
    Flux=Flux/(3*10**8/(Wavelength)**2) # [J/s/Hz/m^2]
    Fluxcgs=Flux/(10**-7*10**4)
    #Calcualte mag for given S\N
    SNmag=-5/2*np.log10(Fluxcgs)-48.6
    print(SNmag)
    #calcualte distance
    dist=10**(1/5*(SNmag-Mag+5)) #[pc]
    return dist    




# Changeable stuff
Sourcemag=14 # mag
Time=0.1 #s
Seeing=10 #arcsec/pixel
Wavelength=308*10**-9 #[m]
Bandwidth=20*10**-9

#########################
# Imaging specifications#
#########################
nrd=3 #e/pixel
Darkcurrent=0.2 #e/s/pixel at -40C
Gain=1 #e/ADU needs changing
FoV=np.sqrt(7*60**2) #arcmin, diameter
Pixelsize=13*10**-6 #m

PS=206256/0.6 #arcsec/m
Platescale= 5 # PS*Pixelsize      #4.64 #arcsec/pixel
QE=0.5


#######################
# Optical transmission#
#######################
lens1 = 0.95
lens2 = 0.95
mirror1 = 0.95
mirror2 = 0.95
lens3 = 0.95
Atmotransmission = 0.4
filter_efficiency = 0.25
Tau=lens1*lens2*lens3*mirror1*mirror2*Atmotransmission*filter_efficiency

# Aperature radius
Apradius=Seeing #arcsec FWHM
Apradiuspix=Apradius/Platescale # pixels
Aparea=np.pi*Apradius**2 #arcsec^2
Apareapix=np.pi*Apradiuspix**2 #pixel^2

################################
# Sky brightness [mag/arcsec^2]#
################################
Topguess = 25
lowguess = 24
Keck = 23.2
Sidingspring = 22.8


def Signal_noise(mag,time, Primary = 0.15, skymag = lowguess, Wavelength = 300e-9, bandwidth = 10e-9, 
					QE = 0.5, Atmo = 0.4, Filter = 0.25, Gain = 1, Darkcurrent = 0.2, 
					Apareapix = 1, Readnoise = 3):
	Secondary=0.45*Primary #m
	Area=4*np.pi*(Primary/2)**2 * 0.7

	lens1 = 0.95
	lens2 = 0.95
	mirror1 = 0.95
	mirror2 = 0.95
	lens3 = 0.95
	
	Tau=lens1*lens2*lens3*mirror1*mirror2*Atmo*Filter

	counts, noise = SN(mag,skymag,time,Wavelength,Bandwidth,Area,QE,Tau,Gain,Darkcurrent,Apareapix,Readnoise)
	return counts, noise
