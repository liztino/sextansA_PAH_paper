# Import smorgasbord
import os
import numpy as np
import scipy.signal
import scipy.interpolate
import scipy.integrate
import astropy.units as u
import specutils
import synphot
import stsynphot as stsyn



# Give data directory, and identify sub-directories we will loop over
data_dir = '/astro/dust_kg/cclark/JWST_M101_PAHs/Spectra/DL07/'
data_subdirs = os.listdir(data_dir)
data_subdirs = [data_subdir for data_subdir in data_subdirs if os.path.isdir(os.path.join(data_dir,data_subdir))]

# Give model types we want for each metallicity
gal_names = ['smc', 'LMC2', 'MW3.1'][::-1]
model_names = [['MW3.1_00','MW3.1_10','MW3.1_20','MW3.1_30','MW3.1_40','MW3.1_50','MW3.1_60'],
               ['LMC2_00','LMC2_05','LMC2_10'],
               ['smc']
               ]

# Define feature & clip wavelengths
features = [3.3, 7.7, 11.2]#, 12.8]
clip_1 = [3.09, 6.9, 10.8]#, 12.2]
clip_2 = [3.52, 9.7, 11.7]#, 13.3]

# Define feature filters and continuum filters
filts_dir = '/astro/dust_kg/cclark/JWST_M101_PAHs/Filters/'
filts_feature = ['NIRCam.F335M', 'MIRI.F770W', 'MIRI.F1130W', 'MIRI.F1280W']
filts_clip_1 = ['NIRCam.F300M', 'MIRI.F560W', 'MIRI.F1000W', 'MIRIF1000W']
filts_clip_2 = ['NIRCam.F360M', 'MIRI.F1000W', 'MIRI.F1500W', 'MIRI.F1500W']

# Provide misc other values
feature_snr_target = 10
snr_quantile = 97.5
jwst_area_sqcm = 331830.72404

# Create dictionary to store feature strengths for each model
out_dict = {}
for g in range(len(gal_names)):
    out_dict[gal_names[g]] = {'model':[],'snr':[]}
data_subdirs = [os.path.join(data_dir,data_subdir) for data_subdir in data_subdirs if os.path.isdir(os.path.join(data_dir,data_subdir))]





# Loop over galaxy model types, then Umin subdirs, then models contained therein
for gal_name in gal_names:
    print('Processing '+gal_name+':')
    for data_subdir in data_subdirs:
        print('Testing '+data_subdir.split('/')[-1]+' models')
        model_names = [model_name for model_name in os.listdir(data_subdir) if ((gal_name in model_name) and ('._' not in model_name))]
        for model_name in model_names:

            # E the very different models where the dust is inside ionized gas in a supernova remnant in an HII region in 30 Dor...
            if 'LMC2_00' in model_name:
                continue

            # Check this model has data; if not, skip to next
            with open(os.path.join(data_subdir, model_name)) as model_file:
                model_lines = len(model_file.readlines())
            if model_lines <= 2:
                continue

            # If model has data, check where actual SED starts
            with open(os.path.join(data_subdir, model_name)) as model_file:
                model_file_count = 0
                for model_line in model_file:
                    model_file_count += 1
                    if model_line.rstrip() == 'lambda    nu*dP/dnu     j_nu':
                        model_file_header = model_file_count + 1
                        break

            # Read in model SED
            data_in = np.genfromtxt(os.path.join(data_subdir, model_name), skip_header=model_file_header)
            data_micron = data_in[:,0]

            # Change units into /kg as oppoed to /H, in order to make numbers work with ETC
            data_in[:,1] /= 1.67E-27
            data_in[:,2] /= 1.67E-27

            # Produce clipped version of wavelength range
            lambda_clip = data_in[:,0].copy()
            for i in range(len(features)):
                where_clip = np.where((data_in[:,0]>clip_1[i]) & (data_in[:,0]<clip_2[i]))
                lambda_clip[where_clip] = np.nan

            # Interpolate over clipped wavelength ranges in log space
            where_not_clip = np.where(np.isnan(lambda_clip)==False)
            interp = scipy.interpolate.interp1d(np.log10(lambda_clip[where_not_clip]),
                                                np.log10(data_in[:,1][where_not_clip]))
            data_cont = 10**interp(np.log10(data_in[:,0]))

            # Subtract clipped spectrum, to get continuum-subtracted line strength
            data_sub = data_in[:,1] - data_cont
            data_sub[np.where(np.isnan(lambda_clip) == False)] = 1.0
            data_sub[np.where(data_sub < 1)] = 1.0




            # Create specutils spectra for each of our models
            spec_in = specutils.Spectrum1D(flux=data_in[:,1]*u.Jy, spectral_axis=data_in[:,0]*u.micron)
            spec_cont = specutils.Spectrum1D(flux=data_cont*u.Jy, spectral_axis=data_in[:,0]*u.micron)
            spec_sub = specutils.Spectrum1D(flux=data_sub*u.Jy, spectral_axis=data_in[:,0]*u.micron)

            # Convert specutils spectra to synphot spectra
            spec_in = synphot.SourceSpectrum.from_spectrum1d(spec_in)
            spec_cont = synphot.SourceSpectrum.from_spectrum1d(spec_cont)
            spec_sub = synphot.SourceSpectrum.from_spectrum1d(spec_sub)

            # Loop over features, to do synthetic photometry
            feature_snr_list = []
            for i in range(len(features)):

                # Perform synphot observations on subtracted and unsubtracted spectra
                filt_feature = synphot.spectrum.SpectralElement.from_file(os.path.join(filts_dir,'JWST_'+filts_feature[i]+'.dat'))
                obs_feature = synphot.Observation(spec_in, filt_feature)
                obs_feature_cont = synphot.Observation(spec_cont, filt_feature)
                obs_feature_flux = obs_feature.countrate(jwst_area_sqcm)
                obs_feature_cont_flux = obs_feature_cont.countrate(jwst_area_sqcm)
                obs_feature_sub_flux = obs_feature_flux - obs_feature_cont_flux

                # Work out relative strength of feature
                obs_feature_rel_filt = 100 * (obs_feature_sub_flux / obs_feature_flux).value

                # Work out S/N required for full spectrum to get S/N>10
                feature_snr = abs(feature_snr_target / (obs_feature_sub_flux / obs_feature_flux).value)
                feature_snr_list.append(feature_snr)

            # Record to results dictionary
            out_dict[gal_name]['model'].append(model_name)
            out_dict[gal_name]['snr'].append(np.array(feature_snr_list))



    # Work out various limit quantile required S/N for each feature
    print(' ')
    snr_array = np.array(out_dict[gal_name]['snr'])
    for i in range(len(features)):
        snr_array_feature = np.array(snr_array[:,i])#SigmaClip(np.array(snr_array[:,i], median=True, sigma_thresh=2))
        snr_quantile_diff = np.abs(np.nanpercentile(snr_array_feature, snr_quantile) - snr_array[:,i])
        snr_quantile_index = np.nanargmin(snr_quantile_diff)
        snr_quantile_snr = snr_array[snr_quantile_index,i]
        snr_quantile_model = out_dict[gal_name]['model'][snr_quantile_index]
        print(gal_name+' '+str(features[i])+'um target S/N='+str(snr_quantile_snr)[:-6]+' with model '+snr_quantile_model)
    print(' ')
    print(' ')
    print(' ')



# Jubilate
print('And furthermore, Carthage must be destroyed!')