#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import fileinput as fileinp
import sys as sys
import constants as const

########################################
##           input_reader.py
########################################

def read_input(parameters,filepath):
    #for line in fileinp.input(sys.argv[1:]):
    for line in fileinp.input(filepath):
        if(line[0]!='#' and line[0]!='\n'):
            key, value = line.split("=")
            parameters[key.strip()] = value.strip()

    ########################################


class parameters:
    def __init__(self,inparams):

        self.lmax = int(inparams['lmax'])
        self.mmax = int(inparams['mmax'])

        self.maxthetapts = int(inparams['maxthetapts'])
        self.maxphipts   = int(inparams['maxphipts'])
        
        self.dEe = float(inparams['dEe'])
        self.Eemaxpts = float(inparams['Eemaxpts'])

        self.dEn = float(inparams['dEn'])
        self.Enmaxpts = float(inparams['Enmaxpts'])

        self.dE  = float(inparams['dE'])
        self.Emax = float(inparams['Emax'])
        
        if(int(inparams['fixed_nuclei']) == 1):
            self.fixed_nuclei = True
        else:
            self.fixed_nuclei = False
            
        # Mass factors
        if(self.fixed_nuclei):
            self.muelec = 1.0
        else:
            self.muelec = ( const.M1 + const.M2 ) / \
                  (const.M1 + const.M2 + 1.0)

        self.munuc = const.M1 * const.M2 / (const.M1 + const.M2)
        self.mufac = 0.5 / self.muelec
        self.Mfac = 0.5 / self.munuc

         # Let's convert the frequency to atomic units
        self.intensity = float(inparams['intensity'])
        self.wavelength = float(inparams['wavelength'])
        self.wavelength0 = self.wavelength / const.aulength_nm
        self.w0 = 2.0 * np.pi * const.speed_light_au / self.wavelength0
        self.period = 2.0 * np.pi / self.w0
        self.no_cycles = float(inparams['no_cycles'])
        self.pulse_duration = self.period * self.no_cycles
        self.pulse_bandwidth = 4.0 * np.pi / self.pulse_duration
        self.quiver = np.sqrt(self.intensity / const.intensity_au) / self.w0 / self.w0
        self.Up = (self.intensity / const.intensity_au) / 4. / self.w0 / self.w0
        self.E0 = np.sqrt(self.intensity / const.intensity_au)
        self.A0 = self.E0 / self.w0 / const.speed_light_au
        self.Ee_ejec_au = float(inparams['Ee_ejec_eV']) / const.energy_au_ev # in au
        self.ke_ejec = np.sqrt(2.0 * self.muelec * self.Ee_ejec_au)
        self.Eg = float(inparams['Eg'])
        self.Ip = self.Eg - 0.5

        self.draw_spherical_amplitude = int(inparams['draw_spherical_amplitude'])
        self.draw_polar_amplitude = int(inparams['draw_polar_amplitude'])
        self.draw_sph_polar_amplitude = int(inparams['draw_sph_polar_amplitude'])
        self.draw_pad = int(inparams['draw_pad'])
        self.draw_polar_pad = int(inparams['draw_polar_pad'])
        self.draw_semipolar_pad = int(inparams['draw_semipolar_pad'])
        self.draw_mes = int(inparams['draw_mes'])
        self.draw_pes = int(inparams['draw_pes'])
        # self.draw_Etotal = int(inparams['draw_Etotal'])
        self.draw_diff_cross = int(inparams['draw_diff_cross'])
        self.draw_total_cross = int(inparams['draw_total_cross'])
        self.draw_sampling_pes = int(inparams['draw_sampling_pes'])
        self.draw_sampling_pad = int(inparams['draw_sampling_pad'])
        self.draw_semipolar_samp_pad = int(inparams['draw_semipolar_samp_pad'])
        self.draw_wavetime = int(inparams['draw_wavetime'])

        self.makeframe = int(inparams['makeframe'])
        self.showframe = int(inparams['showframe'])        
        self.polar_filename = inparams['polar_filename']
        self.spherical_filename = inparams['spherical_filename']        
        self.mes_filename = inparams['mes_filename']
        self.pes_filename = inparams['pes_filename']
        self.samp_pes_filename = inparams['samp_pes_filename']
        self.samp_pad_filename = inparams['samp_pad_filename']
        self.samp_wavetime_filename = inparams['samp_wavetime_filename']
 
