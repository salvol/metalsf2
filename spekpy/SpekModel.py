# For compatibility with Python2 #
from __future__ import print_function, division, absolute_import
##################################
import numpy as np
import spekpy.SpekConstants as Const
from scipy import integrate

## References (Note: Ref. 1-3 describe "legacy" model i.e. SpekCalc)
#[1] Poludniowski G, Evans PM. Calculation of x-ray spectra emerging from an x-ray tube. Part I. electron penetration characteristics in x-ray targets. Med Phys. 2007;34(6):2164-74.
#[2] Poludniowski G. Calculation of x-ray spectra emerging from an x-ray tube. Part II. X-ray production and filtration in x-ray targets. Med Phys. 2007;34(6):2175-86.
#[3] Poludniowski G, Landry G, DeBlois F, Evans PM, Verhaegen F. SpekCalc: a program to calculate photon spectra from tungsten anode x-ray tubes. Phys Med Biol. 2009;54(19):N433-8.
#[4] Seltzer SM, Berger MJ. Bremsstrahlung spectra from electron interactions with screened atomic nuclei and orbital electrons. NIM B 1985;12(1):95-134.

class SpekModel:

    def __init__(self):

        self.norm = 1.0 # Normalization factor to match reference fluence or air kerma

        # Electron attributes
        self.nt = None # Number of depths calculated for electrons in target
        self.nu = None # Max number of electron energies at depth in a target
        self.t = None  # Array with depths in target
        self.dt = None  # Depth increment for depth array
        self.csda_range = None  # cm
     
        # Photon attributes  
        self.number_of_photon_energy_bins = None

        # Self filtration attributes
        self.anode_self_filtration = None

        # Attributes with loaded data tables (functions)
        self.__nist_brem_data = None  # NIST differential bremsstrahlung cross-section data (Ref. [4])
        self.__sig_br = None  # Differential bremstrahlung cross-sections used in model.
        self._pe = None  # Conditional probability function for electrons.
        self._ne = None # Electron frequency at depth
        self._line_data = None # Tabulated data for K and L-lines
       
        # Spectrum arrays
        self.brem_kt = None #2d array (photon energy, depth) for bremsstrahlung emissions
        self.char_kt = None # 2d array (photon energy, depth) for characteristic emissions
        self.brem_k = None # 1d array (photon energy) for bremsstrahlung emissions
        self.char_k = None # 1d array (photon energy) for characteristic emissions
        self.k = None  # 1d array with photon energies

    def get_spectrum_parameters(self, spekpy_obj):

        """
        An internal method to get all of the parameters needed for the spectrum model

        :param Spek spekpy_obj: A spekpy state
        :return Spekmodel self: A Spekmodel object that has been populated with the parameters needed for the model
        """

        self.nt = 100  # Number of depths calculated for electrons in target
        self.nu = 200  # Max number of electron energies at depth in a target       

        E0 = spekpy_obj.state.model_parameters.kvp  # Kinetic energy of incident electrons [keV]
        dk = spekpy_obj.state.model_parameters.dk  # Bin width of photon energies [keV]
        physics = spekpy_obj.state.model_parameters.physics  # String indicating if legacy mode is activated
        mu_data = spekpy_obj.mu_data  # Linear attenuation data

        ntar = Const.avogadros_number * Const.density_tungsten / Const.atomic_weight_tungsten  # Number of Tungsten atoms per cm^3 in target
        constant_factor = ntar * Const.geometry_factor * Const.detour_factor  # Constant normalization factor

        # Load functions for model from numerical tables
        #electron_data = spekpy_obj.electron_data#ElectronData()
        self.__electron_range_csda = spekpy_obj._rg #electron_data.get_rg()
        self.__electron_range_tw = lambda E0: 0.0119 * (E0 ** 1.513) / Const.conversion_g2mg # Thompson-Whidington range g/cm^2 (equation 23 in Ref. [1])
        csda_range = self.__electron_range_csda(E0) / Const.density_tungsten  # [cm]
        tw_range = self.__electron_range_tw(E0) / Const.density_tungsten 
        self._pe = spekpy_obj._pe
        self._ne = spekpy_obj._ne
        self._line_data = spekpy_obj._line_data

        if physics == 'legacy':
            model_param = Const.model_param_legacy
        elif physics == 'default':
            model_param = Const.model_param_default
        else:
            raise Exception('physics mode incorrectly specified!')

        brem_normalization = model_param['nbr']
        L1_line_production = model_param['nL1']
        L2_line_production = model_param['nL2']
        L3_line_production = model_param['nL3']
        K_line_production = model_param['nK']

        if model_param['depth'] == 'csda' and type(model_param['scale']) == type(None):
            E0_ref = E0
            scale = 1.0
            max_emission_depth = csda_range * 0.5
        elif model_param['depth'] == 'csda' and type(model_param['scale']) != type(None):
            E0_ref =  model_param['scale']
            scale = 1.0
            max_emission_depth = csda_range * 0.5
        elif model_param['depth'] == 'tw' and type(model_param['scale']) == type(None):
            E0_ref = E0
            scale = 1.0 
            max_emission_depth = tw_range
        elif model_param['depth'] == 'tw' and type(model_param['scale']) != type(None):
            E0_ref =  model_param['scale']
            scale = (self.__electron_range_tw(E0) / self.__electron_range_tw(E0_ref)) \
                / (self.__electron_range_csda(E0) / self.__electron_range_csda(E0_ref)) 
            max_emission_depth = tw_range
        else:
            raise Exception('depth or scale option incorrectly specified!')

        if model_param['brxs'] == 'nist':
            self.__nist_brem_data = spekpy_obj._nist_brem_data  # Interpolate bremstrahlung cross-sections from NIST
            self.__sig_br = self._sig_br_nist
        elif model_param['brxs'] == 'mebh':
            self.__sig_br = self._sig_br_mebh  # Modified Elwert-Bethe-Heitler bremsstrahlung cross-section (as in SpekCalc)
        else:
            raise Exception('brxs option incorrectly specified!')

        # Calculate depth array and step size
        [self.t, self.dt] = np.linspace(0.0, max_emission_depth, self.nt, retstep=True)
        t_scaled=self.t / (csda_range * scale)

        # Calculate the set of photon energies to calculate based on specified bin width
        self.number_of_photon_energy_bins = int(((E0 - dk * 0.5) - 1.0) / dk) + 1

        self.k = np.linspace(E0 - dk * (self.number_of_photon_energy_bins - 0.5),
                             E0 - dk * 0.5, self.number_of_photon_energy_bins)

        # Calculate spectral contributions 
        # "kt" indicates that these are 2D arrays corresponding to emission energy (k) and emission depth (t)
        self.brem_kt = brem_normalization * constant_factor * self.__brem_kt(E0, E0_ref, t_scaled)

        L1_char_kt = brem_normalization * constant_factor * L1_line_production * self.__char_kt(self._line_data['L_1'], E0, dk, E0_ref, t_scaled)
        L2_char_kt = brem_normalization * constant_factor * L2_line_production * self.__char_kt(self._line_data['L_2'], E0, dk, E0_ref, t_scaled)
        L3_char_kt = brem_normalization * constant_factor * L3_line_production * self.__char_kt(self._line_data['L_3'], E0, dk, E0_ref, t_scaled)
        K_char_kt = brem_normalization * constant_factor * K_line_production * self.__char_kt(self._line_data['K'], E0, dk, E0_ref, t_scaled)
        
        self.char_kt = L1_char_kt + L2_char_kt + L3_char_kt + K_char_kt # Add L and K-line contributions

        self.anode_self_filtration = mu_data.get_mu_over_rho(Const.atomic_number_tungsten, self.k)[:, None] * \
            Const.density_tungsten * self.t[None, :]

        return self

    def __brem_kt(self, E0, E0_ref, t_scaled):
        """
        An internal method to get a 2D array corresponding to bremstrahlung emission at emission energy (k) and emission depth (t)

        :param float E0: Kinetic energy of incident electrons
        :param flocat scale: Factor used when scaling estimations of photon emission
        :param float bfac: electron backscatter factor
        :param float E0_ref: Reference kinetic energy of incident electrons
        :return array brem: 2D array of bremsstrahlung emission with respect to emission energy and emission depth
        """
        brem = np.asarray([self.__brem_t(k, E0, E0_ref, t_scaled) for k in self.k])
        return brem

    def __brem_t(self, k, E0, E0_ref, t_scaled):
        """Internal method to calculate bremstrahlung emissions corresponding to a 1D array with respect to emission depth. 

        :param float k: Emission energy [keV]
        :param float E0: Kinetic energy of incident electrons [keV]
        :param float scale: Factor used when scaling estimations of photon emission 
        :param float E0_ref: Reference kinetic energy of incident electrons
        :return array brem_t: 1D array of bremsstrahlung emissions with respect to emission depth
        """
        n = max(3, int((1.0 - k / E0) * self.nu)) # Number of electron energies (bremsstrahlung emission energy and higher)
        [u, du] = np.linspace(k / E0, 1.0, n, retstep=True)
        Ei = u * E0
        sig_br = self.__sig_br(Ei, k)
        pe = self._pe(E0_ref, t_scaled, u)
        ne = self._ne(E0_ref, t_scaled)
        brem_t = integrate.simps(sig_br[None, :] * pe * ne[:, None], axis=1, dx=du) / self.__norm_t(E0, E0_ref, t_scaled)
        return brem_t

    def __char_kt(self, char_data, E0, dk, E0_ref, t_scaled):
        """
        An internal method to get a 2D array corresponding to Characteristic emission at emission energy (k) and emission depth (t)

        :param float E0: Kinetic energy of incident electrons
        :param flocat scale: Factor used when scaling estimations of photon emission
        :param float bfac: Electron backscatter factor
        :param float E0_ref: Reference kinetic energy of incident electrons
        :return array char_kt: 2D array of bremsstrahlung emission with respect to characteristic emission energy and emission depth
        """

        char_kt = np.zeros([self.number_of_photon_energy_bins, self.nt])
        if self.k[-1] > char_data['Edge']:
          c = self.__char(char_data['Edge'], E0, dk, E0_ref, t_scaled) / dk
          iind = np.floor((char_data['Lines'] - (self.k[0] - dk * 0.5)) / dk).astype(int) # Indices of energy array corresponding to line energies 
          nind = len(iind)
          for i in range(nind):
            char_kt[iind[i], :] = char_kt[iind[i], :] + c * char_data['W'][i] # Puts weighted fraction for each line in correct energy bin
        return char_kt

    def __char(self, edge_energy, E0, dk, E0_ref, t_scaled):
        """
        A method to calculate a 1D array corresponding to characteristic emission with respect to energy.

        :param float edge_energy: Energy of L or K edge
        :param float E0: Kinetic energy of incident electrons
        :param float dk: Width of energy bin 
        :param float scale: Factor used when scaling estimation of photon emission
        :param float bfac: Electron backscatter factor
        :param float E0_ref: Reference kinetic energy of incident electrons. 
        :return array char: 1D array consisting of total bremsstrahlung emitted above edge (differential in depth), at different depths
                            The assumption is that a fixed fraction of these are absorbed and produce characteristic emissions
        """
        nk = max(2, int((E0 - edge_energy) / dk))  # Number of bremsstrahlung energies (from edge energy up)
        [k, dk] = np.linspace(edge_energy, E0, nk, retstep=True)
        brem = [self.__brem_t(k[i], E0, E0_ref, t_scaled) for i in range(nk)]
        char = integrate.simps(brem, axis=0, dx=dk)
        return char

    def __norm_t(self, E0, E0_ref, t_scaled):
        """
        Returns the integral of pe for each interpolated depth
        The pe should integrate to one and dividing by this function compensates for numerical errors

        :param float E0: Kinetic energy of incident electrons
        :param float scale: Factor used when scaling estimation of photon emission
        :param float E0_ref: Reference kinetic energy of incident electrons. 
        :return array norm: 1D array consisting of normalization factors at different depths
        """
        n = self.nu
        [u, du] = np.linspace(0.0, 1.0, n, retstep=True)
        pe = self._pe(E0_ref, t_scaled, u)
        norm = integrate.simps(pe, axis=1, dx=du)
        return norm

    def _sig_br_nist(self,Ei, k):
        """
        Internal method to calculate bremsstrahlung differential x-section using
        interpolation of NIST data (see Ref. [4])
        
        :param array Ei: Electron kinetic energies prior to emission [keV]
        :param array k: Energy of photon emission [keV]
        :return array sig: Array of differential bremstrahlung cross-sections
        """
        with np.errstate(invalid='ignore', divide='ignore'):
            Ti = Ei + Const.electron_mass  # Initial total energy of electron
            Ef = Ei - k  # Final kinetic energy of electron
            gi = Ti / Const.electron_mass
            bi = np.sqrt(1.0 - gi ** -2)
            sigval = [1.0e-27 * (Const.atomic_number_tungsten ** 2 / (k * bi[i] ** 2)) *
                      self.__nist_brem_data(Ei[i], k / Ei[i])[0][0] for i in range(len(Ei))]
            sig = np.array(sigval)
            sig[Ef < 0.0] = 0.0
        return sig

    def _sig_br_mebh(self, Ei, k):
        """
        Static method to calculate the modified semi-relativistic Elwert-Bethe-Heitler (srMEBH) bremsstrahlung 
        differential cross-sections (see Ref. [2]). This method was used in the SpekCalc algorithm.

        :param array Ei: Array of electron kinetic energies prior to emission [keV].
        :param array k: Energy of photon emission [keV].
        :return array sig: Array of differential bremsstrahlung cross-sections.
        """

        with np.errstate(invalid='ignore', divide='ignore'):
            Ti = Ei + Const.electron_mass  # Initial energy of electron
            Ef = Ei - k  # Final kinetic energy of electron
            Tf = Ef + Const.electron_mass  # Final energy of electron equivalent to
            qi = np.sqrt(Ti ** 2 - Const.electron_mass ** 2)  # Equivalent to p_i * c
            qf = np.sqrt(Tf ** 2 - Const.electron_mass ** 2)  # Equivalent to p_f * c
            elwert_factor = qi / qf # This is the modified Elwert factor described in Ref. [2]
            l = 2.0 * np.log((Ti * Tf + qi * qf - Const.electron_mass ** 2.0) / (Const.electron_mass * k))
            phival = (2.0 / 3.0) * (Ei / (k * qi ** 2.0)) * (4.0 * Ti * Tf * l - 7.0 * qi * qf) * elwert_factor
            phi = Const.phi0 * np.where(Ef < 0.0, 0.0,
                                        np.where(Ef < 1.0, (2. / 3.) * (1. + 8. * Const.electron_mass / k), phival))
            sig= phi / Ei # Converting from Heitler's notation to conventional
        return sig



