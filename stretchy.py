#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import time
from matplotlib import pyplot as plt
import stretchy_plot_photon_weights as pspw


CHAND_SLAB_SOL = np.array(
    [[00.00, 1.26938], [18.19, 1.22945], [25.84, 1.18947], [31.79, 1.14943],
     [36.87, 1.10931], [41.41, 1.06911], [45.57, 1.02882], [49.46, 0.98842],
     [53.13, 0.94789], [56.63, 0.90722], [60.00, 0.86637], [63.26, 0.82530],
     [66.42, 0.78938], [69.51, 0.74234], [72.54, 0.70029], [75.52, 0.65770],
     [78.46, 0.61439], [81.37, 0.57001], [84.26, 0.52397], [87.13, 0.47490],
     [90.00, 0.41441]])


# Constants which control the general workings of the simulation, you should not
# need to change these variables unless you want to change the seed.

NEUTRAL_WEIGHT = 1
MIN_WEIGHT = 1e-5 * NEUTRAL_WEIGHT
MAX_WEIGHT = 1e5 * NEUTRAL_WEIGHT
LMAX = 1
randseed = np.random.randint(int(1e9))
SEED = 349903800
print("Random seed in use:", SEED)
N_BINS = CHAND_SLAB_SOL.shape[0]
VERBOSE = False
WRITE_PHOT_INFO = False  # This will probably make things slowwwww

# Parameters which control photon transport and path stretching, changing these
# will alter the performance of the script and algorithms

SCAT_ALBEDO = 1.0      # 1 == pure scattering, 0 == pure absorption
N_PHOTONS = int(1e6)
P_KILL = 0.0
TAU_MAX = 5
ALPHA_FLOOR = 0

ALPHA = None
if len(sys.argv) == 2:
    try:
        ALPHA = float(sys.argv[1])
    except ValueError:
        print("Uh oh, something went wrong trying to convert input to a float: ", sys.argv[1])
        exit(1)

print("")
print("PARAMETERS:")
print("ALPHA", ALPHA)
print("ALPHA_FLOOR", ALPHA_FLOOR)
print("SCAT_ALBEDO", SCAT_ALBEDO)
print("N_PHOTONS", N_PHOTONS)
print("P_KILL", P_KILL)
print("TAU_MAX", TAU_MAX)
print("")


class IntensitySpectrum:
    """
    Contains all the guff required to create the intensity spectrum of intensity
    vs angle.
    """

    def __init__(self, nbins):
        """
        Initialise the spectrum structures.

        nbins   int
                the number of bins for the spectrum
        """

        self.nphot = 0
        self.nbins = nbins
        self.weight = np.zeros(nbins)
        self.theta = np.zeros(nbins)
        self.intensity = np.zeros(nbins)

        return

    def init_bins(self):
        """
        Initialise the theta bins for the spectrum.
        """

        dtheta = 1 / self.nbins
        hwidth = 0.5 * dtheta
        for i in range(self.nbins):
            self.theta[i] = np.arccos(i * dtheta + hwidth)

        return

    def bin_photon_direction(self, pcostheta, pweight):
        """
        Bin a photon weight into the appropriate theta bin.

        pcostheta    float
                     the escape direction of a photon packet
        pweight      float
                     the weight of the photon
        """

        if pweight == 0:
            return

        ibin = np.abs(int(pcostheta * self.nbins))
        self.weight[ibin] += pweight
        self.nphot += 1

        return

    def calculate_intensity(self):
        """
        Calculate the intensity from the hist bins.
        """

        self.intensity = (self.weight * self.nbins) / (2 * N_PHOTONS * np.cos(self.theta))

        return

    def write(self, name):
        """
        Write the spectrum to file.

        name        str
                    the base name of the output file
        """

        fname = name + ".txt"
        f = open(fname, "w")

        f.write("theta\tcostheta\tw\tintensity\n")

        for i in range(self.nbins):
            f.write("{}\t{}\t{}\t{}\n".format(self.theta[i], np.cos(np.deg2rad(self.theta[i])),
                    self.weight[i], self.intensity[i]))

        f.close()

        return


class PhotonPacket:
    """
    Contains all of the functions required to generate and transport a photon.
    """

    def __init__(self, nphot):
        """
        Initialise a photon. Photons are initialised at the origin of the slab,
        with a direction so that they are pointing upwards.
        """

        self.x = 0
        self.y = 0
        self.z = 0
        self.inslab = True
        self.nphot = nphot
        self.weight = NEUTRAL_WEIGHT
        self.costheta = np.sqrt(np.random.rand())
        self.sintheta = np.sqrt(1 - self.costheta ** 2)
        self.cosphi = np.cos(2 * np.pi * np.random.rand())
        self.sinphi = np.sqrt(1 - self.cosphi ** 2)
        self.nscats = 0                # number of scats undergone
        self.nrr = 0                   # number of rr games
        self.tau_path = 0              # optical depth to escape
        self.nstretch = 0              # number of stretched paths
        self.tau_scat = 0              # optical depth till next scatter
        self.tau_depth = TAU_MAX       # optical depth within the slab

        self.weight_history = []       # the history of the photon's weight
        self.tau_scat_history = []     # the history of the sampled tau_scats

        return

    def play_russian_roulette(self):
        """
        Play a game of Russian Roulette with the photon.
        """

        self.nrr += 1
        if np.random.rand() <= P_KILL:
            self.weight = -1
            return
        self.weight *= 1 / (1 - P_KILL)

        return self.weight

    def reflect_bottom_slab(self):
        """
        Reflect the photon if it hits the bottom of the slab.
        """

        tmp = np.deg2rad(90 - np.arccos(self.costheta))
        self.x = self.z * np.tan(tmp)
        self.z = 0
        self.costheta = np.sqrt(np.random.rand())
        self.sintheta = np.sqrt(1 - self.costheta ** 2)
        self.cosphi = np.cos(2 * np.pi * np.random.rand())
        self.sinphi = np.sqrt(1 - self.cosphi ** 2)

        return

    def find_tau_path(self):
        """
        Find the maximum optical depth to escape
        """

        self.tau_path = 0

        # This should stop photons from going underneath the slab
        if self.tau_depth > TAU_MAX:
            self.reflect_bottom_slab()

        if self.costheta > 0:
            self.tau_path = self.tau_depth / self.costheta
        elif self.costheta < 0:
            self.tau_path = (self.tau_depth - TAU_MAX) / self.costheta

        return self.tau_path

    def reweight_photon(self, alpha):
        """
        Reweight a photon after it has taken a stretched path

        alpha   float
                the stretching parameter, has a value 0 <= alpha <= 1
        """

        self.weight *= np.exp(self.tau_scat * (alpha - 1)) / alpha

        return self.weight

    def bias_tau_scat(self, alpha):
        """
        Sample from q(tau) = alpha * e ** (- alpha * tau), where alpha is the
        stretching parameter of the distribution between 0 and 1

        alpha   float
                the stretching parameter, has a value 0 <= alpha <= 1
        """

        if 0 > alpha > 1:
            print("error: expected 0 <= alpha <= 1, but got alpha = {}".format(alpha))
            print("       setting tau_scat form normal distribution")
            self.tau_scat = self.unbiased_tau_scat()
            return self.tau_scat
        self.nstretch += 1
        self.tau_scat = -1.0 * (np.log(1 - np.random.rand()) / alpha)
        self.reweight_photon(alpha)

        return self.tau_scat

    def unbiased_tau_scat(self):
        """
        Sample from p(tau) = e ** -tau.
        """

        self.tau_scat = -1.0 * np.log(1 - np.random.rand())

        return self.tau_scat

    def scatter(self):
        """
        Isotropic scattering of a photon.
        """

        self.nscats += 1
        self.costheta = 2 * np.random.rand() - 1
        self.sintheta = np.sqrt(1 - self.costheta ** 2)
        self.cosphi = np.cos(2 * np.pi * np.random.rand())
        self.sinphi = np.sqrt(1 - self.cosphi ** 2)

        return

    def translate(self, ds):
        """
        Update the position of a photon. This will also update the tau_depth
        variable which indicates how deep the photon is within the slab.

        ds  float
            distance the photon is being moved
        """

        self.x += ds * self.sintheta * self.cosphi
        self.y += ds * self.sintheta * self.sinphi
        self.z += ds * self.costheta
        self.tau_depth -= self.costheta * self.tau_scat

        return


def append_to_phot_file(root, phot):
    """
    Write certain photon properties to file for later analysis.

    Parameters
    ----------
    root                str
                        The root name for the output file
    phot                Photon
                        the photon to write to file
    """

    fname = root + "phot_info.txt"
    with open(fname, "a") as f:
        if len(phot.weight_history) != len(phot.tau_scat_history):
            print("phot {} weight history and tau scat history different length somehow".format(phot.nphot))
            return
        for i in range(len(phot.weight_history)):
            f.write("{} {} {}\n".format(phot.nphot, i, phot.weight_history[i], phot.tau_scat_history[i]))

    return


def trans_phot_single(phot, path_stretch_on, rr_on):
    """
    Transport a single photon packet.

    phot                Photon
                        the photon to transport
    path_stretch_on     bool
                        if True, path stretching will be used
    rr_on               bool
                        if True, Russian Roulette will be used
    """

    phot.weight_history.append(phot.weight)
    phot.tau_scat_history.append(phot.tau_scat)

    while phot.inslab and phot.weight > 0:
        if path_stretch_on:
            if ALPHA:
                alpha = ALPHA
            else:
                tau_path = phot.find_tau_path()
                alpha = 1 / (1 + tau_path)
            tau_scat = phot.bias_tau_scat(alpha)
            if rr_on and phot.weight <= MIN_WEIGHT:
                phot.play_russian_roulette()
                if phot.weight == -1:
                    break
        else:
            tau_scat = phot.unbiased_tau_scat()

        phot.weight_history.append(phot.weight)
        phot.tau_scat_history.append(phot.tau_scat)

        ds = tau_scat / TAU_MAX
        phot.translate(ds)

        # The photon has gone below the slab, reflect photon off the bottom
        if phot.z < 0:
            phot.reflect_bottom_slab()
        # The photon has escaped from the top of the slab
        elif phot.z > LMAX:
            phot.inslab = False
        # Otherwise, the photon is still within the slab and may scatter
        else:
            # The photon scatters
            if np.random.rand() <= SCAT_ALBEDO:
                phot.scatter()
            # The photon is completely absorbed
            else:
                phot.weight = 0

    return phot


def trans_phot(root, path_stretch_on):
    """
    The main loop which will transport a photon through the slab

    root                str
                        the root name for any file output
    path_stretch_on     bool
                        flag for enabling path stretching
    """

    if SEED:
        np.random.seed(SEED)

    # Set if Russian Roulette is turned off or on
    rr_on = False
    if P_KILL > 0:
        rr_on = True

    f = open(root + "phot_info.txt", "w")
    f.close()

    spectrum = IntensitySpectrum(N_BINS)
    spectrum.init_bins()

    for iphot in range(N_PHOTONS):
        phot = PhotonPacket(iphot)
        phot = trans_phot_single(phot, path_stretch_on, rr_on)
        if phot.inslab is False and phot.weight > 0:
            spectrum.bin_photon_direction(phot.costheta, phot.weight)
        append_to_phot_file(root, phot)
        if (iphot + 1) % (N_PHOTONS / 10) == 0:
            print("{:3.0f}% photons transported".format(iphot / N_PHOTONS * 100.0))

    sname = "intens_" + root
    if path_stretch_on:
        sname += "path_str_on"
    spectrum.calculate_intensity()
    spectrum.write(sname)

    return spectrum


def main():
    """
    Main control function.
    """

    if ALPHA is not None:
        root = "tmax_{}_pkill_{}_alpha_{}_".format(TAU_MAX, P_KILL, ALPHA)
    else:
        root = "tmax_{}_pkill_{}_alpha_tpath_".format(TAU_MAX, P_KILL, ALPHA)

    print("Running MCRT with path stretching enabled\n")
    start = time.time()
    pasMCRT = trans_phot(root, path_stretch_on=True)
    end = time.time()
    print("\nRun time: {} seconds".format(int(end - start)))

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    ax.semilogy(np.cos(np.deg2rad(CHAND_SLAB_SOL[:, 0])), CHAND_SLAB_SOL[:, 1], label="Solution")
    ax.semilogy(np.cos(pasMCRT.theta), pasMCRT.intensity, label="Path Stretched MCRT")
    ax.set_xlabel(r"$\mu = \cos(\theta)$")
    ax.set_ylabel(r"Normalised Intensity, I")
    ax.set_xlim(0, 1)
    ax.legend()
    plt.savefig(root + "intens_plot.png")
    plt.close()

    return


if __name__ == "__main__":
    main()
