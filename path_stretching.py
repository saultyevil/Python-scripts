#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import time
from matplotlib import pyplot as plt
import path_stretching_plot_weights as pspw

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
SEED = 23165489
N_BINS = CHAND_SLAB_SOL.shape[0]
VERBOSE = False

# Parameters which control photon transport and path stretching, changing these
# will alter the performance of the script and algorithms

SCAT_ALBEDO = 1      # 1 == pure scattering, 0 == pure absorption
N_PHOTONS = int(1e7)
P_KILL = 0.0
TAU_MAX = 5
SCAT_THRESHOLD = 0
TAU_THRESHOLD = np.sqrt(SCAT_THRESHOLD)


print("")
print("SCAT_ALBEDO", SCAT_ALBEDO)
print("N_PHOTONS", N_PHOTONS)
print("P_KILL", P_KILL)
print("TAU_MAX", TAU_MAX)
print("SCAT_THRESHOLD", SCAT_THRESHOLD)
print("TAU_THRESHOLD", TAU_THRESHOLD)
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
        self.hist = np.zeros(nbins)
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

    def bin_photon_direction(self, costheta, weight):
        """
        Bin a photon weight into the appropriate theta bin.

        costheta    float
                    the escape direction of a photon packet
        weight      float
                    the weight of the photon
        """

        if weight == 0:
            return

        bin = np.abs(int(costheta * self.nbins))
        self.hist[bin] += weight
        self.nphot += 1

        return

    def calculate_intensity(self):
        """
        Calculate the intensity from the hist bins.
        """

        self.intensity = (self.hist * self.nbins) / (2 * N_PHOTONS * np.cos(self.theta))

        return


class PhotonPacket:
    """
    Contains all of the functions required to generate and transport a photon.
    """

    def __init__(self):
        """
        Initialise a photon. Photons are initialised at the origin of the slab.
        """

        self.nrr = 0                # number of rr games
        self.cscat = 0              # number of scatters between stretched paths
        self.tau_path = 0           # optical depth to escape
        self.nstretch = 0           # number of stretched paths
        self.tau_scat = 0           # optical depth till next scatter
        self.tau_depth = TAU_MAX    # optical depth within the slab
        self.weight = NEUTRAL_WEIGHT
        self.x = 0
        self.y = 0
        self.z = 0
        self.costheta = np.sqrt(np.random.rand())
        self.sintheta = np.sqrt(1 - self.costheta ** 2)
        self.cosphi = np.cos(2 * np.pi * np.random.rand())
        self.sinphi = np.sqrt(1 - self.cosphi ** 2)
        self.inslab = True

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

    def find_tau_path(self):
        """
        Find the maximum optical depth to escape
        """

        self.tau_path = 0
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

    def generate_tau_scat(self):
        """
        Sample from p(tau) = e ** -tau.
        """

        self.tau_scat = -1.0 * np.log(1 - np.random.rand())

        return self.tau_scat

    def generate_bias_tau_scat(self, alpha):
        """
        Sample from q(tau) = alpha * e ** (- alpha * tau), where alpha is the
        stretching parameter of the distribution between 0 and 1

        alpha   float
                the stretching parameter, has a value 0 <= alpha <= 1
        """

        if 0 > alpha > 1:
            print("error: expected 0 <= alpha <= 1, but got alpha = {}".format(alpha))
            print("       setting tau_scat form normal distribution")
            self.tau_scat = self.generate_tau_scat()
            return self.tau_scat

        self.nstretch += 1
        self.tau_scat = -1.0 * (np.log(1 - np.random.rand()) / alpha)
        self.reweight_photon(alpha)

        return self.tau_scat

    def isotropic_scatter_photon(self):
        """
        Isotropic scattering of a photon.
        """

        self.costheta = 2 * np.random.rand() - 1
        self.sintheta = np.sqrt(1 - self.costheta ** 2)
        self.cosphi = np.cos(2 * np.pi * np.random.rand())
        self.sinphi = np.sqrt(1 - self.cosphi ** 2)

        return

    def scatter_bottom_slab(self):
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

    def move_photon(self, ds):
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


def write_phots(phots, pas):
    """
    Write out a bunch of data for all of the photons

    phots       list[PhotonPacket]
                a list containing all of the photons
    pas         bool
                if True, append filename with pas indicator
    """

    fname = "photon_statistics"
    if pas:
        fname += "_pas"
    fname += ".txt"
    f = open(fname, "w")

    f.write("nrr\tnstretch\tweight\n")
    for i in range(N_PHOTONS):
        phot = phots[i]
        str = "{}\t{}\t{}".format(phot.nrr, phot.nstretch, phot.weight)
        f.write("{}\n".format(str))

    f.close()

    return


def trans_phot(path_stretch_on):
    """
    The main loop which will transport a photon through the slab

    path_stretch_on     bool
                        flag for enabling path stretching
    """

    rr_on = True
    if P_KILL == 0:
        rr_on = False

    photstore = []

    spectrum = IntensitySpectrum(N_BINS)
    spectrum.init_bins()

    if SEED:
        np.random.seed(SEED)

    nkilled = 0  # Track the number of photons killed by RR

    for iphot in range(N_PHOTONS):

        phot = PhotonPacket()
        photstore.append(phot)

        # Keep iterating whilst the photon is in the slab
        while phot.inslab:
            # Generate the optical depth to the next scattering event, this can
            # either be from the normal distribution or from the biased distribution
            tau_scat = phot.generate_tau_scat()
            if path_stretch_on and phot.cscat >= SCAT_THRESHOLD:
                tau_path = phot.find_tau_path()

                if tau_path >= TAU_THRESHOLD:
                    phot.cscat = 0
                    alpha = 1 / (1 + tau_path)
                    tau_scat = phot.generate_bias_tau_scat(alpha)

                    # Play RR if the photon weight is below the minimum weight
                    if rr_on and phot.weight < MIN_WEIGHT:
                        phot.play_russian_roulette()
                        if phot.weight == -1:
                            phot.weight = 0
                            nkilled += 1
                            break

            # Now that there is a value for tau_scat, move the photon to this
            # location
            ds = tau_scat / TAU_MAX
            phot.move_photon(ds)

            # Check boundary conditions of the slab, and scatter the photon if
            # it is still in the slab
            zphot = phot.z
            if zphot < 0:       # Has gone below the slab, reflect photon
                phot.scatter_bottom_slab()
            elif zphot > LMAX:  # Has escaped the slab from top
                spectrum.bin_photon_direction(phot.costheta, phot.weight)
                phot.inslab = False
            else:               # Else, photon scatters
                if np.random.rand() <= SCAT_ALBEDO:
                    phot.isotropic_scatter_photon()
                    phot.cscat += 1
                else:
                    phot.weight = 0
                    phot.inslab = 0

        if VERBOSE:
            print("photon {} nstretch {:5d} nrr {:5d}\n".format(iphot, phot.nstretch, phot.nrr))

        if (iphot + 1) % (N_PHOTONS / 10) == 0:
            print("{:3.0f}% photons transported".format(iphot / N_PHOTONS * 100.0))

    spectrum.calculate_intensity()

    if nkilled:
        print("{} photons killed by Russian Roulette".format(nkilled))

    write_phots(photstore, path_stretch_on)

    return spectrum


def main():
    """
    Main control function - runs trans_phot with path stretching enabled first
    and then without path stretching. Plots the intensity as a function of angle
    for both simulations against some analytical solution for a slab.
    """

    print("Running MCRT with path stretching enabled")
    start = time.time()
    pasMCRT = trans_phot(path_stretch_on=True)
    end = time.time()
    print("Run time: {} seconds".format(int(end - start)))

    print("\nRunning MCRT without any acceleration")
    start = time.time()
    MCRT = trans_phot(path_stretch_on=False)
    end = time.time()
    print("Run time: {} seconds".format(int(end - start)))

    print("\nOriginal photon weight before transport: {:2e}".format(N_PHOTONS * NEUTRAL_WEIGHT))
    print("\n{} photons contributed to pasMCRT spectrum".format(pasMCRT.nphot))
    print("Total weight in pasMCRT spectrum: {:.2e}".format(np.sum(pasMCRT.hist)))
    print("\n{} photons contributed to MCRT spectrum".format(MCRT.nphot))
    print("Total weight in MCRT spectrum: {:.2e}".format(np.sum(MCRT.hist)))

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    ax.semilogy(np.cos(np.deg2rad(CHAND_SLAB_SOL[:, 0])), CHAND_SLAB_SOL[:, 1], label="Analytic Soution")
    ax.semilogy(np.cos(pasMCRT.theta), pasMCRT.intensity, label="Path Stretching Enabled")
    ax.semilogy(np.cos(MCRT.theta), MCRT.intensity, label="Path Stretching Disabled")
    ax.set_xlabel(r"$\mu = cos(\theta)$")
    ax.set_ylabel(r"Normalised Intensity")
    ax.legend()

    fig.tight_layout()
    plt.savefig("path_stretching_taumax_{}.png".format(TAU_MAX, P_KILL))
    plt.show()

    pspw.plot_weight_file(True)

    return


if __name__ == "__main__":
    main()
