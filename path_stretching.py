#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import time
from matplotlib import pyplot as plt

CHAND_SLAB_SOL = np.array(
        [[00.00, 1.26938],
         [18.19, 1.22945],
         [25.84, 1.18947],
         [31.79, 1.14943],
         [36.87, 1.10931],
         [41.41, 1.06911],
         [45.57, 1.02882],
         [49.46, 0.98842],
         [53.13, 0.94789],
         [56.63, 0.90722],
         [60.00, 0.86637],
         [63.26, 0.82530],
         [66.42, 0.78938],
         [69.51, 0.74234],
         [72.54, 0.70029],
         [75.52, 0.65770],
         [78.46, 0.61439],
         [81.37, 0.57001],
         [84.26, 0.52397],
         [87.13, 0.47490],
         [90.00, 0.41441]])


# Constants which control the general workings of the simulation

NEUTRAL_WEIGHT = 1
MIN_WEIGHT = 1e-5 * NEUTRAL_WEIGHT
MAX_WEIGHT = 1e5 * NEUTRAL_WEIGHT
LMAX = 1
SEED = 23165489
VERBOSE = False
N_BINS = 30

# Parameters which control photon transport and path stretching

SCAT_ALBEDO = 1      # 1 == pure scattering, 0 == pure absorption
N_PHOTONS = int(1e7)
P_KILL = 0.5
TAU_MAX = 5
SCAT_THRESHOLD = 10
TAU_THRESHOLD = np.sqrt(SCAT_THRESHOLD)


class IntensitySpectrum:
    """
    Contains all the guff to create the intensity spectrum
    """

    def __init__(self, nbins):
        """
        Initialise the spectrum structures
        """

        self.nbins = nbins
        self.hist = np.zeros(nbins)
        self.theta = np.zeros(nbins)
        self.intensity = np.zeros(nbins)

        return

    def init_bins(self):
        """
        Initialise the theta bins for the spectrum
        """

        dtheta = 1 / self.nbins
        hwidth = 0.5 * dtheta
        for i in range(self.nbins):
            self.theta[i] = np.arccos(i * dtheta + hwidth)

        return

    def bin_photon_direction(self, costheta, weight):
        """
        Bin a photon weight into the appropriate theta bin
        """

        bin = np.abs(int(costheta * self.nbins))
        self.hist[bin] += weight

        return

    def calculate_intensity(self):
        """
        Calculate the intensity from the hist bins
        """

        self.intensity = (self.hist * self.nbins) / (2 * N_PHOTONS * np.cos(self.theta))

        return


class PhotonPacket:
    """
    Contains all of the functions required to generate and transport a photon.
    """

    def __init__(self):
        """
        Initialise a photon
        """

        self.nrr = 0
        self.cscat = 0
        self.tau_cell = 0
        self.nstretch = 0
        self.tau_scat = 0
        self.weight = NEUTRAL_WEIGHT
        self.coordinates = np.zeros(3)
        self.tau_depth = TAU_MAX
        self.costheta = np.sqrt(np.random.rand())
        self.sintheta = np.sqrt(1 - self.costheta ** 2)
        self.cosphi = np.cos(2 * np.pi * np.random.rand())
        self.sinphi = np.sqrt(1 - self.cosphi ** 2)
        self.in_slab = True

        return

    def play_russian_roulette(self):
        """
        Play a game of Russian Roulette with the photon
        """

        self.nrr += 1
        if np.random.rand() <= P_KILL:
            return -1
        self.weight *= 1 / (1 - P_KILL)

        return self.weight

    def find_tau_cell(self):
        """
        Find the maximum optical depth to escape
        """

        self.tau_cell = 0
        if self.costheta > 0:
            self.tau_cell = self.tau_depth / self.costheta
        elif self.costheta < 0:
            self.tau_cell = (self.tau_depth - TAU_MAX) / self.costheta

        return self.tau_cell

    def reweight_photon(self, alpha):
        """
        Reweight a photon after it has taken a stretched path

        alpha   float
                the stretching parameter, has a value 0 <= alpha <= 1
        """

        self.weight *= np.exp(self.tau_scat * (alpha - 1)) / alpha

        return self.weight

    def generate_bias_tau_scat(self, alpha):
        """
        Sample from q(tau) = alpha * e ** (- alpha * tau), where alpha is the
        stretching parameter of the distribution between 0 and 1

        alpha   float
                the stretching parameter, has a value 0 <= alpha <= 1
        """

        if 0 > alpha > 1:
            return -1
        self.nstretch += 1
        self.tau_scat = -1.0 * (np.log(1 - np.random.rand()) / alpha)
        self.reweight_photon(alpha)

        return self.tau_scat

    def generate_tau_scat(self):
        """
        Sample from p(tau) = e ** -tau
        """

        self.tau_scat = -1.0 * np.log(1 - np.random.rand())

        return self.tau_scat

    def isotropic_scatter_photon(self):
        """
        Isotropic scattering of a photon
        """

        self.costheta = 2 * np.random.rand() - 1
        self.sintheta = np.sqrt(1 - self.costheta ** 2)
        self.cosphi = np.cos(2 * np.pi * np.random.rand())
        self.sinphi = np.sqrt(1 - self.cosphi ** 2)

        return

    def update_photon_position(self, ds):
        """
        Update the position of a photon

        ds  float
            distance the photon is being moved
        """

        self.coordinates[0] += ds * self.sintheta * self.cosphi
        self.coordinates[1] += ds * self.sintheta * self.sinphi
        self.coordinates[2] += ds * self.costheta
        self.tau_depth -= self.costheta * self.tau_scat

        return


def trans_phot(path_stretch_on):
    """
    The main loop which will transport a photon through the slab

    path_stretch_on     bool
                        flag for enabling path stretching
    """

    spectrum = IntensitySpectrum(N_BINS)
    spectrum.init_bins()

    if SEED:
        np.random.seed(SEED)

    n_scats_tot = 0           # Track total number of scatters in simulation
    tot_nrr = 0               # Track total number of times RR played
    n_terminated = 0          # Track how many photons were killed by RR
    n_survived = 0            # Track how many photons survived playing RR
    n_stretch_tot = 0         # Track total number of stretched paths

    for iphot in range(N_PHOTONS + 1):

        phot = PhotonPacket()
        total_nscat = 0

        if VERBOSE:
            str = "    Photon {:5d}    ".format(iphot)
            print("-" * len(str))
            print("{}".format(str))
            print("-" * len(str), "\n")

        while LMAX >= phot.coordinates[2] >= 0 and phot.in_slab is True:

            # Generate a location for tau_scat, re-weight the photon and play RR
            # if required. Finally move the photon to the interaction...

            tau_scat = phot.generate_tau_scat()

            if path_stretch_on and phot.cscat >= SCAT_THRESHOLD:

                # find the optical depth to the edge of the slab
                tau_cell = phot.find_tau_cell()

                # Only do the extra work when photons are far enough away from
                # the edge of the slab
                if tau_cell > TAU_THRESHOLD:
                    # we use tau_cell to set the value of alpha, thus photons
                    # which are closer to the surface will take shorter steps
                    alpha = 1 / (1 + tau_cell)
                    weight_orig = phot.weight
                    tau_scat = phot.generate_bias_tau_scat(alpha)

                    # check for correct exit of generate_tau_scat, -1 indicates that
                    # the function failed probably because of a bad value of alpha
                    if tau_scat == -1:
                        tau_scat = phot.generate_tau_scat()
                        phot.weight = weight_orig

                    # Reset the number of times the photon has scattered in slab
                    phot.cscat = 0

                    # Play RR if the photon weight is below the minimum weight
                    if phot.weight < MIN_WEIGHT:
                        tot_nrr += 1
                        if phot.play_russian_roulette() == -1:
                            if VERBOSE:
                                print("this photon was killed by russian roulette")
                            n_terminated += 1
                            break
                        else:
                            n_survived += 1

            ds = tau_scat / TAU_MAX
            phot.update_photon_position(ds)

            # Check to see if photon if the photon is still in the slab and
            # scatter the photon if it is

            if phot.coordinates[2] < 0:
                phot = PhotonPacket()
                total_nscat = 0
            elif phot.coordinates[2] > LMAX:
                phot.in_slab = False
                spectrum.bin_photon_direction(phot.costheta, phot.weight)
            else:
                if np.random.rand() <= SCAT_ALBEDO:
                    phot.isotropic_scatter_photon()
                    phot.cscat += 1
                    total_nscat += 1
                else:  # photon is absorbed
                    break

        n_scats_tot += total_nscat
        n_stretch_tot += phot.nstretch

        if VERBOSE:
            print("nscat {:5d} nstretch {:5d} rrplays {:5d}\n".format(total_nscat, phot.nstretch, phot.nrr))

        if iphot % (N_PHOTONS / 20) == 0:
            print("{:3.0f}% photons transported".format(iphot / N_PHOTONS * 100.0))

    print("\n----------------------------------------------")
    print("path_stretch_on {}".format(path_stretch_on))
    print("TAU_MAX {}".format(TAU_MAX))
    print("P_KILL {}".format(P_KILL))
    print("SCAT_THRES {}".format(SCAT_THRESHOLD))
    print("TAU_THRES {}".format(TAU_THRESHOLD))
    print("n_scats_tot {}".format(n_scats_tot))
    print("n_stretch_tot {} ({:3.2f}% of scats)".format(n_stretch_tot, n_stretch_tot / n_scats_tot * 100.0))
    if tot_nrr != 0:
        print("{:3.2f}% of entire photon sample were killed".format(n_terminated / N_PHOTONS * 100.0))
        print("{:3.2f}% of photons who played RR survived".format(n_survived / tot_nrr * 100.0))
        print("tot_nrr {}".format(tot_nrr))
    print("----------------------------------------------\n")

    spectrum.calculate_intensity()

    return spectrum, n_scats_tot


def main():
    """
    Main control function - runs trans_phot with path stretching enabled first
    and then without path stretching. Plots the intensity as a function of angle
    for both simulations against some analytical solution for a slab.
    """

    # trans_phot with path stretching
    print("starting trans_phot w/ path stretching")
    start = time.time()
    spec_pas, nscats_pas = trans_phot(True)
    end = time.time()
    print("trans_phot with path stretching: {} seconds".format(int(end - start)))

    # trans_phot without path stretching
    print("\nstarting trans_phot w/o path stretching")
    start = time.time()
    spec, nscats = trans_phot(False)
    end = time.time()
    print("trans_phot without path stretching: {} seconds".format(int(end - start)))

    print("nscats with pas / nscats without pas = {}".format(nscats_pas / nscats))

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    ax.semilogy(np.cos(np.deg2rad(CHAND_SLAB_SOL[:, 0])), CHAND_SLAB_SOL[:, 1], label="actual")
    ax.semilogy(np.cos(spec_pas.theta), spec_pas.intensity, label="with path stretching")
    ax.semilogy(np.cos(spec.theta), spec.intensity, label="without path stretching")
    ax.set_xlabel(r"$cos(\theta)$")
    ax.set_ylabel(r"Intensity")
    ax.legend()
    plt.savefig("path_stretching.png")
    plt.show()

    return


if __name__ == "__main__":

    """
    Here are the global variables which control the general performance and
    magic of the algorithm above. Changing them will change the amount of noise
    introduced and the general performance.
    """

    VERBOSE = False
    N_PHOTONS = int(1e4)
    SCAT_ALBEDO = 1        # 1 == pure scattering, 0 == pure absorption
    P_KILL = 0.5
    TAU_MAX = 5
    SCAT_THRESHOLD = 10
    TAU_THRESHOLD = np.sqrt(SCAT_THRESHOLD)

    main()
