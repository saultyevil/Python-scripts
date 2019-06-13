#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

NEUTRAL_WEIGHT = 1
MIN_WEIGHT = 1e-5 * NEUTRAL_WEIGHT
MAX_WEIGHT = 1e5 * NEUTRAL_WEIGHT
P_KILL = 0.5
TAU_MAX = 5
LMAX = 1
SCAT_ALBEDO = 1
SCAT_THRES = 10
N_PHOTONS = 10
SEED = 23165489


class Photon:
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

        if self.costheta > 0:
            self.tau_cell = self.tau_depth / self.costheta
        elif self.costheta < 0:
            self.tau_cell = (self.tau_depth - TAU_MAX) / self.costheta
        return self.tau_cell

    def reweight_photon(self, alpha):
        """
        Reweight a photon after it has taken a stretched path
        """

        self.weight *= np.exp(self.tau_scat * (alpha - 1)) / alpha
        return self.weight

    def generate_bias_tau_scat(self, alpha):
        """
        Sample from q(tau) = alpha * e ** (- alpha * tau), where alpha is the
        stretching parameter of the distribution between 0 and 1
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
        """

        self.coordinates[0] += ds * self.sintheta * self.cosphi
        self.coordinates[1] += ds * self.sintheta * self.sinphi
        self.coordinates[2] += ds * self.costheta
        self.tau_depth -= self.costheta * self.tau_scat
        return


def trans_phot():
    """
    The main loop which will transport a photon through the slab
    """

    if SEED:
        np.random.seed(SEED)

    n_scats_tot = 0
    n_terminated = 0
    n_stretch_tot = 0

    for i in range(N_PHOTONS + 1):

        phot = Photon()
        phot_nscat = 0

        str = "    Photon {:5d}    ".format(i)
        print("-" * len(str))
        print("{}".format(str))
        print("-" * len(str), "\n")

        while LMAX >= phot.coordinates[2] >= 0 and phot.in_slab is True:

            #
            # Generate a location for tau_scat, re-weight the photon and play RR
            # if required. Finally move the photon to the interaction...
            #

            if phot.cscat < SCAT_THRES:
                tau_scat = phot.generate_tau_scat()
            else:
                tau_cell = phot.find_tau_cell()
                alpha = 1 / (1 + tau_cell)
                tau_scat = phot.generate_bias_tau_scat(alpha)
                phot.cscat = 0

                if phot.weight < MIN_WEIGHT:
                    if phot.play_russian_roulette() == -1:
                        print("killed")
                        n_terminated += 1
                        break

            ds = tau_scat / TAU_MAX
            phot.update_photon_position(ds)

            #
            # Check to see if photon if the photon is still in the slab, otherwise
            # scatter the photon
            #

            if phot.coordinates[2] < 0:
                phot = Photon()
                phot_nscat = 0
            elif phot.coordinates[2] > LMAX:
                phot.in_slab = False
            else:
                if np.random.rand() <= SCAT_ALBEDO:
                    phot.isotropic_scatter_photon()
                    phot.cscat += 1
                    phot_nscat += 1

        n_scats_tot += phot_nscat
        n_stretch_tot += phot.nstretch

        print("nscat {:5d} nstretch {:5d} rrplays {:5d}\n".format(phot_nscat, phot.nstretch, phot.nrr))

    print("TAUMAX {}".format(TAU_MAX))
    print("PKILL {}".format(P_KILL))
    print("{:3.2f}% photons were killed".format(n_terminated / N_PHOTONS * 100.0))

    return


def main():
    """
    Main control function
    """

    trans_phot()

    return


if __name__ == "__main__":
    main()
