#!/usr/bin/env python


from sys import argv
import numpy as np
from matplotlib import pyplot as plt
from PyPython.Constants import PARSEC, C
from PyPython.SpectrumUtils import smooth, read_spec
from PyPython.SpectrumUtils import common_lines, plot_line_ids
from PyPython.PythonUtils import file_len, array_index
from tqdm import tqdm


def construct_spectrum_from_weights(photons: np.ndarray, spectrum: np.ndarray, spec_norm: float,
                                    nphotons: int = None, nspec: int = None, nbins: int = None, nline: int = -1,
                                    dnorm: float = 100) -> np.ndarray:
    """
    Construct a spectrum from photon weights.

    If nphotons, nspec or nbins are not provided, then the function will
    automatically detect what these numbers are.

    There should be no NaN values in the photons array.

    Parameters
    ----------
    photons: np.ndarray (nphotons, 3)
        An array containing the photon frequency, weight and spectrum number.
    spectrum: np.ndarray (nbins, nspec)
        An array containing the frequency bins and empty bins for each
        inclination angle.
    spec_norm: float
        The spectrum normalization amount - usually the number of spectrum
        cycles.
    nphotons: [optional] int
        The number of photons in the photons array.
    nspec: [optional] int
        The number of inclination angles to bin.
    nbins: [optional] int
        The number of frequency bins in the spectrum.
    dnorm: [optional]
        The distance normalization for the flux calculation in parsecs. By
        default this is 100 parsecs.

    Returns
    -------
    spectrum: np.ndarray (nbins, nspec)
        The constructed spectrum in units of F_lambda erg/s/cm/cm/A.
    """

    n = construct_spectrum_from_weights.__name__

    assert(photons.shape[0] != 0)
    assert(photons.shape[1] == 3)
    assert(spectrum.shape[0] != 0)
    assert(spectrum.shape[1] > 1)

    if np.isnan(photons).any():
        print("{}: There are NaN values in photon array. Cannot bin so returning!".format(n))
        return spectrum

    if not nphotons:
        nphotons = photons.shape[0]
    if not nbins:
        nbins = spectrum.shape[0]
    if not nspec:
        nspec = np.max(photons[:, 3]) + 1

    distance_normalisation = 4 * np.pi * (dnorm * PARSEC) ** 2

    photon_freqs = photons[:, 0]
    photon_weights = photons[:, 1]
    spec_index = photons[:, 2].astype(int) + 1  # + 1 as index 0 in spectrum is frequency bins

    # Slowly bin the photons into frequency space
    for p in tqdm(range(nphotons), desc="Binning photons by frequency", unit=" photons", smoothing=0):
        freq_index = array_index(spectrum[:, 0], photon_freqs[p])
        if freq_index == -1:
            continue
        spectrum[freq_index, spec_index[p]] += photon_weights[p]

    # Now convert the binned weights into a flux
    for n in tqdm(range(nspec), desc="Converting weights to flux", unit=" bins", smoothing=0):
        for p in range(nbins - 1):
            freq = spectrum[p, 0]
            dfreq = spectrum[p + 1, 0] - spectrum[p, 0]
            spectrum[p, n + 1] *= (freq ** 2 * 1e-8) / (dfreq * distance_normalisation * C)
        spectrum[:, n + 1] /= spec_norm

    return spectrum


def create_spectra_from_delay_dump(root: str, loc: str = ".", fmin: float = None, fmax: float = None, nbins: int = None,
                                   spec_norm: float = 1, verbose: bool = False) -> np.ndarray:
    """
    Create a spectrum for each inclination angle using the filtered
    root.delay_dump diagnostic file.

    Parameters
    ----------
    root: str
        The root name of the simulation.
    loc: [optional] str
        The directory containing the simulation.
    fmin: [optional] float
        The smallest frequency bin.
    fmax: [optional] float
        The largest frequency bin
    nbins: [optional] int
        The number of frequency bins.
    spec_norm: float
        A normalization constant for the spectrum. Usually the number of photon
        cycles.
    verbose: [optional] bool
        Enable verbose logging.

    Returns
    -------
    filtered_spectrum: np.ndarray
        A 2D array containing the frequency in the first column and the
        fluxes for each inclination angle in the other columns.
    """

    file = "{}/{}.delay_dump".format(loc, root)
    nlines = file_len(file)

    file_columns = {
        "Freq.": 0, "Lambda": 1, "Weight": 2, "Last X": 3, "Last Y": 4,
        "Last Z": 5, "Scat.": 6, "RScat": 7, "Delay": 8, "Spec.": 9,
        "Orig.": 10, "Res.": 11, "LineRes.": 12
    }

    extract = ["Freq.", "Weight", "Spec."]  # , "Scat", "Rscat", "Orig", "Res.", "LineRes."]
    nextract = len(extract)
    dumped_photons = np.zeros((nlines, nextract))
    dumped_photons[:, :] = np.nan   # because some lines will not be photons mark them as NaN

    # Read in photons line by line to extract frequency, weight and nspec
    with open(file, "r") as f:
        for i in tqdm(range(nlines), desc="Reading {}.delay_dump".format(root), unit=" lines", smoothing=0):
            line = f.readline().split()
            if len(line) == 13:
                for j, e in enumerate(extract):
                    dumped_photons[i, j] = float(line[file_columns[e]])

    # Remove the nan lines from the photons array
    dumped_photons = dumped_photons[~np.isnan(dumped_photons).any(axis=1)]
    nphotons = dumped_photons.shape[0]
    nspec = int(np.max(dumped_photons[:, 2])) + 1

    if not fmin:
        fmin = np.min(dumped_photons[:, 0])
    if not fmax:
        fmax = np.max(dumped_photons[:, 0])
    if not nbins:
        nbins = int(1e4)

    if verbose:
        print("\n----------------------------------")
        print(" There are {:2.1e} photons to bin".format(nphotons))
        print(" There are {} spectra to create".format(nspec))
        print(" Freq. min {:2.1e} Hz".format(fmin))
        print(" Freq. max {:2.1e} Hz".format(fmax))
        print("----------------------------------\n")

    # Now bin the photons to create the filtered spectrum
    filtered_spectrum = np.zeros((nbins, 1 + nspec))
    filtered_spectrum[:, 0] = np.logspace(np.log10(fmin), np.log10(fmax), nbins, endpoint=True)
    filtered_spectrum = construct_spectrum_from_weights(dumped_photons, filtered_spectrum, spec_norm, nphotons, nspec,
                                                        nbins)

    # TODO: improve writing out spectrum, e.g. add headers, output lambda etc.
    np.savetxt("{}/{}.delay_dump.spec".format(wd, root), filtered_spectrum)

    return filtered_spectrum


if __name__ == "__main__":

    wd = "."
    root = "tde_all_filter"

    full_spectrum = read_spec("{}/{}.log_spec".format(wd, root))
    fmin = np.min(full_spectrum["Freq."])
    fmax = np.max(full_spectrum["Freq."])
    nbins = len(full_spectrum["Freq."])
    norm = 10

    filtered_spectrum = create_spectra_from_delay_dump(root, wd, fmin, fmax, nbins, norm, True)

    smooth_amount = 5
    inclinations = ["10", "30", "45", "60", "75", "85"]
    print("Plotting spectra for ", inclinations, "degrees")
    for i, inclination in enumerate(inclinations):
        fig, ax = plt.subplots(figsize=(12, 5))

        ax.plot(C * 1e8 / filtered_spectrum[:-1, 0], smooth(filtered_spectrum[:-1, i + 1], smooth_amount),
                label="Filtered spectrum", linewidth=1.4, alpha=0.75)

        ax.plot(full_spectrum["Lambda"].values, smooth(full_spectrum[inclination].values, smooth_amount),
                label="Spectrum", linewidth=1.4, alpha=0.75)

        ax.set_yscale("log")
        ax.set_xscale("log")

        ax.set_xlabel(r"Wavelength [$\AA$]", fontsize=15)
        ax.set_ylabel(r"Flux at 100 Pc [ergs s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]", fontsize=15)
        ax.legend(fontsize=13)
        ax = plot_line_ids(ax, common_lines(), logx=True)

        fig.tight_layout(rect=[0.015, 0.015, 0.985, 0.985])
        fig.savefig("{}/{}_i{}.delay_dump_spectrum.png".format(wd, root, inclination), dpi=300)
        fig.savefig("{}/{}_i{}.delay_dump_spectrum.pdf".format(wd, root, inclination), dpi=300)

        plt.close()
