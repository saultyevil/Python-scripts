#!/usr/bin/env python

from consts import G, MPROT, C, THOMPSON, PI, MSOL, MSOL_PER_YEAR


def mdot_edd(m_bh, efficiency):
    """
    Calculate the Eddington accretion limit for a black hole. Note that the
    accretion rate can be larger than the Eddington accretion rate. See, for
    example, Foundations of High-Energy Astrophysics by Mario Vietri.

    Parameters
    ----------
    m_bh: float
        The mass of the black hole in units of grams.
    efficiency: float
        The efficiency of the accretion process. Less than 1.

    Returns
    -------
    The Eddington accretion rate in units of grams / second.
    """
    assert(efficiency <= 1)
    assert(m_bh > 0)
    return (4 * PI * G * m_bh * MPROT) / (efficiency * C * THOMPSON)


def l_edd(m_bh):
    """
    Calculate the Eddington luminosity for accretion onto a black hole.

    Parameters
    ----------
    m_bh: float
        The mass of the black hole in units of grams.

    Returns
    -------
    The Eddington luminosity for the black hole in units of ergs / second.
    """
    assert(m_bh > 0)
    return (4 * PI * G * m_bh * C * MPROT) / THOMPSON


if __name__ == "__main__":

    mdot = 1e-3
    m_bh = 10 ** 6
    efficiency = 0.1

    print("Efficiency       = {:3.2f}".format(efficiency))
    print("M_bh             = {:3.2e} msol".format(m_bh))
    print("Mdot_disc        = {:3.2e} g/s".format(mdot * MSOL_PER_YEAR))
    print("Mdot_disc        = {:3.2f} msol/yr".format(mdot))
    print("L_edd            = {:3.2e} ergs/s".format(l_edd(m_bh * MSOL)))
    print("Mdot_edd         = {:3.2e} g/s".format(mdot_edd(m_bh * MSOL, efficiency)))
    print("Mdot_edd         = {:3.2f} msol/yr".format(mdot_edd(m_bh * MSOL, efficiency) / MSOL_PER_YEAR))
    print("Mdot / Mdot_edd  = {:3.2f}".format(mdot * MSOL_PER_YEAR / mdot_edd(m_bh * MSOL, efficiency)))
