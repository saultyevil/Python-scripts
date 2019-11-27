from PyPython import Simulation
from PyPython import Utils
from matplotlib import pyplot as plt
import numpy as np

pfs = Utils.find_parameter_files("./")
for i in range(len(pfs)):
    root, wd = Utils.split_root_directory(pfs[i])
    convergence = Simulation.check_convergence_criteria(root, wd)
    n_tr = convergence[0]
    n_te = convergence[1]
    n_te_max = convergence[2]
    n_hc = convergence[3]
    convergence = Simulation.check_convergence(root, wd, return_per_cycle=True)
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    ax.plot(np.arange(0, len(convergence)), convergence, label="Convergence")
    ax.plot(np.arange(0, len(n_te), 1), n_te, "--", label="Electron temperature")
    ax.plot(np.arange(0, len(n_tr), 1), n_tr, "--", label="Radiation temperature")
    ax.plot(np.arange(0, len(n_hc), 1), n_hc, "--", label="Heating/Cooling")
    ax.plot(np.arange(0, len(n_te_max), 1), n_te_max, "--", label="Electron temperature max")
    ax.legend()
    ax.set_xlabel("Cycle")
    ax.set_ylabel("Fraction of cells")
    plt.savefig("{}/{}_convegence_detailed.png".format(wd, root))

