import numpy as np
from matplotlib import pyplot as plt

plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14

total_time = []
avg_phot_trans_ions = []
avg_phot_trans_spec = []

ncores = [
    1, 2, 4, 8, 16, 32, 40, 80, 120, 160, 200
]

model_root = "cv"

root = "/home/saultyevil/PySims/release/{}_parallel/ncores_".format(model_root)

# Loop over each model run with n processes

for n in ncores:

    # For each n process run, loop over each diag file

    print ("Checking nproc =", n)

    phot_trans_times_ions = []
    phot_trans_times_spec = []

    for nn in range(n):
        phot_trans_times = []
        fname = root + str(n) + "/diag_{}/{}_{}.diag".format(model_root, model_root, nn)
        with open(fname, "r") as f:
            lines = f.readlines()
        for l in lines:
            if l.find("photon transport completed in:") != -1:
                phot_trans_times.append(float(l.split()[5]))
            if l.find("At program completion, the elapsed TIME was") != -1 and nn == 0:
                total_time.append(float(l.split()[-1]))

        # Sample from cycle 1->19 for ion cycles, and all cycles for spec cycles

        phot_trans_times_ions.append(np.average(phot_trans_times[1:19]))
        phot_trans_times_spec.append(np.average(phot_trans_times[-10:]))

    avg_phot_trans_ions.append(np.average(phot_trans_times_ions))
    avg_phot_trans_spec.append(np.average(phot_trans_times_spec))

print("\nAverage Ionization Cycle per Nprocs\n", avg_phot_trans_ions)
print("\nAverage Spectrum Cycle per Nprocs\n", avg_phot_trans_spec)
print("\nTotal Time per Nprocs\n", total_time)

ions = np.array(avg_phot_trans_ions)
spec = np.array(avg_phot_trans_spec)
totl = np.array(total_time)

fig, ax = plt.subplots(1, 1, figsize=(13, 8))
ax.loglog(ncores, ncores, "k--", label="Perfect Scaling", linewidth=3)
ax.plot(ncores, ions[0] / ions, "D-", label="Ionization Cycles", linewidth=1.5, markersize=8)
ax.plot(ncores, spec[0] / spec, "o-", label="Spectrum Cycles", linewidth=1.5, markersize=8)
ax.plot(ncores, totl[0] / totl, "s-", label="Total Time", linewidth=1.5, markersize=8)
ax.legend(fontsize=14)
ax.set_xlabel(r"Number of Processors N$_{i}$", fontsize=14)
ax.set_ylabel(r"Parallel Speedup T$_{N = 1}$ / T$_{N = i}$", fontsize=14)
fig.tight_layout()
fig.savefig("parallel_speed.pdf", dpi=300)
fig.savefig("parallel_speed.png", dpi=300)
plt.show()
