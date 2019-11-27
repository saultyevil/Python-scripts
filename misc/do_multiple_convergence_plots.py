from PyPython import Utils
from subprocess import Popen, PIPE

i = [2, 3, 8, 8, 6, 9, 24, 15]
j = [0, 0, 8, 7, 6, 6, 24, 15]
nx = 30
nz = 30

pfs = Utils.find_parameter_files("./")
for ii in range(len(pfs)):
    root, wd = Utils.split_root_directory(pfs[ii])
    print(pfs[ii])
    for jj in range(len(i)):
        cmd = "cd {}; py_plot.py {} tau_spec; py_plot_cell_convergence.py python {} {}; py_plot_cell_spec.py {} {} {} {} {}"\
            .format(wd, root, i[jj], j[jj], root, nx, nz, i[jj], j[jj])
        sh = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        stdout, stderr = sh.communicate()
        if stderr:
            print(stderr.decode("utf-8"))
Utils.remove_data_sym_links("./")
