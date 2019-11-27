from PyPython import Utils
from subprocess import Popen, PIPE
from time import sleep

pfs = Utils.find_parameter_files("./")
print(pfs)

#sleep(5)

for i in range(len(pfs)):
    print(pfs[i])
    root, wd = Utils.split_root_directory(pfs[i])
    sh = "cd {}; rm *.png; py_windsave_plot.py; py_check_run.py {}".format(wd, root)
    cmd = Popen(sh, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    if stderr:
        print(stderr.decode("utf-8"))
    print(stdout.decode("utf-8"))
    sh = "cd {};  rm out.mp4; ffmpeg -framerate 1 -pattern_type glob -i 'python*.png' -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4".format(wd)
    cmd = Popen(sh, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    if stderr:
        print(stderr.decode("utf-8"))
    print(stdout.decode("utf-8"))
