#!/usr/bin/env python3

from subprocess import Popen, PIPE


def main():
    command = "find . -type f -name '*.pf'"
    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    err = stderr.decode("utf-8")
    if err:
        print(err)

    print("Simulations found:")
    print(output)
    with open("simulations", "w") as f:
        for sim in output:
            sim = sim.strip("\r\n")
            f.write("{}\n".format(sim))

    return output


if __name__ == "__main__":
    main()