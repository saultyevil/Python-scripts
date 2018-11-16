#!/usr/bin/env python3

import argparse
from platform import system
from subprocess import Popen, PIPE


def avg_world_ping(world_number, n_tests):
    """
    Returns the average latency of an OldSchool RuneScape world over n_tests.

    Parameters
    ----------
    world_number: int.
        The world to ping.
    n_tests: int.
        The number of times to ping the world.

    Returns
    -------
    avg_latency: float.
       The average latency over all ping tests of the world pinged.
    """

    if type(world_number) != int:
        exit_msg = "world_number is not an int"
        print(exit_msg)
        return 999, exit_msg
    if type(n_tests) != int:
        exit_msg = "n_tests is not an int"
        print(exit_msg)
        return 999, exit_msg
    if n_tests <= 1:
        exit_msg = "bad value n_tests <= 1"
        print(exit_msg)
        return 999, exit_msg

    # Create the ping command to send to the terminal - note that the repeat
    # command is different on Windows and Unix based operating systems
    OS = system().lower()
    hostname = "oldschool{}.runescape.com".format(world_number)
    if OS == "windows":
        repeat = "-n {}".format(n_tests)
    elif OS == "darwin" or OS == "linux":
        repeat = "-c {}".format(n_tests)
    ping_command = "ping {} {}".format(repeat, hostname)

    print("{}".format(ping_command))

    # Now send the ping command to the terminal and capture stdout and stderr,
    # then decode the bytes string into UTF-8 and split the list
    ping = Popen(ping_command, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = ping.communicate()
    stdout_split = stdout.decode("utf-8").split()

    return_code = ping.returncode
    if return_code != 0:
        exit_msg = "ping exited with exit status".format(return_code)
        print("{} exited with non-zero exit status {}:".format(ping_command, return_code))
        if return_code == 1:
            print("Error code 1 means no reply packets received.")
            print("Bad internet or hostname ({})?".format(hostname))
        else:
            print("Check the man page for return code {}".format(return_code))
        return 999, exit_msg

    # Now we need to find where the average latency is output - note that this
    # is different between Windows and Unix based operating systems again
    avg_ping_idx = -1
    if OS == 'windows':
        for i in range(len(stdout_split)):
            out = stdout_split[i]
            if out == "Average":
                avg_ping_idx = i + 2
                break
        if avg_ping_idx == -1:
            raise IndexError("Couldn't find average latency for {}".format(hostname))
        avg_latency = stdout_split[avg_ping_idx]
    elif OS == "darwin" or OS == "linux":
        for i in range(len(stdout_split)):
            out = stdout_split[i]
            if out == "min/avg/max/stddev" or out == "min/avg/max/mdev":
                avg_ping_idx = i + 2
                break
        if avg_ping_idx == -1:
            raise IndexError("Couldn't find average latency for {}".format(hostname))
        avg_latency = stdout_split[avg_ping_idx].replace("/", " ").split()[1]

    exit_msg = "successful exit"

    return avg_latency, exit_msg


def parse_world_tests_args():
    """
    Parse the choice for worlds to ping and the number of repeats for the ping command

    Returns
    -------
    worlds: list of ints
        A list of the worlds which are to be pinged
    n_tests: int
        The number of times each world should be pinged
    """

    # Parse arguments from the command line, note that nargs="?" indicates that n_tests is optional
    parser = argparse.ArgumentParser(description="Input the worlds to be tested and optionally the number of ping repeats.")
    parser.add_argument("world_list", type=str, help="A single or list of worlds to ping.")
    parser.add_argument("n_tests", nargs="?", type=int, help="The number of times to ping tested worlds.")
    args = parser.parse_args()

    n_tests = 5  # default value of n_tests
    world_type = args.world_list.lower()

    # Set up the world choices
    if world_type == "uk":
        print("Pinging all of the non-PvP UK worlds:")
        worlds = [2, 9, 10, 17, 18, 33, 34, 41, 42, 49, 50, 58, 65, 66, 73]
    elif world_type == "allworlds":
        print("Pinging all of the RuneScape worlds:")
        worlds = range(1, 125)
    elif world_type == "custom_list":
        print("List the worlds wanted to be pinged separated by a space, e.g. 1 12 42:")
        custom_worlds = input("World choices > ")
        worlds = [int(i) for i in custom_worlds.split()]
    elif world_type.isdigit() is True:
        worlds = [int(world_type)]
    else:
        world_types_choices = ["uk", "allworlds", "custom_list", "an integer for a single world"]
        raise ValueError("Disallowed values for world choice. Allowed values are:\n\t{}".format(world_types_choices))
    n_worlds = len(worlds)

    # If n_tests has been provided, check for valid values else exit
    if args.n_tests:
        n_tests = args.n_tests
        if n_tests <= 1:
            raise ValueError("Disallowed value for n_tests: n_tests > 1")
    else:
        print("Using default values of n_tests {}".format(n_tests))

    return worlds, n_tests


def main():
    """
    The main control function for the script
    """

    worlds, n_tests = parse_world_tests_args()
    n_worlds = len(worlds)

    print("n_tests: {}".format(n_tests))
    print("Worlds: {}".format(worlds))

    if n_tests > 5 and n_worlds > 5:
        print("ALART: n_worlds {} n_tests {} may take a while to run".format(n_worlds, n_tests))
    print("")

    # Now iterate through each world and call the ping test function
    n_err = 0
    pings = []
    err_worlds = []
    err_exit_msg = []

    for world in worlds:
        world_ping, exit_msg = avg_world_ping(world, n_tests)
        if world_ping == 999:
            n_err += 1
            err_worlds.append(world)
            err_exit_msg.append(exit_msg)
            continue
        pings.append(world)
        pings.append(world_ping)

    # Figure out which world had the lowest ping
    # The latency is index 1, 3, 5, ... etc and 0, 2, 4, .. etc is the
    # world number
    lowest_ping = min(pings[1::2])
    lowest_ping_world_number = pings[pings.index(lowest_ping) - 1]

    if n_err:
        print("Number of errors with ping encountered: {}".format(n_err))
        for i in range(n_err):
            print("World {}: {}".format(err_worlds[i], err_exit_msg[i]))
        if n_err == n_worlds:
            print("n_err = n_worlds, the script should exit now")
            return 1

    print("\n-------------------------------------------------------------------------------\n")
    print("\tLowest ping world:\n\tWorld {} with an average ping of {} ms over {} tests"
          .format(lowest_ping_world_number, lowest_ping, n_tests))
    print("\n-------------------------------------------------------------------------------")

    return


if __name__ == '__main__':
    main()
