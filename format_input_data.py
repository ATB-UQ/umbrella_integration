# This module offers some functions to help reformat input data to be passed to umbrella_integration.py.
# More specifically, it reformats the first two columns of the input_file to be a list of lines of the form:
#  time  position  equilibrium_position  force_constant(kJ mol^-1 distance^-2)
#
# Note: assumes the first 2 columns in input_file are white-space separated and of the form:
# time(ps) reaction_cooridinate_position
from argparse import ArgumentParser
from os.path import basename
import numpy as np
import glob
import gzip
import itertools
import logging

logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s]: %(message)s')

SKIP_LINE_CHARACTERS = ["#", "@"]
UI_INPUT_FILE_HEADER = "#  time  position  equilibrium_position  force_constant(kJ mol^-1 distance^-2)"
UI_INPUT_LINE_TEMPLATE = "{0:>10g}   {1:>10g}      {2:>10g}                  {3:>10g}"

# Reformat the first two columns of the input_file to be a list of lines of the form:
# time  position  equilibrium_position  force_constant(kJ mol^-1 distance^-2)
#
# Note: assumes the first 2 columns in input_file are white-space separated and of the form:
# time reaction_cooridinate_position
#
# transform_from_relative: Some simulation packages write out instantaneous positions relative to the equilibrium
#     position of the harmonic potential, setting transform_from_relative=True will apply a shift to get absolute
#     position along the reaction coordinate.
def generate_ui_input_lines(input_file, eq_pos, fc, transform_from_relative=False):
    time, pos = parse_input_file(input_file)
    if transform_from_relative:
        pos = pos + eq_pos
    n = len(time)
    return [UI_INPUT_LINE_TEMPLATE.format(t, p, eq_p, f) for t, p, eq_p, f in zip(time, pos, [eq_pos]*n, [fc]*n)]

# Reformat and merge data into a single file.
# run_data needs to be a list of the form:
#     [[file_1, equilibrium_position_1, force_constant_1], [file_2, equilibrium_position_2, force_constant_2], ...]
#
# transform_from_relative: Some simulation packages write out instantaneous positions relative to the equilibrium
#     position of the harmonic potential, setting transform_from_relative=True will apply a shift to get absolute
#     position along the reaction coordinate.
def generate_ui_input_lines_merged(run_data, transform_from_relative=False):
    ui_input_lines = []
    # sort data on eq_pos first and then input_file
    for input_file, eq_pos, fc in sorted(run_data, key=lambda x:(x[2], x[0])):
        logging.debug("Processing {0}".format(input_file))
        ui_input_lines.extend(generate_ui_input_lines(input_file, eq_pos, fc, transform_from_relative))
    return ui_input_lines

def write_output_file(output_file, ui_input_lines):
    with open(output_file, "w") as fh:
        fh.write("\n".join([UI_INPUT_FILE_HEADER] + ui_input_lines))

def parse_input_file(input_file, keep_every=1):
    open_function = gzip.open if input_file.endswith(".gz") else open
    fh = open_function(input_file)
    lines = fh.read().splitlines()
    #print "Parsing file, keeping every {0}(st/nd/rd/th) value".format(keep_every)
    line_iterator = itertools.islice(lines, None, None, keep_every)
    lines = np.array([map(float, l.split()[:2]) for l in line_iterator if l.strip()[0] not in SKIP_LINE_CHARACTERS])
    fh.close()
    return lines.T

def parse_commandline():

    parser = ArgumentParser(prog=basename(__file__))
    # python -i <input data file> -r <reaction coordinate equilibrium position> -f <force constant> [-s]")
    # required
    parser.add_argument("-i", "--input_file",
        help="Input data file, white space separated columns of the form: col1=time, col2=position",
        action="store", required=True)

    parser.add_argument("-o", "--output_file",
        help="Output data file: col1=time, col2=position, col5=umbrella_equilibrium_position, col4=force_constant (kJ mol^-1 distance^-2)",
        action="store", required=True)

    parser.add_argument("-f", "--force_constant",
        help="Force constant used for the umbrella potential (kJ mol^-1 distance^-2)",
        action="store", required=True, type=float)

    parser.add_argument("-r", "--reaction_coordinate_position",
        help="Equilibrium position of the umbrella potential along the reaction coordinate",
        action="store", required=True, type=float)

    parser.add_argument("-s", "--transform_from_relative",
        help="Shift all reaction coordinate positions by value set for reaction_coordinate_position (helpful since some packages output positions relative to the umbrella equilibrium position)",
        action="store_true", required=False)

    args = parser.parse_args()

    return args.input_file, args.output_file, args.force_constant, args.reaction_coordinate_position, args.transform_from_relative

def run_umbrella_integration():
    input_file, output_file, fc, eq_pos, transform_from_relative = parse_commandline()
    ui_input_lines = generate_ui_input_lines(input_file, eq_pos, fc, transform_from_relative)
    write_output_file(output_file, ui_input_lines)

if __name__=="__main__":
    run_umbrella_integration()
    # Below is an example of a more advanced use of this module to reformat and merge data from
    # multiple umbrella simulations. Also see examples/example_analysis_script.py
    if False:
        # example script for generating merged umbrella_integration.py input file
        # force constant is the same for all windows
        fc = 500000.0 # J mol-1 nm-2
        fc /= 1000.0 # kJ mol-1 nm-2
        run_data = []
        # use a glob pattern to list relavent files
        for f in glob.glob("Nafc500*.dat"):
            # extract equilibrium position (window_center) from filename
            window_center = float(basename(f).split("_")[3])
            run_data.append( (f, window_center, fc) )
        ui_input_lines = generate_ui_input_lines_merged(run_data, transform_from_relative=False)
        write_output_file("na_ion_merged_ui_input.dat", ui_input_lines)