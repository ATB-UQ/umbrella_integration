# Example script for running umbrella integration
import sys
import glob
from os.path import basename

sys.path.append("../")
from format_input_data import generate_ui_input_lines_merged, write_output_file
from umbrella_integration import run_umbrella_integration

# example of the generation of merged ui input file
fc = 500000.0 # J mol-1 nm-2
fc /= 1000.0 # kJ mol-1 nm-2
run_data = []
umbrella_integration_data_file = "./data/merged_ui_input.dat"
for f in glob.glob("./data/raw_data/Nafc500*.dat"):
    window_center = float(basename(f).split("_")[3])
    run_data.append( (f, window_center, fc) )
ui_input_lines = generate_ui_input_lines_merged(run_data, transform_from_relative=False)
write_output_file(umbrella_integration_data_file, ui_input_lines)

run_umbrella_integration(
    [umbrella_integration_data_file],
    298.15,
    bin_width=0.2,
    minimum_maximum_value=[2.5, 5.1],
    integration_error_reference="left",
    plot_derivatives="derivatives.eps",
    plot_pmf="pmf.eps",
    position_histogram_plot="position_histograms.png",
    output_pmf_file="example.pmf",
    n_blocks=15,
    integration_method="trapezoidal",
    )
