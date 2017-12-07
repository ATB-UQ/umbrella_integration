# File: umbrella_integration.py
# Authors: Martin Stroet, Evelyne Deplazes
# Date: October 2016

# Calculate the potential of mean force (PMF) using the
# umbrella integration method as described in Kaestner and Thiel (2005)
#
# When using this code, please cite the following:
#
# GitHub project using the following DOI:
#    DOI: 10.5281/zenodo.164996
#
# Kaestner J. and Thiel W., J. Chem. Phys. 123, 144104 (2005)
# doi: 10.1063/1.2052648 
# Kaestner J. and Thiel W., J. Chem. Phys. 124, 234106 (2006)
# doi: 10.1063/1.2206775
#
# python umbrella_integration.py -i infile -o outfile -b binwidth
# 
# -i input_file (required)
#    Path(s) to the input file(s).
#    The input file needs to contain the following data in white-space separated column format:
#    col 1 = time (while not used in the PMF calculation directly it is useful when assessing
#         the convergence of umbrella simulations)
#    col 2 = instantaneous value of reaction coordinate
#    col 2 = centre (equilibrium position) of the harmonic potential used in the umbrella simulation
#    col 1 = harmonic force constant used for the umbrella simulation in kJ mol^-1 distance^-2
# -t temperature
#    Temperature at which the calculation is to be performed.
# -b bin_width (required if -n number_bins not provided)
#    The step size between point where the potential will be calculated along the reaction coordinate.
#    This is independent of the spacing between umbrellas/windows used 
#    in the umbrella sampling simulations.
# -n number_bins (required if -b bin_width not provided)
#    Number of positions along the reaction coordinate at which the potential will be calculated.
#    This is independent of the umbrellas/windows used in umbrella sampling simulations.
# -m min_max_value
#    Minimum and maximum reaction coordinate positions at which the potential will be calculated.
#    e.g. -m 0.2 5.5
#    Note that data beyond this range will be included in the calculation.
# -p show_plots
#    Interactively show plots of the derivative (which are numerically integrated) and final PMF
# -r integration_error_reference
#    Integration error estimation reference point, left or right
# -im integration_method
#    Numerical method used to the integrate derivatives. Simpon's method has a lower error
#    for smooth data, while the trapezoidal method is more robust on noisy data.
# -o output_pmf_file
#    Name of output file that contains the final PMF.
#    Format, col 1 = reaction coordinate, col 2 = energy, col 3 = uncertainty
# -d derivatives_file
#    Name of output file that contains the derivatives for each reaction coordinate point.
# -ph position_histrogram_plot
#    File name for the histogram plot.
# -np n_blocks
#    Number of blocks used in error analysis to calculate block averaged error estimate for each window.
# -im integration_method
#    Integration method used to numerically integration derivatives

import numpy as np
from os.path import basename
import os
import logging
import argparse
from functools import partial
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from scipy.integrate import simps
    SCIPY = True
except ImportError as e:
    SCIPY = False
    print "An error occurred while importing scipy.integrate, integration with Simpson's method will be disabled: {0}".format(e)

try:
    from helpers.plot_data import simple_plot, create_figure,\
        add_axis_to_figure, finalise_figure, generate_histogram, plot_with_kde
    CAN_PLOT = True
except ImportError as e:
    CAN_PLOT = False
    print "An error occurred while importing helpers.plot_data, plotting will be disabled: {0}".format(e)

try:
    from integration_errors.calculate_error import trapz_integrate_with_uncertainty
    CAN_USE_TRAPZ_INTEGRATION_ERRORS = True
except ImportError as e:
    CAN_USE_TRAPZ_INTEGRATION_ERRORS = False

INTEGRATION_ERROR_REFERENCE_OPTIONS = ("left", "right")
INTEGRATION_METHODS = ("simpsons", "trapezoidal")

# Recommended parameters for block averaging limits as described in Kaestner and Thiel (2006).
# Note that these are guides only! The block size should be set to the correlation time of
# the slowest fluctuation in the position data. Block averaged error estimation can provide
# a reasonable estimate for this provided the slowest fluctuation has been well samples.
N_BLOCKS_LOWER_LIMIT = 24 # segments
MINIMUM_BLOCK_SIZE = 100 # frames

logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s]: %(message)s')

def run_umbrella_integration(input_files,
    temperature,
    bin_width=None,
    number_bins=None,
    minimum_maximum_value=None,
    plot_derivatives=None,
    plot_pmf=None,
    integration_error_reference=INTEGRATION_ERROR_REFERENCE_OPTIONS[0],
    output_pmf_file=None,
    position_histogram_plot=None,
    derivatives_file=False,
    n_blocks=N_BLOCKS_LOWER_LIMIT,
    integration_method="trapezoidal",
    integration_error_analysis_method="Kaestner",
    ):

    if not input_files:
        logging.error("No input files found!")
        return
    if integration_method == "simpsons" and not SCIPY:
        logging.error("Cannot integration with Simpson's method due to import error. Please install scipy for this feature. Will continue with Trapezoidal method.")
        integration_method = "trapezoidal"
    input_data = parse_input_files(input_files)
    # check distances between windows relative to the force constant used
    window_separation_distance_check(input_data, temperature)

    if position_histogram_plot:
        if CAN_PLOT:
            generate_window_distributions(input_data, position_histogram_plot)
        else:
            logging.error('Cannot plot histogram due to import error. Please check that you have matplotlib installed.')
    calculate_window_statistics(input_data, n_blocks)

    # bin the derivatives (ie combine derivatives from all windows for calculation of final pmf by integration)
    bin_centers = get_bins(sorted(input_data.keys()), bin_width, number_bins, minimum_maximum_value)
    logging.info("N evaluation points: {0}; evaluation point range: {1}-{2}".format(len(bin_centers), bin_centers[0], bin_centers[-1]))
    window_functions = get_window_functions(input_data, temperature)
    bin_derivatives, bin_var = calculate_bin_derivatives(window_functions, bin_centers)

    # check whether force constant used for harmonic umbrella constraint is sufficient
    maximum_second_derivative_check(bin_derivatives, bin_centers, input_data)

    if plot_derivatives:
        if CAN_PLOT:
            simple_plot(bin_centers, bin_derivatives, np.sqrt(bin_var), r"$\xi$", r"$dA/d\xi\ (\mathrm{kJ} \, \mathrm{mol}^{-1})$", fig_name=plot_derivatives, show_auc=False)
        else:
            logging.error('Cannot plot derivatives due to import error. Please check that you have matplotlib installed.')
    if derivatives_file:
        write_derivatives_to_file(bin_centers, bin_derivatives, bin_var, derivatives_file)

    # integrate derivatives to calculate final PMF, incl. propagation of errors and plot data
    reaction_coordinate_positions, integral = integrate_derivatives(bin_centers, bin_derivatives, integration_method)
    mu_sigma_windows = np.mean([window_data["sigma_xi_b_i"] for window_data in input_data.values()])
    integral_point_var = integral_point_variance(bin_centers, bin_derivatives, bin_var, mu_sigma_windows,
        integration_method, integration_error_reference=integration_error_reference, integration_error_analysis_method=integration_error_analysis_method)

    # shift values such that lowest value is at 0
    integral = integral - np.min(integral)
    if plot_pmf:
        if CAN_PLOT:
            simple_plot(reaction_coordinate_positions, integral, np.sqrt(integral_point_var), r"$\xi$", r"$A\ (\mathrm{kJ} \, \mathrm{mol}^{-1})$", fig_name=plot_pmf)
        else:
            logging.error('Cannot plot PMF due to import error. Please check that you have matplotlib installed.')

    if output_pmf_file:
        outfile_pmf(reaction_coordinate_positions, integral, np.sqrt(integral_point_var), output_pmf_file)

def integrate_derivatives(bin_centers, bin_derivatives, integration_method):

    if integration_method == "simpsons":
        xs, integral = zip(*[(bin_centers[i+1], simps(bin_derivatives[:i + 2], bin_centers[:i + 2])) for i in range(len(bin_derivatives) - 2)])
    elif integration_method == "trapezoidal":
        xs, integral = zip(*[(np.mean([bin_centers[i], bin_centers[i+1]]), np.trapz(bin_derivatives[:i + 1], bin_centers[:i + 1])) for i in range(len(bin_derivatives) - 1)])
    else:
        raise Exception("Unrecognised integration method: {0}, please use: {1}".format(integration_method, " or ".join(INTEGRATION_METHODS)))
    return xs, integral

# -------- ROUTINES FOR CALCULATING PROBABILITIES AND DERIVATIVES -----------
def get_window_functions(input_data, temperature):
    RT = 8.31451*temperature/1000.0 #kJ/mol
    BETA = 1.0/RT
    # Eq 5.Kaestner and Thiel 2005
    def probability_b_i(sigma_xi_b_i, mu_xi_b_i, xi):
        return (1.0/(sigma_xi_b_i*np.sqrt(2.0*np.pi)))*np.exp(-0.5*(((xi-mu_xi_b_i)/sigma_xi_b_i)**2))

    # calculate variance of derivative, given in Eq. 5 in Kaestner and Thiel 2006
    def var_dA_b_u_over_d_xi(var_mu_segmented, var_var_segmented, sigma_xi_b_i, mu_xi_b_i, xi):
        return (1./(BETA**2*sigma_xi_b_i**4))*(var_mu_segmented + ((xi - mu_xi_b_i)**2/sigma_xi_b_i**4)*var_var_segmented)

    # Eq 6. in Kaestner and Thiel 2005
    def dA_b_u_over_d_xi(xi_i, sigma_xi_b_i, mu_xi_b_i, k, xi):
        return (1.0/BETA)*(xi-mu_xi_b_i)/sigma_xi_b_i**2 - k*(xi-xi_i)

    pmf_window_functions = {}
    for eq_pos, umbrella_data in sorted(input_data.items()):
        pmf_window_functions[eq_pos] = {
            "p_b_i":partial(probability_b_i, umbrella_data["sigma_xi_b_i"], umbrella_data["mu_xi_b_i"]), #set window parameters for Eq. 5 Kaestner and Thiel 2005
            "dA/dxi":partial(dA_b_u_over_d_xi, eq_pos, umbrella_data["sigma_xi_b_i"], umbrella_data["mu_xi_b_i"], umbrella_data["k"]), #set window parameters for Eq. 6 Kaestner and Thiel 2005
            "var(dA/dxi)":partial(var_dA_b_u_over_d_xi, umbrella_data["var_mu_segmented"], umbrella_data["var_var_segmented"], umbrella_data["sigma_xi_b_i"], umbrella_data["mu_xi_b_i"]), #set window parameters for Eq. 5 in Kaestner and Thiel 2006
            }
        pmf_window_functions[eq_pos]["N"] = umbrella_data["com"].shape[0]
    return pmf_window_functions

# Calculation of the core umbrella integration equations written in a functional style.
# This make for simpler relation to equations in the paper but comes at a slight performance cost.
def calculate_bin_derivatives(window_functions, bin_centers, calculate_var=True):

    # Eq. 7 Kaestner and Thiel 2005
    def dA_dxi_bin(xi_bin, p_i_list, dA_dxi_list):
        return np.sum([p_i(xi_bin)*dA_dxi(xi_bin) for p_i, dA_dxi in zip(p_i_list, dA_dxi_list)])

    # Eq. 9 Kaestner and Thiel 2006
    def var_dA_dxi_bin(xi_bin, p_i_list, var_dA_dxi_list):
        return np.sum([p_i(xi_bin)**2*var_dA_dxi(xi_bin) for p_i, var_dA_dxi in zip(p_i_list, var_dA_dxi_list)])

    # Eq. 8 Kaestner and Thiel 2005
    def probability_i(N_i, p_b_i, norm_factor, xi):
        return N_i*p_b_i(xi)/float(norm_factor)

    # Eq. 8 total_norm_factor Kaestner and Thiel 2005
    def calc_total_norm_factor(N_i_list, p_b_i_list, xi):
        return np.sum([N_i*p_b_i(xi) for N_i, p_b_i in zip(N_i_list, p_b_i_list)])

    p_b_i_list, dA_dxi_list, var_dA_dxi_list, N_i_list = zip(*[(funcs["p_b_i"], funcs["dA/dxi"], funcs["var(dA/dxi)"], funcs["N"]) for _, funcs in sorted(window_functions.items())])
    derivatives = np.zeros(len(bin_centers))
    var_derivatives = np.zeros(len(bin_centers))
    for i, xi_bin in enumerate(bin_centers):
        total_norm_factor = calc_total_norm_factor(N_i_list, p_b_i_list, xi_bin)
        assert total_norm_factor != 0, "Nomalisation factor is equal to zero, this will produce NaN derivatives. "\
        "This occures when the umbrella positions do not appropriately cover the reaction coordinate range that has been specified in this calculation."
        p_i_list = [partial(probability_i, N_i, p_b_i, total_norm_factor) for p_b_i, N_i in zip(p_b_i_list, N_i_list)]
        derivatives[i] = dA_dxi_bin(xi_bin, p_i_list, dA_dxi_list)
        if calculate_var:
            var_derivatives[i] = var_dA_dxi_bin(xi_bin, p_i_list, var_dA_dxi_list)
    return derivatives, var_derivatives

# A declarative implementation of the core umbrella integration equations.
# Slightly more efficient than the functional implementation.
def calculate_bin_derivatives_faster(input_data, bin_centers, temperature):
    RT = 8.31451*temperature/1000.0 #kJ/mol
    BETA = 1.0/RT
    sqrt_2_pi = np.sqrt(2.0*np.pi)
    window_eq_xi = sorted(input_data.keys())
    window_means, window_sigmas, force_constants, Ns = zip(*[(np.mean(window_data["com"]), np.std(window_data["com"]), window_data["k"], len(window_data["com"])) for _, window_data in sorted(input_data.items())])

    n_bins = len(bin_centers)
    da_dxi = np.zeros(n_bins)
    for i, xi_bin in enumerate(bin_centers):
        total_weight = 0
        for xi_i, mu_xi_b_i, sigma_xi_b_i, k, N in zip(window_eq_xi, window_means, window_sigmas, force_constants, Ns):
            xi_bin_minus_mu_xi_b_i = xi_bin-mu_xi_b_i
            var_xi_b_i = sigma_xi_b_i**2
            weight = N*np.exp(-0.5*(xi_bin_minus_mu_xi_b_i**2/var_xi_b_i))/(sigma_xi_b_i*sqrt_2_pi)
            da_dxi[i] += weight*(xi_bin_minus_mu_xi_b_i/(var_xi_b_i*BETA) - k*(xi_bin-xi_i))
            total_weight += weight
        da_dxi[i] = da_dxi[i]/float(total_weight)
    return da_dxi

# get_bins
# routine to combine the derivatives from the different windows
# to a single data set (based on their corresponding com-distances)
def get_bins(window_centers, bin_width, number_bins, minimum_maximum_value, exclude_end_points=False):

    # get evenly spaced bins between min and max bin boundaries using user provided binwidth
    if bin_width is not None and number_bins is not None:
        raise Exception("Only one of the following arguments can be provided: bin_width, number_bins")
    first_bin_center = window_centers[0] if not minimum_maximum_value else minimum_maximum_value[0]
    last_bin_center = window_centers[-1] if not minimum_maximum_value else minimum_maximum_value[1]
    if exclude_end_points:
        bin_width = (last_bin_center - first_bin_center)/float(number_bins)
        start = first_bin_center + 0.5*bin_width
        end = last_bin_center - 0.5*bin_width
    else:
        start = first_bin_center
        end = last_bin_center
    if bin_width is not None:
        bin_centers = np.arange(start, end + bin_width, bin_width)
    if number_bins is not None:
        bin_centers = np.linspace(start, end, number_bins)
    return bin_centers

# -------- ROUTINES FOR ERROR ANALYSIS -----------

# routine to calculate the error in the integration of derivatives when calculating the final PMF 
def integral_point_variance(bin_centers, bin_derivatives, bin_var, mu_sigma_windows,
    integration_method, integration_error_reference="left", integration_error_analysis_method="Kaestner"):

    n_remove = 1 if integration_method == "trapezoidal" else 2 # 2 for simpsons
    interval_indexes = range(len(bin_centers)-n_remove)

    if integration_error_analysis_method == "trapz_analysis" and integration_method == "trapezoidal" and CAN_USE_TRAPZ_INTEGRATION_ERRORS:
        if integration_error_reference == "right":
            return [
                trapz_integration_error_analysis(bin_centers[-i-n_remove:], bin_derivatives[-i-n_remove:], bin_var[-i-n_remove:])
                for i in interval_indexes
                ][::-1]
        else:
            return [
                trapz_integration_error_analysis(bin_centers[:i+n_remove], bin_derivatives[:i+n_remove], bin_var[:i+n_remove])
                for i in interval_indexes
                ]
    elif integration_error_analysis_method == "Kaestner":
        if integration_error_reference == "right":
            return [var_delta_A(bin_var[-i-n_remove:], bin_centers[-i-n_remove:], mu_sigma_windows) for i in interval_indexes][::-1]
        else:
            return [var_delta_A(bin_var[:i+n_remove], bin_centers[:i+n_remove], mu_sigma_windows) for i in interval_indexes]
    else:
        if integration_error_analysis_method == "trapz_analysis":
            raise Exception("To use trapezoidal integration error analysis: (1) select trapezoidal integration method, " \
            "(2) ensure the integration_errors (<github address here>) module is in your python path.")
        else:
            raise Exception("Unknown integration method: {0}".format(integration_error_analysis_method))

def trapz_integration_error_analysis(xs, ys, es):
    if len(xs) == 1:
        return 0.0
    return trapz_integrate_with_uncertainty(xs, ys, np.sqrt(es))[1]

# Eq. 15 Kaestner and Thiel 2006
def var_delta_A(var_derivatives, xis, mu_sigma_windows):
    if xis[-1] == xis[0]:
        return 0.0
    mu_var_derivatives = np.mean(var_derivatives)
    delta_xi = abs(xis[-1] - xis[0])
    return mu_var_derivatives*(delta_xi*mu_sigma_windows*np.sqrt(2*np.pi) - 2.0*mu_sigma_windows**2)

def calculate_window_statistics(input_data, n_blocks):
    for _, window_data in sorted(input_data.items()):
        block_size = len(window_data["com"])/n_blocks
        if block_size < MINIMUM_BLOCK_SIZE:
            logging.warning("The number of frames in each block {0} is below the recommended minimum of {1}".format(block_size, MINIMUM_BLOCK_SIZE))
        # calculate mean and standard deviation for segments of window (based on block size from convergence analysis)
        mu_segmented, var_segmented = get_segmented_window_statistics(window_data["com"], block_size)
        window_data["mu_xi_b_i"] = np.mean(window_data["com"])
        window_data["sigma_xi_b_i"] = np.std(window_data["com"])
        # block averaged estimate of the variance of the mean
        window_data['var_mu_segmented'] = np.var(mu_segmented, ddof=1)/float(mu_segmented.shape[0])
        # block averaged estimate of the variance of the variance
        window_data['var_var_segmented'] = np.var(var_segmented, ddof=1)/float(mu_segmented.shape[0])

def get_segmented_window_statistics(com_distances, block_size):
    block_boundary_indexes = range(0, len(com_distances), block_size)
    block_means = np.zeros(len(block_boundary_indexes))
    block_vars = np.zeros(len(block_boundary_indexes))
    for i, block_index in enumerate(block_boundary_indexes):
        block_data = com_distances[block_index:block_index+block_size+1]
        block_means[i] = np.mean(block_data)
        block_vars[i] = np.var(block_data, ddof=1)
    return np.array(block_means), np.array(block_vars)

def maximum_second_derivative_check(derivatives, bin_centers, input_data):
    caveat = "Johannes Kaestner and Walter Thiel 2006 (DOI: 10.1063/1.2206775). NOTE: this implementation uses numerical second derivatives to estimate kappa which can be very sensitive to "\
        "noise in the data and thus the result should be treated with caution."
    ddf_ddx = np.diff(derivatives)/np.diff(bin_centers)
    kappa, bin_center = sorted(zip(ddf_ddx, bin_centers[:len(ddf_ddx)]))[0]
    # take the nearest umbrella to the bin_center to get the force constant corresponding to that second derivative
    umbrella_eq_pos = input_data.keys()
    nearest_eq_pos = umbrella_eq_pos[(np.abs(np.array(umbrella_eq_pos)-bin_center)).argmin()]
    fc = input_data[nearest_eq_pos]["k"]
    kappa = -kappa
    if fc < 6*kappa:
        logging.info("Lowest second derivative was: {0:.3g} at {1}".format(-kappa, bin_center))
        logging.warning("The harmonic force constant used ({0:.0g} kJ mol-1 distance-2) is less than the recommended cutoff of 6*kappa={1:.0g} kJ mol-1 distance-2: {2}".format(fc, 6*kappa, caveat))
    else:
        logging.info("Lowest second derivative was: {0:.3g} at {1} distance (recommended harmonic force constant ~{2:.0g}): {3}".format(-kappa, bin_center, 6*kappa, caveat))

# routine to check whether the distances between the windows are
# sufficiently small, as given by eq. 18 in Kaestner and Thiel 2006
def window_separation_distance_check(data, temperature):
    RT = 8.31451*temperature/1000.0 #kJ/mol
    BETA = 1.0/RT
    window_centers =  sorted(data.keys())
    displayed_messages = []
    for i in range(len(window_centers)-1):
        window_delta = window_centers[i+1] - window_centers[i]
        # maximum force constant used for a given window separation distance
        max_fc = np.max([np.max(data[window_centers[i]]["k"]), np.max(data[window_centers[i+1]]["k"])])
        window_size_cutoff = 3./np.sqrt(BETA*max_fc)
        if window_delta > window_size_cutoff:
            warn_msg = "Distance between umbrellas {0} is greater than the recommended cutoff of 3/sqrt(BETA*K_fc) = {1:.3g}: "\
                "Johannes Kaestner and Walter Thiel 2006 (DOI: 10.1063/1.2206775)".format(window_delta, window_size_cutoff)
            if warn_msg not in displayed_messages:
                logging.warning(warn_msg)
                logging.debug("Foce constant was: {0:.3g} kJ mol^-1 distance^-2".format(max_fc))
                logging.debug("Umbrella equilibrium points used: {0} - {1}".format(window_centers[i], window_centers[i+1]))
                displayed_messages.append(warn_msg)
        else:
            info_msg = "Distance between umbrellas {0} is within the recommended cutoff of 3/sqrt(BETA*K_fc) = {1:.3g}: "\
                "Johannes Kaestner and Walter Thiel 2006 (DOI: 10.1063/1.2206775)".format(window_delta, window_size_cutoff)
            if info_msg not in displayed_messages:
                logging.info(info_msg)
                logging.debug("\tFoce constant was: {0:.3g} kJ mol^-1 distance^-2".format(max_fc))
                logging.debug("\tUmbrella equilibrium points used: {0} - {1}".format(window_centers[i], window_centers[i+1]))
                displayed_messages.append(info_msg)

# -------- INPUT / OUTPUT ROUTINES-----------
def parse_input_files(input_files):
    # read data from input file
    # format of input file
    #col 1 = time (that corresponds to the instantaneous value of the reaction coordinate in col 4)
    #col 2 = instantaneous value of reaction coordinate at time t in col 3, produced with umbrella a center col 2 and
    #col 3 = center of the umbrella potential used in the simulation for that window (reference position)
    #col 4 = force constant (fc) used in the simulation for that window (in kJ mol-1 distance-2)
    #     force constant col 1)

    data = {}
    for file_name in input_files:
        logging.info("reading: {0}".format(file_name))
        with open(file_name) as fh:
            lines = fh.read().splitlines()
        for l in lines:
            if l.strip()[0] == "#":
                continue
            t, com, eq_pos, fc = map(float, l.split()[:4])
            data.setdefault(eq_pos, {"time":[], "com":[], "fc":[]})
            data[eq_pos]["time"].append(t)
            data[eq_pos]["com"].append(com)
            data[eq_pos]["fc"].append(fc)
    # sort data on time
    for eq_pos, umbrella_data in data.items():
        t, com, fc = umbrella_data["time"], umbrella_data["com"], umbrella_data["fc"]
        t, com, fc = zip(*sorted(zip(t, com, fc)))
        assert np.unique(fc).shape[0] == 1, "All force constants for this window were not equal: {0}".format(eq_pos)
        data[eq_pos] = {"time":np.array(t), "com":np.array(com), "k":fc[0]}
    return data

def generate_window_distributions(data, position_histogram_plot):
    create_dir(position_histogram_plot)
    fig_distributions = create_figure(figsize=(6,4))
    ax_distributions = add_axis_to_figure(fig_distributions)
    for _, window_data in sorted(data.items()):
        centers, his, kde = generate_histogram(window_data["com"])
        plot_with_kde(ax_distributions, centers, his, kde=kde)

    fig_name="{0}.png".format(position_histogram_plot) if not "." in position_histogram_plot else position_histogram_plot
    finalise_figure(fig_distributions, ax_distributions, xlabel=r"$\xi$", ylabel="occurrence", fig_name=fig_name)

def create_dir(filename):
    if os.path.dirname(filename) and not os.path.exists(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename))

# routine to write derivatives for each bin centre to a file 
def write_derivatives_to_file(bin_centers, bin_derivatives, bin_vars, output_file):
    derivatives = "\n".join(["{0:>10f} {1:>10f} {2:>10f}".format(bin_center, derivative, bin_var) for bin_center, derivative, bin_var in zip(bin_centers, bin_derivatives, bin_vars)])
    with open(output_file, "w") as fh:
        fh.write(derivatives)

# routine to write final pmf data to output file  
def outfile_pmf(bin_centers, integral, integral_point_std, pmf_output_file):
    with open(pmf_output_file, "w") as outfile:
        outfile.write("\n".join(["# reaction_coordinate  energy(kJ/mol)  error(kJ/mol)"] + ["{0:5f} {1:5f} {2:5f}".format(bc, v, e) for bc, v, e in zip(bin_centers, integral, integral_point_std)]))

def parse_commandline():

    parser = argparse.ArgumentParser(prog=basename(__file__))
    parser.add_argument("-i", "--input_files",
        help="List of input data files (e.g. data/*.ui_dat); input file format: col1=time, col2=position, col3=umbrella potential centre, col4=harmonic force constant (kJ mol^-1 distance^-2).",
        action='store', nargs='+', required=True)
    parser.add_argument("-t", "--temperature",
        help="Temperature at which the analysis is to be performed.",
        action='store', required=True, type=float)
    bin_group = parser.add_mutually_exclusive_group(required=True)
    bin_group.add_argument("-b", "--bin_width",
        help="Bin width along the reaction coordinate (independent of umbrellas used in the sampling).", action="store", type=float, required=False)
    bin_group.add_argument("-n", "--number_bins",
        help="Number of evaluation points to be used along the reaction coordinate (independent of umbrellas used in the sampling).", action="store", type=int, required=False)
    parser.add_argument("-m", "--minimum_maximum_value",
        help="Minimum and maximum reaction coordinate positions to be considered e.g. -m 0.2 5.5", action="store", type=float, required=False, nargs=2)
    parser.add_argument("-pd", "--plot_derivatives",
        help="Plot the derivatives at each evaluation point along the reaction coordinate. File name ending will determine image format; default format is png.",
         action="store", required=False, default=None)
    parser.add_argument("-pp", "--plot_pmf",
        help="Plot the PMF. File name ending will determine image format; default format is png.",
         action="store", required=False, default=None)
    parser.add_argument("-r", "--integration_error_reference",
        help="Integration error estimation reference point: ({0})".format("|".join(INTEGRATION_ERROR_REFERENCE_OPTIONS)), action="store", required=False)
    parser.add_argument("-ph", "--position_histogram_plot",
        help="Position histogram plot. File name ending will determine image format; default format is png.", action="store", required=False)
    parser.add_argument("-o", "--output_pmf_file",
        help="PMF output data file.", action="store", required=False)
    parser.add_argument("-d", "--derivatives_file",
        help="Write derivative at each evalutation position.", action="store", required=False, default=False)
    parser.add_argument("-nb", "--n_blocks",
        help="Number of blocks used in error analysis to calculate block averaged error estimate for each window.", action="store", type=int, required=False, default=N_BLOCKS_LOWER_LIMIT)
    parser.add_argument("-im", "--integration_method",
        help="Integration method used to numerically integration derivatives.", action="store", required=False, default="simpsons")

    args = parser.parse_args()

    return (args.input_files,
        args.temperature,
        args.bin_width,
        args.number_bins,
        args.minimum_maximum_value,
        args.plot_derivatives,
        args.plot_pmf,
        args.integration_error_reference,
        args.output_pmf_file,
        args.position_histogram_plot,
        args.derivatives_file,
        args.n_blocks,
        args.integration_method,
        )

if __name__=="__main__":
    (input_files,
    temperature,
    bin_width,
    number_bins,
    minimum_maximum_value,
    plot_derivatives,
    plot_pmf,
    integration_error_reference,
    output_pmf_file,
    position_histogram_plot,
    derivatives_file,
    n_blocks,
    integration_method) = parse_commandline()

    run_umbrella_integration(input_files,
        temperature,
        bin_width,
        number_bins,
        minimum_maximum_value,
        integration_error_reference=integration_error_reference,
        plot_derivatives=plot_derivatives,
        plot_pmf=plot_pmf,
        output_pmf_file=output_pmf_file,
        position_histogram_plot=position_histogram_plot,
        derivatives_file=derivatives_file,
        n_blocks=n_blocks,
        integration_method=integration_method)
