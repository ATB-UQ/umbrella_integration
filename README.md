-----------------------------------
Umbrella integration
-----------------------------------

Python implementation of the umbrella integration method for potential of mean force (PMF) calculations

Authors: Martin Stroet (University of Queensland),
         Evelyn Deplazes (Curtin University, University of Queensland)

Please also cite the GitHub project using the following DOI:

    DOI: 10.5281/zenodo.164996

The calculation of the derivatives, the free energy profile (potential of mean
force, PMF) is based on the equations in the following paper by Kästner & Thiel:

    Kästner J. and Thiel W., J. Chem. Phys. 123, 144104 (2005)
    doi: 10.1063/1.2052648

The estimates of uncertainties for the PMF are based on equations in
the equations in the following paper by Kästner & Thiel:

    Kästner J. and Thiel W., J. Chem. Phys. 124, 234106 (2006)
    doi: 10.1063/1.2206775


--------------------
About the program
--------------------

umbrella_integration implements the umbrella integration method for the calculation of a potential of mean force (PMF) as described by Kästner and Thiel (2005, J. Chem. Phys.). Also implemented is the error estimation method described in Kästner and Thiel (2006, J. Chem. Phys.). The program is independent of the MD package used for the umbrella sampling simulations. All the information required as input can easily be extracted and/or calculated from the trajectories or other simulation output files.

The project contains the following two modules:
- umbrella_integration.py: calculates the PMF along with error estimates
- format_input_data.py: a tool to format umbrella simulation data to be compatible with umbrella_integration.py

** IMPORTANT **
This project doesn't aim to be a PMF black-box. Some degree of Python programming will most likely be required.

# Reaction coordinate
No assumptions are made about the nature of the reaction coordinate, only that it must be one-dimensional and that the restraining potential is harmonic.

#### One-dimensional only
The program only works for a one-dimensional PMF i.e. a single reaction
coordinate. The following paper by Kästner describes an extension of umbrella
integration for two or more reaction coordinates:
Kästner J., J. Chem. Phys. 131, 034109 (2009) http://dx.doi.org/10.1063/1.3175798

#### Restraining potential
The program currently only handles harmonic bias potentials. However extending the program to handle
other restraining potentials should be relatively simple.

# Energy units
Currently only kJ/mol is supported, however it will not be difficult to extend this if required.

---------------------
Requirements
----------------------

    Python 2.7 or Python 2.6
    numpy

---------------------
Optional
----------------------

    matplotlib (show plots)
    scipy (Simpson's integration method, plot histogram kde's)

--------------------
Running the program
--------------------
E.g. python umbrella_integration.py -i input_file -t 298.15
Other input parameters include:

     -i input_file (required)
        Path(s) to the input file(s).
        The input file needs to contain the following data in white-space separated column format:
        col 1 = time (while not used in the PMF calculation directly it is useful when assessing
             the convergence of umbrella simulations)
        col 2 = instantaneous value of reaction coordinate
        col 2 = centre (equilibrium position) of the harmonic potential used in the umbrella simulation
        col 1 = harmonic force constant used for the umbrella simulation in kJ mol^-1 distance^-2
     -t temperature
        Temperature at which the calculation is to be performed.
     -b bin_width (required if -n number_bins not provided)
        The step size between point where the potential will be calculated along the reaction coordinate.
        This is independent of the spacing between umbrellas/windows used
        in the umbrella sampling simulations.
     -n number_bins (required if -b bin_width not provided)
        Number of positions along the reaction coordinate at which the potential will be calculated.
        This is independent of the umbrellas/windows used in umbrella sampling simulations.
     -m min_max_value
        Minimum and maximum reaction coordinate positions at which the potential will be calculated.
        e.g. -m 0.2 5.5
        Note that data beyond this range will be included in the calculation.
     -p show_plots
        Interactively show plots of the derivative (which are numerically integrated) and final PMF
     -r integration_error_reference
        Integration error estimation reference point, left or right
     -im integration_method
        Numerical method used to the integrate derivatives. Simpon's method has a lower error
        for smooth data, while the trapezoidal method is more robust on noisy data.
     -o output_pmf_file
        Name of output file that contains the final PMF.
        Format, col 1 = reaction coordinate, col 2 = energy, col 3 = uncertainty
     -d derivatives_file
        Name of output file that contains the derivatives for each reaction coordinate point.
     -ph position_histrogram_plot
        File name for the histogram plot.
     -np n_blocks
        Number of blocks used in error analysis to calculate block averaged error estimate for each window.


----------------------
Input file format
----------------------

The input file(s) needs to contain the following data in the format:

    col 1 = time (while not used in the PMF calculation directly it is useful when assessing the convergence of umbrella simulations)
    col 2 = instantaneous value of reaction coordinate
    col 3 = centre (equilibrium position) of the harmonic potential used in the umbrella simulation
    col 4 = harmonic force constant used for the umbrella simulation in kJ mol^-1 distance^-2

For an example see the input file merged_ui_input.dat in the example/data folder.

----------------------------
Generating input files
----------------------------

The script generate_input_data.py contains a number of routines as well as a commandline interface to format umbrella simulation data to be compatible with umbrella_integration.py. While not required, data from multiple umbrella simulations can be combined into a single file. The translation of relative reaction coordinate positions to absolute one is also included.

----------------------------
Example
----------------------------

A comprehensive example (from raw data to PMF) illustrating the programs use can be found in the ./example folder. The example includes both bash and python scripts to run the same calculation (albeit with slightly different input parameters).
