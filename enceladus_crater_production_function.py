#!/usr/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt

"""
Crater Production Function (CPF) Implementation and Plotting 

This script provides a reproducible implementation of the Crater Production Function (CPF) of Enceladus, published in Wong et al (updated years at the end). It models the expected number of craters on Enceladus as a function of impactors' size, impactor velocities and other parameters.

Key features:
- Implements the polynomial function of the CPF.
- Provides functions to compute crater densities for arbitrary impactor diameters.
- Includes comparisons with other published impactor and crater Size-Frequency Distributions (SFDs).

Users can reproduce Figure 1 from Wong et al (year) and adapt the CPF to updated satellite-specific input parameters, such as surface gravity, impact velocities, and impact probability.

Definitions:
- CPF: Crater Production Function, a mathematical model describing the number of craters formed over time.
- SFD: Size-Frequency Distribution, a statistical representation of the number of objects (impactors or craters) as a function of their size.

Usage:
To run the script and reproduce results, execute:
    python cpf.py [crater_diameters_in_km]

You may modify input parameters within the script to adapt the CPF for different icy satellites or scenarios.
"""


def usage():
    print("Usage: python cpf.py [crater_diameter_in_km [crater_diameter_in_km ...]]")
    print("Example: python cpf.py 1.0 5.0 10.0")
    print("- Allowed crater diameter range: 0.8 km to 44 km")
    print(
        "- If no arguments are provided, 1 km will be used as the default crater diameters to calculate the expected crater density."
    )


# Evaluate the normalized CPF
def generate_CPF():
    """
    Generate a normalised crater production function (CPF) as a 10th-degree polynomial.
    Coefficients are in increasing order (a0, a1, ..., a10).
    """
    # === 1) Coefficient of the 10-degree polynomial crater production function ===
    cpf_coefficients = np.array(
        [
            -0.065237099,  # a0
            -2.35027511,  # a1
            -5.55729243,  # a2
            24.95085739,  # a3
            -77.26792659,  # a4
            212.99933586,  # a5
            -406.68327268,  # a6
            464.57706063,  # a7
            -303.36628993,  # a8
            104.74296626,  # a9
            -14.84990662,  # a10
        ]
    )
    # np.poly1d expects coefficients from highest to lowest degree, so reverse the input
    crater_production_function = np.poly1d(cpf_coefficients[::-1])
    return crater_production_function


def wong_crater_sfd(crater_diameter):
    """
    Calculate the expected cumulative crater density (crater/km^2) for given crater diameter(s) in
    km
    e.g., wong_crater_sfd([1, 2, 3])
    """
    # Our crater production function (CPF) for Enceladus
    ##  Create the normalised CPF function within the crater diameter ranges
    crater_production_function = generate_CPF()
    return 10 ** (crater_production_function(np.log10(crater_diameter)))


def define_impact_parameter_enceladus():
    return {
        # Enceladus' surface gravity in m/s2
        "surface_gravity": 0.1134,
        # Impactor's density to satellite's density ratio
        #   Bulk density of impactor ~ 400 kg/m3, satellite's icy crust ~ 1000 kg/m3
        "density_ratio": 0.4,
        # Average impact velocity in km/s from Genga N-body simulations (Wong et al 2021, 2023)
        "average_impact_velocity": 22.174,  # 21.7
        # Simple-to-complex crater transition diameter in km (Schenk 1991)
        "transition_diameter": 15,
        # average impact probability of Enceladus (Wong et al 2023)
        "impact_probability": 2.4854e-8,  # 1.355e-8
        # global spherical area of Enceladus in km2
        "spherical_area": 798_648,
    }


def define_impact_parameter_tethys():
    return {
        # Enceladus' surface gravity in m/s2
        "surface_gravity": 0.1462,
        # Impactor's density to satellite's density ratio
        #   Bulk density of impactor ~ 400 kg/m3, satellite's icy crust ~ 1000 kg/m3
        "density_ratio": 0.4,
        # Average impact velocity in km/s from Genga N-body simulations (Wong et al 2021, 2023)
        "average_impact_velocity": 19.997,
        # Simple-to-complex crater transition diameter in km (Schenk 1991)
        "transition_diameter": 15,
        # average impact probability of Enceladus (Wong et al 2023)
        "impact_probability": 4.74e-8,  # TODO check for updated value
        # global spherical area of Enceladus in km2
        "spherical_area": 3.543e6,
    }


def unpack_moon(moon_impact_parameters):
    return (
        moon_impact_parameters["surface_gravity"],
        moon_impact_parameters["density_ratio"],
        moon_impact_parameters["average_impact_velocity"],
        moon_impact_parameters["transition_diameter"],
        moon_impact_parameters["impact_probability"],
        moon_impact_parameters["spherical_area"],
    )


# === 2) Crater scaling law (Zahnle et al 2003) ===
def impactor_to_crater(impactor_diameter, moon_impact_parameters):
    """
    Convert impactor diameter (km) to resulting crater diameter (km) on Enceladus.
        impactor_to_crater(impactor_diameter) = crater_diameter
    """
    surface_gravity, density_ratio, average_impact_velocity, transition_diameter = (
        unpack_moon(moon_impact_parameters)[:4]
    )
    crater_diameter = (
        impactor_diameter
        * (0.1431**-1)
        * (surface_gravity**-0.282)
        * (density_ratio**0.427)
        * (average_impact_velocity**0.564)
        * (transition_diameter**-0.192)
    ) ** (1 / 1.0897)

    return crater_diameter


def crater_to_impactor(crater_diameter, moon_impact_parameters):
    """
    Convert crater diameter (km) on Enceluads to typical impactor diameter (km), i.e., heliocentric
    comet.
        crater_to_impactor(crater_diameter) = impactor_diameter
    """
    surface_gravity, density_ratio, average_impact_velocity, transition_diameter = (
        unpack_moon(moon_impact_parameters)[:4]
    )
    impactor_diameter = (
        0.1431
        * (crater_diameter**1.0897)
        * (surface_gravity**0.282)
        * (density_ratio**-0.427)
        * (average_impact_velocity**-0.564)
        * (transition_diameter**0.192)
    )
    return impactor_diameter


# === 3) Calculate the number of impactor in the outer Solar System ===


# Singer et al. (2019) size-frequency distribution of impactors
#   Studied the Pluto-Charon crater population
def singer_impactor_sfd(impactor_diameter):
    return np.where(
        impactor_diameter >= 1.0,
        (1.0 / impactor_diameter) ** 2.0,
        (1.0 / impactor_diameter) ** 0.7,
    )


# Cumulative number of Scattered disk object larger than certain diameter
def number_of_impactor(impactor_diameter, impactor_sfd):
    """
    Calculate the cumulative number of scattered disk objects larger than a given diameter,
    normalized to the reference size (10 km) using Singer's size-frequency distribution.

    Parameters:
        impactor_diameter (float): Diameter of the impactor in km.

    Returns:
        float: Scaled number of impactors larger than the specified diameter.
    """
    # Known current number of outer Solar System small bodies of > 10 km (reference)
    reference_impactor_diameter = 10.0  # km
    current_number_of_10km_impactor = 2e7
    reference_impactor_count = impactor_sfd(reference_impactor_diameter)
    scaled_impactor_count = impactor_sfd(impactor_diameter)
    return (
        current_number_of_10km_impactor
        * scaled_impactor_count
        / reference_impactor_count
    )


# Expected impact density (estimated from current number of 10 km object)
def expected_impact_density(impactor_diameter, moon_impact_parameters, impactor_sfd):
    impact_probability, spherical_area = unpack_moon(moon_impact_parameters)[4:6]
    return (
        number_of_impactor(impactor_diameter, impactor_sfd)
        / spherical_area
        * impact_probability
    )


# === 4) Creating the cumulative crater density function with ...
# our crater production function and
# different impactor or crater size-frequency distribution ===


# Define the other size-frequency distributions (SFDs) for comparison
# (already included Singer et al 2019 above)


# Zahnle et al 2003 case A impactor size-frequency distribution
#   Studied the Jovican satellites crater
def zahnle_impactor_sfd(impactor_diameter):
    return np.where(
        impactor_diameter >= 5.0,
        0.129 * (impactor_diameter / 5) ** (-2.5),
        np.where(
            impactor_diameter >= 1.5,
            (impactor_diameter / 1.5) ** (-1.7),
            (impactor_diameter / 1.5) ** (-1.0),
        ),
    )


# Kirchoff and Schenk 2009 for Mid-latitude crater plains on Enceladus
def kirchoff_crater_sfd(crater_diameter):
    return np.where(
        crater_diameter >= 7.0,
        (2.3671 / crater_diameter) ** 2.99,
        np.where(
            crater_diameter >= 4.0,
            (1.7203 / crater_diameter) ** 2.310,
            (1.0 / crater_diameter) ** 1.406,
        ),
    )


def print_crater_density_table(crater_diameter_list):
    print("\nExpected crater density on Enceladus using the CPF of this work:")
    print(f"{'Crater Diameter (km)':>20} | {'Expected Density (crater/km^2)':>30}")
    print(f"{'-' * 20} | {'-' * 30}")
    for crater_diameter in crater_diameter_list:
        expected_crater_density = wong_crater_sfd(crater_diameter)
        print(f"{crater_diameter:20.2g} | {expected_crater_density:30.3g}")


def calculate_impactor_density(
    moon_crater_diameter,
    moon_impact_parameters,
    impactor_sfd,
    minimum_impactor_diameter,
    maximum_impactor_diameter,
):
    # Expected number of 1 km impactor from the impactor size-frequency distribution
    one_km_impactor = impactor_sfd(1.0)
    expected_one_km_impactor = expected_impact_density(
        1.0, moon_impact_parameters, singer_impactor_sfd
    )

    # Creating function with the impactor diameter range plotted/discussed in the cited publication
    reliable_impactor_diameter = 10 ** (
        np.linspace(
            np.log10(minimum_impactor_diameter),
            np.log10(maximum_impactor_diameter),
            1000,
        )
    )
    # Calculating the impactor density using the 'impactor_sfd'
    impactor_density = impactor_sfd(reliable_impactor_diameter)
    # Normalising the impactor density to the expected number of 1 km impactor onto Enceladus
    normalised_density = impactor_density / one_km_impactor * expected_one_km_impactor

    # Creating function with the entire observed crater ranges in this work
    all_impactor_diameter = crater_to_impactor(
        moon_crater_diameter, moon_impact_parameters
    )
    impactor_density_all_impactor_diameter = impactor_sfd(all_impactor_diameter)
    normalised_density_all_impactor_diameter = (
        impactor_density_all_impactor_diameter
        / one_km_impactor
        * expected_one_km_impactor
    )
    return (
        reliable_impactor_diameter,
        normalised_density,
        all_impactor_diameter,
        normalised_density_all_impactor_diameter,
    )


def calculate_crater_density(
    enceladus_crater_diameter,
    moon_impact_parameters,
    crater_sfd,
    minimum_enceladus_crater_diameter,
    maximum_enceladus_crater_diameter,
):
    # Gather the impact parameters of Enceladus
    enceladus_impact_parameters = define_impact_parameter_enceladus()
    # From Kirchoff and Schenk (2009) *crater* size-frequency distribution
    expected_enceladus_crater_size_from_one_km_impactor = impactor_to_crater(
        1.0, enceladus_impact_parameters
    )
    one_km_impactor_on_enceldaus = crater_sfd(
        expected_enceladus_crater_size_from_one_km_impactor
    )
    one_km_impactor_on_moon = expected_impact_density(
        1.0, moon_impact_parameters, singer_impactor_sfd
    )

    # Creating function with the crater diameter range plotted/discussed in the cited publication
    enceladus_reliable_crater_diameter = calculate_crater_diameter(
        minimum_enceladus_crater_diameter, maximum_enceladus_crater_diameter
    )
    moon_reliable_crater_diameter = impactor_to_crater(
        crater_to_impactor(
            enceladus_reliable_crater_diameter, enceladus_impact_parameters
        ),
        moon_impact_parameters,
    )
    density = crater_sfd(enceladus_reliable_crater_diameter)
    normalised_density = (
        density / one_km_impactor_on_enceldaus * one_km_impactor_on_moon
    )

    # Creating function with the entire observed crater ranges in this work
    density_all_crater_diameter = crater_sfd(enceladus_crater_diameter)
    normalised_density_all_crater_diameter = (
        density_all_crater_diameter
        / one_km_impactor_on_enceldaus
        * one_km_impactor_on_moon
    )
    return (
        moon_reliable_crater_diameter,  # list of all crater sizes
        normalised_density,  # list of how many craters of each size (within the reliable range) exist
        normalised_density_all_crater_diameter,  # list of how many craters of each size exist (extrapolated)
    )


def calculate_crater_diameter(min, max):
    # Creating function with the crater diameter range plotted/discussed in the cited publication
    return 10 ** (np.linspace(np.log10(min), np.log10(max), 1000))


def main():
    for x in sys.argv[1:]:
        try:
            input_crater_diameter = float(x)
            if input_crater_diameter < 0.8 or input_crater_diameter > 44:
                usage()
                exit()
        except ValueError:
            usage()
            exit()

    # Get the user input (list of) crater diameter(s) or use 1 km is none provided
    crater_diameter_list = [float(arg) for arg in sys.argv[1:]] or [1.0]
    # print the expected crater density for the input crater diameter(s)
    # if the user input is not provided, the expected crater density for 1 km crater on Enceladus is 0.8605238288 crater per km2 with our CPF
    print_crater_density_table(crater_diameter_list)

    # Gather the impact parameters of Enceladus
    enceladus_impact_parameters = define_impact_parameter_enceladus()
    # Gather the impact parameters of the target moon
    # moon_impact_parameters = define_impact_parameter_tethys()
    moon_impact_parameters = define_impact_parameter_enceladus()

    ##  Define the crater diameter ranges for the plot.
    ##  same as the diameter ranges observed on the surface images of Enceladus
    enceladus_crater_diameter = calculate_crater_diameter(min=0.8, max=44.0)
    moon_crater_diameter = impactor_to_crater(
        crater_to_impactor(enceladus_crater_diameter, enceladus_impact_parameters),
        moon_impact_parameters,
    )

    # Expected number of 1 km impactor from Singer et al. (2019) impactor size-frequency distribution
    (
        singer_reliable_impactor_diameter,
        singer_normalised_density,
        singer_all_impactor_diameter,
        singer_normalised_density_all_impactor_diameter,
    ) = calculate_impactor_density(
        moon_crater_diameter,
        moon_impact_parameters,
        singer_impactor_sfd,
        minimum_impactor_diameter=0.1,
        maximum_impactor_diameter=15.0,
    )

    # From Zahnle et al. (2003) impactor size-frequency distribution
    (
        zahnle_reliable_impactor_diameter,
        zahnle_normalised_density,
        zahnle_all_impactor_diameter,
        zahnle_normalised_density_all_impactor_diameter,
    ) = calculate_impactor_density(
        moon_crater_diameter,
        moon_impact_parameters,
        zahnle_impactor_sfd,
        minimum_impactor_diameter=0.1,
        maximum_impactor_diameter=20.0,
    )

    # From kirchoff and Schenk (2009) *crater* size-frequency distribution
    (
        kirchoff_reliable_crater_diameter,
        kirchoff_normalised_density,
        kirchoff_normalised_density_all_crater_diameter,
    ) = calculate_crater_density(
        enceladus_crater_diameter,
        moon_impact_parameters,
        kirchoff_crater_sfd,
        minimum_enceladus_crater_diameter=1.0,
        maximum_enceladus_crater_diameter=30.0,
    )

    # Creating function with the crater diameter range plotted/discussed in the cited publication
    (
        wong_reliable_crater_diameter,
        wong_normalised_density,
        _wong_normalised_density_all_crater_diameter,
    ) = calculate_crater_density(
        enceladus_crater_diameter,
        moon_impact_parameters,
        wong_crater_sfd,
        minimum_enceladus_crater_diameter=0.8,
        maximum_enceladus_crater_diameter=44.0,
    )

    all_data = [
        # Plotting the crater production function of Enceladus of this work
        [
            wong_reliable_crater_diameter,
            wong_normalised_density,
            "Enceladus CPF (this work)",
            "black",
            "solid",
        ],
        # Plotting Kirchoff and Schenk 2009, crater size-frequency distribution of mid-latitude crater plains on Encealdus
        ## Plot reliable crater ranges in solid line
        [
            kirchoff_reliable_crater_diameter,
            kirchoff_normalised_density,
            "Enceladus cratered plains",  # Kirchoff 2009",
            "#144e62",
            "solid",
        ],
        ## Plot crater ranges outside the reliable ranges in dotted line
        [
            moon_crater_diameter,
            kirchoff_normalised_density_all_crater_diameter,
            None,
            "#144e62",  # BLUE
            "dotted",
        ],
        # Plotting Zahnle et al 2003 Case A impactor size-frequency distribution
        ## Plot reliable crater ranges in solid line
        [
            impactor_to_crater(
                zahnle_reliable_impactor_diameter, moon_impact_parameters
            ),
            zahnle_normalised_density,
            "Jupiter family comet",  # Zahnle 2003",
            "#D29343",  # ORANGE
            "solid",
        ],
        ## Plot crater ranges outside the reliable ranges in dotted line
        [
            impactor_to_crater(zahnle_all_impactor_diameter, moon_impact_parameters),
            zahnle_normalised_density_all_impactor_diameter,
            None,
            "#D29343",  # ORANGE
            "dotted",
        ],
        # Plotting Singer et al 2019 impactor size-frequency distribution
        ## Plot reliable crater ranges in solid line
        [
            impactor_to_crater(
                singer_reliable_impactor_diameter, moon_impact_parameters
            ),
            singer_normalised_density,
            "Kuiper belt objects",  # Singer 2019",
            "#fdb7bc",  # PINK
            "solid",
        ],
        ## Plot crater ranges outside the reliable ranges in dotted line
        [
            impactor_to_crater(singer_all_impactor_diameter, moon_impact_parameters),
            singer_normalised_density_all_impactor_diameter,
            None,
            "#fdb7bc",  # PINK
            "dotted",
        ],
    ]

    # === 5) Plot the normalised functions
    plot_figure(all_data)


def plot_figure(all_data):
    # loop through all the crater size-frequency distribution (including the Enceladus CPF in this
    # study) and overplot them in one figure
    plt.figure(figsize=(7, 10))

    for data in all_data:
        crater_diameter, size_function, label, color, linestyle = data
        plt.plot(
            crater_diameter,
            size_function,
            label=label,
            color=color,
            linestyle=linestyle,
        )

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Crater diameter on Enceladus (km)", fontsize=14)
    plt.ylabel("Cumulative crater density   (km$^{-2}$)", fontsize=14)
    plt.xlim([8e-1, 4e1])
    plt.ylim([1e-6, 2e0])
    plt.legend()
    plt.grid(True, which="both", color="#aeaeae", linewidth=0.5)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.savefig("./fig/Wongetal2025_fig1.png")
    plt.show()


if __name__ == "__main__":
    main()
