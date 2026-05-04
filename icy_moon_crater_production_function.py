#!/usr/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt

"""
Crater Production Function (CPF) Implementation and Plotting 

This script provides a reproducible implementation of the Crater Production Function (CPF) of Enceladus, published in Wong et al (updated years at the end). It models the expected number of craters on Enceladus as a function of impactors' size, impactor velocities and other parameters.

Key features:
- Provide the polynomial function of the CPF.
- Compute crater densities for arbitrary impactor diameters.
- Compare with other published impactor and crater Size-Frequency Distributions (SFDs), i.e.,
    - crater SFD of Enceladus mid-latitude cratered plains from Kirchoff and Schenk (2009),
    - jupiter family comet impactor SFD from Case A of Zahnle et al. (2003),
    - triton crater SFD from Case B of Zahnle et al. (2003),
    - kuiper belt object impactor SFD from Singer et al. (2019).

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


# ==================================================================================
# === Icy Moon Class Definition ===
class IcyMoon:
    """
    Represents an icy moon with impact parameters for crater scaling calculations.
    """

    def __init__(
        self,
        name: str,
        surface_gravity: float,
        density_ratio: float,
        average_impact_velocity: float,
        transition_diameter: float,
        impact_probability: float,
        spherical_area: float,
    ) -> None:

        # Name of the icy momon
        self.name: str = name
        # Surface gravity in m/s²
        self.surface_gravity: float = surface_gravity
        # Impactor's density to satellite's density ratio
        # Bulk density of impactor ~ 400 kg/m3, satellite's icy crust ~ 1000 kg/m3
        self.density_ratio: float = density_ratio
        # Average impact velocity in km/s from Genga N-body simulations (Wong et al 2021, 2023)
        self.average_impact_velocity: float = average_impact_velocity
        # Simple-to-complex crater transition diameter in km (Schenk 1991)
        self.transition_diameter: float = transition_diameter
        # average impact probability of Enceladus (Wong et al 2023)
        self.impact_probability: float = impact_probability
        # global spherical area of Enceladus in km^{2}
        self.spherical_area: float = spherical_area

        self.crater_production_functions: Dict[str, Any] = {}

    # === Crater scaling law (Zahnle et al 2003) ===
    def impactor_to_crater(self, impactor_diameter: float) -> float:
        """
        Convert impactor diameter (km) to resulting crater diameter (km).
        Uses the moon's specific impact parameters.
        """
        crater_diameter = (
            impactor_diameter
            * (0.1431**-1)
            * (self.surface_gravity**-0.282)
            * (self.density_ratio**0.427)
            * (self.average_impact_velocity**0.564)
            * (self.transition_diameter**-0.192)
        ) ** (1 / 1.0897)
        return crater_diameter

    # === Inverse crater scaling law (Zahnle et al 2003) ===
    def crater_to_impactor(self, crater_diameter: float) -> float:
        """
        Convert crater diameter (km) to typical impactor diameter (km).
        Uses the moon's specific impact parameters.
        """
        impactor_diameter = (
            0.1431
            * (crater_diameter**1.0897)
            * (self.surface_gravity**0.282)
            * (self.density_ratio**-0.427)
            * (self.average_impact_velocity**-0.564)
            * (self.transition_diameter**0.192)
        )
        return impactor_diameter

    def __repr__(self):
        return f"IcyMoon('{self.name}', g={self.surface_gravity} m/s²)"


# Create specific icy moon instances
# with the impact parameters of Enceladus
def create_enceladus() -> IcyMoon:
    return IcyMoon(
        name="Enceladus",
        surface_gravity=0.1134,
        density_ratio=0.4,
        average_impact_velocity=22.174,
        transition_diameter=15,
        impact_probability=2.4854e-8,
        spherical_area=798_648,
    )

    # with the impact parameters of Tethys


def create_tethys() -> IcyMoon:
    return IcyMoon(
        name="Tethys",
        surface_gravity=0.1462,
        density_ratio=0.4,
        average_impact_velocity=19.997,
        transition_diameter=15,
        impact_probability=4.74e-8,
        spherical_area=3.543e6,
    )


# ==================================================================================
# === 1) Evaluate the normalized CPF ===
def generate_CPF() -> np.poly1d:
    """
    Generate a normalised crater production function (CPF) as a 10th-degree polynomial.
    Coefficients are in increasing order (a0, a1, ..., a10).
    """
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
    crater_production_function = np.poly1d(cpf_coefficients[::-1])
    return crater_production_function


def wong_crater_sfd(crater_diameter: float | np.ndarray) -> float | np.ndarray:
    """
    Calculate the expected cumulative crater density (crater/km^2) for given crater diameter(s) in km
    """
    crater_production_function = generate_CPF()
    return 10 ** (crater_production_function(np.log10(crater_diameter)))


# ==================================================================================
# === 2) Calculate the number of impactor in the outer Solar System ===
def singer_impactor_sfd(impactor_diameter: float | np.ndarray) -> float | np.ndarray:
    """Singer et al. (2019) size-frequency distribution of impactors"""
    return np.where(
        impactor_diameter >= 1.0,
        (1.0 / impactor_diameter) ** 2.0,
        (1.0 / impactor_diameter) ** 0.7,
    )


def number_of_impactor(
    impactor_diameter: float | np.ndarray,
    impactor_sfd: callable,
) -> float | np.ndarray:
    """
    Calculate the cumulative number of scattered disk objects larger than a given diameter,
    normalized to the reference size (10 km) using Singer's size-frequency distribution.
    """
    reference_impactor_diameter = 10.0  # km
    current_number_of_10km_impactor = 2e7
    reference_impactor_count = impactor_sfd(reference_impactor_diameter)
    scaled_impactor_count = impactor_sfd(impactor_diameter)
    return (
        current_number_of_10km_impactor
        * scaled_impactor_count
        / reference_impactor_count
    )


def expected_impact_density(
    impactor_diameter: float | np.ndarray,
    IcyMoon,
    impactor_sfd: callable,
) -> float | np.ndarray:
    """Expected impact density estimated from current number of 10 km object"""
    return (
        number_of_impactor(impactor_diameter, impactor_sfd)
        / IcyMoon.spherical_area
        * IcyMoon.impact_probability
    )


# ==================================================================================
# === 3) Size-Frequency Distribution Functions ===
def zahnle_jfc_impactor_sfd(
    impactor_diameter: float | np.ndarray,
) -> float | np.ndarray:
    """Zahnle et al 2003 case A impactor size-frequency distribution (Jovian-family comets)"""
    return np.where(
        impactor_diameter >= 5.0,
        (1.8749 / impactor_diameter) ** 2.5,
        np.where(
            impactor_diameter >= 1.5,
            (1.1817 / impactor_diameter) ** 1.7,
            (1.0 / impactor_diameter) ** 1.0,
        ),
    )


def zahnle_tc_impactor_sfd(impactor_diameter: float | np.ndarray) -> float | np.ndarray:
    """Zahnle et al 2003 case B impactor size-frequency distribution (Triton's crater)"""
    return np.where(
        impactor_diameter >= 1.5,
        (1.1385 / impactor_diameter) ** 2.5,
        (1.0 / impactor_diameter) ** 1.7,
    )


def kirchoff_crater_sfd(crater_diameter: float | np.ndarray) -> float | np.ndarray:
    """Kirchoff and Schenk 2009 for Mid-latitude crater plains on Enceladus"""
    return np.where(
        crater_diameter >= 7.0,
        (2.3671 / crater_diameter) ** 2.99,
        np.where(
            crater_diameter >= 4.0,
            (1.7203 / crater_diameter) ** 2.310,
            (1.0 / crater_diameter) ** 1.406,
        ),
    )


def tethys_crater_sfd(crater_diameter: float | np.ndarray) -> float | np.ndarray:
    """Kirchoff and Schenk 2009 for crater plains on Tethys"""
    return np.where(
        crater_diameter >= 10.0,
        (1.7131 / crater_diameter) ** 2.22,
        (1.0 / crater_diameter) ** 1.701,
    )


# ==================================================================================
# === 4) Helper Functions ===
def print_crater_density_table(crater_diameter_list: list[float]) -> None:
    print("\nExpected crater density on Enceladus using the CPF of this work:")
    print(f"{'Crater Diameter (km)':>20} | {'Expected Density (crater/km^2)':>30}")
    print(f"{'-' * 20} | {'-' * 30}")
    for crater_diameter in crater_diameter_list:
        expected_crater_density = wong_crater_sfd(crater_diameter)
        print(f"{crater_diameter:20.2g} | {expected_crater_density:30.3g}")


def make_diameter_range(min: float, max: float) -> np.ndarray:
    """Creating function with the crater or impactor diameter range plotted/discussed in the cited publication"""
    return 10 ** (np.linspace(np.log10(min), np.log10(max), 1000))


def calculate_impactor_density(
    crater_diameter: float | np.ndarray,
    IcyMoon,
    impactor_sfd: callable,
    minimum_impactor_diameter: float,
    maximum_impactor_diameter: float,
) -> tuple[np.ndarray, np.ndarray, float | np.ndarray, float | np.ndarray]:
    """Calculate impactor density using the impactor size-frequency distribution"""
    # Expected number of 1 km impactor from the impactor size-frequency distribution
    one_km_impactor = impactor_sfd(1.0)
    expected_one_km_impactor = expected_impact_density(
        1.0, IcyMoon, singer_impactor_sfd
    )

    # Creating function with the impactor diameter range plotted/discussed in the cited publication
    reliable_impactor_diameter = make_diameter_range(
        minimum_impactor_diameter, maximum_impactor_diameter
    )
    # Calculating the impactor density using the 'impactor_sfd'
    impactor_density = impactor_sfd(reliable_impactor_diameter)
    normalised_density_reliable_impactor_diameter = (
        impactor_density / one_km_impactor * expected_one_km_impactor
    )

    # Creating function with the entire observed crater ranges in this work
    all_impactor_diameter = IcyMoon.crater_to_impactor(crater_diameter)
    impactor_density_all_impactor_diameter = impactor_sfd(all_impactor_diameter)
    normalised_density_all_impactor_diameter = (
        impactor_density_all_impactor_diameter
        / one_km_impactor
        * expected_one_km_impactor
    )
    return (
        reliable_impactor_diameter,
        normalised_density_reliable_impactor_diameter,
        all_impactor_diameter,
        normalised_density_all_impactor_diameter,
    )


def calculate_crater_density(
    crater_diameter: float | np.ndarray,
    IcyMoon,
    crater_sfd: callable,
    minimum_crater_diameter: float,
    maximum_crater_diameter: float,
) -> tuple[np.ndarray, np.ndarray, float | np.ndarray]:
    """Calculate crater density using the crater size-frequency distribution"""
    # Get the Enceladus instance to convert between impactor and crater diameters
    expected_crater_size_from_one_km_impactor = IcyMoon.impactor_to_crater(1.0)
    # Calculate expected 1 km impactor density from the provided *crater* size-frequency distribution
    unscaled_one_km_impact_density_on_moon = crater_sfd(
        expected_crater_size_from_one_km_impactor
    )
    expected_one_km_impact_density_on_moon = expected_impact_density(
        1.0, IcyMoon, singer_impactor_sfd
    )

    # Creating function with the crater diameter range plotted/discussed in the cited publication
    reliable_crater_diameter = make_diameter_range(
        minimum_crater_diameter, maximum_crater_diameter
    )

    # Calculating the crater density using "crater_sfd"
    unscaled_crater_density = crater_sfd(reliable_crater_diameter)
    normalised_density_reliable_crater_diameter = (
        unscaled_crater_density
        / unscaled_one_km_impact_density_on_moon
        * expected_one_km_impact_density_on_moon
    )

    unscaled_density_all_crater_diameter = crater_sfd(crater_diameter)
    normalised_density_all_crater_diameter = (
        unscaled_density_all_crater_diameter
        / unscaled_one_km_impact_density_on_moon
        * expected_one_km_impact_density_on_moon
    )
    return (
        reliable_crater_diameter,
        normalised_density_reliable_crater_diameter,
        normalised_density_all_crater_diameter,
    )


# ==================================================================================
# === 5) Plotting ===
def plot_figure(all_data: list[tuple[np.ndarray, np.ndarray, str, str, str]]) -> None:
    """Plot the normalised functions"""
    plt.figure(figsize=(7, 10))

    # loop through all the crater size-frequency distribution (including the Enceladus CPF in this
    # study) and overplot them in one figure
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
    plt.ylim([1e-6, 3e0])
    plt.legend()
    plt.grid(True, which="both", color="#aeaeae", linewidth=0.5)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.savefig(
        "./Tethys_trial_ZahnleSFD_fro_241015_current_rate_calculation_excelfile_lumo.png"
    )
    plt.show()


def main() -> None:
    for x in sys.argv[1:]:
        try:
            input_crater_diameter = float(x)
            if input_crater_diameter < 0.8 or input_crater_diameter > 44:
                usage()
                exit()
        except ValueError:
            usage()
            exit()

    crater_diameter_list = [float(arg) for arg in sys.argv[1:]] or [1.0]
    print_crater_density_table(crater_diameter_list)

    # enceladus = create_enceladus()
    moon = create_enceladus()  # Change to create_tethys() for Tethys

    moon_crater_diameter = make_diameter_range(min=0.8, max=44.0)

    # Expected number of 1 km impactor from Singer et al. (2019) impactor size-frequency distribution
    (
        singer_reliable_impactor_diameter,
        singer_normalised_density_reliable_impactor_diameter,
        singer_all_impactor_diameter,
        singer_normalised_density_all_impactor_diameter,
    ) = calculate_impactor_density(
        moon_crater_diameter, moon, singer_impactor_sfd, 0.1, 15.0
    )

    # From Zahnle et al. (2003) Case A JFC impactor size-frequency distribution
    (
        zahnle_jfc_reliable_impactor_diameter,
        zahnle_jfc_normalised_density_reliable_impactor_diameter,
        zahnle_jfc_all_impactor_diameter,
        zahnle_jfc_normalised_density_all_impactor_diameter,
    ) = calculate_impactor_density(
        moon_crater_diameter, moon, zahnle_jfc_impactor_sfd, 0.1, 20.0
    )

    # From Zahnle et al. (2003) Case B Triton's impactor size-frequency distribution
    (
        zahnle_tc_reliable_impactor_diameter,
        zahnle_tc_normalised_density_reliable_impactor_diameter,
        zahnle_tc_all_impactor_diameter,
        zahnle_tc_normalised_density_all_impactor_diameter,
    ) = calculate_impactor_density(
        moon_crater_diameter, moon, zahnle_tc_impactor_sfd, 0.1, 30.0
    )

    # From kirchoff and Schenk (2009) *crater* size-frequency distribution
    (
        kirchoff_reliable_crater_diameter,
        kirchoff_normalised_density_reliable_crater_diameter,
        kirchoff_normalised_density_all_crater_diameter,
    ) = calculate_crater_density(
        moon_crater_diameter, moon, kirchoff_crater_sfd, 1.0, 30.0
    )

    # From kirchoff and Schenk (2009) Tethys' *crater* size-frequency distribution
    (
        tethys_reliable_crater_diameter,
        tethys_normalised_density_reliable_crater_diameter,
        tethys_normalised_density_all_crater_diameter,
    ) = calculate_crater_density(
        moon_crater_diameter, moon, tethys_crater_sfd, 1.0, 30.0
    )

    # Creating function with the crater diameter range plotted/discussed in the cited publication
    (
        wong_reliable_crater_diameter,
        wong_normalised_density_reliable_crater_diameter,
        _wong_normalised_density_all_crater_diameter,
    ) = calculate_crater_density(moon_crater_diameter, moon, wong_crater_sfd, 0.8, 44.0)

    all_data = [
        # Plotting the crater production function of Enceladus of this work
        [
            wong_reliable_crater_diameter,
            wong_normalised_density_reliable_crater_diameter,
            "Enceladus CPF (this work)",
            "black",
            "solid",
        ],
        # Plotting Kirchoff and Schenk 2009, crater size-frequency distribution of mid-latitude crater plains on Encealdus
        ## Plot reliable crater ranges in solid line
        [
            kirchoff_reliable_crater_diameter,
            kirchoff_normalised_density_reliable_crater_diameter,
            "Enceladus cratered plains (Kirchoff and Schenk, 2009)",
            "#144e62",
            "solid",
        ],
        ## Plot crater ranges outside the reliable ranges in dotted line
        [
            moon_crater_diameter,
            kirchoff_normalised_density_all_crater_diameter,
            None,
            "#144e62",
            "dotted",
        ],
        # Plotting Zahnle et al 2003 Case A JFC impactor size-frequency distribution
        ## Plot reliable crater ranges in solid line
        [
            moon.impactor_to_crater(zahnle_jfc_reliable_impactor_diameter),
            zahnle_jfc_normalised_density_reliable_impactor_diameter,
            "Jupiter-family comet (Zahnle et al, 2003)",
            "#D29343",
            "solid",
        ],
        ## Plot crater ranges outside the reliable ranges in dotted line
        [
            moon.impactor_to_crater(zahnle_jfc_all_impactor_diameter),
            zahnle_jfc_normalised_density_all_impactor_diameter,
            None,
            "#D29343",
            "dotted",
        ],
        # Plotting Zahnle et al 2003 Case B Triton Crater impactor size-frequency distribution
        ## Plot reliable crater ranges in solid line
        [
            moon.impactor_to_crater(zahnle_tc_reliable_impactor_diameter),
            zahnle_tc_normalised_density_reliable_impactor_diameter,
            "Triton's crater (Zahnle et al, 2003)",
            "#687B3E",
            "solid",
        ],
        ## Plot crater ranges outside the reliable ranges in dotted line
        [
            moon.impactor_to_crater(zahnle_tc_all_impactor_diameter),
            zahnle_tc_normalised_density_all_impactor_diameter,
            None,
            "#687B3E",
            "dotted",
        ],
        # Plotting Singer et al 2019 impactor size-frequency distribution
        ## Plot reliable crater ranges in solid line
        [
            moon.impactor_to_crater(singer_reliable_impactor_diameter),
            singer_normalised_density_reliable_impactor_diameter,
            "Kuiper belt objects (Singer et al, 2009)",
            "#fdb7bc",
            "solid",
        ],
        ## Plot crater ranges outside the reliable ranges in dotted line
        [
            moon.impactor_to_crater(singer_all_impactor_diameter),
            singer_normalised_density_all_impactor_diameter,
            None,
            "#fdb7bc",
            "dotted",
        ],
        # Plotting Kirchoff and Schenk 2009, crater size-frequency distribution of mid-latitude crater plains on Encealdus
        ## Plot reliable crater ranges in solid line
        [
            tethys_reliable_crater_diameter,
            tethys_normalised_density_reliable_crater_diameter,
            "Tethys cratered plains (Kirchoff and Schenk, 2009)",
            "#00ffff",
            "solid",
        ],
    ]

    plot_figure(all_data)


if __name__ == "__main__":
    main()

# TODO: 
# computational: 
# - function oriented --> object oriented: 
#   create class for the size-frequency distribution
# - create input CSV for reading in the attributes of the different IcyMoon instant and add more IcyMoon (Europe, Ganymede, Callisto, Titan, etc.)
# divided into module files (modularise the code)
# set up file for test case
# set up configuration file as user interface for input parameters
#
"""
Jotting for future steps:

separated into different module file
** with simple test case

add "attribute" host star

write the SFD of cratre and impactor as object 
- impactor or crater
- the function 
- for which icy moon, prompt if not correct
- default color line
- reliable ranges
** with simple test case

have configurataion file for the following one case: 
- diff size for one icy moon
icy moon: could be list (repeated plot)
SFD: could be list
color: could be list (or none use the default color)
line style: 
save_files
** with simple test case

have configruation file for the following case:
- SFD for different icy moon
icy moon: list
SFD one to one correspondance
color:
linetsyle
save_files
** with simple test case

add the cratre density overplotting and three plot overplotting
crater data (only one file, later expand to list of file? )
SFD (only one)

import the crater_data libarray for the other SFD

===========================

inside a class: method (function), attribute (variable)
python file (*.py): module
collection of module (*.py): package
"""

