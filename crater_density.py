import numpy as np
import matplotlib.pyplot as plt

output_png = "crater_density.png"

## === 1) Original coefficient of the 10-degree polynomial_1km (dimensionless) ===
coefficient = np.array(
    [
        6.30769779,
        -2.35027511,
        -5.55729243,
        24.95085739,
        -77.26792659,
        212.99933586,
        -406.68327268,
        464.57706063,
        -303.36628993,
        104.74296626,
        -14.84990662,
    ]
)


# === 2) Crater scaling law (Zahnle et al 2003) ===
def impactor_to_crater(impactor_diameter):
    # impactor_to_crater(impactor_diameter) = crater_diameter

    # Enceladus' surface gravity in m/s2
    surface_gravity = 0.1134  
    # impactor_density-to-satellites_density ratio
    # bulk density of impactor ~ 400 kg/m3, icy crust ~ 1500 kg/m3
    density_ratio = 0.4  
    average_impact_velocity = 22.174  # in km/s (Wong et al 2021, 2023)
    # simple-to-complex crater transition diameter in km (Schenk 1991)
    transition_diameter = 15

    crater_diameter = (
        impactor_diameter
        * (0.1431**-1)
        * (surface_gravity**-0.282)
        * (density_ratio**0.427)
        * (average_impact_velocity**0.564)
        * (transition_diameter**-0.192)
    ) ** (1 / 1.0897)

    return crater_diameter


def crater_to_impactor(crater_diameter):
    # crater_to_impactor(crater_diameter) = impactor_diameter
    #
    # Enceladus' surface gravity in m/s2
    surface_gravity = 0.1134  
    # impactor_density-to-satellites_density ratio
    # bulk density of impactor ~ 400 kg/m3, icy crust ~ 1500 kg/m3
    density_ratio = 0.4  
    average_impact_velocity = 22.174  # in km/s (Wong et al 2021, 2023)
    # simple-to-complex crater transition diameter in km (Schenk 1991)
    transition_diameter = 15

    impactor_diameter = (
        0.1431
        * (crater_diameter**1.0897)
        * (surface_gravity**0.282)
        * (density_ratio**-0.427)
        * (average_impact_velocity**-0.564)
        * (transition_diameter**0.192)
    )
    return impactor_diameter


# === 3) Calculate the number of impactor ===

# Singer et al. (2019) size-frequency distribution of impactors
# study form the Pluto-Charon crater population
def singer_impactor_sfd(impactor_diameter):
    return np.where(
        impactor_diameter >= 1.0,
        (1.0 / impactor_diameter) ** 2.0,
        (1.0 / impactor_diameter) ** 0.7,
    )


# Cumulative number of Scattered disk object larger than certain diameter
def number_of_impactor(impactor_diameter):
    current_number_of_10km_impactor = 2e7
    reference_impactor_diameter = 10.0 #km
    number_reference_impactor = singer_impactor_sfd(reference_impactor_diameter)
    number_xkm_impactor = singer_impactor_sfd(impactor_diameter)
    #    mean_removal_rate = 320
    #    return mean_removal_rate * current_number_of_10km_impactor * distribution_xkm / distribution_10km
    return current_number_of_10km_impactor * number_xkm_impactor / number_reference_impactor


# Expected impact density (estimated from current number of 10 km object)
def expected_impact_density(impactor_diameter):
    impact_probabliity = 2.4854e-8  # Wong et al 2023 
    spherical_area = 798648
    return number_of_impactor(impactor_diameter) / spherical_area * impact_probabliity


# === 4) Normalised the CPF and provide a physical dimension as crater density ===
reference_crater_diameter = impactor_to_crater(1.0)
## expected output: 22.538654175960747 km
print(f"  Enceladus crater size created by 1 km impactor = {reference_crater_diameter:.2f} km")
log10_crater_diameter = np.log10(reference_crater_diameter)

# Polynomial evaluation at 1 km impactor
polynomial_1km = sum(
    coefficient[i] * log10_crater_diameter**i for i in range(len(coefficient))
)
## expected output: 2.1670057703366865
print(f"   - Polynomial evaluated at 1 km impactor = {polynomial_1km:.4f}")

# Estimated impact density for 1 km impactor
modelled_density_1km_impactor = expected_impact_density(1.0)

# Normalisation constant
normalisation_constant = np.log10(modelled_density_1km_impactor) - polynomial_1km
## expected output: -6.372934889096243
print(f"   - Normalisation constant = {normalisation_constant:.4f}")

# Adjust the zero-degree coefficient
zero_degree_coefficient = coefficient[0] + normalisation_constant
## expected output: -0.06523709909624298
print(f"   - New 0th-degree coefficient = {zero_degree_coefficient:.4f}") 

coefficient_norm = coefficient.copy()
coefficient_norm[0] = zero_degree_coefficient
print(f"coefficient_norm = {coefficient_norm}")


# Evaluate the normalized CPF
def cpf_normalised(crater_diameter, coefficient):
    log_diameter = np.log10(crater_diameter)
    logS = sum(coefficient[i] * (log_diameter**i) for i in range(len(coefficient)))
    return 10**logS


# Define the other size-frequency distributions (SFDs) for comparison 
# (already included Singer et al 2019 above)

# Zahnle et al 2003 case A impactor size-frequency distribution
# studied for the Jovican satellites crater
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


# d) Kirchoff and Schenk 2009 for Mid-latitude crater plains
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


# === 5) Creating the cumulative crater density function with 
# our crater production function and 
# different impactor or crater size-frequency distribution ===


# Define the crater diameter ranges for the plot
crater_range = 10 ** (np.linspace(np.log10(0.8), np.log10(44), 100))


# Our crater production function (CPF) for Enceladus
cpf_normalised_function = cpf_normalised(crater_range, coefficient_norm)
crater_density_1km = cpf_normalised(1, coefficient_norm)
## expected output: 0.8605238288
print(f"  Crater density for 1 km crater on Enceladus: {crater_density_1km:.2f} crater per km^2")


# From Singer et al. (2019) impactor size-frequency distribution
singer_1km_impactor = singer_impactor_sfd(1.0)

# Creating function with the range plotted/discussed in the cited publication
reliable_impactor_range_singer = 10 ** (np.linspace(np.log10(0.1), np.log10(15), 1000))
singer_impactor_function = singer_impactor_sfd(reliable_impactor_range_singer)
singer_normalised_function = (
    singer_impactor_function / singer_1km_impactor * modelled_density_1km_impactor
)

# Creating function with the entire crater ranges in this work
singer_all_impactor_range = crater_to_impactor(crater_range)
singer_impactor_function_all_impactor_range = singer_impactor_sfd(singer_all_impactor_range)
singer_normalised_function_all_impactor_range = (
    singer_impactor_function_all_impactor_range
    / singer_1km_impactor
    * modelled_density_1km_impactor
)

# From Zahnle et al. (2003) impactor size-frequency distribution
zahnle_1km_impactor = zahnle_impactor_sfd(1.0)

# Creating function with the range plotted/discussed in the cited publication
reliable_impactor_range_zahnle = 10 ** (np.linspace(np.log10(0.1), np.log10(20), 1000))
zahnle_impactor_function = zahnle_impactor_sfd(reliable_impactor_range_zahnle)
zahnle_normalised_function = (
    zahnle_impactor_function / zahnle_1km_impactor * modelled_density_1km_impactor
)

# Creating function with the entire crater ranges in this work
zahnle_all_impactor_range = crater_to_impactor(crater_range)
zahnle_impactor_function_all_impactor_range = zahnle_impactor_sfd(zahnle_all_impactor_range)
zahnle_normalised_function_all_impactor_range = (
    zahnle_impactor_function_all_impactor_range
    / zahnle_1km_impactor
    * modelled_density_1km_impactor
)


# From kirchoff and Schenk (2009) *crater* size-frequency distribution
kirchoff_1km_impactor = kirchoff_crater_sfd(impactor_to_crater(1.0))
reliable_crater_range_kirchoff = 10 ** (np.linspace(np.log10(1.0), np.log10(30), 1000)) #np.linspace(1, 30, 1000)

# Creating function with the range plotted/discussed in the cited publication
kirchoff_crater_function = kirchoff_crater_sfd(reliable_crater_range_kirchoff)
kirchoff_normalised_function = (
    kirchoff_crater_function / kirchoff_1km_impactor * modelled_density_1km_impactor
)

# Creating function with the entire crater ranges in this work
kirchoff_crater_function_all_crater_range = kirchoff_crater_sfd(crater_range)
kirchoff_normalised_function_all_crater_range = (
    kirchoff_crater_function_all_crater_range
    / kirchoff_1km_impactor
    * modelled_density_1km_impactor
)


# === 7) Plot the normalised functions
plt.figure(figsize=(7, 10))

# Plotting the crater production function of Enceladus of this work
plt.plot(crater_range, cpf_normalised_function, label="This work", color="black")

# Plotting Kirchoff and Schenk 2009, crater size-frequency distribution of mid-latitude crater plains on Encealdus 
## Plot reliable crater ranges in solid line
plt.plot(
    reliable_crater_range_kirchoff,
    kirchoff_normalised_function,
    label="Kirchoff 2009",
    color="#144e62",
)
## Plot crater ranges outside the reliable ranges in dotted line
plt.plot(
    crater_range,
    kirchoff_normalised_function_all_crater_range,
    color="#144e62",
    linestyle="dotted",
)

# Plotting Zahnle et al 2003 Case A impactor size-frequency distribution
## Plot reliable crater ranges in solid line
plt.plot(
    impactor_to_crater(reliable_impactor_range_zahnle),
    zahnle_normalised_function,
    label="Zahnle 2003",
    color="#D29343",
)
## Plot crater ranges outside the reliable ranges in dotted line
plt.plot(
    crater_range,
    zahnle_normalised_function_all_impactor_range,
    color="#D29343",
    linestyle="dotted",
)

# Plotting Singer et al 2019 impactor size-frequency distribution
## Plot reliable crater ranges in solid line
plt.plot(
    impactor_to_crater(reliable_impactor_range_singer),
    singer_normalised_function,
    label="Singer 2019",
    color="#fdb7bc",
)
## Plot crater ranges outside the reliable ranges in dotted line
plt.plot(
    crater_range,
    singer_normalised_function_all_impactor_range,
    color="#fdb7bc",
    linestyle="dotted",
)


plt.xscale("log")
plt.yscale("log")
plt.xlabel("Crater diameter on Enceladus", fontsize=14)
plt.ylabel("Cumulative crater density   (km$^{-2}$)", fontsize=14)
plt.xlim([8e-1, 4e1])
plt.ylim([1e-6, 2e0])
plt.legend()
plt.grid(True, which="both", color="#aeaeae", linewidth=0.5)  # , linestyle='dotted')
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.savefig(output_png)  # Save the plot as an image
plt.show()

# coefficient = [ -14.84990662, 104.74296626, -303.36628993, 464.57706063
#              , -406.68327268, 212.99933586, -77.26792659, 24.95085739
#              , -5.55729243, -2.35027511, 6.30769779]
#
# def cpf(crater_diameter, coefficient):
#    poly = np.poly1d(coefficient)
#    return 10**(poly(np.log10(crater_diameter)))
