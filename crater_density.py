import numpy as np
import matplotlib.pyplot as plt

output_png = "crater_density.png"

## === 1) Original coefficient of the 10-degree polynomial (dimension) ===
coefficient = np.array([
    6.30769779, -2.35027511, -5.55729243, 24.95085739, -77.26792659,
    212.99933586, -406.68327268, 464.57706063, -303.36628993, 104.74296626, 
    -14.84990662])

# === 2) Crater scaling law (Zahnle et al 2003) ===
def impactor_to_crater(impactor_diameter):
# impactor_to_crater(impactor_diameter) = crater_diameter
    surface_gravity = 0.1134 # Emi: give unit 
    # impactor density-to-satellites bulk density ratio
    density_ratio = 0.4 #(Emi: perhaps need to explain more) 
    average_impact_velocity = 22.174 # in km/s (Wong et al 2021, 2023)
    # simple-to-complex crater transition diameter (Schenk 1991)
    transition_diameter = 15 

    crater_diameter = (
        impactor_diameter
        * (0.1431 ** -1)
        * (surface_gravity ** -0.282)
        * (density_ratio ** 0.427)
        * (average_impact_velocity ** 0.564)
        * (transition_diameter ** -0.192)
    ) ** (1 / 1.0897)

    return crater_diameter

# Example: 1-km impactor
D_cr_1km = impactor_to_crater(1.0)
print(f"Typical crater diameter from 1-km impactor on Enceladus: {D_cr_1km:.2f} km")

# === 3) Calculate the number of impactor ===

# Singer et al. (2019) size-frequency distribution of impactors
def impactor_size_frequency_distribution(impactor_diameter):
    if impactor_diameter > 1.0:
        return (1.0 / impactor_diameter) ** 2.0
    else:
        return (1.0 / impactor_diameter) ** 0.7

# Cumulative number of Scattered disk object larger than certain diameter
def number_of_impactor(impactor_diameter):
    current_number_of_10km_impactor = 2e7
    reference_impactor_diameter = 10.0
    reference_number = impactor_size_frequency_distribution(reference_impactor_diameter)
    x_number = impactor_size_frequency_distribution(impactor_diameter)
#    mean_removal_rate = 320
#    return mean_removal_rate * current_number_of_10km_impactor * distribution_xkm / distribution_10km
    return current_number_of_10km_impactor * x_number / reference_number

# Expected impact density (estimated from current number of 10 km object) 
def expected_impact_density(impactor_diameter):
    impact_probabliity = 2.4854E-8 # reference
    spherical_area = 798648
    return number_of_impactor(impactor_diameter) / spherical_area * impact_probabliity 

# === 4) Normalised the CPF and provide a physical dimension ===
reference_crater_diameter = impactor_to_crater(1.0)
print("Enceladus crater created by 1 km impactor =", reference_crater_diameter)

log10_crater_diameter = np.log10(reference_crater_diameter)  
# Polynomial evaluation at 1 km impactor
polynomial = sum(coefficient[i] * log10_crater_diameter**i for i in range(len(coefficient)))  
print("Polynomial evaluated at 1 km impactor=", polynomial)
#
# Expected impact density for 1 km impactor
expected_density_1km_impactor = expected_impact_density(1.0)  
#
# Normalisation constant
normalisation_constant = np.log10(expected_density_1km_impactor) - polynomial  
print("Normalisation constant =", normalisation_constant)

# Adjust the zero-degree coefficient
zero_degree_coefficient = coefficient[0] + normalisation_constant  
print("New 0th-degree coefficient =", zero_degree_coefficient)  # -0.06523709909624298

coefficient_norm = coefficient.copy()
coefficient_norm[0]  = zero_degree_coefficient 

# === Calulating the expected impact density or crater density ===
# Define the x values for the plot
x_values = 10**(np.linspace(np.log10(0.8), np.log10(44), 100))  # Adjust as needed

# === Evaluate the normalized CPF (dimensionless) ===
def cpf_new(crater_diameter, coefficient):
    log_diameter = np.log10(crater_diameter)
    logS = sum(coefficient[i] * (log_diameter**i) for i in range(len(coefficient)))
    return 10**logS  # dimensionless normalized CPF

#======================================================================

# Calculate the normalised values at x = 22.5709
#normalised = cpf(x_values, coefficient) / cpf(impactor_to_crater(1.0), coefficient) * expected_impact_density(1.0)
#
normalised = cpf_new(x_values, coefficient_norm)

#crater_density_1km = cpf(1, coefficient) / cpf(impactor_to_crater(1.0), coefficient) * expected_impact_density(1.0)

crtaer_density_1km = cpf_new(1, coefficient_norm)

#print(f"crater density for 1 km crater on Enceladus: {crater_density_1km:.2f} km$^2$")
print(f"crater density for 1 km crater on Enceladus: {crtaer_density_1km:.2f} km$^2$")



# Plot the normalised functions
plt.figure(figsize=(8, 10))
#plt.plot(x_values, normalised, label='Enceladus CPF', color='black')
plt.plot(x_values, normalised, label='Enceladus CPF_new', color='black')

plt.xscale('log')
plt.yscale('log')
plt.xlabel("Crater diameter on Enceladus", fontsize=15)
plt.ylabel("Cumulative crater density   (km$^{2}$)", fontsize=15)
plt.xlim([8e-1, 4e1])
plt.ylim([1e-6, 2e0])
plt.legend()
plt.grid(True, which='both', color='#aeaeae', linewidth=0.5) #, linestyle='dotted')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig(output_png)  # Save the plot as an image
plt.show()
#coefficient = [ -14.84990662, 104.74296626, -303.36628993, 464.57706063
#              , -406.68327268, 212.99933586, -77.26792659, 24.95085739
#              , -5.55729243, -2.35027511, 6.30769779]
#
#def cpf(crater_diameter, coefficient):
#    poly = np.poly1d(coefficient)
#    return 10**(poly(np.log10(crater_diameter)))


# crater_to_impactor(crater_diameter) = impactor_diameter
#dreference_ef crater_to_impactor(crater_diameter):
# 
#    surface_gravity = 0.1134 # Emi: give unit 
#    # impactor density-to-satellites bulk density ratio
#    density_ratio = 0.4 #(Emi: perhaps need to explain more) 
#    average_impact_velocity = 22.174 # in km/s (Wong et al 2021, 2023)
#    # simple-to-complex crater transition diameter (Schenk 1991)
#    transition_diameter = 15 
#
#    impactor_diameter = (
#        0.1431
#        * (crater_diameter ** 1.0897)
#        * (surface_gravity ** 0.282)
#        * (density_ratio ** -0.427)
#        * (average_impact_velocity ** -0.564)
#        * (transition_diameter ** 0.192)
#    )
#    return impactor_diameter


## === Normlised cumulative probability function (CPF) ===
#def cpf_normalised(crater_diameter):
#    """
#    Normalised CPF: cumulative crater density on Enceladus [km^-2],
#    anchored at the density expected from 1-km impactor.
#    """
#    # crater diameter corresponding to 1 km impactor
#    crater_by_1km_impactor = impactor_to_crater(1.0)
#
#    # normalisation constant
#    norm_factor = expected_impact_density(1.0) / cpf(crater_by_1km_impactor, coefficient)
#
#    # return the normalised CPF
#    return norm_factor * cpf(crater_diameter, coefficient)



