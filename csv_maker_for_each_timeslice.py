# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 13:00:33 2025

@author: User
"""

import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap
import pandas as pd
import math
from scipy.ndimage import zoom
from matplotlib.animation import FuncAnimation
import os


# Create a folder to store output files
output_folder = "biomass_outputs"
os.makedirs(output_folder, exist_ok=True)

# Define geological period names
geological_periods = {
    12: "Triassic",
    13: "Late_Triassic",
    14: "Early_Jurassic",
    15: "Jurassic",
    16: "Late_Jurassic",
    17: "Cretaceous",
    18: "Late_Cretaceous",
    19: "Eocene",
    20: "Oligocene",
    21: "Miocene",
    22: "Present_Day"
}




# User inputs the timeslice range
start_timeslice = int(input("Enter start timeslice: "))
end_timeslice = int(input("Enter end timeslice: "))
num_runs = 10  # Number of model runs per timeslice

lat = pd.read_csv('lat.csv', header=None).values.flatten()
biomass_results = []

simulation_time_steps = 10
a = 400.17  # Scaling constant
# Constants for the calculations
b = 0.80    # Scaling exponent
T_opt = 8  # Optimal temperature (°C)
O2_opt = 20.72  # Optimal oxygen level (for normalization)
sigma = 25 # Width of Gaussian curve
C = 1       # Scaling constant for metabolic demand (optional for adjustment)
theta = 0.05
r = 100




# Define the temperature factor function (hat shape)
def temp_factor_function(T, k=0.05):
    """
    Models the temperature dependency using a hat-shaped function.
    Parameters:
    T : float or array-like
        Temperature(s) in degrees Celsius.
    k : float
        Steepness of the drop-off (default = 0.05).
    """
    # Left drop-off for T < -10
    left = np.exp(-k * (T + 1)**2) * (T < -1)
    # Plateau for -10 <= T <= 40
    middle = 1.0 * ((T >= -1) & (T <= 40))
    # Right drop-off for T > 40
    right = np.exp(-k * (T - 40)**2) * (T > 40)
    #Gaussian drop-off for T < -1
      
    # Combine all regions
    return left + middle + right

# Updated function to calculate metabolic rate R_d
def R_d(biomass):
    return a * (biomass ** b) * C

# Function to calculate daily intake I_d
def daily_intake(biomass):
    return R_d(biomass)

# Function to calculate carrying capacity K
def K(NPP, biomass_density):
    if biomass_density == 0:
        return NPP  # No competition, carrying capacity is fully available
    else:
        return NPP / biomass_density  # Decrease carrying capacity with higher biomass density

# Updated survival probability function
def survival_probability(K, daily_intake, T, pO2, k=0.05):
    """
    Calculate the survival probability.
    Parameters:
    K : float
        Carrying capacity.
    daily_intake : float
        Daily intake of energy or resources.
    T : float
        Temperature in degrees Celsius.
    pO2 : float
        Partial pressure of oxygen.
    k : float
        Steepness parameter for the temperature factor.
    """
    temp_factor = temp_factor_function(T, k)  # Use the hat-shaped temperature factor
    oxygen_factor = pO2 / O2_opt
    return min(1, max(0, (K * temp_factor * oxygen_factor) / (K + daily_intake)))




for timeslice in range(start_timeslice, end_timeslice + 1):
    biomass_per_run = []
    
    for run in range(num_runs):
        # Load data for current timeslice
        scale_factor = 1
        land_values_lowres_df = pd.read_csv(f'Land_TimeSlice{timeslice}.csv', header=None)
        land = zoom(np.matrix.transpose(land_values_lowres_df.values), zoom=(scale_factor, scale_factor))
        
        temperature_data = pd.read_csv(f'Temperature_TimeSlice{timeslice}.csv', header=None)
        T = zoom(np.nan_to_num(np.matrix.transpose(temperature_data.values), nan=999), zoom=(scale_factor, scale_factor))
        
        resource_values_lowres_df = pd.read_csv(f'NPP_TimeSlice{timeslice}.csv', header=None)
        resource_transpose = np.matrix.transpose(resource_values_lowres_df.values)
        resource_scaled = resource_transpose * 0.1
        resource_scaled = np.nan_to_num(resource_scaled, nan=0)
        
        
        po2_values = {
            22: 20.72, 21: 21.57, 20: 21.41, 19: 25.12, 18: 26.65, 17: 25.86,
            16: 25.86, 15: 25.13, 14: 28.6, 13: 32.38, 12: 30.77, 11: 31.05,
            10: 27.25, 9: 25.52, 8: 24.76, 7: 16.22, 6: 16.76, 5: 14.94,
            4: 5.05, 3: 3.74, 2: 5.14, 1: 1.9
        }
        
        #for time_slice in range(12, 20):  # Looping from 1 to 22
        pO2 = po2_values[timeslice]
        
        
        # Apply zoom (if scaling factor > 1)
        if scale_factor > 1:
            NPP = zoom(resource_scaled, zoom=(scale_factor, scale_factor))
        else:
            NPP = resource_scaled


        lat = pd.read_csv('lat.csv', header=None)
        #lat = zoom(np.squeeze(lat_df.values), zoom=(scale_factor))

        lon= pd.read_csv('lon.csv', header=None)
        #lon = zoom(np.squeeze(lon_df.values), zoom=(scale_factor))
        
        # Initialize biomass density
        biomass_density = np.zeros(land.shape)
        
        # Run the model (simplified logic)
        grid_size_lon, grid_size_lat = land.shape

        # Initialize population
        population = []
        for _ in range(1000):
            while True:
                x = np.random.randint(0, grid_size_lon)
                y = np.random.randint(0, grid_size_lat)
                if land[x, y] == 1 and NPP[x, y] > 0:
                    break
            population.append({'x': x, 'y': y, 'biomass': np.random.uniform(1, 500), 'age': 0})

        # Prepare figure and animation function
        fig, axs = plt.subplots(1, 4, figsize=(40, 20))
        time_text = fig.text(0.5, 0.95, '', ha='center', va='center', fontsize=16, color='blue')

        total_biomass_over_time = []
        total_population_size_over_time = []

        # Simulation loop
        # Simulation loop
        for t in range(simulation_time_steps):
            new_population = []  # To store new offspring
            surviving_population = []  # To store individuals that survive
            num_births = 0
            num_deaths = 0
            num_moves = 0
            
            # Reset biomass and population density for the current time step
            biomass_density = np.zeros((grid_size_lon, grid_size_lat))
            population_density = np.zeros((grid_size_lon, grid_size_lat))
            
            # Reset survival probability distribution for the current time step
            survival_prob_distribution = np.zeros((grid_size_lon, grid_size_lat))
            
            for individual in population:
                x, y, biomass = individual['x'], individual['y'], individual['biomass']
                local_T = T[x, y]
                R_d_value = R_d(biomass)
                I_d = daily_intake(biomass)
                local_NPP = NPP[x, y]
                
                # Calculate the total biomass in the current cell
                cell_biomass_density = biomass_density[x, y]
                
                # Update carrying capacity considering biomass density
                K_value = K(local_NPP, cell_biomass_density)
                
                # Update survival probability based on new K function
                survival_probability_value = survival_probability(K_value, I_d, local_T, pO2)
                individual['survival_probability'] = survival_probability_value

                # Save survival probability in the matrix
                #survival_prob_distribution[x, y] += survival_probability_value

                # Default new_x and new_y to the current position
                new_x, new_y = x, y

                
                    
                if random.random() < (0.05 if survival_probability_value > 0.75 else 0.1 if survival_probability_value > 0.25 else 0.25):
                    # List of all 8 possible moves: vertical, horizontal, and diagonal
                    possible_moves = [(0, 1), (1, 0), (0, -1), (-1, 0), (1, 1), (1, -1), (-1, 1), (-1, -1)]
                    dx, dy = random.choice(possible_moves)
                    new_x, new_y = (x + dx) % grid_size_lon, (y + dy) % grid_size_lat  # Ensure boundary wrap
                    individual['x'], individual['y'] = new_x, new_y
                    num_moves += 1
                    
                    


                # Calculate new survival after movement
                individual['R_d'] = R_d(biomass)
                I_d = daily_intake(biomass)
                cell_biomass_density = biomass_density[new_x, new_y]
                
                # Recalculate carrying capacity and survival probability after movement
                K_value = K(NPP[new_x, new_y], cell_biomass_density)
                individual['survival_probability'] = survival_probability(K_value, I_d, T[new_x, new_y], pO2)

                
                # Handle death based on the new D formula
                D = 1 - survival_probability_value*r*2
                if random.random() < D:  # Death happens based on D
                    num_deaths += 1
                    biomass_density[x, y] -= biomass
                    population_density[x, y] -= 1
                    # Do not add the individual to the surviving_population
                else:
                    # Individual survives, add to surviving_population
                    surviving_population.append(individual)

                    # Handle reproduction (birth) using the new B formula
                    B = r * (1 - cell_biomass_density / K_value)
                    if random.random() < B:  # Reproduction happens based on B
                        offspring_biomass = individual['biomass']   # Offspring biomass as % of parent's biomass
                        new_population.append({'x': x, 'y': y, 'biomass': offspring_biomass, 'age': 0})
                        num_births += 1

                # Update biomass and population density in grid cell
                biomass_density[new_x, new_y] += biomass
                population_density[new_x, new_y] += 1

            population = surviving_population + new_population

            total_biomass = sum(ind['biomass'] for ind in population)
            total_population_size = len(population)
            total_biomass_over_time.append(total_biomass)
            total_population_size_over_time.append(total_population_size)
            
            #np.savetxt(f'biomass_density_timeslice{timeslice}_run{run}.csv', biomass_density, delimiter=",")
            # Get geological period name
            # Get geological period name or use default "Timeslice_X" if not in dictionary
            period_name = geological_periods.get(timeslice, f"Timeslice_{timeslice}")
            
            # Create a folder for the period if it doesn't exist
            period_folder = f"biomass_outputs/{period_name}"
            os.makedirs(period_folder, exist_ok=True)
            
            # Save the CSV file inside the corresponding folder
            filename = f"{period_folder}/{period_name}_{run}.csv"
            np.savetxt(filename, biomass_density, delimiter=",")


            
            #print(f"Time Step {t}: Births = {num_births}, Deaths = {num_deaths}, Moves = {num_moves}, Biomass = {total_biomass}, Population Size = {len(population)}")
            print(f"Time Step {t}: Biomass = {total_biomass}, Population Size = {len(population)}")
        
#         biomass_per_run.append(np.sum(biomass_density, axis=0))  # Sum over longitude
    
#     biomass_results.append(np.array(biomass_per_run))

# # Compute mean and std deviation across runs
# biomass_means = [np.mean(b, axis=0) for b in biomass_results]
# biomass_stds = [np.std(b, axis=0) for b in biomass_results]

# # Plot results
# plt.figure(figsize=(10, 6))
# lat2 = np.squeeze(lat)
# for i, timeslice in enumerate(range(start_timeslice, end_timeslice + 1)):
#     plt.errorbar(lat2, biomass_means[i], yerr=biomass_stds[i], label=f'Timeslice {timeslice}', alpha=0.7)

# plt.xlabel("Latitude")
# plt.ylabel("Mean Biomass Density")
# plt.title("Biomass Density vs. Latitude Across Timeslices")
# plt.legend()
# plt.grid(True)
# plt.show()
