import simulation_functions as sim

#############################
###          INIT         ###
#############################

## physical constants ##
height = 3  # height of the simulation space, m
T_0 = 25.0  # Temperature at t = 0, °C
T_air = 25.0  # Temperature of the surrounding air, °C
h_conv = 10.0  # Convective heat transfer coefficient, W/(m^2 K)
Q = 100  # Total heat generation rate, W, artificially high for troubleshooting
loading = 1e-6  # mass fraction of CB in PDMS, g/g

PDMS_thermal_conductivity_WpmK = 0.2 + loading
PDMS_density_gpmL = 1.02
PDMS_heat_capacity_JpgK = 1.67
PDMS_thermal_diffusivity_m2ps = PDMS_thermal_conductivity_WpmK / (PDMS_density_gpmL * PDMS_heat_capacity_JpgK)

beam_radius_m = 0.5
abs_coeff = 1 + (loading * 2300) # abs lerps between 1 and ~230 over the loading range of 0% to 10%

## simulation parameters ##
Nx = Ny = 50
dx = dy = height / (Nx - 1)
dt = dx**2 / (PDMS_thermal_diffusivity_m2ps * 4)  # time step, s 
x = y = sim.np.linspace(0, height, Nx)
X, Y = sim.np.meshgrid(x, y)

beam_mask = X <= beam_radius_m
q = sim.np.zeros((Nx, Ny))
Q_abs = Q / (sim.np.pi * beam_radius_m**2 * height)  # before implementing exponential decay
q[beam_mask] = Q / Q_abs / dx**3
q = sim.set_up_volumetric_power_distribution(q, beam_mask, abs_coeff, dx, dy, height)
transmittance, absorbed_power = sim.Preview_Decay(q, dx, dy, Q, height)
output_times = [0, 1, 3, 5]

#############################
### MAIN CODE STARTS HERE ###
#############################

print(f'dt: {dt:2g}, transmittance: {transmittance:3g}, absorbed power: {absorbed_power:3g}')
output_temperatures_RK = sim.Compute_T(output_times, Nx, Ny, T_0, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q, h_conv, T_air)
sim.Plot_T_Slices(X, Y, output_temperatures_RK, output_times, height)
# save_output_temperatures(output_temperatures_RK, "output_temperatures.npz")