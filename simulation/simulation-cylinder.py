import numpy as np
import matplotlib.pyplot as plt

def is_within_beam(x, y, r):
    """Check if a point is within the beam radius."""
    return (x - r)**2 + (y - r)**2 <= r**2

def compute_next_temperature(T, q_xyz, delta_x, delta_tau, alpha, h, k, T_inf, r):
    """Compute the temperature for the next time step using finite difference."""
    T_next = T.copy()
    M = (delta_x**2) / (alpha * delta_tau)
    
    for i in range(1, T.shape[0] - 1):
        for j in range(1, T.shape[1] - 1):
            T_next[i, j] = (alpha * delta_tau / delta_x**2) * (T[i+1, j] + T[i-1, j] + T[i, j+1] + T[i, j-1]) \
                         + (1 - 4*alpha*delta_tau/delta_x**2) * T[i, j]
            
            if is_within_beam(i*delta_x, j*delta_y, r):
                T_next[i, j] += (q_xyz[i, j] / (c_PDMS * rho_PDMS)) * delta_tau
    
    # Adiabatic boundary conditions
    T_next[0, :] = T_next[1, :]
    T_next[-1, :] = T_next[-2, :]
    T_next[:, 0] = T_next[:, 1]
    T_next[:, -1] = T_next[:, -2]
    
    return T_next

# Material properties and parameters from README
r = 1.5e-2
alpha_PDMS = 0.133e-6
k_PDMS = 0.2
rho_PDMS = 1.03e3
c_PDMS = 1.46e3
h = 10
T_inf = 25
alpha_abs = 0.1e3
P_high = 600  # Updated to the higher power
Q0_new = P_high / (np.pi * r**2)

# Simulation parameters
delta_x = r/20
delta_y = delta_x
delta_tau = 0.01
end_time = 2
time_steps = int(end_time / delta_tau)

# Initialize 2D grid
x_points = int(2*r / delta_x)
y_points = x_points
T = np.full((x_points, y_points), T_inf)

# Compute q_xyz for each node with the new power
q_xyz_new = np.zeros((x_points, y_points))
for i in range(x_points):
    for j in range(y_points):
        y_start = j * delta_y
        y_end = y_start + delta_y
        integral_value = Q0_new * (np.exp(-alpha_abs * y_start) - np.exp(-alpha_abs * y_end)) / alpha_abs
        q_xyz_new[i, j] = integral_value * delta_x * delta_z / (np.pi * r**2)

# Time-stepping
T_snapshots_new = [T.copy()]
for t in range(time_steps):
    T = compute_next_temperature(T, q_xyz_new, delta_x, delta_tau, alpha_PDMS, h, k_PDMS, T_inf, r)
    if (t+1) * delta_tau in [1, 2]:
        T_snapshots_new.append(T.copy())

# Plotting the results
fig, axs = plt.subplots(1, 3, figsize=(15, 5))
times = [0, 1, 2]
for i, ax in enumerate(axs):
    c = ax.imshow(T_snapshots_new[i], cmap='hot', origin='lower', extent=[0, 2*r, 0, 2*r])
    ax.set_title(f"Temperature Distribution at t = {times[i]} s (Higher Power)")
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    fig.colorbar(c, ax=ax, label='Temperature (°C)')
plt.tight_layout()
plt.show()
