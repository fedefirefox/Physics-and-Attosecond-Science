import numpy as np
import matplotlib.pyplot as plt

# Parameters
N = 512          # number of spatial grid points
x_min = -10.0
x_max = 10.0
x = np.linspace(x_min, x_max, N)

dx = x[1] - x[0]

dt = 0.01        # time step
steps = 500      # number of time steps

hbar = 1.0
mass = 1.0
omega = 1.0

# Harmonic potential V(x) = 0.5*m*omega^2*x^2
V = 0.5 * mass * omega**2 * x**2

# Initial wavefunction: gaussian packet
sigma = 1.0
psi = (1/(np.pi*sigma**2)**0.25) * np.exp(-x**2/(2*sigma**2))

# Precompute constants for the Crank-Nicolson scheme
# interior points exclude boundaries assuming psi=0 at edges
x_int = x[1:-1]
V_int = V[1:-1]
Nint = N - 2

lap_const = 1.0/dx**2

H_diag = lap_const + V_int
H_off = -0.5 * lap_const

A_diag = 1.0 + 0.5j * dt/hbar * H_diag
A_off = 0.5j * dt/hbar * H_off
B_diag = 1.0 - 0.5j * dt/hbar * H_diag
B_off = -0.5j * dt/hbar * H_off

# Construct banded matrices
A = np.diag(A_diag) + np.diag(A_off*np.ones(Nint-1), 1) + np.diag(A_off*np.ones(Nint-1), -1)
B = np.diag(B_diag) + np.diag(B_off*np.ones(Nint-1), 1) + np.diag(B_off*np.ones(Nint-1), -1)

# Time propagation
for _ in range(steps):
    psi_int = psi[1:-1]
    rhs = B @ psi_int
    psi_next = np.linalg.solve(A, rhs)
    psi[1:-1] = psi_next
    psi[0] = psi[-1] = 0.0

# Probability density
prob_density = np.abs(psi)**2

plt.plot(x, prob_density)
plt.xlabel('x')
plt.ylabel('Probability density')
plt.title('1D TDSE with harmonic potential after {} steps'.format(steps))
plt.show()
