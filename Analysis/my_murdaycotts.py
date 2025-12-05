import numpy as np
import math 
import matplotlib.pyplot as plt

def my_murdaycotts(Delta, delta, r, D0, bval):
    """
    Function that calculates diffusion attenuation, mlnS = - ln(S/S0), 
    inside a perfectly reflecting sphere according to Murday and Cotts, JCP 1968
    
    Reference value: g = 0.01070 for 40 mT/m

    Here, bardelta = delta/td, parameter of applicability of Neuman's
    approximation:
    for bardelta >> 1, Neuman's limit mlnSneuman can be used, independent of Delta.
    for bardelta << 1, narrow pulse limit mlnSnp can be used.

    (c) Dmitry Novikov, June 2021
    https://github.com/palombom/SANDI-Matlab-Toolbox-Latest-Release/blob/main/functions/support_functions/my_murdaycotts.m

    Args:
        Delta (float): distance between fronts of pulses, in s
        delta (float): pulse width, in s 
        r     (float): radius of the sphere, in m
        D0    (float): free diffusion coefficient, in m²/s
        bval  (float): b-values, in rad² * s/m²

    Returns:
        mlnS       (float): - ln(S/S0), without units
        mlnSneuman (float): - ln(S/S0) with Neuman's approximation (bardelta >> 1), without units
        mlnSnp     (float): - ln(S/S0) with narrow pulse approximation (bardelta << 1), without units
        bardelta   (float): delta / td, without units. td is the diffusion time, in s

    """

    # In s
    if Delta.size > 1:
        Delta[Delta==0] = 1e-9
    elif Delta == 0:
        Delta = 1e-9

    # In s
    if delta.size > 1:
        delta[delta==0] = 1e-9
    elif delta == 0:
        delta = 1e-9


    GAMMA = 2.675987E8 # in rad/(s*T)

    ginput = np.sqrt(bval / (Delta - delta/3)) / (GAMMA * delta) # in T/m 

    g = ginput * GAMMA # in 1/m*s

    td = r**2 / D0 # in s
    bardelta = delta / td # without unit
    barDelta = Delta / td # without unit

    # precompute beta_{1,k}, zeros of x dJ_{3/2}(x)/dx = 1/2 * J_{3/2}(x)
    N = 20; # max # of terms in the sum
    # dJ = @(x) besselj(3/2,x) - x.*(besselj(1/2,x)-besselj(5/2,x));
    # beta = @[k] fzero(dJ, [(k-1)*pi+eps, k*pi]); 
    # b = zeros(1,N); for k=1:N, b[k] = beta[k]; end 
    # One can now tabulate b[k] and never calculate them again: 
    b = [2.0816, 5.9404, 9.2058, 12.4044, 15.5792, 18.7426, 21.8997, 25.0528, 28.2034, 31.3521, 
         34.4995, 37.6460, 40.7917, 43.9368, 47.0814, 50.2257, 53.3696, 56.5133, 59.6567, 62.8000]

    mlnS = 0 
    for k in range(N):
        if bardelta.size == 1:
            mlnS = mlnS + (2 / (b[k]**6 * (b[k]**2 - 2))) * (-2 + 2*b[k]**2 * bardelta + \
                   2 * (math.exp(-b[k]**2 * bardelta) + math.exp(-b[k]**2 * barDelta)) - \
                   math.exp(-b[k]**2 * (bardelta + barDelta)) - math.exp(-b[k]**2 * (barDelta - bardelta)))        
        else:
            mlnS = mlnS + (2 / (b[k]**6 * (b[k]**2 - 2)))*(-2 + 2 * b[k]**2 * bardelta[k] + \
                   2 * (math.exp(-b[k]**2 * bardelta[k]) + math.exp(-b[k]**2 * barDelta[k])) - \
                    math.exp(-b[k]**2 * (bardelta[k] + barDelta[k])) - math.exp(-b[k]**2 * (barDelta[k] - bardelta[k])))
                

    mlnS = mlnS * D0 * g**2 * td**3

    mlnSneuman = 16 / 175 * g**2 * delta * r**4 / D0

    mlnSnp = (g * delta)**2 * r**2 / 5

    return mlnS, mlnSneuman, mlnSnp, bardelta, b


G         = np.array([0, 0.015, 0.034, 0.048, 0.059]) # in T/m
Delta     = np.array([0.05]*G.size)  # in s
delta     = np.array([0.0165]*G.size)# in s
TE        = np.array([0.067]*G.size) # in s
r_soma    = 10e-6  # m
r_neurite = 0.5e-6 # m
D0        = 2.5e-9 # m²/s
gamma     = 2.6751525e8 # rad/(s*T)
bval      = gamma**2 * G**2 * delta**2 * (Delta - delta/3) # rad² * s/m²

# Calculate the volume fraction of neurites & soma, needed for the signal weighted sum
nb_neurites     = 20
l_neurite       = 245e-6 # m
volume_neurites = nb_neurites * np.pi*r_neurite**2*l_neurite # in m³
volume_soma     = 4/3 * np.pi * r_soma**3 # in m³
volume_neuron   = volume_neurites + volume_soma
neurite_fraction= volume_neurites / volume_neuron
soma_fraction   = volume_soma / volume_neuron

# Calculate the signal from the soma only, neurites only or combination of both at different b-values
soma_signal     = []
neurites_signal = []
both_signal     = []
for i in range(bval.size):
    # Soma signal:
    mlnS, mlnSneuman, mlnSnp, bardelta, b = my_murdaycotts(Delta[i], delta[i], r_soma, D0, bval[i])
    # Store S/S0
    soma_signal.append(math.exp(-mlnS))

    # Signal intra sticks:
    # If b = 0, put the signal to 1 (no diffusion => no attenuation)
    if bval[i] == 0:
        neurites_signal.append(1)
        both_signal.append(1)
    else:
        # Palombo et al, 2020 (SANDI)
        neurite_signal = np.sqrt(np.pi/(4 * bval[i] * D0)) * math.erf(np.sqrt(bval[i] * D0))
        neurites_signal.append(neurite_signal)

        # The weighted sum of soma & neurites signals
        both_signal.append(neurite_fraction * neurite_signal + soma_fraction * math.exp(-mlnS))


# plt.plot(bval, soma_signal, label='sphere')
# plt.plot(bval, neurites_signal, label='neurites')
# plt.plot(bval, both_signal, label='sphere & neurites')
# plt.ylabel('S/S0')
# plt.xlabel('b [m/s]')
# plt.title('Analytical solution: simple neuron')
# plt.legend()
# plt.show()



