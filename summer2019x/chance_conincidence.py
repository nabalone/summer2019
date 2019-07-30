import numpy as np

def calculate_coincidence(separation, size, magnitude):
    '''
    Calculate the chance that a galaxy of size R_h and magnitude M falls
    within a separation R of a transient. The galaxies with the lowest
    chance probability will be selected as the best candidate hosts.

    Parameters
    ---------------
    separation: Separation between the host and transient [Arcseconds]
    size: Half light radius of the galaxy [Arcseconds]
    Magnitude: Magnitude of the galaxy

    Output
    ---------------
    P_cc = Probability of chance coincidence
    '''

    # Observed number density of galaxies brighter than magnitude M (From Berger 2010)
    sigma = 10 ** (0.33 * (magnitude - 24) - 2.44) / (0.33 * np.log(10))

    # Effective radius
    #R_effective = np.sqrt(4 * separation ** 2 + size ** 2)
    R_effective = np.sqrt(np.abs(separation) + 4 * np.abs(size) ** 2)

    # Probability of chance coincidence
    chance_coincidence = 1 - np.exp(-np.pi * R_effective ** 2 * sigma)

    return chance_coincidence
