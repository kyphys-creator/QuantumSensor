import numpy as np

GeV = 1e9
GRAM = 5.62e23 * GeV
CM = (1.98e-14 * GeV) ** (-1)
SEC = (6.58e-25 * GeV) ** (-1)
YR = np.pi * 1e7 * SEC

BINS = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1])
M_DETECTOR = 1.0  # kg
