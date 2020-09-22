# Matthew Jackson
# PHYS 513
# September 20, 2020

import numpy as np
import scipy as sp
from matplotlib import pyplot


def HW4_BiotSavart(p=None, XYZ=None, I=None, closed=True):
    """
    This code will calculate the magnetic field from the current, I,
     going through a current with geometry XYZ, at point p. XYZ and p
     must be 3 dimensional.

    Biot-Savart law is defined by B = mu * I / ( 4 * pi ) * integral{ dl CROSS R / ( R ^ 3 ) }

    :param p: The point to evaluate the B field at. This must be a three dimensional point
    :param XYZ: The segments of wire that are being integrated over. These must be three dimensional points
    :param I: The current in the wire given by XYZ
    :return: The components to the magnetic field that are being evaluated at point p
    """
    assert p.size % 3 == 0  # Normally I shouldn't do assert, but whatever
    assert XYZ.size % 3 == 0

    # I am normalizing my constants here. They can be backed out later if I want them to be
    mu_over_4pi = 1

    if closed:
        iterationStrategy = np.empty((XYZ.shape[0], 2, XYZ.shape[1]))
        iterationStrategy[:, 0] = XYZ
        iterationStrategy[:-1, 1] = XYZ[1:, :]
        iterationStrategy[-1, 1] = XYZ[0, :]  # This is closing the loop
    else:
        iterationStrategy = np.empty((XYZ.shape[0] - 1, 2, XYZ.shape[1]))
        iterationStrategy[:, 0] = XYZ[:-1, :]
        iterationStrategy[:, 1] = XYZ[1:, :]

    if p.ndim == 1:
        p = p[np.newaxis, :]

    B = np.empty_like(p, dtype=np.float64)

    for iii, point in enumerate(p):
        B_point = np.zeros_like(point)
        for dL1, dL2 in iterationStrategy:
            '''
            My implementation is working with the midpoints of the segments provided to create the dL 
            segments I need to get a cross product with R. As the segments get smaller, that approximation
            gets more accurate. This was the only approach I could think of that would allow me to 
            approximate the dL vector with any accuracy
            '''
            midpoint = (dL1 + dL2) / 2
            dL = dL2 - dL1
            R = point - midpoint
            R_mag = np.power(np.power(R, 2).sum(), 0.5)
            dL_cross_R = np.cross(dL, R)  # This going from d1 to d2
            B_point = B_point + dL_cross_R / np.power(R_mag, 3)

        B[iii] = mu_over_4pi * B_point

    return B


if __name__ == '__main__':
    '''
    This is a code check that was laid out in the Homework
    '''
    p = np.array((0, 0, 0))
    XYZ = np.array([1, 1, 0, -1, 1, 0, -1, -1, 0, 1, -1, 0]).reshape(4, 3)
    I = 1
    B = HW4_BiotSavart(p, XYZ, I)

