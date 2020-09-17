## This is my calculation attempt for homework 1 assignment 2
# Matthew Jackson
# August 31, 2020
# Phys 513

import numpy as np
import matplotlib.pyplot as pyplot

"""
The exact solution to an electric field from a line charge of length 2L about the origin is

E = (2 * lam * L) /  
    (4 * pi * eps * y (y**2 + L**2)**(1/2))
Source is from example 2.2 on page 64 of Griffths E&M
    
Griffiths, D. (2017). Electrostatics. In Introduction to Electrodynamics (pp. 59-112).
 Cambridge: Cambridge University Press. doi:10.1017/9781108333511.004 
"""

# I am going to normalize all my constants here
lam = 0.5  # This is the charge per unit length
# I set this to 0.5 because q = lam * L_tot where L_tot = 2L, so I am setting lam = 0.5 to normalize q
y = 1  # This is the range from the wire
L = 1  # This is the length of the wire
k = 1  # This would be 1 / (4 * pi * epsilon_naught)
# This E field is only in the y direction. For that reason, I am going to assert that my E_x values all sum to near 0
E_exact = (k * 2 * lam * L) / (y * (y ** 2 + L ** 2) ** (1 / 2))


# Define a point charge calculation
def point_E_Field(q, charge_loc, field_loc):
    """
    This is calculating the electric field from a given point charge. This is done so that I can calculate the E field
    contribution from each point charge. I am doing all my calculations in the x-y plane

    This is the equation I am calculating

    E = k * q * ((x - x_prime) * x_hat + (y - y_prime) * y_hat) /  \
             (((x - x_prime) ** 2 + (y - y_prime) ** 2) ** (3 / 2))

    :param pos:        this is the fraction of the total charge at a given point
    :param charge_loc: this is the location of my charge. These are the primed values in my equation above. The inputs
                       for this are normalized, so they should be defined between -1 and 1
    :param field_loc:  this is the location of my field. These are the regular values in my equation above. These values
                       are always going to be [0, 1]
    :return:           This is the x and y components of the E field
    """
    k = 1
    E = np.empty((2,))
    x_prime, y_prime = charge_loc
    x, y = field_loc
    # assert x == 0
    # assert y == 1
    E_x = k * q * (x - x_prime) / \
          (((x - x_prime) ** 2 + (y - y_prime) ** 2) ** (3 / 2))
    E_y = k * q * (y - y_prime) / \
          (((x - x_prime) ** 2 + (y - y_prime) ** 2) ** (3 / 2))
    E[0] = E_x
    E[1] = E_y
    return E


# I am going to start iterating my solution from two point charges of distance 2 * L
num_charges = 1
error = []
error_2 = []
converged = False

ii_count = 0

while not converged and ii_count < 25:
    charge_locs = np.zeros((num_charges, 2))  # Initialize my points
    charge_locs_2 = np.zeros((num_charges, 2))  # I am changing the point scheme that I am working with
                                                # this was done as a recommendation to smooth my graph
    if num_charges == 1:
        charge_locs[:, 1] = 0
        charge_locs_2[:, 1] = 0
    else:
        charge_locs[:, 0] = np.linspace(-1, 1, num_charges)
        charge_locs_2[:, 0] = np.arange(-num_charges+1, num_charges, 2) * (1 / (num_charges + 1))
    E = np.zeros((2,))
    E_2 = np.zeros((2,))
    for loc, loc_2 in zip(charge_locs, charge_locs_2):
        E += point_E_Field(q=(1 / num_charges), charge_loc=loc, field_loc=np.array([0, 1]))
        E_2 += point_E_Field(q=(1 / num_charges), charge_loc=loc_2, field_loc=np.array([0, 1]))
    # First, I want to make sure E_x is close to 0
    # assert np.abs(E[0]) < 1e-5
    error.append((np.abs(E_exact - E[1]) / E_exact) * 100)
    error_2.append((np.abs(E_exact - E_2[1]) / E_exact) * 100)  # For some reason I am getting half the value here
    converged = True if error[-1] < 0.01 else False # I want to 1% to make the plot look better
    ii_count += 1
    num_charges += 1

# Good that you checked both options! There is an ambiguity that will matter
# for small N.
# Great plots.
# For a plot like this, you may want to show dots and a line.
print(num_charges-1) # Shown answer .... (I'd put in legend.)
error = np.asarray(error)
pyplot.figure(1)
pyplot.loglog(np.arange(1, error.size + 1), error)
e_ind = (error > 10).nonzero()[0].size
pyplot.loglog(np.arange(1, e_ind+2), error[:e_ind+1])
pyplot.loglog(e_ind+1, error[e_ind], '*')
pyplot.title('Error Between Exact and Approximation\nUsing fixed end points')
pyplot.ylabel('Error (%)')
pyplot.xlabel('Number of charges')
pyplot.legend(['>1% Error Threshold', '>10% Error Threshold', 'Error {:.2f}'.format(error[e_ind])])
error1 = error

error = np.asarray(error_2)
pyplot.figure(2)
pyplot.loglog(np.arange(1, error.size + 1), error)
e_ind = (error > 10).nonzero()[0].size
pyplot.loglog(np.arange(1, e_ind+2), error[:e_ind+1])
pyplot.loglog(e_ind+1, error[e_ind], '*')
pyplot.grid()
pyplot.title('Error Between Exact and Approximation\nUsing deltas from end points')
pyplot.ylabel('Error (%)')
pyplot.xlabel('Number of charges')
pyplot.legend(['>1% Error Threshold', '>10% Error Threshold', 'Error {:.2f}'.format(error[e_ind])])
error2 = error
