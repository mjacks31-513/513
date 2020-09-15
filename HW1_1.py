## This is my plotting function homework 1 assignment 1
# Matthew Jackson
# August 31, 2020
# Phys 513

# As a side note, this was a really fun little challenge. I have never really thought about doing a solver that operates
# on the position of something, since my background is in finite difference on static grids

import numpy as np
import scipy.integrate as spint
import matplotlib.pyplot as pyplot

# I am defining this on a one dimensional scale
pos_charge = np.array([1,0])
neg_charge = np.array([-1,0])

# I am assuming that k and q are normalized to 1
def E_Field(pos, charge_loc, field_loc):
    """
    This is calculating the electric field from a given charge. This is done so that I can do two calculations: one
    from the positive charge and one from the negative charge. I am doing all my calculations in the x-y plane, so I
    will assume that z_prime == z for all my calculations.

    This is the equation I am calculating

    E = k*q* ((x-x_prime)*x_hat + (y-y_prime)*y_hat)/  \
             (((x-x_prime)**2 + (y-y_prime)**2)**(3/2))

    :param pos:        this is a simple true or false for positive or negative
    :param charge_loc: this is the location of my charge. These are the primed values in my equation above. The inputs
                       for this are normalized, so they should always be np.array([1,0]) and np.array([-1,0])
    :param field_loc:  this is the location of my field. These are the regular values in my equation above. This is all
                       normalized to charge locations
    :return:           This is the x and y components of the E field

    """
    # all values are normalized
    q = 1 if pos else -1
    k = 1
    E = np.empty((2,))
    x_prime, y_prime = charge_loc
    x, y = field_loc
    E_x = k*q*(x-x_prime) /  \
             (((x-x_prime)**2 + (y-y_prime)**2)**(3/2))
    E_y = k*q*(y-y_prime) /  \
             (((x-x_prime)**2 + (y-y_prime)**2)**(3/2))
    E[0] = E_x
    E[1] = E_y
    return E


def propagate(s,X):
    """
    This

    :param s: This is the range that I am changing my field by. I am not using this currently
    :param X: This the position of my field that I am calculating
    :return: I am returning E_x_hat and E_y_hat
    """
    # s is not used in this calculation. Might need to check this
    E_1 = E_Field(pos=True, charge_loc=pos_charge, field_loc=X)
    E_2 = E_Field(pos=False, charge_loc=neg_charge, field_loc=X)
    E = E_1 + E_2
    return E/np.linalg.norm(E)


ds = 0.05
t_end = 5


# I didn't want to kill my computation, so I am only taking the angles that
# are facing towards my negative charge in
# increments of 15 degrees
angles = np.arange(90, 271, 15)

for angle in angles:
    x_0 = np.cos(np.deg2rad(angle)) / 20 + 1
    y_0 = np.sin(np.deg2rad(angle)) / 20

    x = []
    y = []
    x.append(x_0)
    y.append(y_0)

    pyplot.figure(1)
    ii_count = 0
    # I know where I want my solution to converge
    while np.linalg.norm([x[-1] + 1, y[-1]]) > 0.05 and ii_count < 100:
        # The next generalization is to stop when line steps out of box of a
        # given size.
        dxyds = propagate(0, [x[-1], y[-1]])
        x.append(x[-1] + ds * dxyds[0])
        y.append(y[-1] + ds * dxyds[1])
        ii_count += 1
    pyplot.plot(x, y)

    pyplot.figure(2)
    # FYI: solve_ivp has "events" function that stop integration when an 
    # event is true.
    output = (spint.solve_ivp(propagate, [0,t_end], [x_0, y_0], t_eval=np.arange(0, t_end, ds)))
    pyplot.plot(output.y[0, :], output.y[1, :])

pyplot.figure(1)
pyplot.title('Forward-Euler solver')
pyplot.xlabel('x/d') # More common syntax
pyplot.ylabel('y/d')
pyplot.ylim([-2, 2])
pyplot.xlim([-2, 2])

pyplot.figure(2)
pyplot.title('Runge-Kutta 45 solver')
pyplot.xlabel('x/d')
pyplot.ylabel('y/d')
pyplot.ylim([-2, 2])
pyplot.xlim([-2, 2])

pyplot.show()
