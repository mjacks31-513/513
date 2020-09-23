import numpy as np

from HW4_BiotSavart import HW4_BiotSavart
import numpy as np
from matplotlib import pyplot

if __name__ == '__main__':
    '''
    I am trying to get the flux through a ring
    '''
    # I am initializing at 0 and then filling in x and y
    lineStart = 0
    lineEnd = 2
    lineDiff = 0.01
    size = int(lineEnd / lineDiff)
    dims = 3
    x_ind = 0
    y_ind = 1
    z_ind = 2

    x = np.cos(np.arange(lineStart, lineEnd, lineDiff) * np.pi)
    y = np.sin(np.arange(lineStart, lineEnd, lineDiff) * np.pi)
    XYZ_line_1 = np.zeros((size, dims))
    XYZ_line_1[:, x_ind] = x
    XYZ_line_1[:, y_ind] = y

    # I am initializing at 2 and then filling in x and y
    ringHt = 2
    XYZ_line_2 = np.full_like(XYZ_line_1, ringHt)
    XYZ_line_2[:, x_ind] = x
    XYZ_line_2[:, y_ind] = y

    I_line_1 = I_line_2 = 1

    '''
    Magnetic flux strategy:
    Since both of my current loops are symmetric on the Z axis, I can use an azimuthal symmetry to
    remove the need to discretize the space inside my current loop into little squares (my first thought)
    I want might field components to have an equal area so I am going to use the fact that my area 
    scales with r^2 to find the correct radii for my points. I am going to build the area sections first 
    and let that inform where I should calculate the B field
    '''
    # I am setting a static area of 0.04
    startNumSections = 2
    endNumSections = 11
    areas = 1 / (np.arange(startNumSections, endNumSections) ** 2)
    magneticFluxes = []
    for sectionArea in areas:
        endRadiusSquared = 1
        startRadiusSquared = 0
        numRadii = int(endRadiusSquared / sectionArea)
        edgeRadii = np.sqrt(np.linspace(startRadiusSquared, endRadiusSquared, numRadii + 1))
        '''
        I want my x values to be in the "middle" of the radii values declared above. The problem is, using 
        the midpoint is actually going to give me incorrect values because the flux is the strength of the 
        B field times the size of the area that I am getting that field over, so I want the "middle area"
        radius, so I am going to take the average area, and then get the b field at that point
        '''
        cumulativeRadiiSquared = np.zeros(numRadii + 2)
        cumulativeRadiiSquared[1:] = np.cumsum(np.linspace(startRadiusSquared, endRadiusSquared, numRadii + 1))
        x = np.sqrt((cumulativeRadiiSquared[2:] - cumulativeRadiiSquared[:-2]) / 2)  # This is a moving average
        '''
        I am testing out other radius schemes
        '''
        # cumulativeRadii = np.zeros(numRadii + 2)
        # cumulativeRadii[1:] = edgeRadii
        # x = (cumulativeRadiiSquared[2:] - cumulativeRadiiSquared[:-2]) / 2  # This is a moving average
        #
        # x = np.sqrt(edgeRadii[1:]*edgeRadii[:-1])

        y = np.zeros_like(numRadii)
        points = np.full((numRadii, dims), ringHt, dtype=np.float64)
        points[:, x_ind] = x  # This is not working?
        points[:, y_ind] = y

        B_line_1 = HW4_BiotSavart(points, XYZ_line_1, I_line_1)
        # I reread the assignment, and I don't need this line
        # B_line_2 = HW4_BiotSavart(points, XYZ_line_2, I_line_2)

        B_total = B_line_1  # + B_line_2

        '''
        This is where the magic happens. I am going to calculate the area of each nested circle and then 
        dot that into the total B field. Since the normal vector for area is only in the z direction, I
        can ignore the x and y components to the B field. 
        '''
        z = np.full_like(x, sectionArea)
        areaVector = np.zeros_like(B_total)
        areaVector[:, z_ind] = z

        magneticFluxes.append(np.pi * (areaVector * B_total).sum())

    pyplot.plot(np.arange(startNumSections, endNumSections), magneticFluxes)
