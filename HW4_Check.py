# Matthew Jackson
# PHYS 513
# September 20, 2020

from matplotlib import pyplot
import numpy as np
from HW4_BiotSavart import HW4_BiotSavart

if __name__ == "__main__":
    """
    I am going to evaluate the mag field from a line charge. This check seems to work well
    """
    point = np.array((0, 1, 0))
    allB_1 = []
    lengths = np.arange(0.2, 10.2, 0.2)
    normalizedCurrent = 1
    for halfLength in lengths:
        start = -1 * halfLength
        end = halfLength
        segments = int(halfLength * 10 + 1)
        line = np.zeros((segments, 3))
        line[:, 0] = np.linspace(start, end, segments)
        allB_1.append(HW4_BiotSavart(point, line, normalizedCurrent, False))
    B_1 = np.asarray(allB_1).squeeze()
    pyplot.figure(1)
    pyplot.plot(2 * lengths, B_1[:, 2], 'ro-')
    # Analytical solution is (mu * I) / (2 * pi * R)
    pyplot.plot([0, 20], [2] * 2, 'k--')
    pyplot.title('Magnetic Field From a Long Straight Wire')
    pyplot.xlabel('Ratio of Wire Length to Distance')
    pyplot.ylabel('Magnetic Field Strength\n(normalized to mu/(4*pi))')
    pyplot.legend(['Solver', 'Analytical'])

    points = np.zeros((36, 3))
    startAngle = 0.0
    endAngle = 2.0
    angleDiff = endAngle / 12.0
    for iii in range(3):
        points[(12 * iii):((iii + 1) * 12), 1] = (iii + 1) * np.cos(np.arange(startAngle, endAngle,
                                                                              angleDiff) * np.pi)
        points[(12 * iii):((iii + 1) * 12), 2] = (iii + 1) * np.sin(np.arange(startAngle, endAngle,
                                                                              angleDiff) * np.pi)
    B_3D_1 = HW4_BiotSavart(points, line, normalizedCurrent)
    B_mag = np.sqrt(B_3D_1[:, 0] ** 2 + B_3D_1[:, 1] ** 2 + B_3D_1[:, 2] ** 2)
    fig = pyplot.figure(2)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(*line.T)
    ax.quiver(*points.T, *B_3D_1.T, length=0.3, pivot='middle', normalize=True,
              linewidths=3 * (B_mag / B_mag.max()))
    pyplot.title('Magnetic Field Around a Wire')

    lineStart = 0
    lineEnd = 1
    numberOfSegments = np.array(list(range(3, 48)))
    point = np.array((0, 0, 0))
    allB_2 = []
    for size in numberOfSegments:
        lineDiff = 1 / size
        XYZ = np.zeros((size, 3))
        XYZ[:, 0] = np.cos(np.arange(lineStart, lineEnd, lineDiff) * 2 * np.pi)
        XYZ[:, 1] = np.sin(np.arange(lineStart, lineEnd, lineDiff) * 2 * np.pi)
        allB_2.append(HW4_BiotSavart(point, XYZ, normalizedCurrent))
    B_2 = np.asarray(allB_2).squeeze()
    pyplot.figure(3)
    pyplot.plot(numberOfSegments, B_2[:, 2], 'ro-')
    # Analytical solution is (mu * I) / (2 * R)
    pyplot.plot([0, 50], [2 * np.pi] * 2, 'k--')
    pyplot.title('Magnetic Field From a Loop of Current')
    pyplot.xlabel('Number of Line Segments')
    pyplot.ylabel('Magnetic Field Strength\n(normalized to mu/(4*pi))')
    pyplot.legend(['Solver', 'Analytical'])
