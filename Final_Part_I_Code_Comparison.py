import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

global SCALING_FACTOR_ON
SCALING_FACTOR_ON = False

def ladderCircuit(w, L, C, N, Z_L, V_0):
    Z_total = np.empty(N, dtype=complex)  # I am going to test this
    I_total = np.empty(N, dtype=complex)
    V_total = np.empty(N, dtype=complex)

    Z_total[-1] = Z_L
    for n in range(N-2, -1, -1):
        Z_total[n] = (1j*w*C + (Z_total[n+1]+1j*w*L)**-1)**-1

    # I am not including V_0 in my V_total.
    # The reason I am doing this is because it makes my loop look nicer.
    # I am going to scale my V_0
    if SCALING_FACTOR_ON:
        rho = (Z_total[0]-1)/(Z_total[0]+1)
        coeff = np.exp(2j*np.pi * N * (w*np.sqrt(L*C))**-1)
        V_0 = V_0 * coeff * (1+rho)
    I_total[0] = V_0/Z_total[0]
    V_total[0] = V_0 - 1j*w*L*I_total[0]
    for n in range(0, N-1):
        I_total[n+1] = I_total[n]-1j*w*C*V_total[n]
        V_total[n+1] = V_total[n]-1j*w*L*I_total[n+1]

    return {'V': V_total, 'I': I_total, 'Z': Z_total}


if __name__ == '__main__':
    ## These are the inputs for my ladder circuit
    V_0 = (1+0j)  # V
    L = 1  # H
    C = 1  # F
    w = 0.005  # s**-1
    N = 1000  #
    Z_L = (5*(L/C) ** 0.5 + 0j)

    ## These are inputs for my analytical solution
    x = np.arange(N)
    lam = 2 * np.pi / w
    # Z_0 = (L/C)**5
    Z_0 = ((L/C) - (w ** 2 * L ** 2 / 4)) ** 0.5
    rho_0 = (Z_L - Z_0)/(Z_L + Z_0)


    ## animate 1 = animation window
    ## animate 2 = 2 plots at 0 and Pi/2
    ## animate 3 = 1000 iterations VSWR
    animate = 1

    # Z_L = np.exp(1j * np.pi / 4)
    output = ladderCircuit(w, L, C, N, Z_L, V_0)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1)

    global domainData
    domainData = np.empty((N, 3))
    for jjj in range(3):
        domainData[:, jjj] = np.arange(N)

    # x = np.arange(0, 2 * np.pi, 0.01)
    ladderData = np.empty((N, 3))
    ladderData[:, 0] = np.real(output['V'])
    ladderData[:, 1] = np.real(output['I'])
    ladderData[:, 2] = np.real(output['Z'])
    line_V, = ax1.plot(np.real(output['V']), label='Voltage')
    line_I, = ax1.plot(np.real(output['I']), label='Current')
    line_Z, = ax1.plot(np.real(output['Z']), label='Impedance')
    # ladderPlots, = ax1.plot([], [],'-')
    # ax1.legend(['Voltage', 'Current', 'Impedance'])
    ax1.legend()
    plt.autoscale()

    global analyticalData
    analyticalData = np.zeros((N, 3), dtype=complex)
    V_i, I_i, Z_i = analyticalData.T
    V_plot, = ax2.plot(np.real(V_i), label='Current')
    I_plot, = ax2.plot(np.real(I_i), label='Current')
    Z_plot, = ax2.plot(np.real(Z_i), label='Impedance')
    # analyticalPlots, = ax2.plot([], [],'-')
    # plt.autoscale()
    # ax2.legend(['Voltage', 'Current', 'Impedance'])
    ax2.legend()

    differenceData = ladderData - np.real(analyticalData)
    V_diff, = ax3.plot(np.real(output['V'] - V_i), label='V_{Ladder}-V_{Analytical}')
    I_diff, = ax3.plot(np.real(output['I'] - I_i), label='I_{Ladder}-I_{Analytical}')
    Z_diff, = ax3.plot(np.real(output['Z'] - Z_i), label='Z_{Ladder}-Z_{Analytical}')
    # differencePlots, = ax3.plot([], [], '-')
    # ax3.legend(['V Diff', 'I Diff', 'Z Diff'])
    ax3.legend()

    ax1.set_title('Ladder')
    ax2.set_title('Analytical')
    ax3.set_title('Difference')

    # lines = [ladderPlots, analyticalPlots, differencePlots]

    def animateFunc(iii, V_t=None):
        if V_t is None:
            # I am interested in giving this a try, but I am not sure it'll work
            waveMod = 0.05 / w  # This is a modifier to increase speed
            V_t = np.exp(1j * w * iii * waveMod)  #The extra 50 is to make the solution animate faster
        newYdata = ladderCircuit(w, L, C, N, Z_L, V_t)
        ladderData[:, 0] = np.real(newYdata['V'])
        ladderData[:, 1] = np.real(newYdata['I'])
        ladderData[:, 2] = np.real(newYdata['Z'])

        rho = rho_0 * np.exp(4j * np.pi * (x - N + 1) / lam)

        # I need to maintain the boundary condition
        V_c = V_t / ((1 + rho[0]))
        V_i = V_c * np.exp(-2j * np.pi * x / lam) * (1 + rho)
        I_i = V_c / Z_0 * np.exp(-2j * np.pi * x / lam) * (1 - rho)
        Z_i = V_i / I_i

        analyticalData[:, 0] = V_i
        analyticalData[:, 1] = I_i
        analyticalData[:, 2] = Z_i

        differenceData = ladderData - np.real(analyticalData)

        if animate == 1:
            line_V.set_ydata(np.real(newYdata['V']))  # update the data.
            line_I.set_ydata(np.real(newYdata['I']))  # update the data.
            line_Z.set_ydata(np.real(newYdata['Z']))  # update the data.

            V_plot.set_ydata(np.real(V_i))
            I_plot.set_ydata(np.real(I_i))
            Z_plot.set_ydata(np.real(Z_i))

            V_diff.set_ydata(np.real(newYdata['V'] - V_i))  # update the data.
            I_diff.set_ydata(np.real(newYdata['I'] - I_i))  # update the data.
            Z_diff.set_ydata(np.real(newYdata['Z'] - Z_i))  # update the data.

            # ladderPlots.set_data(domainData, ladderData)
            # analyticalPlots.set_data(domainData, analyticalData)
            # differencePlots.set_data(domainData, differenceData)

            if iii == 2:
                ax1.relim()
                ax1.autoscale_view()

                ax2.relim()
                ax2.autoscale_view()

                ax3.relim()
                ax3.autoscale_view()
                plt.draw()

            # lines = ladderPlots#, analyticalPlots, differencePlots]
            emptyPlot = ax1.plot([])

            return emptyPlot

        # if animate == 1:
        #     if iii == 0:
        #         ax1.relim()
        #         ax1.autoscale_view()
        #         plt.draw()
        #     return line_V, line_I, line_Z


    if animate == 1:
        ani = animation.FuncAnimation(
            fig, animateFunc, interval=200, blit=True, save_count=10)
    elif animate == 2:
        plt.close(fig)
        fig2, ((ax21, ax22), (ax23, ax24)) = plt.subplots(2, 2)
        fig2.suptitle('Ladder Circuit')

        newYdata1 = ladderCircuit(w, L, C, N, Z_L, V_0)
        animateFunc(0,1)
        V_i, I_i, Z_i = analyticalData.T
        V_plot1, = ax21.plot(np.real(newYdata1['V']), label='Voltage')
        I_plot1, = ax21.plot(np.real(newYdata1['I']), label='Current')
        Z_plot1, = ax21.plot(np.real(newYdata1['Z']), label='Impedance')
        ax21.set_title('Values at 0')
        ax21.legend()
        plt.autoscale()

        V_plot1, = ax23.plot(np.real(V_i), label='Voltage')
        I_plot1, = ax23.plot(np.real(I_i), label='Current')
        Z_plot1, = ax23.plot(np.real(Z_i), label='Impedance')
        ax23.legend()
        plt.autoscale()

        V_1 = np.exp(1j*np.pi/2)
        animateFunc(0, V_1)
        V_i, I_i, Z_i = analyticalData.T
        newYdata2 = ladderCircuit(w, L, C, N, Z_L, V_1)
        V_plot2, = ax22.plot(np.real(newYdata2['V']), label='Voltage')
        I_plot2, = ax22.plot(np.real(newYdata2['I']), label='Current')
        Z_plot2, = ax22.plot(np.real(newYdata2['Z']), label='Impedance')
        ax22.legend()
        ax22.set_title('Values at $\pi / (2 \omega)$')
        plt.autoscale()

        V_plot2, = ax24.plot(np.real(V_i), label='Voltage')
        I_plot2, = ax24.plot(np.real(I_i), label='Current')
        Z_plot2, = ax24.plot(np.real(Z_i), label='Impedance')
        ax24.legend()
        plt.autoscale()

    else:
        V_max = np.empty(1000)
        for iii in range(1000):
            V_t = np.exp(1j * 0.005 * iii * 2)  # The extra 50 is to make the solution animate faster
            newYdata = ladderCircuit(w, L, C, N, Z_L, V_t)
            V_max[iii] = np.abs(np.real(newYdata['V'])).max()
        ax1.clear()
        line_V, = ax1.plot(newYdata['V'], label='Voltage')
        Line_I, = ax1.plot(newYdata['I'], label='Current')
        line_Z, = ax1.plot(newYdata['Z'], label='Impedance')
        ax1.legend()
