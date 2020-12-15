import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

global SCALING_FACTOR_ON
SCALING_FACTOR_ON = True

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
    V_0 = (1+0j)  # V
    L = 1  # H
    C = 1  # F
    w = 5E-3  # s**-1
    N = 1000  #
    Z_L = (3*(L/C) ** 0.5 + 0j)
    ## animate 1 = animation window
    ## animate 2 = 2 plots at 0 and Pi/2
    ## animate 3 = 1000 iterations VSWR
    animate = 1
    # Z_L = np.exp(1j * np.pi / 4)
    output = ladderCircuit(w, L, C, N, Z_L, V_0)

    fig, ax = plt.subplots()

    # x = np.arange(0, 2 * np.pi, 0.01)
    line_V, = ax.plot(np.arange(N), np.real(output['V']), label='Voltage')
    line_I, = ax.plot(np.arange(N), np.real(output['I']), label='Current')
    line_Z, = ax.plot(np.arange(N), np.real(output['Z']), label='Impedance')
    ax.legend()
    plt.autoscale()


    def animateFunc(i):
        # I am interested in giving this a try, but I am not sure it'll work
        V_t = np.exp(1j * 0.005 * i * 50)  #The extra 50 is to make the solution animate faster
        newYdata = ladderCircuit(w, L, C, N, Z_L, V_t)
        if animate == 1:
            line_V.set_ydata(np.real(newYdata['V']))  # update the data.
            line_I.set_ydata(np.real(newYdata['I']))  # update the data.
            line_Z.set_ydata(np.real(newYdata['Z']))  # update the data.
            if i == 0:
                ax.relim()
                ax.autoscale_view()
                plt.draw()
            return line_V, line_I, line_Z

    if animate == 1:
        ani = animation.FuncAnimation(
            fig, animateFunc, interval=200, blit=True, save_count=10)
    elif animate == 2:
        plt.close(fig)
        fig2, (ax21, ax22) = plt.subplots(2, 1)
        fig2.suptitle('Ladder Circuit')

        newYdata1 = ladderCircuit(w, L, C, N, Z_L, V_0)
        V_plot1, = ax21.plot(np.real(newYdata1['V']), label='Voltage')
        I_plot1, = ax21.plot(np.real(newYdata1['I']), label='Current')
        Z_plot1, = ax21.plot(np.real(newYdata1['Z']), label='Impedance')
        ax21.set_title('Values at 0')
        ax21.legend()
        plt.autoscale()

        V_1 = np.exp(1j*np.pi/2)
        newYdata2 = ladderCircuit(w, L, C, N, Z_L, V_1)
        V_plot2, = ax22.plot(np.real(newYdata2['V']), label='Voltage')
        I_plot2, = ax22.plot(np.real(newYdata2['I']), label='Current')
        Z_plot2, = ax22.plot(np.real(newYdata2['Z']), label='Impedance')
        ax22.legend()
        ax22.set_title('Values at $\pi / (2 \omega)$')
        plt.autoscale()
    else:
        V_max = np.empty(1000)
        for iii in range(1000):
            V_t = np.exp(1j * 0.005 * iii * 2)  # The extra 50 is to make the solution animate faster
            newYdata = ladderCircuit(w, L, C, N, Z_L, V_t)
            V_max[iii] = np.abs(np.real(newYdata['V'])).max()
        ax.clear()
        line_V, = ax.plot(newYdata['V'], label='Voltage')
        Line_I, = ax.plot(newYdata['I'], label='Current')
        line_Z, = ax.plot(newYdata['Z'], label='Impedance')
        ax.legend()
