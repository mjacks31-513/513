import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def ladderCircuit(w, L, C, N, Z_L, V_0):
    Z_total = np.empty(N, dtype=complex)  # I am going to test this
    I_total = np.empty(N, dtype=complex)
    V_total = np.empty(N, dtype=complex)

    Z_total[-1] = Z_L
    for n in range(N-2, -1, -1):
        Z_total[n] = (1j*w*C + (Z_total[n+1]+1j*w*L)**-1)**-1

    # I am not including V_0 in my V_total.
    # The reason I am doing this is because it makes my loop look nicer.
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
    Z_L = (5*(L/C) ** 0.5 + 0j)
    # Z_L = np.exp(1j * np.pi / 4)
    output = ladderCircuit(w, L, C, N, Z_L, V_0)

    fig, ax = plt.subplots()

    # x = np.arange(0, 2 * np.pi, 0.01)
    line_V, = ax.plot(np.arange(N), np.real(output['V']), label='Voltage')
    line_I, = ax.plot(np.arange(N), np.real(output['I']), label='Current')
    line_Z, = ax.plot(np.arange(N), np.real(output['Z']), label='Impedance')
    ax.legend()
    ax.set_ylim([-3, 3])

    def animate(i):
        # I am interested in giving this a try, but I am not sure it'll work
        V_t = np.exp(1j * 0.005 * i * 50)  #The extra 50 is to make the solution animate faster
        newYdata = ladderCircuit(w, L, C, N, Z_L, V_t)
        line_V.set_ydata(np.real(newYdata['V']))  # update the data.
        line_I.set_ydata(np.real(newYdata['I']))  # update the data.
        line_Z.set_ydata(np.real(newYdata['Z']))  # update the data.
        return line_V, line_I, line_Z

    ani = animation.FuncAnimation(
        fig, animate, interval=200, blit=True, save_count=10)
    # V_max = np.empty(1000)
    # for iii in range(1000):
    #     V_t = np.exp(1j * 0.005 * iii * 2)  # The extra 50 is to make the solution animate faster
    #     newYdata = ladderCircuit(w, L, C, N, Z_L, V_t)
    #     V_max[iii] = np.abs(np.real(newYdata['V'])).max()
    #
    # line_V, = ax.plot(newYdata['V'], label='Voltage')
    # Line_I, = ax.plot(newYdata['I'], label='Current')
    # line_Z, = ax.plot(newYdata['Z'], label='Impedance')
