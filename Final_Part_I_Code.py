import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def ladderCircuit(w, L, C, N, Z_L, V_0):
    Z_total = np.empty(N)  # I am going to test this
    I_total = np.empty(N)
    V_total = np.empty(N)

    Z_total[-1] = Z_L
    for n in range(N-2, -1, -1):
        Z_total[n] = (1j*w*C + (Z_total[n+1]+1j*w*L)**-1)**-1

    # I am not including V_0 in my V_total.
    # The reason I am doing this is because it makes my loop look nicer.
    I_total[0] = V_0/Z[0]
    V_total[0] = V_0 - 1j*w*L*I[0]
    for n in range(0, N-1):
        I_total[n+1] = I_total[n]-1j*w*C*V_total[n]
        V_total[n+1] = V_total[n]-1j*w*L*I_total[n+1]

    return {'V':V_total, 'I':I_total, 'Z':Z_total}


if __name__ == '__main__':
    V_0 = 1  # V
    L = 1  # H
    C = 1  # F
    w = 5E-3  # s**-1
    N = 1000  #
    L = 1  # H
    Z_L = (L/C)**(0.5)
    t = np.arange()
    output = ladderCircuit(w, L, C, N, Z_L, V_0)

    fig, ax = plt.subplots()

    x = np.arange(0, 2 * np.pi, 0.01)
    line, = ax.plot(x, np.sin(x))

    # I am interested in giving this a try, but
    def animate(i):
        line.set_ydata(np.sin(x + i / 50))  # update the data.
        return line,

    ani = animation.FuncAnimation(
        fig, animate, interval=20, blit=True, save_count=50)