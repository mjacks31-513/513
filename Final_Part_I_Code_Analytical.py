import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


if __name__=="__main__":
    # V_total = np.zeros_like(x, dtype=complex)
    Mx = Nx = 1000
    x = np.arange(Nx)
    w = 0.005
    L = 1
    C = 1
    lam = 2*np.pi/w
    ## animate 1 = animation window
    ## animate 2 = 2 plots at 0 and Pi/2
    ## animate 3 = 1000 iterations VSWR check
    animate = 2

    # Z_0 is the static impedance in the transmission line
    Z_0 = 1 + 0j  # This is defined as np.sqrt(L/C)
    # Z_L is the static impedance in the load line
    Z_L = 3+0j
    # This is my constant rho value
    rho_0 = (Z_L - Z_0)/(Z_L + Z_0)
    # global Vi
    # global V_total
    # global V_reflected
    # global V_transmit
    # Vi = np.zeros_like(x, dtype=complex)
    # V_reflected = np.zeros_like(x, dtype=complex)
    # V_transmit = np.zeros_like(x, dtype=complex)
    global V_i
    global I_i
    global Z_i
    V_i = np.zeros_like(x, dtype=complex)
    I_i = np.zeros_like(x, dtype=complex)
    Z_i = np.zeros_like(x, dtype=complex)

    fig, ax1 = plt.subplots(1, 1)
    # V_plot, = ax1.plot(np.real(Vi), label='V initial')
    # I_plot, = ax1.plot(np.real(V_total), label='V Total')
    # Z_plot, = ax1.plot(np.real(V_reflected), label='V Reflected')
    V_plot, = ax1.plot(np.real(V_i), label='Voltage')
    I_plot, = ax1.plot(np.real(I_i), label='Current')
    Z_plot, = ax1.plot(np.real(Z_i), label='Impedance')
    plt.autoscale()
    ax1.legend()


    def animateFunc(iii, V_0=None):
        # e^(wt)
        if V_0 is None:
            V_0 = np.exp(0.005j * iii * 50)
        # if complex(Z_L) == complex(Z_0):
        #     rho = np.full_like(x, rho_0, dtype=complex)
        # else:
        global rho
        rho = rho_0 * np.exp(4j * np.pi * (x - Nx + 1) / lam)

        # I need to maintain the source voltage
        V_c = V_0 / (1 + rho[0])
        V_i[:] = V_c * np.exp(-2j*np.pi*x/lam)*(1 + rho)
        I_i[:] = V_c/Z_0 * np.exp(-2j*np.pi*x/lam) * (1-rho)
        Z_i[:] = V_i/I_i

        if animate == 1:
            V_plot.set_ydata(np.real(V_i))
            I_plot.set_ydata(np.real(I_i))
            Z_plot.set_ydata(np.real(Z_i))
            if iii == 0:
                ax1.relim()
                ax1.autoscale_view()
                plt.draw()
            return V_plot, I_plot, Z_plot


    if animate == 1:
        ax1.set_title('Analytical')
        ani = animation.FuncAnimation(
            fig, animateFunc, interval=200, blit=True, save_count=10)
        ax1.relim()
        ax1.autoscale_view()
    elif animate == 2:
        plt.close(fig)
        fig2, (ax21, ax22) = plt.subplots(2, 1)
        fig2.suptitle('Analytical')

        V_0 = 1 + 0j
        animateFunc(0, 1)
        V_plot1, = ax21.plot(np.real(V_i), label='Voltage')
        I_plot1, = ax21.plot(np.real(I_i), label='Current')
        Z_plot1, = ax21.plot(np.real(Z_i), label='Impedance')
        ax21.set_title('Values at 0')
        ax21.legend()
        plt.autoscale()

        V_1 = np.exp(1j*np.pi/2)
        animateFunc(0, V_1)
        V_plot2, = ax22.plot(np.real(V_i), label='Voltage')
        I_plot2, = ax22.plot(np.real(I_i), label='Current')
        Z_plot2, = ax22.plot(np.real(Z_i), label='Impedance')
        ax22.legend()
        ax22.set_title('Values at $\pi / (2 \omega)$')
        plt.autoscale()
    else:
        V_max = np.empty(1000)
        for iii in range(1000):
            animateFunc(iii)
            V_max[iii] = np.abs(np.real(V_i)).max()
        ax1.clear()
        V_plot, = ax1.plot(np.real(V_i), label='Voltage')
        I_plot, = ax1.plot(np.real(I_i), label='Current')
        Z_plot, = ax1.plot(np.real(Z_i), label='Impedance')
        ax1.legend()