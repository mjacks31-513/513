import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


if __name__=="__main__":
# V_total = np.zeros_like(x, dtype=complex)
    Mx = Nx = 1000
    x = np.arange(Nx)
    lam = 1300

    # Z_0 is the static impedance in the transmission line
    Z_0 = 1 + 0j  # This is defined as np.sqrt(L/C)
    # Z_L is the static impedance in the load line
    Z_L = 5+0j
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
    ax1.set_ylim([-3, 3])
    ax1.legend()
    animate = True


    def animateFunc(iii):
        # e^(wt)
        V_0 = np.exp(0.005j*iii*50)
        # if complex(Z_L) == complex(Z_0):
        #     rho = np.full_like(x, rho_0, dtype=complex)
        # else:
        rho = rho_0 * np.exp(4j * np.pi * (x - Nx + 1) / lam)

        V_i[:] = V_0 * np.exp(-2j*np.pi*(x-Nx+1)/lam)*(1 + rho)
        I_i[:] = V_0/Z_0 * np.exp(-2j*np.pi*(x-Nx+1)/lam) * (1-rho)
        Z_i[:] = V_i/I_i

        if animate:
            V_plot.set_ydata(np.real(V_i))
            I_plot.set_ydata(np.real(I_i))
            Z_plot.set_ydata(np.real(Z_i))
            return V_plot, I_plot, Z_plot


    if animate:
        ani = animation.FuncAnimation(
            fig, animateFunc, interval=200, blit=True, save_count=10)
    else:
        V_max = np.empty(1000)
        for iii in range(1000):
            animateFunc(iii)
            V_max[iii] = np.abs(np.real(V_i)).max()

        V_plot, = ax1.plot(np.real(V_i), label='Voltage')
        I_plot, = ax1.plot(np.real(I_i), label='Current')
        Z_plot, = ax1.plot(np.real(Z_i), label='Impedance')
