import motorcalc.dcmotor as dcm
import motorcalc.load_dynamics as ldyn
import physics.inertia as inertia
import physics.constants as phconst
import matplotlib.pyplot as plt
import numpy as np

def plot_overview(t, w, a, E, I):
    fig = plt.figure(figsize=(10,10))
    axes=[]
    for n in range(4):
        axes.append(fig.add_subplot(2,2,n+1))
    ldyn.single_plot(axes[0], x=t, y=w/2/np.pi, xlabel='time (s)', ylabel='freq (Hz)', title='Vibration Frequency')
    ldyn.single_plot(axes[1], x=t, y=dcm.angle_to_number_of_rotations(a), xlabel='time (s)', ylabel='number of rotations', title='Motor Rotations')
    ldyn.single_plot(axes[2], x=t, y=E, xlabel='time (s)', ylabel='kinetic energy (J)', title='Kinetic Energy')
    ldyn.single_plot(axes[3], x=t, y=I, xlabel='time (s)', ylabel='current (A)', title='Motor Current')

    plt.show()



def calc_VC2039():
    m=dcm.CDCMotor(U_N=12, I_0=0.03, k_M=0.015, R=5)
    # m.list_spec_table()
    # m.plotCurves()

    dt=1.0E-6               # time step for integration [s]
    w_0=0.0                 # initial speed
    alpha_0=0.0             # inital angular position
    theta_imbalance=inertia.ring_segment(ri=2E-3, ro=6.6E-3, h=8.5E-3, phi=160/180*np.pi, rho=phconst.brass.density)
    theta = theta_imbalance+2E-8
    loss_torque = 0.0

    t, w, a =  ldyn.integrate_omega_alpha(m=m, w_0=w_0, alpha_0=alpha_0, theta=theta, \
                loss_torque=loss_torque, \
                t_start=0.0, t_stop=0.2, dt=dt)

    E = ldyn.E_rot(theta=theta, w_act=w)
    M = m.calc_M_from_omega(w)
    I = m.calc_I_from_M(M)

    plot_overview(t, w, a, E, I)

    
if __name__=="__main__":
    """
    Hier werden alle Projekte eingehängt, und diese Datei wird für die Berechnung ausgeführt
    """
    calc_VC2039()