import motorcalc.dcmotor as dcm
import motorcalc.load_dynamics as ldyn
import physics.inertia as inertia
import physics.constants as phconst
import matplotlib.pyplot as plt
import numpy as np

def plot_overview(t, w, a, F, I):
    fig = plt.figure(figsize=(10,10))
    axes=[]
    for n in range(4):
        axes.append(fig.add_subplot(2,2,n+1))
    ldyn.single_plot(axes[0], x=t, y=dcm.omega_to_speed_rpm(w), xlabel='time (s)', ylabel='n (RPM)', title='Motor Speed')
    ldyn.single_plot(axes[1], x=t, y=dcm.angle_to_number_of_rotations(a), xlabel='time (s)', ylabel='number of rotations', title='Motor Rotations')
    ldyn.single_plot(axes[2], x=t, y=F, xlabel='time (s)', ylabel='radial force (N)', title='Radial Force')
    ldyn.single_plot(axes[3], x=t, y=I, xlabel='time (s)', ylabel='current (A)', title='Motor Current')

    plt.show()


def calc_mass_half_cylinder(r: float, h: float, rho: float):
    return 1/2*np.pi*(r**2)*h*rho

def calc_r_cog_half_cylinder(r: float):
    return 4/(3*np.pi)*r

def calc_imbalance_force(m: float, omega: float, r_cog: float):
    return m * omega**2 * r_cog


def calc_VC2039():
    I_S = 1.5
    U_nom = 12.5
    R = U_nom/I_S
    m=dcm.CDCMotor(U_N=4, I_0=0.050, k_M=0.013, R=R)
    # m.list_spec_table()
    # m.plotCurves()

    dt=1.0E-6               # time step for integration [s]
    w_0=0.0                 # initial speed
    alpha_0=0.0             # inital angular position
    ro_imbalance = 6.6E-3   # outer radius imbalance, 6.6mm on Constar Drawing
    theta_imbalance=inertia.ring_segment(ri=2E-3, ro=ro_imbalance, h=8.5E-3, phi=160/180*np.pi, rho=phconst.brass.density) \
        + inertia.ring_segment(ri=0,ro=2E-3, h=8.5E-3, phi=2*np.pi, rho=phconst.brass.density)
    theta = theta_imbalance+2E-8

    imbalance_mass = calc_mass_half_cylinder(ro_imbalance, 8.5E-3, phconst.brass.density)
    # print(imbalance_mass)
    r_cog = calc_r_cog_half_cylinder(ro_imbalance)

    loss_torque = 0.003

    t, w, a =  ldyn.integrate_omega_alpha(m=m, w_0=w_0, alpha_0=alpha_0, theta=theta, \
                loss_torque=loss_torque, \
                t_start=0.0, t_stop=0.2, dt=dt)

    # E = ldyn.E_rot(theta=theta, w_act=w)
    F = calc_imbalance_force(imbalance_mass, w, r_cog)
    M = m.calc_M_from_omega(w)
    I = m.calc_I_from_M(M)

    plot_overview(t, w, a, F, I)

    
if __name__=="__main__":
    """
    Hier werden alle Projekte eingehängt, und diese Datei wird für die Berechnung ausgeführt
    """
    calc_VC2039()