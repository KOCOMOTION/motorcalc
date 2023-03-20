import motorcalc.dcmotor as dcm
import numpy as np
import matplotlib.pyplot as plt

"""
Module to perform simple integration of the acceleration curve of a DC motor accelerating a given load defined by its moment of inertia.
"""

def integration_step(func:callable, old_val:float = 0.0, dt:float = 1.0E-3, kwargs: dict = None):
    """Calculate a simple first order integration of function func with keyword args kwargs"""
    return old_val + func(**kwargs)*dt
    
def integrate_omega_alpha(
    m: dcm.CDCMotor, 
    w_0: float=0.0, 
    alpha_0: float=0.0, 
    theta: float=1.0,
    loss_torque: float=0.0,
    t_start: float=0.0, 
    t_stop: float=3.0, 
    dt: float=0.1
):
    """Perfom the integration to calculate omega and alpha. Return value omega refers to the motor"""
    t=np.arange(t_start, t_stop, dt, dtype=np.float32)
    w=np.zeros(t.shape,dtype=np.float32)
    a=np.zeros(t.shape,dtype=np.float32)
    # d\alpha/dt = omega -> use a lamda function that takes omega as input und returns omega
    dalpha_dt = lambda omega:omega
    for ix,_ in enumerate(t):
        if ix==0:
            ww = integration_step(func=dw_dt, old_val=w_0, dt=dt ,kwargs={"motor":m, "theta":theta, "w_act":w_0, "loss_torque":loss_torque})
            aa = integration_step(func=dalpha_dt, old_val=alpha_0, dt=dt, kwargs={"omega":w_0})
        else:
            ww = integration_step(func=dw_dt, old_val=w[ix-1], dt=dt ,kwargs={"motor":m, "theta":theta, "w_act":w[ix-1], "loss_torque":loss_torque})
            # d\alpha/dt = omega -> use a lamda function that takes omega as input und returns omega
            aa = integration_step(func=dalpha_dt, old_val=a[ix-1], dt=dt, kwargs={"omega":w[ix-1]})
        w[ix]=ww
        a[ix]=aa
    return t, w, a


def plot_n_over_t(ax: plt.Axes, t:np.array, n:np.array):
    ax.plot(t, n)
    ax.set_xlabel(r'time (s)')
    ax.set_ylabel(r'n (rpm)')
    ax.grid(True)
    ax.set_title(r'motor speed', loc='Left')

def plot_rotations_over_t(ax: plt.Axes, t:np.array, number_of_rotations:np.array):
    ax.plot(t, number_of_rotations)
    ax.set_xlabel(r'time (s)')
    ax.set_ylabel(r'number of rotations')
    ax.grid(True)
    ax.set_title(r'number of rotations', loc='Left')

def plot_E_rot_over_t(ax: plt.Axes, t:np.array, E_rot:np.array):
    ax.plot(t, E_rot)
    ax.set_xlabel(r'time (s)')
    ax.set_ylabel(r'kinetic energy (J)')
    ax.grid(True)
    ax.set_title(r'kinetic energy', loc='Left')

def single_plot(ax:plt.Axes, x: float, y: float, xlabel: str = None, ylabel: str = None, title: str = None):
    ax.plot(x,y)
    ax.grid(True)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, loc='left')


def plot_data(t: np.array, w:np.array, a:np.array, E:np.array, I: np.array):
    fig = plt.figure(figsize=(10,10))
    axes=[]
    for n in range(4):
        axes.append(fig.add_subplot(2,2,n+1))
    single_plot(axes[0], x=t, y=dcm.omega_to_speed_rpm(w), xlabel='time (s)', ylabel='n (rpm)', title='Motor Speed')
    single_plot(axes[1], x=t, y=dcm.angle_to_number_of_rotations(a), xlabel='time (s)', ylabel='number of rotations', title='Motor Rotations')
    single_plot(axes[2], x=t, y=E, xlabel='time (s)', ylabel='kinetic energy (J)', title='Kinetic Energy')
    single_plot(axes[3], x=t, y=I, xlabel='time (s)', ylabel='current (A)', title='Motor Current')

    plt.show()

def dw_dt(
        motor:dcm.CDCMotor = None,
        theta:float = None, 
        w_act: float = None,
        loss_torque: float = 0.0,
) -> float:
    """Calculate the angular acceleration d\omega/dt [rad/s^2]"""
    if not motor or theta==None or w_act==None:
        return None
    res = (motor.calc_M_from_omega(omega=w_act) - loss_torque)/ theta
    if res<0:
        res=0
    return res

def E_rot(
        theta: float = None,
        w_act: float = None,
) -> float:
    """Calucate the rotational kinetic energy of an object with angular speed w_act and moment of inertia theta"""
    return 0.5*theta*w_act**2

if __name__ == '__main__':
    pass