from motorcalc.dcmotor import CDCMotor
import numpy as np
import matplotlib.pyplot as plt

def calc_motor_schulung():
    # m=CDCMotor(U_N=24, I_0=1, k_M=0.0643, R=0.085)
    # m.list_spec_table()
    # m.plotCurves()

    m=CDCMotor(U_N=24, I_0=0.080, k_M=0.55, R=7)
    # m.list_spec_table()
    # m.plotCurves()

    dt=1.0E-4       # time step for integration [s]
    w_0=0.0         # initial speed
    alpha_0=0.0     # inital angular position
    t=np.arange(0.0, 3.0, dt, dtype=np.float32)
    w=np.zeros(t.shape,dtype=np.float32)
    a=np.zeros(t.shape,dtype=np.float32)
    theta=0.6
    reduction_ratio=42
    loss_torque = 0.2
    for ix,_ in enumerate(t):
        if ix==0:
            ww, aa = integration_step(motor=m, theta=theta/reduction_ratio, w_act=w_0, alpha_act=alpha_0, loss_torque=loss_torque, dt=dt)
        else:
            ww, aa = integration_step(motor=m, theta=theta/reduction_ratio, w_act=w[ix-1], alpha_act=a[ix-1], loss_torque=loss_torque, dt=dt)
        w[ix]=ww
        a[ix]=aa
    E = E_rot(theta=theta, w_act=w/reduction_ratio)
    M = m.calc_M_from_omega(w)
    I = m.calc_I_from_M(M)

    plot_data(t, w, a, E, I, reduction_ratio=42)

def plot_data(t: np.array, w:np.array, a:np.array, E:np.array, I: np.array, reduction_ratio: float = 42):
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "Helvetica",        
    })
    fig = plt.figure(figsize=(10,10))
    ax1=fig.add_subplot(2,2,1)
    ax1.plot(t, w*30/np.pi)
    ax1.set_xlabel(r'time (s)')
    ax1.set_ylabel(r'n (rpm)')
    ax1.grid(True)
    ax1.set_title(r'motor speed', loc='Left')
    ax2=fig.add_subplot(2,2,2)
    ax2.set_xlabel(r'time (s)')
    ax2.set_ylabel(r'$\alpha$ (°)')
    ax2.grid(True)
    ax2.set_title(r'angle of gate', loc='Left')
    ax2.plot(t,a/reduction_ratio/np.pi*180)
    ax3=fig.add_subplot(2,2,3)
    ax3.set_xlabel(r'time (s)')
    ax3.set_ylabel(r'kinetic energy (J)')
    ax3.grid(True)
    ax3.set_title(r'kinetic energy of gate', loc='Left')
    ax3.plot(t,E)
    ax4=fig.add_subplot(2,2,4)
    ax4.set_xlabel(r'time (s)')
    ax4.set_ylabel(r'current (A)')
    ax4.grid(True)
    ax4.set_title(r'motor current', loc='Left')
    ax4.plot(t,I)

    plt.show()

def dw_dt(
        motor:CDCMotor = None,
        theta:float = None, 
        w_act: float = None,
        loss_torque: float = 0.0,

) -> float:
    """Calculate the angular acceleration d\omega/dt [rad/^2]"""
    if not motor or theta==None or w_act==None:
        return None
    res = (motor.calc_M_from_n(n=w_act*30.0/np.pi) - loss_torque)/ theta
    if res<0:
        res=0
    return res

def E_rot(
        theta: float = None,
        w_act: float = None,
) -> float:
    """Calucate the rotational kinetic energy of an object with angular speed w_act and moment of inertia theta"""
    return 0.5*theta*w_act**2

def integration_step(
        motor: CDCMotor = None, 
        theta: float = None, 
        w_act: float = None,
        alpha_act: float = None,
        loss_torque: float = 0.0,
        dt: float = 1.0E-3,
    ):
    if not motor:
        return None
    if theta==None:
        return None
    if w_act==None:
        return None
    if alpha_act==None:
        return None
    w_new = w_act + dw_dt(motor=motor, theta=theta, w_act=w_act, loss_torque=loss_torque)*dt
    alpha_new = alpha_act + w_act*dt
    return (w_new, alpha_new)
    
    



    
    
if __name__=="__main__":
    """
    Hier werden alle Projekte eingehängt, und diese Datei wird für die Berechnung ausgeführt
    """
    calc_motor_schulung()