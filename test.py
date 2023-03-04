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
    t=np.arange(0.0, 2.0, dt, dtype=np.float32)
    w=np.zeros(t.shape,dtype=np.float32)
    a=np.zeros(t.shape,dtype=np.float32)
    for ix,_ in enumerate(t):
        if ix==0:
            ww, aa = integration_step(motor=m, theta=0.6/42.0, w_act=w_0, alpha_act=alpha_0, dt=dt)
        else:
            ww, aa = integration_step(motor=m, theta=0.6/42.0, w_act=w[ix-1], alpha_act=a[ix-1], dt=dt)
        w[ix]=ww
        a[ix]=aa
    M = m.calc_M_from_omega(w)
    I = m.calc_I_from_M(M)

    plot_data(t, w, a, M, I, reduction_ratio=42)

def plot_data(t: np.array, w:np.array, a:np.array, M:np.array, I: np.array, reduction_ratio: float = 42):
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "Helvetica",        
    })
    fig = plt.figure(figsize=(10,8))
    ax1=fig.add_subplot(2,2,1)
    ax1.plot(t, w*30/np.pi)
    ax1.set_xlabel('time (s)')
    ax1.set_ylabel('n (rpm)')
    ax1.grid(True)
    ax2=fig.add_subplot(2,2,2)
    ax2.set_xlabel('time (s)')
    ax2.set_ylabel('$\alpha$ (°)')
    ax2.grid(True)
    ax2.plot(t,a/reduction_ratio/np.pi*180)
    ax3=fig.add_subplot(2,2,3)
    ax3.set_xlabel('time (s)')
    ax3.set_ylabel('Torque (Nm)')
    ax3.grid(True)
    ax3.plot(t,M)
    ax4=fig.add_subplot(2,2,4)
    ax4.set_xlabel('time (s)')
    ax4.set_ylabel('Current (A)')
    ax4.grid(True)
    ax4.plot(t,I)

    plt.show()

def dw_dt(
        motor:CDCMotor = None,
        theta:float = None, 
        w_act: float = None,
) -> float:
    """Calculate the angular acceleration d\omega/dt [rad/^2]"""
    if not motor or theta==None or w_act==None:
        return None
    res = motor.calc_M_from_n(n=w_act*30.0/np.pi) / theta
    return res

def E_rot(
        theta: float = None,
        w_act: float = None,
) -> float:
    """Calucate the rotational kinetic energy of an object with angular speed w_act and moment of inertia theta"""
    if not theta or not w_act:
        return None
    return 0.5*theta*w_act**2

def integration_step(
        motor: CDCMotor = None, 
        theta: float = None, 
        w_act: float = None,
        alpha_act: float = None,
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
    w_new = w_act + dw_dt(motor=motor, theta=theta, w_act=w_act)*dt
    alpha_new = alpha_act + w_act*dt
    return (w_new, alpha_new)
    
    



    
    
if __name__=="__main__":
    """
    Hier werden alle Projekte eingehängt, und diese Datei wird für die Berechnung ausgeführt
    """
    calc_motor_schulung()