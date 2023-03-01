from motorcalc.dcmotor import CDCMotor

def calc_motor_schulung():
    m=CDCMotor(U_N=13.5, I_0=0.06, k_M=0.016969, R=2.88)
    m.plotCurves()

if __name__=="__main__":
    """
    Hier werden alle Projekte eingehängt, und diese Datei wird für die Berechnung ausgeführt
    """
    calc_motor_schulung()