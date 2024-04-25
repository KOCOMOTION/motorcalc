import numpy as np

def rotational_energy(omega: float, theta: float) -> float:
    """
    Calculate the rotational energy of an object with given moment of inertia theta and angular speed omega.

    Parameters:
    -----------
    omega: float
        angular speed (unit: radians/s)
    theta: float
        moment of inertia (unit: kgm^2)

    Returns:
    --------
    float
        rotational energy of the object (unit: kgm^2/s^2 = J)
    """
    return 0.5*theta*omega**2


def Main():
    pass

if __name__ =="__main__":
    Main()