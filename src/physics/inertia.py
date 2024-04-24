import numpy as np

def ring_segment(ri: float=0.640, ro: float=0.660, h: float = 0.050, phi: float=2*np.pi, rho: float=2698.9 ) -> float:
    """
    Calculate the moment of inertia of a ring segment with rectangular cross-section.
    The rotation axis is in the ring center perpendicular to the ring (symmetry axis).
    
    Parameters:
    ----------
    ri: float
        inner radius of the ring segment (unit: meter)
    ro: float
        outer radius of the ring segment (unit: meter)
    h: float
        height the ring segment (unit: meter)
    phi: float
        opening angle of the segment (unit: radian)
    rho: float
        density of the material (unit: kg/m^3)
    
    Returns:
    --------
    float
        The moment of inertia of the ring segment (unit: kgm^2)
    """
    return 1/4*rho*h*phi*(ro**4-ri**4) # kg/m^3 * m * m^4 -> kgm^2


def hollow_ring_segment(ri: float=0.640, ro: float=0.660, h: float = 0.050, d: float = 0.005, phi: float=2*np.pi, rho: float=2698.9) -> float:
    """
    Calculate the moment of inertia of a hollow ring segment with rectangular cross-section and a given wall thickness.
    The rotation axis is in the ring center perpendicular to the ring (symmetry axis).
    
    Parameters:
    -----------
    ri: float
        inner radius of the ring segment (unit: meter)
    ro: float
        outer radius of the ring segment (unit: meter)
    h:  float
        height the ring segment (unit: meter)
    d:  float
        wall thickness (unit: meter)
    phi: float
        opening angle of the segment (unit: radian)
    rho: float 
        density of the material (unit: kg/m^3)
    
    Returns:
    --------
    float
        The moment of inertia of the ring segment (unit: kgm^2)
    """
    return ring_segment(ri, ro, h, phi, rho) - ring_segment(ri+d, ro-d, h-2*d, phi, rho)

def torus_segment(r: float=0.03, R: float=0.65, phi: float=2*np.pi, rho: float=2689.9)->float:
    """
    Calculate the moment of inertia of a torus segment (with circular cross-section).
    The rotation axis is in the torus center perpendicular to the torus (symmetry axis).
    
    Parameter:
    ----------
    r: float
        radius of the circular cross section (unit: meter)
    R: float 
        radius of the torus segment (unit: meter)
    phi: float
        opening angle of the segment (unit: radian)
    rho: float
        density of the material (unit: kg/m^3)
    
    Returns:
    --------
    float
        The moment of inertia of the torus segment (unit: kgm^2)
    """
    return phi*np.pi*rho*R*r**2*(0.75*r**2+R**2)


def hollow_torus_segment(r: float=0.015, R: float=0.65, d: float = 0.003, phi: float=2*np.pi, rho: float=2689.9)->float:
    """
    Calculate the moment of inertia of a hollow torus segment (with circular cross-section) and given wall thickness.
    The rotation axis is in the torus center perpendicular to the torus (symmetry axis).
    
    Parameter:
    ----------
    r: float
        radius of the circular cross section (unit: meter)
    R: float 
        radius of the torus segment (unit: meter)
    d: float
        wall thickness
    phi: float
        opening angle of the segment (unit: radian)
    rho: float
        density of the material (unit: kg/m^3)
    
    Returns:
    --------
    float
        The moment of inertia of the hollow torus segment (unit: kgm^2)
    """
    return torus_segment(r=r, R=R, phi=phi, rho=rho)-torus_segment(r=r-d, R=R, phi=phi, rho=rho)

def swinggate(r_torus: float, R_torus: float, h:float, d:float, phi:float, rho:float) -> float:
    """
    Calculate the moment of inertia of a swinggate barrier segment and given wall thickness.
    The rotation axis is in the barrier center perpendicular to the barrier (symmetry axis).
    
    Parameter:
    ----------
    r_torus: float
        outer radius of the circular cross-section of the half torus on the top and bottom of the element (unit: meter)
    R_torus: float 
        radius of the torus segment (unit: meter)
    h: float 
        height of the center piece between the two half-toruses (unit: meter)
    d: float
        wall thickness
    phi: float
        opening angle of the segment (unit: radian)
    rho: float
        density of the material (unit: kg/m^3)
    
    Returns:
    --------
    float
        The moment of inertia of the swinggate segment (unit: kgm^2)
    """
    ri=R_torus-r_torus
    ro=R_torus+r_torus
    I_torus = hollow_torus_segment(r=r_torus, R=R_torus, d=d, phi=phi, rho=rho)
    I_inner_wall = ring_segment(ri=ri, ro=ri+d, h=h, phi=phi, rho=rho)
    I_outer_wall = ring_segment(ri=ro-d, ro=ro, h=h, phi=phi, rho=rho)
    return I_torus+I_inner_wall+I_outer_wall


def solid_swinggate(r_torus: float, R_torus: float, h:float, phi:float, rho:float) -> float:
    """
    Calculate the moment of inertia of a solid swinggate barrier segment.
    The rotation axis is in the barrier center perpendicular to the barrier (symmetry axis).
    
    Parameter:
    ----------
    r_torus: float
        outer radius of the circular cross-section of the half torus on the top and bottom of the element (unit: meter)
    R_torus: float 
        radius of the torus segment (unit: meter)
    h: float 
        height of the center piece between the two half-toruses (unit: meter)
    phi: float
        opening angle of the segment (unit: radian)
    rho: float
        density of the material (unit: kg/m^3)
    
    Returns:
    --------
    float
        The moment of inertia of the swinggate segment (unit: kgm^2)
    """
    ri=R_torus-r_torus
    ro=R_torus+r_torus
    I_torus = torus_segment(r=r_torus, R=R_torus, phi=phi, rho=rho)
    I_ring = ring_segment(ri=ri, ro=ro, h=h, phi=phi, rho=rho)
    return I_torus+I_ring


def solid_swinggate_sandwich(ri: float, ro: float, h1:float, h2:float, phi:float, rho1: float, rho2: float) -> float:
    """
    Calculate the moment of inertia of a solid swinggate barrier segment.
    The rotation axis is in the barrier center perpendicular to the barrier (symmetry axis).
    
    Parameter:
    ----------
    ri: float
        inner radius of the ring segment (unit: meter)
    ro: float 
        outer radius of the ring segment (unit: meter)
    h1: float 
        combined height of of sandwich top and bottem (unit: meter)
    h2: float 
        height of of sandwich center (unit: meter)
    phi: float
        opening angle of the segment (unit: radian)
    rho1: float
        density of the material top and bottom (unit: kg/m^3)
    rho2: float
        density of the material in the center (unit: kg/m^3)
    
    Returns:
    --------
    float
        The moment of inertia of the swinggate segment (unit: kgm^2)
    """
    return ring_segment(ri=ri,ro=ro,h=h1,phi=phi,rho=rho1) + ring_segment(ri=ri,ro=ro,h=h2,phi=phi,rho=rho2)

def Main():
    print(hollow_ring_segment(ri=0.635, ro=0.665, h=0.06, d=0.003, phi=110/180*np.pi))
    print(hollow_torus_segment(r=0.015, R=0.650, d=0.003, phi=110/180*np.pi))
    print(swinggate(r_torus=0.015, R_torus=0.650, h=0.03, d=0.003, phi=110/180*np.pi, rho=2689.9))

    print(ring_segment(ri=0.635, ro=0.665, h=0.06, phi=110/180*np.pi))
    print(torus_segment(r=0.015, R=0.650, phi=110/180*np.pi))
    print(solid_swinggate(r_torus=0.015, R_torus=0.650, h=0.03, phi=110/180*np.pi, rho=2689.9))




if __name__=='__main__':
    Main()