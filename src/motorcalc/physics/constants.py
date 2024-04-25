from dataclasses import dataclass
from enum import Enum
import numpy as np

class SwingGateType(Enum):
    O_shape = 1
    O_shape_hollow = 2
    sandwich = 3
    rect_shape = 4
    rect_shape_hollow = 5

@dataclass
class Material():
    """Class to store all SWINGgate relevant material constants"""
    name: str
    density: float = None #kg/m^3
    emodule: float = None #kN/mm^2 (GPa)

@dataclass
class SwingGate():
    """Class to store the different SWINGgate geometries and materials"""
    name: str
    sg_type: int
    material: Material
    material_2: Material = None
    ri: float = 0.5                 #inner radius [ri] = m
    ro: float = 0.53                #outer radius [ro] = m
    r_top_bottom: float = 0.0       #radius of the 0-shape top and bottom
    h: float = 0.03                 #height of material, in case of the 0-shape, the height of the straight walls without the radius of the top and bottom curve
    h_2: float = 0.0                #height of material_2 in case of sandwich geometry
    d: float = 0.003                #thickness of the walls in case of hollow material
    phi: float = 110/180*np.pi      #opening angle of the ring segment
    delta: float = 90/180*np.pi     #angle of movement between the two positions


almg3 = Material(name='Aluminum', density=2680, emodule=70)
pmma = Material(name='Polymethlymetacrylat', density=1190)
brass = Material(name='Brass', density=8500, emodule=100)

sandwich = SwingGate(
    name='Sandwich Alu-PMMA-Alu', 
    sg_type=SwingGateType.sandwich, 
    material=almg3, 
    material_2=pmma, 
    ri=0.5, 
    ro=0.53, 
    h=0.016, 
    h_2=0.03
)


def Main():
    print(almg3)
    print(pmma)
    print(sandwich)
    pass

if __name__=='__main__':
    Main()