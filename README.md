# motorcalc - A package to calculate and plot DC - motors

Gerrit Kocherscheidt
KOCO MOTION / KOCO automotive, March 2023

## Intention
With this package a simple calculation of DC motor parameters based on the input voltage, no-load current, terminal resistance and torque constant is performed. Motor performance curves can be plotted.

## Installation

### Installation using `pip`
`pip install motorcalc`

### Installation of sources from Github
`git clone https://github.com/KOCOMOTION/motorcalc.git`

## Basic usage
To create a DC-motor object use the following code
```python
import motorcalc.dcmotor as dcm

m = dcm.CDCmotor(U_N=12.0, I_0=0.02, k_M=0.014, R=25.0)
```
Here the parameter are 
* U_N: Driving voltage [V]
* I_0: No-load current [A]
* k_M: Torque constant [Nm/A]
* R: Terminal resistance [Ohm]

Note that if not stated explicitely, parameter in this package are considered to be in SI-units without magnitude scaling, i.e., I use ampere instead of milli-ampere etc.

### Get information about the object
To get information about the object the print function can be used
```python
print(m)

INPUT PARAMETER:
------------------------------------------------------------------------------------------------------
parameter	voltage		term. resist.	no-load cur.	no-load speed	torque const.
------------------------------------------------------------------------------------------------------
unit		Volt		Ohm		    Ampere		RPM		    Nm/A
value		12.0		25.00		0.020		7844		0.014

PERFORMANCE DATA:
------------------------------------------------------------------------------------------------------
parameter	unit	no-load		@max eff.	@max power	stall		@working point
------------------------------------------------------------------------------------------------------
speed		RPM	7844		6514		3922		0		    7844
current		A	0.020		0.098		0.250		0.480		0.020
torque		Nm	0.000		0.001		0.003		0.006		0.000
power		W	0.00		0.74		1.32		0.00		0.00
eff.		%	0.0		    63.3		44.1		0.0		    0.0
```

In addition parameter can be plotted with the ```list_spec_table()```-function of the CDCMotorExport-Module.

```python
import motorcalc.dcmotorexport as dcmexp

dcmexp.CDCMotorExport(m).list_spec_table()

+----+-----------------------+---------+----------+
| No | Parameter             | Unit    | Value    |
+----+-----------------------+---------+----------+
| 1  | Voltage               | V       | 12       |
+----+-----------------------+---------+----------+
| 2  | Terminal resistance   | Ω       | 25       |
+----+-----------------------+---------+----------+
| 3  | Terminal inductance   | mH      | 0        |
+----+-----------------------+---------+----------+
| 4  | No-load speed         | rpm     | 7844     |
+----+-----------------------+---------+----------+
| 5  | No-load current       | A       | 0.020    |
+----+-----------------------+---------+----------+
| 6  | Nominal torque        | mNm     | 0        |
+----+-----------------------+---------+----------+
| 7  | Nominal speed         | rpm     | 7844     |
+----+-----------------------+---------+----------+
| 8  | Nominal current       | A       | 0.020    |
+----+-----------------------+---------+----------+
| 9  | Max. output power     | W       | 1.323    |
+----+-----------------------+---------+----------+
| 10 | Max. efficiency       | %       | 63.342   |
+----+-----------------------+---------+----------+
| 11 | Back-EMF constant     | mV/rpm  | 133.690  |
+----+-----------------------+---------+----------+
| 12 | Torque constant       | mNm/A   | 14       |
+----+-----------------------+---------+----------+
| 13 | Speed/torque gradient | rpm/mNm | 1218.023 |
+----+-----------------------+---------+----------+
| 14 | Rotor inertia         | gcm^2   | 0        |
+----+-----------------------+---------+----------+
```
### Additional parameter

Missing or additional parameter can be added either on construction or by directly accessing the member variables of the CDCMotor-class.

```python
m = dcm.CDCMotor(U_N=12.0, I_0=0.002, k_M=0.014, R=25.0, H=1.0E-3 Theta=1.0E-8, nPoints=100, n_WP=6000.0, M_WP=0.0015, motor_name:str='test-motor' application:str='test-application')
```
    Parameters
    ----------
    U_N : float
        nominal Voltage in V
    I_0 : float
        noload current in A
    k_M : float
        torque constant in Nm/A
    R : float
        terminal resistance in Ohm
    H : float
        terminal inductance in mH (for information only)
    Theta : float
        rotor moment of inertia in gcm^2 (for information only)
    nPoints : int
        number of points to be plotted in graph
    n_WP : float
        required speed at working point
    M_WP : float
        required torque at working point
    motor_name : str
        name of the motor used for graph plotting
    application : str
        short description of application for graph plotting

## Plotting
To plot performance curves the ```CDCMotorPlot```-class from the dcmotorplot-module is needed. 
```python
import motorcalc.dcmotorplot as dcmplt
dcmplt.CDCMotorPlot(m).plot_curves()
```
With this command the following Matplotlib-graph is generated:

![Test graph](https://github.com/KOCOMOTION/motorcalc/blob/master/doc/img/performance_curve.png?raw=true)
