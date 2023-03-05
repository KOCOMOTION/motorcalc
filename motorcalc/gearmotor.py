## Calculation of a DC gear motor based on dcmotor.py
## Gerrit Kocherscheidt, KOCO automotive GmbH
## Date 01-Mar-2023

APP_NAME = "gearmotor.py"
APP_VERSION = "1.0"


import motorcalc.dcmotor as dcmotor
import numpy as np
import matplotlib.pyplot as plt
from openpyxl import load_workbook


class CDCMotorWithGearbox(dcmotor.CDCMotor):
    def __init__(
            self,
            U_N:float=0.0, 
            I_0:float=0.0, 
            k_M:float=0.0, 
            R:float=0.0, 
            H:float=0.0, 
            Theta:float=0.0,
            nPoints:int=100,
            n_WP:float=0.0, 
            M_WP:float=0.0, 
            motor_name:str='', 
            application:str='', 
            file_name:str='output.xlsx',
            GB_ratio:float=1.0,
            GB_eta:float=0.7,
            GB_name:str='Gearbox',
            n_WP_system:float=0,
            M_WP_system:float=0
        ):
        super(CDCMotorWithGearbox, self).__init__(U_N=U_N, I_0=I_0, k_M=k_M, R=R, H=H, Theta=Theta,
            nPoints=nPoints, n_WP=n_WP, M_WP=M_WP, motor_name=motor_name, application=application, file_name=file_name)
        self.M_WP_system = M_WP_system
        self.n_WP_system = n_WP_system
        self.GB_name = GB_name
        self.GB_eta = GB_eta
        self.GB_ratio = GB_ratio
        self.calc_system_values()

    def calc_system_values(self):
        self.M_system = self.M*self.GB_ratio*self.GB_eta
        self.M_0_system = self.M_0*self.GB_ratio*self.GB_eta
        self.M_maxpower_system = self.M_maxpower*self.GB_ratio*self.GB_eta
        self.M_S_system = self.M_S*self.GB_ratio*self.GB_eta
        self.M_meff_system = self.M_meff*self.GB_ratio*self.GB_eta
    
        self.n_system = self.n/self.GB_ratio
        self.n_0_system = self.n_system[0]
        self.n_meff_system = self.calc_n_from_M(self.M_meff)/self.GB_ratio   
        self.n_maxpower_system = self.calc_n_from_M(self.M_maxpower)/self.GB_ratio  

        self.P_mech_system = self.M_system * self.n_system *np.pi/30
        self.P_meff_system = self.M_meff_system * self.n_meff_system *np.pi/30
        self.P_maxpower_system = self.M_maxpower_system * self.n_maxpower_system *np.pi/30
        self.eta_system = self.eta*self.GB_eta
        self.eta_max_system = self.eta_max*self.GB_eta
        self.eta_maxpower_system = self.calc_eta_from_M(self.M_maxpower)*self.GB_eta

    def tune_voltage_to_working_point(self):
        self.n_WP=self.n_WP_system*self.GB_ratio
        self.M_WP=self.M_WP_system/self.GB_ratio/self.GB_eta
        super().tune_voltage_to_working_point()
        self.calc_system_values()

    def print_parameter(self):
        super().print_parameter()
        print('')
        print('gearbox input parameter')
        print('parameter\tred. ratio\tgb efficiency')
        print('unit\t\t-\t\t%')
        print('value\t\t{:0.1f}:1\t\t{:0.1f}'.format(self.GB_ratio,self.GB_eta*100.0))
        print('')

        print('system performance data:')
        print('parameter\tunit\tno-load\t\t@max eff.\t@max power\tstall')
        print('speed\t\tRPM\t{:0.0f}\t\t{:0.0f}\t\t{:0.0f}\t\t{:0.0f}'.format(self.n_0_system, self.n_meff_system, self.n_maxpower_system, 0))
        print('current\t\tA\t{:0.3f}\t\t{:0.3f}\t\t{:0.3f}\t\t{:0.3f}'.format(self.I_0\
            , self.I_meff, self.calc_I_from_M(self.M_maxpower),self.I_S))
        print('torque\t\tNm\t{:0.3f}\t\t{:0.3f}\t\t{:0.3f}\t\t{:0.3f}'.format(self.M_0_system, self.M_meff_system, self.M_maxpower_system, self.M_S_system))
        print('power\t\tW\t{:0.2f}\t\t{:0.2f}\t\t{:0.2f}\t\t{:0.2f}'.format(0, self.P_meff_system, self.P_maxpower_system,0))
        print('eff.\t\t%\t{:0.1f}\t\t{:0.1f}\t\t{:0.1f}\t\t{:0.1f}'.format(0, self.eta_max_system*100.0,\
            self.eta_maxpower_system,0))


    def export_to_excel(self):
        row, col = super().export_to_excel()
        wb = load_workbook(self.file_name)
        ws = wb.active
        row += 2
        col = 1
        titleStr = "System performance data"
        headerStr = ['', 'unit', '@no-load', '@max. eff.', '@max. power', 'stall']
        data = [\
            ['speed', 'RPM', self.n_0_system, self.n_meff_system, self.n_maxpower_system, 0],\
            ['current', 'A', self.I_0*self.GB_eta, self.I_meff*self.GB_eta, self.calc_I_from_M(self.M_maxpower)*self.GB_eta, self.I_S*self.GB_eta],\
            ['torque', 'Nm', self.M_0_system, self.M_meff_system, self.M_maxpower_system, self.M_S_system],\
            ['power', 'W', 0, self.P_meff_system, self.P_maxpower_system,0],\
            ['eff.', '%', 0, self.eta_max_system*100.0, self.eta_maxpower_system, 0]\
            ]
        self._excel_add_table_to_worksheet(ws, row, col, titleStr, headerStr, data)
        wb.save(self.file_name)
        

    def plotSystemCurves(self):
        plt.figure(figsize=(12,8))
        host = plt.subplot(111)
        host.set_title(self.application)
        #generate host plot
        plt.subplots_adjust(right=0.75)
        plt.subplots_adjust(left=0.1)
        plt.subplots_adjust(bottom=0.10)
        host.plot(self.M_system*1000.0,self.I,color="red")
        host.plot([1000.0*self.M_meff_system,1000.0*self.M_meff_system],[0.0, 1.2*self.I_S],":",color="black")
        host.plot([1000.0*self.M_maxpower_system,1000.0*self.M_maxpower_system],[0.0, 1.2*self.I_S],":",color="black")
        if self.M_WP_system != 0:
            host.plot([1000.0*self.M_WP_system,1000.0*self.M_WP_system],[0.0, 1.2*self.I_S],"-",color="cyan")
        host.set_xlabel("torque (mNm)")
        host.set_ylabel("current (A)")
        host.yaxis.label.set_color("red")
        host.spines["left"].set_color("red")
        host.tick_params(axis="y", colors="red")
        host.grid(True)
        host.set_xlim(0,1000.0*self.M_S_system)
        host.set_ylim(0,1.2*self.I_S)

        
        ax_power=host.twinx()
        ax_power.plot(self.M_system*1000.0,self.P_mech_system,".-",color="black")
        ax_power.plot(1000.0*self.M_meff_system,self.P_meff_system,"d",color="red")
        ax_power.plot(1000.0*self.M_maxpower_system,self.P_maxpower_system,"d",color="red")
        if self.M_WP_system != 0 and self.n_WP_system != 0:
            P_mech_WP_system = self.M_WP_system * self.n_WP_system * np.pi / 30.0
            ax_power.plot(1000.0*self.M_WP_system, P_mech_WP_system,"o",markerfacecolor="white", markeredgecolor="black", markersize=7)
        ax_power.spines["right"].set_position(("outward",0))
        ax_power.spines["left"].set_color("none")
        ax_power.spines["top"].set_color("none")
        ax_power.spines["bottom"].set_color("none")
        ax_power.yaxis.label.set_color("black")
        ax_power.set_ylabel("power (W)")
        ax_power.set_ylim(0,)

        ax_eta=host.twinx()
        ax_eta.plot(self.M_system*1000.0,self.eta_system*100.0,color="green")
        ax_eta.plot(self.M_meff_system*1000.0,self.eta_max_system*100.0,"d",color="red")
        ax_eta.spines["right"].set_position(("outward",120))
        ax_eta.spines["right"].set_color("green")
        ax_eta.spines["left"].set_color("none")
        ax_eta.spines["top"].set_color("none")
        ax_eta.spines["bottom"].set_color("none")
        ax_eta.yaxis.label.set_color("green")
        ax_eta.tick_params(axis="y", colors="green")
        ax_eta.set_ylabel("$\eta$ (%)")
        ax_eta.set_ylim(0,100)
        
        ax_speed=host.twinx()
        ax_speed.plot(self.M_system*1000.0,self.n_system,color="blue")
        if self.n_WP_system != 0 and self.M_WP_system != 0:
            ax_speed.plot(1000.0*self.M_WP_system,self.n_WP_system,"o",markerfacecolor="white",markeredgecolor="blue",markersize=7)
        ax_speed.spines["right"].set_position(("outward",60))
        ax_speed.spines["right"].set_color("blue")
        ax_speed.spines["left"].set_color("none")
        ax_speed.spines["top"].set_color("none")
        ax_speed.spines["bottom"].set_color("none")
        ax_speed.yaxis.label.set_color("blue")
        ax_speed.set_ylabel("speed (rpm)")
        ax_speed.yaxis.label.set_color("blue")
        ax_speed.tick_params(axis="y", colors="blue")
        ax_speed.set_ylim(0,)
        tstr = 'motor: {}\ngearbox: {}'.format(self.motor_name, self.GB_name)

        th=plt.text(0.65,0.95, tstr,
                    horizontalalignment='left',
                    verticalalignment='center',
                    transform=host.transAxes,
                    bbox=dict(fc='white',ec='black')
                    )
        plt.show()
        

def Main():
    # R=1.44 (-25°C)
    # R=2.00 (125°C)
    
    c = 4.7
    U_N = 10
    T = -40
    m_T = (2.0-1.44)/(125+25)
    T_0 = 1.44 - m_T*(-25)
    R = m_T*T+T_0*c
    I_S = U_N / R
    k_M = 0.005985*np.sqrt(c)
    M_S = 0.017099
    I_0 = 0.13*6/U_N
    n_WP = 300
    M_WP = 0.005
    dcm=dcmotor.CDCMotor(U_N=U_N, I_0=I_0, k_M=k_M, R=R, H=0.61, Theta=6.7, n_WP=n_WP, M_WP=M_WP, application="166_A / LiDAR", motor_name="BO2015_Version 10V")
    dcm.print_parameter()
    dcm.tune_voltage_to_working_point()
    dcm.list_spec_table()
 
    # GB_eta=0.8
    # GB_ratio=100
    # GB_name='100:1'
    # gbmotor=CDCMotorWithGearbox(U_N=U_N, I_0=I_0, k_M=k_M, R=R, GB_ratio=GB_ratio, GB_eta=GB_eta, \
    #     n_WP_system=n_WP/GB_ratio, M_WP_system=M_WP*GB_ratio*GB_eta, \
    #     application="Test with GB", motor_name="BO2015_Version 10V", GB_name=GB_name)
    # gbmotor.print_parameter()
    # gbmotor.tune_voltage_to_working_point()
    # gbmotor.print_parameter()
    # gbmotor.plotSystemCurves()


if __name__ == "__main__":
    Main()