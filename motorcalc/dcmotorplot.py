import numpy as np
import matplotlib.pyplot as plt
import motorcalc.dcmotor as dcm

class CDCMotorPlot():
    def __init__(self, dcm:dcm.CDCMotor):
        self.dcm = dcm
        self.dcm.calc_performance_curves()

    def plot_curves(self, addVoltagesSpeed:list[float]=None):
        """
        Plots system values into a matplotlib graph

        Parameters:
        -----------
        addVoltagesSpeed : list of floats
            List of additional voltages that are used to plot additional n-over-M-curves
        """
        fig=plt.figure(figsize=(12,8))
        fig.patch.set_facecolor('white')
        host = plt.subplot(111)
        host.set_title(self.dcm.application)
        host.patch.set_facecolor('white')

        #generate host plot
        plt.subplots_adjust(right=0.75)
        plt.subplots_adjust(left=0.1)
        plt.subplots_adjust(bottom=0.10)
        host.plot(self.dcm.M*1000.0,self.dcm.I,color="red")
        host.plot([1000.0*self.dcm.M_meff,1000.0*self.dcm.M_meff],[0.0, 1.2*self.dcm.I_S],":",color="black")
        host.plot([1000.0*self.dcm.M_maxpower,1000.0*self.dcm.M_maxpower],[0.0, 1.2*self.dcm.I_S],":",color="black")
        if self.dcm.M_WP != 0:
            host.plot([1000.0*self.dcm.M_WP,1000.0*self.dcm.M_WP],[0.0, 1.2*self.dcm.I_S],"-",color="cyan")
        host.set_xlabel("torque (mNm)")
        host.set_ylabel("current (A)")
        host.yaxis.label.set_color("red")
        host.spines["left"].set_color("red")
        host.tick_params(axis="y", colors="red")
        host.grid(True)
        host.set_xlim(0,1000.0*self.dcm.M_S)
        host.set_ylim(0,1.2*self.dcm.I_S)
        

        ax_power=host.twinx()
        ax_power.patch.set_facecolor('white')
        ax_power.plot(self.dcm.M*1000.0,self.dcm.P_mech,".-",color="black")
        ax_power.plot(1000.0*self.dcm.M_meff,self.dcm.P_meff,"d",color="red")
        ax_power.plot(1000.0*self.dcm.M_maxpower,self.dcm.P_maxpower,"d",color="red")
        if self.dcm.M_WP != 0 and self.dcm.n_WP != 0:
            P_mech_WP = self.dcm.M_WP * self.dcm.n_WP * np.pi / 30.0
            ax_power.plot(1000.0*self.dcm.M_WP, P_mech_WP,"o",markerfacecolor="white", markeredgecolor="black", markersize=7)
        ax_power.spines["right"].set_position(("outward",0))
        ax_power.spines["left"].set_color("none")
        ax_power.spines["top"].set_color("none")
        ax_power.spines["bottom"].set_color("none")
        ax_power.yaxis.label.set_color("black")
        ax_power.set_ylabel("power (W)")
        ax_power.set_ylim(0,)

        ax_eta=host.twinx()
        ax_eta.patch.set_facecolor('white')
        ax_eta.plot(self.dcm.M*1000.0,self.dcm.eta*100.0,color="green")
        ax_eta.plot(self.dcm.M_meff*1000.0,self.dcm.eta_max*100.0,"d",color="red")
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
        ax_speed.patch.set_facecolor('white')
        ax_speed.plot(self.dcm.M*1000.0,self.dcm.n,color="blue")
        if self.dcm.n_WP != 0 and self.dcm.M_WP != 0:
            ax_speed.plot(1000.0*self.dcm.M_WP,self.dcm.n_WP,"o",markerfacecolor="white",markeredgecolor="blue",markersize=7)
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
        ax_speed.annotate('U={:0.1f}V'.format(self.dcm.U_N),(self.dcm.M[80]*1000*1.01,self.dcm.n[80]*1.01), color='blue', fontweight='bold')
        tstr = 'motor: {}'.format(self.dcm.motor_name)

        if addVoltagesSpeed:
            for voltage in addVoltagesSpeed:
                m_tmp = dcm.CDCMotor(U_N=voltage, R=self.dcm.R, I_0=self.dcm.I_0, k_M=self.dcm.k_M)
                m_tmp.calc_performance_curves()
                ax_speed.plot(m_tmp.M*1000.0, m_tmp.n,color="blue", linestyle=":")
                ax_speed.annotate('U={:0.1f}V'.format(m_tmp.U_N),(m_tmp.M[80]*1000*1.01,m_tmp.n[80]*1.01), color='blue', fontweight='normal')
            print('voltage constant = {:0.0f}RPM/V'.format((self.dcm.n_0-m_tmp.n_0)/(self.dcm.U_N-m_tmp.U_N)))


        th=plt.text(0.65,0.95, tstr,
                    horizontalalignment='left',
                    verticalalignment='center',
                    transform=host.transAxes,
                    bbox=dict(fc='white',ec='black')
                    )
        plt.show()



def Main():
    pass


if __name__ == "__main__":
    Main()