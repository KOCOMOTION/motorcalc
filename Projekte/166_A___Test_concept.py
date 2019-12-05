import matplotlib.pyplot as plt
import numpy as np

class CVisualizeTestConcept():
    def __init__(self):
        self.x = np.linspace(-1, 1, num=201, endpoint=True, dtype=float)

    def _gauss(self, sigma=1, x0=0):
        sigma_sqr=sigma**2
        return(1/(np.sqrt(2*np.pi*sigma_sqr))*np.exp(-np.power(self.x-x0,2)/(2*sigma_sqr)))

    def _zeta_plot(self, sigma=[0.15, 0.1, 0.2], x0=[-0.2, 0.0, 0.35], label=('n vs. golden', 'golden vs. golden', 'm vs. golden')):
        fig = plt.figure(figsize=(10,8))
        fig.suptitle=("Distribution of measurement values")
        ax = fig.add_subplot(111)
        ph=[]
        col = ('red', 'black', 'green')
        for ix, s_val in enumerate(sigma):
            ph.append(ax.plot(VTC.x, VTC._gauss(sigma=s_val, x0=x0[ix]), color=col[ix], label=label[ix]))
        ax.legend()
        ax.grid(True)
        ax.set_xticks(x0)
        ax.set_xticklabels(("$\overline{\zeta}_{n}$", "$\overline{\zeta}_{g}$", "$\overline{\zeta}_{m}$"))
        ax.set_yticklabels(())
        ax.set_xlabel('quantified parameter $\zeta$')
        ax.set_ylabel('probability of occurence P')
        plt.show()


if __name__ == "__main__":
    VTC = CVisualizeTestConcept()
    VTC._zeta_plot()
    VTC._zeta_plot(sigma = [0.35, 0.2, 0.4])
