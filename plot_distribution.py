import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from mult_gauss import MultiGauss
from gauss import Gauss

# Plot between -10 and 10 with .001 steps.
#x_axis = np.arange(-10, 10, 0.001)
# Mean = 0, SD = 2.
#plt.plot(x_axis, norm.pdf(x_axis,0,2), linewidth=10)
frame1 = plt.gca()
frame1.axes.xaxis.set_ticklabels([])
frame1.axes.yaxis.set_ticklabels([])
frame1.axes.xaxis.set_ticks([])
frame1.axes.yaxis.set_ticks([])

mult_gauss1 = MultiGauss([0.1, 0.5],[Gauss(25,4),Gauss(80,4)])
#mult_gauss1.plot_mult_gauss()
mult_gauss1 = mult_gauss1.normalise_gauss()
#mult_gauss1.plot_mult_gauss()
#mult_gauss1.truncate_gauss(0.05)
mult_gauss1.plot_mult_gauss(x=np.arange(0, 100, 0.001))
plt.show()
