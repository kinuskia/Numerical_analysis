import numpy as np 
import matplotlib.pyplot as plt 

x_in, y_in = np.loadtxt("data_points.txt", unpack=True)
x_out, y_out = np.loadtxt("spline.txt", unpack=True)


plt.scatter(x_in, y_in)
plt.plot(x_out, y_out, label="spline")
plt.plot(x_out, np.sin(x_out) +0.4*np.cos(2.*x_out), label="exact")
#plt.plot(x_out, np.cos(x_out), label="exact")
plt.legend(loc="best")
plt.savefig("spline.pdf", format="pdf", bbox_inches="tight")