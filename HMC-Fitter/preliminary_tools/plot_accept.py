import numpy as np 
import matplotlib.pyplot as plt 

accept = np.loadtxt("acceptrates.txt", unpack = True)


plt.figure(1)
plt.hist(accept*100, 50)
plt.xlabel("accept_rate")
plt.ylabel("#")
plt.savefig("accept_rate.pdf", format = "pdf", bbox_inches = "tight")





