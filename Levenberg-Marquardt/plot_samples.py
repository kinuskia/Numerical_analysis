import matplotlib.pyplot as plt 
import numpy as np 

a, b, c = np.loadtxt("fit_samples.txt", unpack = True)
n_bins = 64

plt.figure(1)
plt.hist(a, n_bins, normed = True) 
plt.xlabel("$a$")
plt.ylabel("density")
plt.savefig("a.pdf", format = "pdf", bbox_inches = "tight")
plt.close(1)

plt.figure(2)
plt.hist(b, n_bins, normed = True) 
plt.xlabel("$b$")
plt.ylabel("density")
plt.savefig("b.pdf", format = "pdf", bbox_inches = "tight")
plt.close(2)

plt.figure(3)
plt.hist(c, n_bins, normed = True) 
plt.xlabel("$c$")
plt.ylabel("density")
plt.savefig("c.pdf", format = "pdf", bbox_inches = "tight")
plt.close(3)




