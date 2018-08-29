import matplotlib.pyplot as plt 
import numpy as np 

a, b, c, d, e, chi2 = np.loadtxt("fit_samples.txt", unpack = True)
n_bins = 100

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

plt.figure(4)
plt.hist(d, n_bins, normed = True) 
plt.xlabel("$d$")
plt.ylabel("density")
plt.savefig("d.pdf", format = "pdf", bbox_inches = "tight")
plt.close(4)

plt.figure(5)
plt.hist(e, n_bins, normed = True) 
plt.xlabel("$e$")
plt.ylabel("density")
plt.savefig("e.pdf", format = "pdf", bbox_inches = "tight")
plt.close(5)

plt.figure(6)
plt.hist(chi2, n_bins, normed = True) 
plt.xlabel("$\\chi^2$")
plt.ylabel("density")
plt.savefig("chi2.pdf", format = "pdf", bbox_inches = "tight")
plt.close(6)




