import matplotlib.pyplot as plt 
import numpy as np 

x1, x2, y = np.loadtxt("outfile.txt", unpack = True)
n_bins = 100

plt.figure(1)
plt.hist(x1, n_bins, normed = True) 
plt.xlabel("$x_1$")
plt.ylabel("density")
plt.savefig("x1.pdf", format = "pdf", bbox_inches = "tight")
plt.close(1)

plt.figure(2)
plt.hist(x2, n_bins, normed = True) 
plt.xlabel("$x_2$")
plt.ylabel("density")
plt.savefig("x2.pdf", format = "pdf", bbox_inches = "tight")
plt.close(2)

# plt.figure(3)
# plt.hist(x3, n_bins, normed = True) 
# plt.xlabel("$x_3$")
# plt.ylabel("density")
# plt.savefig("x3.pdf", format = "pdf", bbox_inches = "tight")
# plt.close(3)



plt.figure(4)
plt.hist(y, n_bins, normed = True) 
plt.xlabel("$y$")
plt.ylabel("density")
plt.savefig("y.pdf", format = "pdf", bbox_inches = "tight")
plt.close(4)

print(np.mean(x1))
print(np.mean(x2))
#print(np.mean(x3))
print(np.mean(y))





