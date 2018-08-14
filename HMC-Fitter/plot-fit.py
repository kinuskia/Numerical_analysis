import numpy as np 
import matplotlib.pyplot as plt 

n, a0, a1, a2, U, accept = np.loadtxt("data.txt", unpack = True)

keep = (U < 15)

parameters = [a0, a1, a2]

T = 1e0
p = 3
n_bins = 128


plt.figure(1)
plt.hist(U[keep], n_bins, normed = True) 
plt.xlabel("$\\chi^2_\\mathrm{red}$")
plt.ylabel("density")
plt.savefig("plots/hist_U.pdf", format = "pdf", bbox_inches = "tight")
plt.close(1)
plt.figure(2)
plt.plot(n[keep], U[keep])
plt.xlabel("$n$")
plt.ylabel("$\\chi^2_\\mathrm{red}$")
plt.savefig("plots/U.pdf", format = "pdf", bbox_inches = "tight")
plt.close(1)

plt.figure(3)
plt.plot(n[keep], accept[keep])
plt.xlabel("$n$")
plt.ylabel("acceptance rate")
plt.savefig("plots/acceptance.pdf", format = "pdf", bbox_inches = "tight")
plt.close(3)

counter_param = 0
counter_fig = 4
n = n[keep]

for item in parameters:
	item = item[keep]
	plt.figure(counter_fig)
	plt.plot(n, item)
	plt.xlabel("n")
	lab1 = "popt[" + str(counter_param) + "]"
	plt.ylabel(lab1)
	filename1 = "plots/popt" + str(counter_param) + ".pdf"
	plt.savefig(filename1, format = "pdf", bbox_inches = "tight")
	plt.close(counter_fig)
	counter_fig = counter_fig + 1
	
	plt.figure(counter_fig)
	plt.hist(item, n_bins, normed = True)
	plt.xlabel(lab1)
	plt.ylabel("density")
	filename2 = "plots/hist_popt" + str(counter_param) + ".pdf"
	plt.savefig(filename2, format = "pdf", bbox_inches = "tight")
	plt.close(counter_fig)
	counter_fig = counter_fig + 1
	counter_param = counter_param + 1

mean = np.mean(U[keep])
chi2redmin = np.min(U[keep])
diff_theo = p*T/2
print((mean-chi2redmin)/diff_theo )














