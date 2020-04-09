import numpy as np 
import matplotlib.pyplot as plt 


data = np.loadtxt("output/lattice_dynamics.txt", unpack =True)

range_min = 0
range_max = 1000

thermalization = 20 

plt.figure(1)
plt.plot(data[0,range_min:range_max], data[2,range_min:range_max])
plt.savefig("plots/a.pdf", format="pdf", bbox_inches="tight")
plt.close(1)


plt.figure(4)
plt.plot(data[0,range_min:range_max], data[-2,range_min:range_max])
plt.savefig("plots/S.pdf", format="pdf", bbox_inches="tight")
plt.close(4)

plt.figure(5)
plt.plot(data[0,range_min:range_max], data[-1,range_min:range_max])
plt.savefig("plots/accept.pdf", format="pdf", bbox_inches="tight")
plt.close(5)

# times = np.loadtxt("output/autocorrel_times.txt", unpack=True)
# L_max = 199

# plt.figure(6)
# for i in range(1, len(times[:,0])-1):
# 	plt.plot(times[-1,0:L_max], times[i, 0:L_max])
# #plt.legend(loc="best")
# plt.savefig("plots/autocorr.pdf", format="pdf", bbox_inches="tight")
# plt.close(6)


t, G, G_err = np.loadtxt("output/correlator.txt", unpack=True)
def exponential(t, a, b):
	return a*np.exp(-b*t)
from scipy.optimize import curve_fit

popt, pcov =  curve_fit(exponential, t[:int(len(t)/2)], G[:int(len(t)/2)], sigma=G_err[:int(len(t)/2)], absolute_sigma=True, p0=[1,1])
t_values = np.linspace(0, len(t)/2, 100)
chi2_red = np.sum(((G[:int(len(t)/2)]-exponential(t[:int(len(t)/2)],*popt))/G_err[:int(len(t)/2)])**2)/(len(t[:int(len(t)/2)])-2)
print("chi2_red", chi2_red)
a = 0.1
plt.figure(7)
plt.scatter(t, G, s=10)
plt.plot(t_values, exponential(t_values, *popt), label="$E_1-E_0 = $ {0} +- {1}".format(round(popt[1]/a,4), np.round(np.sqrt(pcov[1][1])/a,4)))
plt.legend(loc="best")
#plt.yscale("log")
#plt.ylim(1e-4,1e2)
plt.errorbar(t, G, yerr=G_err, linestyle = 'none', elinewidth=2, capsize = 6, capthick = 2)
plt.savefig("plots/correlator.pdf", format="pdf", bbox_inches="tight")
plt.close(7)

# plt.figure(8)
# plt.errorbar(t[:-1], np.log(G[:-1]/G[1:])/a, yerr=np.sqrt((G_err[:-1]/G[:-1])**2+(G_err[1:]/G[1:])**2)/a, linestyle = 'none', elinewidth=2, capsize = 6, capthick = 2)
# plt.xlim(0,10)
# plt.ylim(0,2)
# plt.savefig("plots/correlator_dE.pdf", format="pdf", bbox_inches="tight")
# plt.close(8)
