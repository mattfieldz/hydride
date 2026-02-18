

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import math

# log normal distribution of oxide depth
# mean = exp(3) \approx 20 nm
# ad-hoc standard deviation - gives some likely points down to a few nm
mu, sigma = np.log(18), 0.24


# usual log normal PDF
def logNormalPDF(x):
    return np.exp(-((np.log(x) - mu) ** 2) / (2 * sigma**2)) / (
        x * sigma * np.sqrt(2 * np.pi)
    )


# usual log normal CDF
def logNormalCDF(x):
    return 0.5 * (1 + sp.special.erf((np.log(x) - mu) / (math.sqrt(2.0) * sigma)))


# seeding makes it repeatable
np.random.seed(0)

# Coupon size is 1.2cm by 1.2cm : Bazley et al. 2012 SSI paper
# 1.2x1.2 = 1.44 cm^2 coupon
# each "grain" is a sample with a different Loxide depth
# ad hoc choice of grain size: 30umx30um = 900 um^2 area grain
# see Jones et al. 2013 SSI paper fig 5.
# 1.44cm^2 / 900um^2 grain size means 160000 grains in a coupon
Npts = 160000

# for each grain, sample from a log normal distribution
Loxide = np.random.lognormal(mu, sigma, Npts)
# sort the list of oxide depths into increasing order
Loxide.sort()
# Nsub = number of points taken from the left-tail of the overall distribution
# these are the points we expect to initiate first given their smallest depth
Nsub = 1000
print(f"The left-tail {Nsub} grains/points with least oxide depth give (nm):")
print(Loxide[:Nsub])

# Sanity check:
# generate a histogram of data for the Loxide list of samples
# compare it to the expected log-normal PDF
count, bins, ignored = plt.hist(Loxide, 100, density=True, align="mid")
x = np.linspace(1, 50, 400)  # plot from 1 to 50nm
plt.plot(x, logNormalPDF(x), linewidth=2, color="r", label="expected PDF")
plt.title("Ensemble of samples from lognormal - sanity check")
plt.legend()
# plt.xscale("log")
plt.xlabel("Oxide thickness (nm)")
plt.ylabel("Prob. density")
plt.xlim((1, 50))  # up to 50nm shown
plt.savefig("formal_work/DATA_COMP/histogram.png", dpi=200)
plt.clf()

# plot the oxide depth Loxide vs the point index for the first Nsub points
# repeat this process for 10 different virtual 'coupons'
for i in range(0, 10):
    Loxide = np.random.lognormal(mu, sigma, Npts)
    Loxide.sort()
    plt.scatter(np.linspace(1, Nsub, Nsub), Loxide[:Nsub], alpha=0.3)
# this should reproduce the CDF for a log-normal distribution
x = np.linspace(0.1, 10, 400)  # 0.1 to 10nm oxide depths
plt.plot(logNormalCDF(x) * Npts, x, color="r", label="expected from CDF")
plt.title(f"The smallest {Nsub} sampled oxide depths from 10 coupons")
plt.legend()
plt.xlim((-5, 1.2 * Nsub))
plt.ylim((0, 1.2 * Loxide[Nsub]))
plt.xlabel("Ordered index number of sample point")
plt.ylabel("Oxide thickness (nm)")
plt.savefig("formal_work/DATA_COMP/Loxide_vs_nindex.png", dpi=200)
plt.clf()

# Glascott values
m=1

Dm = 2.19e-9  # cm^2/s
Do = 3.86e-10  # cm^2/s
c0 = 1e-7
c0ck = 1e-5

CHS = 3.86e-9

alpha_c = 0.98
N2star = 8.03e-2
reactK = 8.11e3
eps = c0 / N2star

CHSnon = CHS/c0  # terminal solubility relative to surface forcing
# construct a fake initiation time for each Loxide value using Glascott's expression
for i in range(0, 10):
    Loxide = np.random.lognormal(mu, sigma, Npts)
    Loxide.sort()
    tg = (math.pi * Dm / 4) * (CHSnon * Loxide[:Nsub] * 1.0e-7 / Do) ** 2
    # t2a = 3*(1e-7 * Loxide[:Nsub] * 1000 * 1e-7) / (Do * eps * (1-CHSnon))
    t2a = 5000 * (Loxide[:Nsub]/5)
    t3 = 5000 * 1e-7 / (np.sqrt((2*Dm * reactK * (1-CHS/c0ck)**(m+1) * c0ck**(m+1)) / (3*N2star**2 * (1-alpha_c)**2 * (m+1)) ))
    print(tg[0],t2a[0],t3)
    tt = tg + t2a + t3
    plt.scatter(tt, np.linspace(1, Nsub, Nsub), alpha=0.2)
# add t^3.3 line
tvalues = np.linspace(1, 1000000, 1000)
plt.plot(tvalues, 1e-16 * np.power(tvalues, 3.3), color="r", label="t^3.3")
plt.title(f"Time (Glascott) to initiation for the first {Nsub} points from 10 coupons")
plt.xlabel("Time to initiation (sec)")
plt.ylabel("Total number of initiated points")
plt.xscale("log")
plt.yscale("log")
plt.ylim((0.5, 10000))
# predict the initiation time using the CDF instead of a t^g power
depths = np.linspace(1, 20, 100)
cdfpredictions = logNormalCDF(depths)
tnpredictions = (math.pi * Dm / 4) * (CHSnon * depths * 1.0e-7 / Do) ** 2
# plt.plot(tnpredictions, cdfpredictions * Npts, color="b", label="Prediction using CDF")

print('lol',np.mean(Loxide),np.std(Loxide))

plt.legend()
plt.savefig("formal_work/DATA_COMP/tn.png", dpi=200)
