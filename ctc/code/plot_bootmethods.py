#! /usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt

# Files to read
boot_bl = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/PSA64_FRF_RA.5_8.6_CHAN95_115_SEP0,1_IDENTITYMULTWEIGHT_WEIGHTI/pspec_final_sep0,1.npz') # bootstrap over baseline only
#boot_bl = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/PSA64_FRF_RA.5_8.6_CHAN95_115_SEP0,1_RANGEOFPROJECTEDMODES/project_3_modes/pspec_final_sep0,1.npz') # bootstrap over baseline only
boot_time = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/PSA64_FRF_RA.5_8.6_CHAN95_115_SEP0,1_IDENTITYMULTWEIGHT_WEIGHTI/pspec_final_sep0,1_oldboot.npz') # bootstrap over time and bl


noise = 4436767.36822 # XXX

# Plot
plt.plot(boot_bl['kpl'], boot_bl['pIn_err']*2, color='black', label='Bootstrap baselines only')
plt.plot(boot_time['kpl'], boot_time['pIn_err']*2, color='0.5', label='Bootstrap baselines and times')
plt.axhline(y=noise*2, color='green', linestyle='-', label='Analytic')
plt.legend(numpoints=1,prop={'size':14},loc='best')
plt.grid()
plt.yscale('log')
#plt.title('2$\sigma$ errors for PAPER-64 using different bootstrapping methods')
plt.xlabel('$k_{\\parallel}$ [$h$ Mpc$^{-1}$]',fontsize=18)
plt.ylabel('$P(k)$ $[mK^{2}(h^{-1} Mpc)^{3}]$',fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=14)
#plt.ylim(1e6,7e9) # for data
plt.ylim(5e5,1e8)
plt.show()


