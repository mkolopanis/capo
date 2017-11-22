#! /usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt

# Files to read
#blue = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/PSA64_FRF_RA.5_8.6_CHAN95_115_SEP0,1_IDENTITYMULTWEIGHT_SUBOPTIMALBLGROUPING/pspec_final_sep0,1_final.npz') # avg over time; suboptimal bl bootstrapping
blue = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/PSA64_FRF_RA.5_8.6_CHAN95_115_SEP0,1_IDENTITYMULTWEIGHT_SUBOPTIMALBLGROUPING/inject_sep0,1_0.01/pspec_pk_k3pk.npz') # avg over time; suboptimal bl bootstrapping
#black = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/PSA64_FRF_RA.5_8.6_CHAN95_115_SEP0,1_IDENTITYMULTWEIGHT/pspec_final_sep0,1_final.npz') # avg over time; bootstraps last baseline slot only
black = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/PSA64_FRF_RA.5_8.6_CHAN95_115_SEP0,1_IDENTITYMULTWEIGHT/inject_sep0,1_0.01/pspec_pk_k3pk.npz') # avg over time; bootstraps last baseline slot only
#green = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/PSA64_FRF_RA.5_8.6_CHAN95_115_SEP0,1_IDENTITYMULTWEIGHT_SUBOPTIMALBLGROUPING/pspec_final_sep0,1_oldboot.npz') # bootstrap over time and bl; suboptimal bl bootstrapping
green = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/PSA64_FRF_RA.5_8.6_CHAN95_115_SEP0,1_IDENTITYMULTWEIGHT_SUBOPTIMALBLGROUPING/inject_sep0,1_0.01/pspec_pk_k3pk_oldboot.npz') # bootstrap over time and bl; suboptimal bl bootstrapping
#red = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/PSA64_FRF_RA.5_8.6_CHAN95_115_SEP0,1_IDENTITYMULTWEIGHT/pspec_final_sep0,1_oldboot.npz') # bootstrap over time and bl; bootstraps last baseline slot only
red = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/PSA64_FRF_RA.5_8.6_CHAN95_115_SEP0,1_IDENTITYMULTWEIGHT/inject_sep0,1_0.01/pspec_pk_k3pk_oldboot.npz') # bootstrap over time and bl; bootstraps last baseline slot only


noise = 4436767.36822 # XXX

# Plot
plt.plot(blue['kpl'], blue['pIv_err']*2, color='blue', label='Bootstrap baselines only; suboptimal baseline sampling')
plt.plot(black['kpl'], black['pIv_err']*2, color='black', label='Bootstrap baselines only; optimal baseline sampling')
plt.plot(red['kpl'], red['pIv_err']*2, color='red', label='Bootstrap baselines and times; optimal baseline sampling')
plt.axhline(y=noise*2, color='green', linestyle='-', label='Analytic')
plt.legend(numpoints=1,prop={'size':12},loc='best')
plt.grid()
plt.yscale('log')
plt.title('2$\sigma$ errors for PAPER-64 using different bootstrapping methods')
plt.xlabel('$k_{\\parallel}$ [$h$ Mpc$^{-1}$]')
plt.ylabel('$P(k)$ $[mK^{2}(h^{-1} Mpc)^{3}]$')
plt.ylim(1e6,5e9)
plt.show()


