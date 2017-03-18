import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from numpy.random import random
sns.set_context("paper", font_scale=2)

#FILE = 'corr_res.csv'
FILES = ['HERA_243_pm.csv', 'HERA_128_pm.csv', 'PAPER_128_pm.csv']
LABELS = ['HERA243', 'HERA128', 'PAPER128']

def gen_color(l=1):
	colors = []
	for i in range(l): colors.append((random(),random(),random()))
	return np.array(colors)
COLORS = gen_color(len(FILES))

def pairplot(Theta_min=0):
	dflist = []
	for i, file in enumerate(FILES):
		df = pd.read_csv(file)
		df['label'] = LABELS[i]
		df['peak'] /= np.amax(df['peak'])
		df['rho0'] = 0.001*40/df['bl1']
		df['Theta'] = np.sqrt(df['mult'])*df['peak']/np.sqrt(1+df['rho0']*2*np.sqrt(df['mult']))
		df['Theta'] /= np.amax(df['Theta'])
		df = df.loc[df['Theta']>Theta_min]
		dflist.append(df)

	df = pd.concat(dflist)

	g = sns.pairplot(df,hue='label',vars=['dT','peak','Theta','bl1','bl2'],
		plot_kws={'alpha':0.2, "s":30})
	# g = sns.PairGrid(df)
	# g = g.map_diag(sns.kdeplot, lw=3)
	# g = g.map_offdiag(sns.kdeplot, lw=1)

def get_imp(df, Theta_min=0.0):
	dft = df.loc[df['Theta']>Theta_min]
	dfeq = dft.loc[dft['sep']==dft['sep2']]
	dfnq = dft.loc[dft['sep']!=dft['sep2']]
	
	totalsumsq = np.sum(dft['Theta']**2)
	eqsumsq = np.sum(dfeq['Theta']**2)
	totalsum = np.sum(dft['Theta'])
	eqsum = np.sum(dfeq['Theta'])
	totalsens = totalsumsq/totalsum*np.sqrt(len(dft.index))
	eqsens = eqsumsq/eqsum*np.sqrt(len(dfeq.index))
	improve = (totalsens-eqsens)/eqsens
	return totalsens, eqsens

def sensplot():
	print "========= Statistics of sensitibity contribution =========="
	plt.figure()
	for i, file in enumerate(FILES):
		df = pd.read_csv(file)
		df['rho0'] = 0.001*40/df['bl1']
		df['Theta'] = np.sqrt(df['mult'])*df['peak']/np.sqrt(1+df['rho0']*2*np.sqrt(df['mult']))
		df['Theta'] /= np.amax(df['Theta'])
		#df['Theta'] /= 1000

		TL = np.arange(0,1,0.01)
		 
		imp = np.array([get_imp(df, l) for l in TL])
		totalsensL, eqsensL = imp.T
		plt.plot(TL, totalsensL, label=LABELS[i], color=COLORS[i], linewidth=3)
		plt.plot(TL, eqsensL,'--', color=COLORS[i], linewidth=3)
	plt.legend()


if __name__=="__main__":
	sensplot()
	pairplot(0.1)
	plt.show()
	#import IPython; IPython.embed()

