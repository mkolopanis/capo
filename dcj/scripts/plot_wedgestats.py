#! /usr/bin/env python

#plot the output of wedgestats.py
#columns
#obsid, wedge_x,wedge_y,window_x,window_y,residual_wedge_x,residual_wedge_y,galactic_x,galactic_y,point_source_x,point_src_y
import numpy as n
from pylab import *
import sys,ephem,optparse
from astropy import time
from matplotlib.ticker import MultipleLocator,FuncFormatter


o = optparse.OptionParser()
o.set_usage('plot_wedgestats.py [options] statsfile.txt')
o.set_description(__doc__)
o.add_option('--flags',type=str,
        help='text file listing obsid,flagx,flagy where 0=no flag, 1=flag')
o.add_option('--pointingfile',type=str,
        help='text file listing obsid, jd, and pointing index')
o.add_option('--subset',type=str,
        help='file with a list obsids to plot')
o.add_option('--scale',type=float,default=0.5,
        help='yscale scale. how big in the vertical? defaut=0.5')
opts, args = o.parse_args(sys.argv[1:])

majorlocator = MultipleLocator(1)
minorlocator = MultipleLocator(10/60.)
def formathours(lst,pos):
    return '{h}:{m}'.format(h=int(lst),m=n.modf(lst)*60)
def fitsmoothpoly(x,y,windowpoints=10,order=3):
    #fit a smooth polynomial of the given order over a window windowpoints wide
    assert(len(x)==len(y))
    assert(windowpoints>order)
    Npoints = len(x)
    model=[]
    #first handle the window edge where we don't have enough data to do the fit. just insert the input values
    for i in range(0,order):
        model.append([x[i],y[i],0])
    #next handle where we have enough points to fit but not enough to to do the full window
    for i in range(order,windowpoints/2):
        x_window,y_window = x[i:i+windowpoints/2],y[i:i+windowpoints/2]
        P = n.poly1d(n.polyfit(x_window,y_window,order))
        res = n.sqrt(n.var(y_window-P(x[i])))
        model.append([x[i],P(x[i]),res])
    #fit over the majority of the data
    for i in range(windowpoints/2,Npoints-windowpoints/2):
        x_window,y_window = x[i-windowpoints/2:i+windowpoints/2],y[i-windowpoints/2:i+windowpoints/2]
        P = n.poly1d(n.polyfit(x_window,y_window,order))
        res = n.sqrt(n.var(y_window-P(x[i])))
        model.append([x[i],P(x[i]),res])
    #close out the fit on the end
    for i in range(Npoints-windowpoints/2,Npoints-order):
        x_window,y_window = x[i:i+windowpoints/2],y[i:i+windowpoints/2]
        P = n.poly1d(n.polyfit(x_window,y_window,order))
        res = n.sqrt(n.var(y_window-P(x[i])))
        model.append([x[i],P(x[i]),res])
    for i in range(Npoints-order,Npoints):
        model.append([x[i],y[i],0])
    assert(len(model)==len(x))#always return a matching x axis
    return n.array(model)
if not opts.flags is None:#load the txt flag file generated by filter_wedgestats.py
    lines = open(opts.flags).readlines()
    flags = {}
    for l in lines:
        if l.startswith('#'):continue   
        obsid = int(l.split()[0].strip()),
        flags[obsid] = map(int,l.split()[1:])
else:flags={}
minorformatter = FuncFormatter(formathours)
O = ephem.Observer()
O.lat = '-26:42:11.95'
O.lon = '116:40:14.92'
def t2lst_locs(t):
    lsts = []
    for myt in t:
        O.date = myt.jd-2415020 #convert to dublin jd via wikipedia
        lsts.append(float(O.sidereal_time()))
    lsts = n.array(lsts)
    lst_locs = (lsts*12/n.pi + 12)%24 -12 
    return lst_locs
def t2lst(t):
    #input a astropy.time.Time
    #output lst in hours
    #vectorized
    lsts = []
    for myt in t:
        O.date = myt.jd-2415020 #convert to dublin jd via wikipedia
        lsts.append(float(O.sidereal_time()))
    lsts = n.array(lsts)*12/n.pi
    return lsts

for filename in args:
    print "loading ",filename
    D = n.loadtxt(filename)
    #use the raw data to set the scale
    scale = n.std(D[:,3:5])*opts.scale
    MIN = D[:,3:5].min()
    MAX = D[:,3:5].max()
    MEAN = n.mean(D[:,3:5])
    TRANGE = [time.Time(D[:,0].min(),format='gps'),time.Time(D[:,0].max(),format='gps')]
    LSTLOCRANGE = t2lst_locs(TRANGE)
    LSTLOCRANGE = [n.floor(LSTLOCRANGE[0]),n.ceil(LSTLOCRANGE[1])]
    #keep only those in the input subset
    if not opts.subset is None:
        subset = n.loadtxt(opts.subset)
        data = []
        for i,obsid in enumerate(D[:,0]):
            if obsid in subset:
                data.append(D[i,:])
        D = n.array(data)
    t = time.Time(D[:,0],format='gps')

    lst_locs = t2lst_locs(t)#lst_locs are lsts that go from -12 to 12. good for plotting across the divide
    dt = time.TimeDelta(0.5,format='jd')
    nights = n.argwhere(n.abs(t[1:].jd-t[:-1].jd)>0.5).squeeze()+1


    if True:
        figure(1)
        clf()
        #plot the window power vs time
        title('window_power')
        plot_date(t.plot_date,D[:,3],'.',label='x')
        plot_date(t.plot_date,D[:,4],'.',label='y')
        #xwindowmodel = fitsmoothpoly(t.jd,D[:,3],windowpoints=15,order=2)
        #ywindowmodel = fitsmoothpoly(t.jd,D[:,4]rwindowpoints=15,order=2)
        xlabel('observation date')
        ylabel('arbitrary power')
        #vlines(t[nights].plot_date,0,D[:,3-5].max())



    #turn the lsts inside out cause we're never going to cross 12 hrs in eor work
    figure(2)
    clf()
    title("window power")
    axx = subplot(121)
    axx.set_ylim([MEAN-len(nights)*scale,MAX])
    axy = subplot(122)
    axy.set_ylim([MEAN-len(nights)*scale,MAX])
    for i in xrange(len(nights)):
        if i==0:
            inds=range(0,nights[0])
        else:inds = range(nights[i-1],nights[i])
        axx.plot(lst_locs[inds],D[inds,3]-i*scale,'k',lw=0.5)
        axx.plot(lst_locs[inds],D[inds,3]-i*scale,'.k',ms=3)
        axy.plot(lst_locs[inds],D[inds,4]-i*scale,'k',lw=0.5)
        axy.plot(lst_locs[inds],D[inds,4]-i*scale,'.k',ms=3)
    for f in flags:
        print "plotting flags",f
        myt = n.argwhere(f==D[:,0]).squeeze()
        if flags[f][0]==1:
            print myt,'x'
            plot(lst_locs[myt],D[myt,3]-i*scale,'ob',mfc='None')
        if flags[f][1]==1:
            print myt,'y'
            plot(lst_locs[myt],D[myt,4]-i*scale,'og',mfc='None')
    #Put in lines at the pointing changes    
    if not opts.pointingfile is None:
        #load a list of pointings for some obsids supposedly in this file. 
        # go through and find out if all my obsids have a pointing 
        pointingdata = n.loadtxt(opts.pointingfile)
        pointingchanges = n.argwhere(n.logical_and(n.diff(pointingdata[:,2]),n.diff(pointingdata[:,1])<200)).squeeze()
        lstbins = n.linspace(0,24,num=50)
        pointingchange_lsts = t2lst(time.Time(pointingdata[pointingchanges,0],format='gps'))
        pointingchange_lst_locs = (pointingchange_lsts+12)%24-12
        axx.vlines(pointingchange_lst_locs,MIN-len(nights)*scale,MAX,alpha=0.3)
        axy.vlines(pointingchange_lst_locs,MIN-len(nights)*scale,MAX,alpha=0.3)



    #if there are more than two hours
    #set the formating
    axx.set_xlim(LSTLOCRANGE)
    axy.set_xlim(LSTLOCRANGE)
    axx.grid()
    axy.grid()
#    axy.yaxis.set_ticklabels([])
    axx.xaxis.set_major_locator(majorlocator)
    axx.xaxis.set_minor_locator(minorlocator)
    axy.xaxis.set_major_locator(majorlocator)
    axy.xaxis.set_minor_locator(minorlocator)
    #axy.set_xticks(axy.get_xticks()[1:]) #remove an overlapping number
    lst_ticks = axx.get_xticks()
    lst_ticks[lst_ticks<0] += 24
    axx.set_xticklabels(lst_ticks.astype(int))
    axy.set_xticklabels(lst_ticks.astype(int))

    axx.set_xlabel('lst [hours]')
    axy.set_xlabel('lst [hours]')
    axx.set_ylabel('Jy$^2$')
    axx.set_title('x')
    axy.set_title('y')
    #subplots_adjust(wspace=0.001)
    tight_layout()



    show()
    sys.exit()
    figure(3)
    clf()
    ax = subplot(111)
    title('window/galaxy')
    plot(lst_locs,D[:,3]/D[:,7],'.')
    plot(lst_locs,D[:,4]/D[:,8],'.')
    ax.xaxis.set_major_locator(majorlocator)
    ax.xaxis.set_minor_locator(minorlocator)
    lst_ticks = ax.get_xticks()
    lst_ticks[lst_ticks<0] += 24
    ax.set_xticklabels(lst_ticks.astype(int))
    ylabel('window leakage')
    xlabel('lst [hours]')
    lst_ticks = ax.get_xticks()
    lst_ticks[lst_ticks<0] += 24
    #if there are more than two hours
    draw()
    show()
