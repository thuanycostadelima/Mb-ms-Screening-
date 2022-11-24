#!/usr/bin

import numpy as np
import matplotlib.pyplot as plt
from obspy import read
from obspy.core import Trace, Stream
import matplotlib as mpl
import glob
import os, sys
from obspy.io.sac import SACTrace
import time
from mpl_toolkits.basemap import Basemap
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
from matplotlib import rcParams
start_time = time.time()

###############################################################
### This code plots the magnitude of seismic events according to the \
### screening criterion in subsection 7.3.2 "Ms:mb screening" \
### on page 143 of the CTBTO IDC/OPS manual

### Created on the 26th September 2022
### E-mail: Thuany.CostadeLima@ga.gov.au

###############################################################
## FIGURE LAYOUT

mpl.rcParams['font.size'] = 12
fig = plt.figure(figsize=(12.5, 6.5), layout="tight")
gs = fig.add_gridspec(4,4)

ax = [fig.add_subplot(gs[:,0]), fig.add_subplot(gs[0:2,1]),\
              fig.add_subplot(gs[0:2,2]), fig.add_subplot(gs[2:4,1], projection=ccrs.Robinson()), \
              fig.add_subplot(gs[2:4,2], projection=ccrs.Robinson()), fig.add_subplot(gs[:,3])]


###############################################################
## GLOBAL PARAMETERS

sigmab=0.34
sigmas=0.23


###############################################################
## Reading the information from an input file:

data_file = 'iran_mbms.lst_cut2' 
#data_file = 'iran_mbms.lst_OUTPUT' 

odate = np.loadtxt(data_file, usecols=0, skiprows=4, dtype=str) # reading the Event Origin Time (first collumn)
otime = np.loadtxt(data_file, usecols=1, skiprows=4, dtype=str) # reading the Event Origin Time (second collumn)
ORID=np.loadtxt(data_file, usecols=2, skiprows=4, dtype=float)
evla=np.loadtxt(data_file, usecols=3, skiprows=4, dtype=float)
evlo=np.loadtxt(data_file, usecols=4, skiprows=4, dtype=float)
evdp=np.loadtxt(data_file, usecols=5, skiprows=4,dtype=float)
mb=np.loadtxt(data_file, usecols=6, skiprows=4,dtype=float)
ms=np.loadtxt(data_file, usecols=7, skiprows=4,dtype=float)
ml=np.loadtxt(data_file, usecols=8, skiprows=4,dtype=float)
stnm=np.loadtxt(data_file, usecols=9, skiprows=4, dtype=str)


ms_list = np.array(ms)
mb_list = np.array(mb)
evla_list = np.array(evla)
evlo_list = np.array(evlo)
evdp_list = np.array(evdp)

###############################################################
## Reading the magnitude information of each event based 
## on their ID:
## @@@ Please note the coditions in place:
## 1) at least two Ms stations for each event
## and 2) Mb >=3.5

ev_uniq = list()
for ev in np.unique(ORID):
	ev_uniq.append(ev)
print('Total events:', len(ev_uniq))
print("Startting the screening ...")

score_list = list()
mb_mean_condition_list = list()
ms_mean_condition_list = list()
mb_mean_xcondition_list = list()
ms_mean_xcondition_list = list()
evla_sout = list()
evlo_sout = list()
evla_nsout = list()
evlo_nsout = list()
for index, ev_id in enumerate(ev_uniq):
	mb_ntwk_list = list()
	ms_ntwk_list = list()
	evla_stnm_list = list()
	evlo_stnm_list = list()
	evdp_stnm_list = list()

	## Estimating the average 'mb' and 'Ms', body and surface-wave magnitudes, respectively.
	for n in range(0, len(ms_list)):
		if ORID[n] == ev_uniq[index]:
			
			ms_ntwk_list.append(ms_list[n])
			mb_ntwk_list.append(mb_list[n])

			evla_stnm_list.append(evla_list[n])
			evlo_stnm_list.append(evlo_list[n])
			evdp_stnm_list.append(evdp_list[n])

		else:
			continue
	evla_uniq = np.array(evla_stnm_list)[0]
	evlo_uniq = np.array(evlo_stnm_list)[0]
	evdp_uniq = np.array(evdp_stnm_list)[0]
	#print(evdp_uniq)

	#print("ms and mb measures by each event", ev_id, len(ms_ntwk_list))
	ms_mean = np.mean(np.array(ms_ntwk_list))
	mb_mean = np.mean(np.array(mb_ntwk_list))

	## Conditions for Mb >=3.5, and Ns >=2, and event depth=10km.
	depth_th1 = 0
	depth_th2 = 10
	if mb_mean >= 3.5 and ms_mean != -999 and evdp_uniq >= depth_th1 and evdp_uniq < depth_th2:

		Nb = len(np.array(mb_ntwk_list)) # Number of stations that recorded 'mb'
		Ns = len(np.array(ms_ntwk_list)) # Number of stations that recorded 'Ms'
		
		#print("MS NETWORK by event %s" % ev_id, ms_mean, Ns)
		if Ns > 2: 

			## Defining the Ms:mb screening score:
			sigmam = np.sqrt( np.square(1.25) * np.square(sigmab)/Nb + \
						 np.square(sigmas)/Ns)

			var = mb_mean - ms_mean + 1.96 * sigmam
			
			## Score equation
			SCORE = ((0.64 - mb_mean + ms_mean)\
					/  ( 1.96 * sigmam) ) -1
			#print(ev_id, SCORE)


			## Score equation
			#SCORE = ((2.2 - 1.25*mb_mean + ms_mean)\
			#		/  ( 2 *sigmam ) ) - 1
			
			## Now, the condition below needs to be satisfied for each event:
			if SCORE > 0: #or float(var) < float(0.64):
				
				#print("The condition is only satisfied for event #",index, ev_id, ORID[index]) #,\
					 #sigmam, var, mb_mean, ms_mean, Nb, Ns, SCORE )
				mb_mean_condition_list.append(mb_mean)
				ms_mean_condition_list.append(ms_mean)
				score_list.append(SCORE)
				evla_sout.append(evla_uniq)
				evlo_sout.append(evlo_uniq)
				

			else:
				# for events Not Screened Out 
				mb_mean_xcondition_list.append(mb_mean)
				ms_mean_xcondition_list.append(ms_mean)
				evla_nsout.append(evla_uniq)
				evlo_nsout.append(evlo_uniq)
				
	else:
		pass #print('Mb', mb_mean)

# mean values of magnitude for the events that are screened out
mb_mean = np.array(mb_mean_condition_list)
ms_mean = np.array(ms_mean_condition_list)

# mean values of magnitude for the events that are not screend out
mb_xmean = np.array(mb_mean_xcondition_list)
ms_xmean = np.array(ms_mean_xcondition_list)



###############################################################
## Lets select the list of top 50 most used 
## stations.


print("Selecting top 50 most used stations..")

def my_sort(line):
	line_fields = line.strip().split()
	amount = float(line_fields[0])
	return amount

st_uniq = list()
st_uniq_counts = list()

f = open("STNM_mslist.dat", "w")

for st in np.unique(stnm):
	st_uniq.append(st)
for index, st_id in enumerate(st_uniq):

	st_ms_list = list()
	for n in range(0, len(ms_list)):
		if stnm[n] == st_uniq[index]: st_ms_list.append(ms_list[n])
		else: continue
	st_uniq_counts.append(len(st_ms_list))
	#print("Counts of ms for STNM:", st_uniq[index], len(st_ms_list))
	f.write("%s counts of ms for station %s\n" % (len(st_ms_list), st_uniq[index]))

#print("Median of ms counts:", np.median(st_uniq_counts))
#print("Mean of ms counts:", np.mean(st_uniq_counts))
f.close()

# opening the file to retrieve content into a list
contents =  open("STNM_mslist.dat", 'r').readlines()

# sorting 
contents.sort(key=my_sort, reverse=True)

stnm_top50 = list()
mscounts_top50 = list()
# printing the sorting contents to stdout
for line in contents[:50]:
	#print(line.strip())
	stnm_top50.append(line.strip().split()[6])
	mscounts_top50.append(int(line.strip().split()[0]))
#print(stnm_top50)
#print(mscounts_top50)



###############################################################
# !!!!!!!!!!!!! Make Plots from here !!!!!!!!!!!!!

print("Ready to plot..")
print("Number of events considered",int(len(ms_mean)+len(ms_xmean)))

################################################
### Plotting the mb:ms figure

ax[0].set_title('mb:Ms \n %s events \n EVDP range: %s to %s km' % (int(len(ms_mean)+len(ms_xmean)), depth_th1, depth_th2), fontweight='bold', fontsize=10)

## Defining the range of axis
ymin, ymax, ystep = int(min(ms_mean)-1), int(max(ms_mean)+2), 1
xmin, xmax, xstep = int(min(mb_mean)-1), int(max(mb_mean)+2), 1

#yticks = np.arange(ymin, ymax, ystep).astype(int)
#xticks = np.arange(xmin, xmax, xstep).astype(int)
yticks = np.arange(2, 7, 1).astype(int)
xticks = np.arange(2, 7, 1).astype(int)

ax[0].set_yticks(yticks)
ax[0].set_yticklabels(yticks, rotation='horizontal', va='center')
ax[0].set_ylim(2.0, 7.0) #ax[0].set_ylim(ymin, ymax)
ax[0].set_ylabel('Ms')

ax[0].set_xticks(xticks)
ax[0].set_xticklabels(xticks, rotation='horizontal', va='center')
ax[0].set_xlim(2.0, 7.0) #ax[0].set_xlim(xmin, xmax) 
ax[0].set_xlabel('mb')

for x, y in zip(mb_mean, ms_mean):
	ax[0].plot(x,y, '.', color='black', markerfacecolor='cornflowerblue', markeredgewidth=0.2)
ax[0].plot(mb_mean[0], ms_mean[0],'.',color='black', markerfacecolor='cornflowerblue', markeredgewidth=0.2,label="%s Screened Out" % len(mb_mean)) # label

for x, y in zip(mb_xmean, ms_xmean):
	ax[0].plot(x,y, '.', color='black', markerfacecolor='goldenrod', markeredgewidth=0.2)
ax[0].plot(mb_xmean[0], ms_xmean[0],'.',color='black',markerfacecolor='goldenrod', markeredgewidth=0.2, label="%s Not Screened Out"% len(mb_xmean)) # label



################################################
### Plotting the screening line ms=mb-0.64:

i_list = list()
for i in np.arange(0, 10, 0.1):
	i_list.append(i)
i_list=np.array(i_list)

ax[0].plot(i_list, i_list-0.64, '--', color='black', linewidth=0.8,label="Ms = mb - 0.64")
#ax[0].plot(i_list, i_list-0.64 + 1.96*sigmam, '--', color='black', label="Ms = mb - 0.64 + 1.96*sigmam")

#ax[0].plot(i_list, 1.25*i_list-2.2, '--', color='red', linewidth=0.8, label="Ms = 1.25mb - 2.2")


################################################
#### Making the figure prettier

ax[0].grid(ls='--', color='k', lw=0.5, alpha=0.5)        
#ax[0].annotate('(%c)' % (ord('a')), xy=(-0.05, .98), xycoords='axes fraction', ha='right', fontweight='bold', fontsize=12)
ax[0].set(facecolor = "whitesmoke")
ax[0].legend(loc="upper left", fontsize=10)
ax[1].grid(ls='--', color='k', lw=0.5, alpha=0.5) 
ax[2].grid(ls='--', color='k', lw=0.5, alpha=0.5) 

      
################################################
#### Plotting the histograms of Ms, and Mb
bins=np.linspace(2, 7, 100)

ax[1].hist(ms_mean, bins, alpha=0.5, label='Screened Out')
ax[1].set(facecolor = "whitesmoke")

ax[1].hist(ms_xmean, bins, alpha=0.5, label='Not Screened Out')
ax[1].set(facecolor = "whitesmoke")

ax[1].set_title('Histogram of Ms', fontweight='bold', fontsize=10) # (Screened Out %s events)' % len(ms_mean)) 
ax[1].set_ylabel('Number of events')
ax[1].legend(loc='upper right', fontsize=10)
ax[1].set_xlabel('Ms')

ax[2].hist(mb_mean, bins, alpha=0.5, label='Screened Out')
ax[2].set(facecolor = "whitesmoke")

ax[2].hist(mb_xmean, bins, alpha=0.5, label='Not Screened Out')
ax[2].set(facecolor = "whitesmoke")

ax[2].set_title('Histogram of mb', fontweight='bold', fontsize=10) # (Screened Out %s events)' % len(mb_mean), fontweight='normal', fontsize=12) 
ax[2].set_ylabel('Number of events')
ax[2].legend(loc='upper right', fontsize=10)
ax[2].set_xlabel('mb')

y_pos = np.arange(len(stnm_top50))

ax[5].barh(y_pos, mscounts_top50[::-1], align='center', alpha=0.5, color='blue')
#print(min(np.array(mscounts_top50[::-1])))
#print(max(np.array(mscounts_top50[::-1])))
var=np.arange(min(np.array(mscounts_top50[::-1])), max(np.array(mscounts_top50[::-1])), 1100)
ax[5].set_xticks(var, labels=var)

ax[5].set_yticks(y_pos, labels=stnm_top50[::-1])
ax[5].invert_yaxis()  # labels read top-to-bottom

ax[5].set(facecolor = "whitesmoke")
ax[5].set_ylabel('STNM')
ax[5].set_title("50 most used stations", fontweight='bold', fontsize=10)

ax[5].tick_params(axis='y', labelsize=6)
ax[5].tick_params(axis='x', labelsize=6, rotation=90)

ax[5].grid(ls='--', color='k', lw=0.5, alpha=0.5)


################################################
#### Plotting the map of event location of events 
#### screened out and not screened out

## Option1

ax[3].set_title('Screened out', fontweight='bold', fontsize=10)

pts = set()
for x, y in zip(evla_sout, evlo_sout): # events Screened Out
	print(x,y)
	pts.add((x,y))

for x, y in pts:#
	ax[3].plot(y, x, '*', color='black', alpha=0.6,markerfacecolor='red', markeredgewidth=0.3, transform=ccrs.Geodetic())
ax[3].plot(evlo_sout[0], evla_sout[0], '*', color='black',alpha=0.6, markerfacecolor='red', markeredgewidth=0.3, transform=ccrs.Geodetic(), label='Screened Out')

pts_ = set()

for x_, y_ in zip(evla_nsout, evlo_nsout): # events Nots screened Out
	print(x_,y_)
	pts_.add((x_,y_))


for x, y in pts_:#
	ax[4].plot(y, x, '*', color='black', alpha=0.6, markerfacecolor='Yellow', markeredgewidth=0.3, transform=ccrs.Geodetic(), zorder=100)
ax[4].plot(evlo_nsout[0], evla_nsout[0], '*', color='black', alpha=0.6, markerfacecolor='Yellow', markeredgewidth=0.3, \
	transform=ccrs.Geodetic(), label='Not Screened Out',  zorder=1000)

'''
## Option2 #### Adjusting size of markers of earthquake figure based on magnitude
ax[4].set_title('Not screened out', fontweight='bold', fontsize=10)
x_ev = list()
y_ev = list()
z_ev = list()

#pts_ = set()
for x_, y_, z_ in zip(evla_nsout, evlo_nsout, mb_xmean): # events Nots screened Out
	x_ev.append(x_)
	y_ev.append(y_)
	z_ev.append(z_)

x_ev = np.array(x_ev)
y_ev = np.array(y_ev)
z_ev = np.array(z_ev)


print('min mb:', min(z_ev), 'max mb:',max(z_ev))

for i,z in enumerate(z_ev):#
	if z_ev[i] > 6.0:
		#print("range3:", z_ev[i])
		ax[4].plot(y_ev[i], x_ev[i], '*', markersize=10, \
			color='black', alpha=0.6, markerfacecolor='red', \
			markeredgewidth=0.3, transform=ccrs.Geodetic(),  zorder=1)

	if z_ev[i] <= 6.0 and z_ev[i] > 5.0:
		#print("range3:", z_ev[i])
		ax[4].plot(y_ev[i], x_ev[i], '*', markersize=10, \
			color='black', alpha=0.4, markerfacecolor='red', \
			markeredgewidth=0.3, transform=ccrs.Geodetic(),  zorder=2)

	if z_ev[i] <= 5.0 and z_ev[i] > 4.0:
		#print("range2:", z_ev[i])
		ax[4].plot(y_ev[i], x_ev[i], '*', markersize=10, \
			color='black', alpha=0.6, markerfacecolor='goldenrod', \
			markeredgewidth=0.3, transform=ccrs.Geodetic(),  zorder=3)

	if z_ev[i] <= 4.0:
		#print("range1:", z_ev[i])
		ax[4].plot(y_ev[i], x_ev[i], '*', markersize=10, \
			color='black', alpha=0.6, markerfacecolor='Yellow', \
			markeredgewidth=0.3, transform=ccrs.Geodetic(),  zorder=4)

ax[4].plot(evlo_nsout[0], evla_nsout[0], '*', color='black', alpha=0.6, markerfacecolor='red', markeredgewidth=0.3, \
		transform=ccrs.Geodetic(), label='Ms > 6',  zorder=0)
ax[4].plot(evlo_nsout[0], evla_nsout[0], '*', color='black', alpha=0.4, markerfacecolor='red', markeredgewidth=0.3, \
		transform=ccrs.Geodetic(), label='5 < Ms <= 6 ',  zorder=0)
ax[4].plot(evlo_nsout[0], evla_nsout[0], '*', color='black', alpha=0.6, markerfacecolor='goldenrod', markeredgewidth=0.3, \
		transform=ccrs.Geodetic(), label='4 < Ms <= 5',  zorder=0)
ax[4].plot(evlo_nsout[0], evla_nsout[0], '*', color='black', alpha=0.6, markerfacecolor='Yellow', markeredgewidth=0.3, \
		transform=ccrs.Geodetic(), label='Ms <= 4',  zorder=0)



ax[3].set_title('Screened out', fontweight='bold', fontsize=10)
x_ev = list()
y_ev = list()
z_ev = list()

#pts_ = set()
for x_, y_, z_ in zip(evla_sout, evlo_sout, mb_mean): # events Nots screened Out
	x_ev.append(x_)
	y_ev.append(y_)
	z_ev.append(z_)

x_ev = np.array(x_ev)
y_ev = np.array(y_ev)
z_ev = np.array(z_ev)

print('min mb:', min(z_ev), 'max mb:',max(z_ev))

for i,z in enumerate(z_ev):#
	if z_ev[i] > 6.0:
		#print("range3:", z_ev[i])
		ax[3].plot(y_ev[i], x_ev[i], '*', markersize=10, \
			color='black', alpha=0.6, markerfacecolor='red', \
			markeredgewidth=0.3, transform=ccrs.Geodetic(),  zorder=10)

	if z_ev[i] <= 6.0 and z_ev[i] > 5.0:
		#print("range3:", z_ev[i])
		ax[3].plot(y_ev[i], x_ev[i], '*', markersize=10, \
			color='black', alpha=0.4, markerfacecolor='red', \
			markeredgewidth=0.3, transform=ccrs.Geodetic(),  zorder=12)

	if z_ev[i] <= 5.0 and z_ev[i] > 4.0:
		#print("range2:", z_ev[i])
		ax[3].plot(y_ev[i], x_ev[i], '*', markersize=10, \
			color='black', alpha=0.6, markerfacecolor='goldenrod', \
			markeredgewidth=0.3, transform=ccrs.Geodetic(),  zorder=13)

	if z_ev[i] <= 4.0:
		#print("range1:", z_ev[i])
		ax[3].plot(y_ev[i], x_ev[i], '*', markersize=10, \
			color='black', alpha=0.6, markerfacecolor='Yellow', \
			markeredgewidth=0.3, transform=ccrs.Geodetic(),  zorder=14)

ax[3].plot(evlo_sout[0], evla_sout[0], '*', color='black', alpha=0.6, markerfacecolor='red', markeredgewidth=0.3, \
		transform=ccrs.Geodetic(), label='Ms > 6',  zorder=0)
ax[3].plot(evlo_sout[0], evla_sout[0], '*', color='black', alpha=0.4, markerfacecolor='red', markeredgewidth=0.3, \
		transform=ccrs.Geodetic(), label='5 < Ms <= 6 ',  zorder=0)
ax[3].plot(evlo_sout[0], evla_sout[0], '*', color='black', alpha=0.6, markerfacecolor='goldenrod', markeredgewidth=0.3, \
		transform=ccrs.Geodetic(), label='4 < Ms <= 5',  zorder=0)
ax[3].plot(evlo_sout[0], evla_sout[0], '*', color='black', alpha=0.6, markerfacecolor='Yellow', markeredgewidth=0.3, \
		transform=ccrs.Geodetic(), label='Ms <= 4',  zorder=0)
## ENd of option2 ####
'''
####
for i in np.array([3, 4]):
	ax[i].legend(fontsize=10)
	#ax[3].set_extent([20, 80, 10, 55])
	ax[i].set_extent([40, 65, 20, 49])
	#ax[i].set_title('Map of earthquakes', fontweight='bold', fontsize=10)
	gl = ax[i].gridlines(draw_labels=True, xlocs=[40, 60], linestyle='--') #xlocs=[30, 90]
	gl.top_labels = False
	gl.right_labels = False
	ax[i].stock_img()
	ax[i].legend(loc='upper right', fontsize=6, bbox_to_anchor=(1.21, 1.03), fancybox=True, shadow=True, facecolor="white")


## adding shape boundaries
import shapefile
from mapping_tools import drawshapepoly, labelpolygon
wd=os.getcwd()

shpfile = wd+'/tectonicplates-master/PB2002_plates.shp'
sf = shapefile.Reader(shpfile)
for shape in sf.shapeRecords():
   xy = [i for i in shape.shape.points[:]]
   x = [i[0] for i in xy]
   y = [i[1] for i in xy]
   ax[3].plot(x, y, markeredgewidth=0.3, transform=ccrs.Geodetic())


shpfile2 = wd+'/tectonicplates-master/PB2002_plates.shp'
sf2 = shapefile.Reader(shpfile2)
for shape in sf2.shapeRecords():
   xy = [i for i in shape.shape.points[:]]
   x = [i[0] for i in xy]
   y = [i[1] for i in xy]
   ax[4].plot(x, y, markeredgewidth=0.3, transform=ccrs.Geodetic())

# change colour, add condition to mb:ms for reagion; add label of tectonic domains;

################################################
#### Labels of the figures
axLabels = ('(a)','(b)','(c)','(d)', '(e)', '(f)')
for i, ax in enumerate(fig.get_axes()):
	ax.text(-0.25, 1, axLabels[i], transform=ax.transAxes, size=14, weight='bold')

###############################################################
####  Saviing the output and show the result.
plt.savefig('/Users/thuanycostadelima/Desktop/mbms_evdpRange_%s_%sKM.pdf' % (depth_th1, depth_th2), dpi=300)
plt.show()
