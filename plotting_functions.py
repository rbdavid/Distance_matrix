#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# USAGE:
# from fn_plotting.py import *

# PREAMBLE:

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
from matplotlib.ticker import NullFormatter

stdev = np.std
sqrt = np.sqrt
nullfmt = NullFormatter()

# ----------------------------------------
# PLOTTING SUBROUTINES

def make_colormap(seq):
	"""Return a LinearSegmentedColormap
	seq: a sequence of floats and RGB-tuples. The floats should be increasing
	and in the interval (0,1).
	"""
	seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
	cdict = {'red': [], 'green': [], 'blue': []}
	for i, item in enumerate(seq):
		if isinstance(item, float):
			r1, g1, b1 = seq[i - 1]
			r2, g2, b2 = seq[i + 1]
			cdict['red'].append([item, r1, r2])
			cdict['green'].append([item, g1, g2])
			cdict['blue'].append([item, b1, b2])
	return mcolors.LinearSegmentedColormap('CustomMap', cdict)

def plot_1d(xdata, ydata, color, x_axis, y_axis, system, analysis, average = False, t0 = 0, **kwargs):
	""" Creates a 1D scatter/line plot:

	Usage: plot_1d(xdata, ydata, color, x_axis, y_axis, system, analysis, average = [False|True], t0 = 0)
	
	Arguments:
	xdata, ydata: self-explanatory
	color: color to be used to plot data
	x_axis, y_axis: strings to be used for the axis label
	system: descriptor for the system that produced the data
	analysis: descriptor for the analysis that produced the data
	average: [False|True]; Default is False; if set to True, the function will calc the average, standard dev, and standard dev of mean of the y-data
	t0: index to begin averaging from; Default is 0
	
	kwargs:
		xunits, yunits: string with correct math text describing the units for the x/y data
		x_lim, y_lim: list w/ two elements, setting the limits of the x/y ranges of plot
		plt_title: string to be added as the plot title
		draw_line: int value that determines the line style to be drawn; giving myself space to add more line styles if I decide I need them

	"""
	# INITIATING THE PLOT...
	plt.plot(xdata, ydata, '%s' %(color))

	# READING IN KWARG DICTIONARY INTO SPECIFIC VARIABLES
	for name, value in kwargs.items():
		if name == 'xunits':
			x_units = value
			x_axis = '%s (%s)' %(x_axis, value)
		elif name == 'yunits':
			y_units = value
			y_axis = '%s (%s)' %(y_axis, value)
		elif name == 'x_lim':
			plt.xlim(value)
		elif name == 'y_lim':
			plt.ylim(value)
		elif name == 'plt_title':
			plt.title(r'%s' %(value), size='14')
		elif name == 'draw_line':
			draw_line = value
			if draw_line == 1:
				plt.plot([0,max(ydata)],[0,max(ydata)],'r-',linewidth=2)
			else:
				print 'draw_line = %s has not been defined in plotting_functions script' %(line_value)
	
	plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
	plt.xlabel(r'%s' %(x_axis), size=12)
	plt.ylabel(r'%s' %(y_axis), size=12)

	# CALCULATING THE AVERAGE/SD/SDOM OF THE Y-DATA
	if average != False:
		avg = np.sum(ydata[t0:])/len(ydata[t0:])
		SD = stdev(ydata[t0:])
		SDOM = SD/sqrt(len(ydata[t0:]))

		plt.axhline(avg, xmin=0.0, xmax=1.0, c='r')
		plt.figtext(0.680, 0.780, '%s\n%6.4f $\\pm$ %6.4f %s \nSD = %4.3f %s' %(analysis, avg, SDOM, y_units, SD, y_units), bbox=dict(boxstyle='square', ec='r', fc='w'), fontsize=12)

	plt.savefig('%s.%s.plot1d.png' %(system,analysis),dpi=300)
	plt.close()


def hist1d(data, x_axis, system, analysis, num_b = 100, norm = False, average = False, t0 = 0, **kwargs):
	""" Creates a 1D histogram:

	Usage: hist1d(data, x_axis, num_b, system, analysis, norm)
	
	Arguments:
	data: self-explanatory
	x_axis: string to be used for the axis label
	system: descriptor for the system analyzed
	analysis: descriptor for the analysis performed and plotted
	num_b: number of bins to be used when binning the data; Default is 100
	norm = [False][True]; Default is False; if False, plotting a frequency of data; if True, plotting a probability density
	average: [False|True]; Default is False; if set to True, the function will calc the average, standard dev, and standard dev of mean of the y-data
	t0: index to begin averaging from; Default is 0

	kwargs:
		xunits: string with correct math text describing the units for the x data
		x_lim, y_lim: list w/ two elements, setting the limits of the x/y ranges of plot
		plt_title: string to be added as the plot title

	"""
	
	# INITIATING THE PLOT...
	events, edges, patches = plt.hist(data, bins=num_b, histtype = 'bar', normed=norm)
	
	# READING IN KWARG DICTIONARY INTO SPECIFIC VARIABLES
	for name, value in kwargs.items():
		if name == 'xunits':
			x_units = value
			x_axis = '%s (%s)' %(x_axis, value)
		elif name == 'x_lim':
			plt.xlim(value)
		elif name == 'y_lim':
			plt.ylim(value)
		elif name == 'plt_title':
			plt.title(r'%s' %(value), size='14')
	
	plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
	plt.xlabel(r'%s' %(x_axis), size=12)

	# CALCULATING THE AVERAGE/SD/SDOM OF THE Y-DATA
	if average != False:
		avg = np.sum(data[t0:])/len(data[t0:])
		SD = stdev(data[t0:])
		SDOM = SD/sqrt(len(data[t0:]))

		plt.axvline(avg, ymin=0.0, ymax=1.0, c='r')
		plt.figtext(0.680, 0.780, '%s\n%6.4f $\\pm$ %6.4f %s \nSD = %4.3f %s' %(analysis, avg, SDOM, x_units, SD, x_units), bbox=dict(boxstyle='square', ec='r', fc='w'), fontsize=12)
	
	
	if norm == True:
		plt.ylabel('Probability Density')
		plt.savefig('%s.%s.prob1d.png' %(system,analysis),dpi=300)
		nf = open('%s.%s.prob1d.dat' %(system,analysis),'w')
	else:
		plt.ylabel('Frequency', size=12)
		plt.savefig('%s.%s.hist1d.png' %(system,analysis),dpi=300)
		nf = open('%s.%s.hist1d.dat' %(system,analysis), 'w')

	for i in range(len(events)):
		nf.write('%10.1f      %10.4f\n' %(events[i], edges[i]))
	
	plt.close()
	nf.close()
	events = []
	edges = []
	patches = []


def scat_hist(xdata, ydata, color, x_axis, y_axis, system, analysis, num_b = 100, average = False, t0 = 0, **kwargs):
	""" Creates 1D scatter plot w/ a 1D histogram

	Usage: scat_hist(xdata, ydata, color, x_axis, y_axis, system, analysis, num_b)
	
	Arguments:
	xdata, ydata: self-explanatory
	color: color to be used to plot data
	x_axis, y_axis: strings to be printed on the axi labels
	system: descriptor for the system analyzed
	analysis: descriptor for the analysis performed and plotted
	num_b: number of bins to be used when binning the data; Default is 100
	average: [False|True]; Default is False; if set to True, the function will calc the average, standard dev, and standard dev of mean of the y-data
	t0: index to begin averaging from; Default is 0
	
	kwargs:
		xunits, yunits: string with correct math text describing the units for the x/y data
		x_lim, y_lim: list w/ two elements, setting the limits of the x/y ranges of plot
		plt_title: string to be added as the plot title

	"""
	# INITIATING THE PLOT SIZES
	left, width = 0.1, 0.65
	bottom, height = 0.1, 0.8
	bottom_h = left_h = left+width+0.01
	rect_scatter = [left, bottom, width, height]
	rect_histy = [left_h, bottom, 0.2, height]
	
	# INITIATING THE PLOT...
	plt.figure(1, figsize=(10,8))
	axScatter =plt.axes(rect_scatter)
	axScatter.plot(xdata, ydata, '%s.' %(color))
	
	# READING IN KWARG DICTIONARY INTO SPECIFIC VARIABLES
	for name, value in kwargs.items():
		if name == 'xunits':
			x_units = value
			x_axis = '%s (%s)' %(x_axis, value)
		elif name == 'yunits':
			y_units = value
			y_axis = '%s (%s)' %(y_axis, value)
		elif name == 'x_lim':
			plt.xlim(value)
		elif name == 'y_lim':
			plt.ylim(value)
		elif name == 'plt_title':
			plt.title(r'%s' %(value), size='14')
	
	plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
#	plt.xlim((0,500))
	plt.ylabel(r'%s' %(y_axis),size=12)
	plt.xlabel(r'%s' %(x_axis),size=12)

	if average != False:
		avg = np.sum(ydata[t0:])/len(ydata[t0:])
		SD = stdev(ydata[t0:])
		SDOM = SD/sqrt(len(ydata[t0:]))
		plt.axhline(avg, xmin=0.0, xmax=1.0, c='r')

	axHisty = plt.axes(rect_histy)
	axHisty.yaxis.set_major_formatter(nullfmt)
	axHisty.xaxis.set_major_formatter(nullfmt)
	plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
	axHisty.hist(ydata, bins=num_b, orientation='horizontal', color = ['gray'])
	axHisty.set_ylim(axScatter.get_ylim())
	
	# CALCULATING THE AVERAGE/SD/SDOM OF THE Y-DATA
	if average != False:
		plt.axhline(avg, xmin=0.0, xmax=1.0, c='r')
		plt.figtext(0.775, 0.810, '%s\n%6.4f $\\pm$ %6.4f %s \nSD = %4.3f %s' %(analysis, avg, SDOM, y_units, SD, y_units), bbox=dict(boxstyle='square', ec='r', fc='w'), fontsize=12)
	
	plt.savefig('%s.%s.scat_hist.png' %(system, analysis),dpi=300)
	plt.close()


def bar(xdata, ydata, x_axis, y_axis, system, analysis, **kwargs): 
	""" Creates a bar graph

	Usage: bar(xdata, ydata, x_axis, y_axis, **kwarg)

	Arguments:
	xdata, ydata: self-explanatory
	x_axis, y_axis: strings to be printed on the axi labels
	system: descriptor for the system analyzed
	analysis: descriptor for the analysis performed and plotted

	kwargs:
	xunits, yunits: string with correct math text describing the units for the x/y data
	x_lim, y_lim: list (or tuple) w/ two elements, setting the limits of the x/y ranges of plot
	plt_title: string to be added as the plot title

	"""
	
	# INITIATING THE PLOT...
	plt.bar(xdata,ydata)

	# READING IN KWARG DICTIONARY INTO SPECIFIC VARIABLES
	for name, value in kwargs.items():
		if name == 'xunits':
			x_units = value
			x_axis = '%s (%s)' %(x_axis, value)
		elif name == 'yunits':
			y_units = value
			y_axis = '%s (%s)' %(y_axis, value)
		elif name == 'x_lim':
			plt.xlim(value)
		elif name == 'y_lim':
			plt.ylim(value)
		elif name == 'plt_title':
			plt.title(r'%s' %(value), size='16')
	
	plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
	plt.ylabel(r'%s' %(y_axis),size=12)
	plt.xlabel(r'%s' %(x_axis),size=12)

	plt.savefig('%s.%s.bar.png' %(system,analysis),dpi=300)
	plt.close()

def hist2d(xdata, ydata, x_axis, y_axis, num_b, system, analysis, norm):
	""" Creates a 2D histogram (heat map)
	
	Usage: hist2d(xdata, ydata, x_axis, y_axis, num_b, system, analysis, norm)
	
	Arguments:
	xdata, ydata: self-explanatory
	x_axis, y_axis: strings to be printed on the axi labels
	num_b: number of bins to be used when binning the data
	system: descriptor for the system analyzed
	analysis: descriptor for the analysis performed and plotted
	norm = [False][True]; if False, plotting a frequency of data; if True, plotting a probability density
	"""

	my_cmap = plt.cm.get_cmap('jet')
	my_cmap.set_under('w')
	counts, xedges, yedges, image = plt.hist2d(xdata, ydata, bins=num_b, normed=norm, cmap=my_cmap, vmin=0.001)#, cmap=plt.get_cmap('jet')) # cmap: jet (blue to red), blues (white to blue), ...
	cb1 = plt.colorbar()
	if norm == True:
		cb1.set_label('Prob. Density', size=12)
	else:
		cb1.set_label('Frequency')
#	plt.title('Distribution of Base Pair interactions - %s-%s' %(base_a, base_b))
#	plt.xlim((0,8))
#	plt.ylim((0,8))
	plt.xlabel(r'%s' %(x_axis), size=12)
	plt.ylabel(r'%s' %(y_axis), size=12)
	plt.savefig('%s.%s.hist2d.png' %(system, analysis),dpi=300)
	plt.close()
	counts = []
	xedges = []
	yedges = []
	image = []

def matrix2d(matrix, x_axis, y_axis, cb_axis, system, analysis, **kwargs):
	""" Creates a 2D matrix image

	Usage: matrix2d(matrix,x_axis,y_axis,system,analysis)

	Arguments:
	matrix: the data matrix to be plotted (should have shape of MxN, but can have MxNx3 or MxNx4)
	x_axis, y_axis: strings to be printed on the axi labels
	system: descriptor for the system analyzed
	analysis: descriptor for the analysis performed and plotted
        
	kwargs:
		vmin, vmax: floats that define the limits for the color bar; if below vmin, data will be colored white; if above vmax, data will be colored red (might want to change this for aesthetics)
		plt_title: string to be added as the plot title
		cb_units: sting to be added to the color bar label to indicate the units of the color bar 
	"""

	vmin =0.001
	vmax = None

#	c = mcolors.ColorConverter().to_rgb
#	bgr = make_colormap([c('blue'),c('lime'),0.50,c('lime'),c('red'),1.00,c('red')])
#	bgr = make_colormap([c('red'),c('lime'),0.50,c('lime'),c('blue'),1.00,c('blue')])
#	bgr.set_under('k')
#	bgr.set_over('r')
#	bgr.set_over('w')
#	my_cmap = bgr

	my_cmap = plt.cm.get_cmap('jet')

#	my_cmap = plt.cm.get_cmap('gray')

	# READING IN KWARG DICTIONARY INTO SPECIFIC VARIABLES
	for name, value in kwargs.items():
		if name == 'vmin':
			vmin = value
		elif name == 'vmax':
			vmax = value
		elif name == 'cb_units':
			cb_units = value
			cb_axis = '%s (%s)' %(cb_axis, value)
		elif name == 'plt_title':
			plt.title(r'%s' %(value), size='14')
	
	plt.imshow(matrix,cmap=my_cmap,vmin=vmin,vmax=vmax,interpolation='none',origin='lower')
	cb1 = plt.colorbar(extend='max',cmap=my_cmap)
	cb1.set_label(r'%s' %(cb_axis), size=12)
	plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
	plt.xlabel(r'%s' %(x_axis), size=12)
	plt.ylabel(r'%s' %(y_axis), size=12)
	plt.savefig('%s.%s.matrix2d.png' %(system, analysis),dpi=300)
	plt.close()

