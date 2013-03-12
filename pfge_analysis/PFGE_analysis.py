#!/usr/bin/env python2.7
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import spline
# from http://www.lfd.uci.edu/~gohlke/
try:
  import tifffile
except ImportError:
  print '''
The Tifffile module by Christoph Gohlke is required 
for gel image analysis. Download from 
http://www.lfd.uci.edu/~gohlke/, save in src/ 
(same folder as this script) and try again 
(but save the rendered contents of the html 
file tifffile.py.html as tifffile.py).

Warnings about additional modules can be 
ignored after Tifffile installation.
'''
  import sys
  sys.exit()
def getLanes(imagefile,lanes,top,bottom):
    '''Load pixel data summed by pixel rows down length of specified region corresponding to a gel lane'''
    tif = tifffile.TiffFile(imagefile)
    image = tif.asarray()
    lanes1d = []
    for l in lanes:
        lanes1d += [np.sum(image[top:bottom,l[0]:l[1]],axis=1)]
    return(lanes1d)

def getPerLaneLadder(ladder1x,ladder1peaks,ladder2x,ladder2peaks,samplelanes,sampleoffsets):
    '''Linearly interpolate between corresponding bands in ladders on either side of samples for "per lane ladders"'''
    ldiff = ladder2x - ladder1x
    # take min peaks counting from small to big
    nladderpeaks = min(len(ladder1peaks),len(ladder2peaks))
    ladder1peaks = ladder1peaks[-nladderpeaks:]
    ladder2peaks = ladder2peaks[-nladderpeaks:]
    ladderbylanes = []
    for n1,l in enumerate(samplelanes):
        thesepeaks = []
        soffset = sampleoffsets[n1] - ladder1x
        for p in range(nladderpeaks):
            peakdiff = ladder1peaks[p]-ladder2peaks[p]
            thesepeaks += [ladder1peaks[p]-peakdiff*(soffset/float(ldiff))]
        ladderbylanes += [np.array(thesepeaks)]
    return(ladderbylanes)

def Quantify(ladderbylanes,laddersizes,allpeaks,samplenames,doplot,outfolder=None,prefix=None,desired_range=False,gelimagefile=None):
  '''Interpolate ladder sizes versus migration distance with second order splines and quantify samples'''
  sizes = {}
  for n,ladder_peak_dists in enumerate(ladderbylanes):
    sizes[samplenames[n]] = []
    # interpolate and smooth between ladder peaks (which may not be linear or log-linear etc)
    interpolated_lsizes = np.linspace(laddersizes.min(),laddersizes.max(),1000)
    smoothed_lpeaks_dists = spline(laddersizes[::-1],ladder_peak_dists[::-1],interpolated_lsizes[::-1],order=2)[::-1]
    if doplot:
      fig = plt.figure(figsize=(12, 9))
      axes = fig.add_subplot(111)
      axes.plot(interpolated_lsizes,smoothed_lpeaks_dists,'-r',lw=0.7)
      axes.plot(laddersizes,ladder_peak_dists,'og')
      axes.set_xlabel('Fragment size (bp)')
      axes.set_ylabel('Band migration (pixels)')
      if desired_range:
        drange = 'Desired range: %s bp, s'%desired_range
      else:
        drange = 'S'
      #axes.set_title('Band fragment sizes by interpolation of ladder migration distances\n%sample lane %s, %s' % (drange,n+1,samplenames[n]))
      plt.title('Band fragment sizes by interpolation of ladder migration distances\n%sample lane %s, %s, %s' % (drange,n+1,samplenames[n],gelimagefile))
      for b,band in enumerate(ladder_peak_dists):
        axes.annotate('%s bp' % laddersizes[b], xy=(laddersizes[b],band), xytext=(laddersizes[b]+5000,band), arrowprops=None)
    # find where sample peaks intersect ladder spline and get corresponding DNA fragment size
    xmin,xmax = axes.get_xlim()
    ymin,ymax = axes.get_ylim()
    for peak in allpeaks[n+1]:
      s = [i for i,a in enumerate(smoothed_lpeaks_dists[1:]) if a < peak and smoothed_lpeaks_dists[i-1] >= peak][0]
      sizes[samplenames[n]] += [int(round(interpolated_lsizes[s]))]
      if doplot:
        axes.axvline(x=interpolated_lsizes[s],ymax=(peak-ymin)/(ymax-ymin),color = 'purple',lw=1)
        axes.axhline(y=peak,xmax=(interpolated_lsizes[s]-xmin)/(xmax-xmin),color = 'purple',lw=1)
    if doplot:
      axes.annotate('%s' % '\n'.join(['Quantification:']+[str(s)+' bp' for s in sorted(sizes[samplenames[n]])]), xy=((xmax-xmin)*0.6+xmin,(ymax-ymin)*0.8+ymin), xytext=((xmax-xmin)*0.6+xmin,(ymax-ymin)*0.8+ymin), arrowprops=None)
      plt.savefig(outfolder+os.sep+prefix+'quantification_sample_lane_%s_%s%ssvg' % (n+1,samplenames[n],os.extsep), format='svg')
  return(sizes)

def AnalyseGel(gelimagefile,samplenames,lane_x_ranges,top,bottom,allpeaks,laddersizes,outfolder,prefix,doplot,desired_range):
  # get lane data from image file
  lanes1d = getLanes(gelimagefile,lane_x_ranges,top,bottom)
  if not os.path.exists(outfolder): os.mkdir(outfolder)
  if doplot:
    fig = plt.figure(figsize=(12, 9))
    # plot ladder peak intensities
    axes = fig.add_subplot(111)
    axes.plot(lanes1d[0])
    axes.set_xlabel('Distance along lane (pixels)')
    axes.set_ylabel('Pixel row intensity')
    plt.title('Lane intensities summed across pixel rows versus distance along ladder lane 1\n%s' % gelimagefile)
    [axes.axvline(x=a,color = 'red') for a in allpeaks[0]]
    plt.savefig(outfolder+os.sep+prefix+os.extsep.join(['ladder1','svg']), format='svg')
    axes.clear()
    axes.plot(lanes1d[-1])
    axes.set_xlabel('Distance along lane (pixels)')
    axes.set_ylabel('Pixel row intensity')
    plt.title('Lane intensities summed across pixel rows versus distance along ladder lane 2\n%s' % gelimagefile)
    [axes.axvline(x=a,color = 'red') for a in allpeaks[-1]]
    plt.savefig(outfolder+os.sep+prefix+os.extsep.join(['ladder2','svg']), format='svg')
    # plot samples 
    for n,lane in enumerate(lanes1d[1:-1]):
      plt.cla()
      axes.plot(lane)
      axes.set_xlabel('Distance along lane (pixels)')
      axes.set_ylabel('Pixel row intensity')
      plt.title('Lane intensities summed across pixel rows versus distance along\nsample lane %s, %s, %s' % (n+1,samplenames[n],gelimagefile))
      [axes.axvline(x=a,color = 'red') for a in allpeaks[n+1]]
      plt.savefig(outfolder+os.sep+prefix+'intensities_sample_lane_%s_%s.svg' % ((n+1),samplenames[n]), format='svg')
  # linear regress between ladder bands to create per lane ladders
  ladder1x = sum(lane_x_ranges[0])/2
  ladder2x = sum(lane_x_ranges[-1])/2
  ladder1peaks = allpeaks[0]
  ladder2peaks = allpeaks[-1]
  samplepeaks = [np.array(a) for a in allpeaks[1:-1]]
  sampleoffsets = [sum(n)/2 for n in lane_x_ranges[1:-1]]
  ladderbylanes = getPerLaneLadder(ladder1x,ladder1peaks,ladder2x,ladder2peaks,samplepeaks,sampleoffsets)
  # interpolate with a smoothed spline
  sizes = Quantify(ladderbylanes,laddersizes,allpeaks,samplenames,doplot,outfolder,prefix,desired_range,gelimagefile)
  # print sizes
  for sample,thesesizes in sizes.items():
    print 'Sample %s: %s' % (sample,', '.join([str(s) for s in sorted(thesesizes)]))
  
  return(sizes)


# low range ladder
# https://www.neb.com/products/n0350-low-range-pfg-marker
lowrangeladder = [23100,9420,6550,4360,2320,2030]
# includes lambda ladder

# lambda ladder:
# http://www.neb.com/nebecomm/products/productn0340.asp
lambdaladder = [1018500,970000,921500,873000,824500,776000,727500, 679000, 630500, 582000, 533500, 485000, 436500, 388000, 339500, 291000, 242500, 194000, 145500, 97000, 48500]

# H. wingei chromosomes ladder:
# http://www.bio-rad.com/prd/en/US/adirect/biorad?cmd=catProductDetail&vertical=LSR&country=US&lang=en&productID=170-3667
Hwingeiladder = [3130000,2700000,2350000,1810000,1660000,1370000,1050000]

doplot = True

samples = ['E_1_37','A_3_34','C_4_22','D_4_27','B_4_28','K12']
outfolder = os.sep.join([os.pardir,'image','analysis'])

sizes = {}
def main():

  ## Small fragment ICeuI PFGE gel
  desired_range = '< 145,500'
  gelimagefile = os.sep.join([os.pardir,'image','gels',os.extsep.join(['ICeuI_small','tif'])])
  laddersizes = np.array(lambdaladder[-4:]+lowrangeladder[:1])
  ## samples in left lanes
  samplenames = samples[:3]
  # x pixel ranges of lanes
  lane_x_ranges = ((176,274),(334,434),(494,584),(642,726),(798,904))
  # upper and lower pixels
  (top,bottom) = (10,1128)
  # band peak intensity distance down gel image -top
  allpeaks = [np.array([  25,  237,  549,  856, 1030]),
   np.array([286, 415, 930]),
   np.array([280, 428, 912]),
   np.array([399, 512, 909]),
   np.array([  33,  243,  554,  862, 1032])]

  prefix = 'ICeuI_small_123_'
  print('\nAnalysing %s, plotting to %s%s%s*.svg\n' % (gelimagefile,outfolder,os.sep,prefix))
  sizes[prefix] =  AnalyseGel(gelimagefile,samplenames,lane_x_ranges,top,bottom,allpeaks,laddersizes,outfolder,prefix,doplot,desired_range)

  ## samples in right lanes
  samplenames = samples[3:]
  # x pixel ranges of lanes
  lane_x_ranges = ((798,904),(962,1082),(1114,1224),(1264,1368),(1422,1530))
  # band peak intensity distance down gel image -top
  allpeaks = [np.array([  33,  243,  554,  862, 1032]),
   np.array([ 76, 408, 914]),
   np.array([304, 575, 913]),
   np.array([331, 588, 927]),
   np.array([  36,  255,  568,  879, 1053])]

  prefix = 'ICeuI_small_456_'
  print('\nAnalysing %s, plotting to %s%s%s*.svg\n' % (gelimagefile,outfolder,os.sep,prefix))
  sizes[prefix] = AnalyseGel(gelimagefile,samplenames,lane_x_ranges,top,bottom,allpeaks,laddersizes,outfolder,prefix,doplot,desired_range)

  ## Mid fragment ICeuI PFGE gel
  desired_range = '> 145,500 & < 1,000,000'
  gelimagefile = os.sep.join([os.pardir,'image','gels',os.extsep.join(['ICeuI_medium','tif'])])
  laddersizes = np.array(lambdaladder[-18:-9])
  samplenames = samples
  # x pixel ranges of lanes
  lane_x_ranges = ((23,63),(96,148),(173,230),(255,305),(333,389),(415,466),(493,546),(584,622))
  # upper and lower pixels
  (top,bottom) = (0,260)
  # band peak intensity distance down gel image -top
  allpeaks = [np.array([  8,  31,  55,  79, 106, 135, 167, 202, 240]),
   np.array([ 38,  90, 193]),
   np.array([105, 107, 210]),
   np.array([ 79, 108, 207]),
   np.array([ 32,  92, 202]),
   np.array([ 88, 131, 202]),
   np.array([ 99, 121, 212]),
   np.array([ 17,  39,  62,  87, 113, 143, 175, 209, 247])]

  prefix = 'ICeuI_medium_123456_'
  print('\nAnalysing %s, plotting to %s%s%s*.svg\n' % (gelimagefile,outfolder,os.sep,prefix))
  sizes[prefix] = AnalyseGel(gelimagefile,samplenames,lane_x_ranges,top,bottom,allpeaks,laddersizes,outfolder,prefix,doplot,desired_range)

  ## Large fragment ICeuI PFGE gel
  desired_range = '> 2,000,000'
  gelimagefile = os.sep.join([os.pardir,'image','gels',os.extsep.join(['ICeuI_large','tif'])])
  laddersizes = np.array(Hwingeiladder[:3])
  # upper and lower pixels
  (top,bottom) = (130,482)

  ## sample in left lane
  samplenames = samples[:1]
  lane_x_ranges = ((162,204),(316,386),(472,514))
  # band peak intensity distance down gel image -top
  allpeaks = [np.array([ 84, 165, 263]), np.array([119]), np.array([ 87, 158, 254])]
  prefix = 'ICeuI_large_1_'
  print('\nAnalysing %s, plotting to %s%s%s*.svg\n' % (gelimagefile,outfolder,os.sep,prefix))
  sizes[prefix] = AnalyseGel(gelimagefile,samplenames,lane_x_ranges,top,bottom,allpeaks,laddersizes,outfolder,prefix,doplot,desired_range)

  ## samples in middle lane
  samplenames = samples[1:3]
  lane_x_ranges = ((472,514),(598,660),(728,802),(878,922))
  # band peak intensity distance down gel image -top
  allpeaks = [np.array([ 87, 158, 248]), np.array([163]), np.array([116]), np.array([ 90, 163, 251])]
  prefix = 'ICeuI_large_23_'
  print('\nAnalysing %s, plotting to %s%s%s*.svg\n' % (gelimagefile,outfolder,os.sep,prefix))
  sizes[prefix] = AnalyseGel(gelimagefile,samplenames,lane_x_ranges,top,bottom,allpeaks,laddersizes,outfolder,prefix,doplot,desired_range)

  ## samples in right lanes
  samplenames = samples[3:]
  lane_x_ranges = ((878,922),(1020,1076),(1160,1220),(1292,1362),(1436,1508))
  # band peak intensity distance down gel image -top
  allpeaks = [np.array([ 90, 163, 251]), np.array([107]), np.array([130]), np.array([145]), np.array([83, 166, 256])]
  prefix = 'ICeuI_large_456_'
  print('\nAnalysing %s, plotting to %s%s%s*.svg\n' % (gelimagefile,outfolder,os.sep,prefix))
  sizes[prefix] = AnalyseGel(gelimagefile,samplenames,lane_x_ranges,top,bottom,allpeaks,laddersizes,outfolder,prefix,doplot,desired_range)
  allsizes = {}

  for sz in sizes.values():
    for sm in samples:
      if sm in sz:
        if sm in allsizes:
          allsizes[sm] += sz[sm]
        else:
          allsizes[sm] = sz[sm]

  print('\n')

  for sm,szs in sorted(allsizes.items()):
    print('%s: %s' % (sm,sum(szs)))


  print('\nSee ../image/analysis/ for analysis plots\n')

if __name__ == '__main__':
  main()
