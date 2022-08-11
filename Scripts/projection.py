#-----------------------------------------------------------------------------#
# projection.py
#
# NPS Night Skies Program
#
# Last updated: 2021/02/22
#
# This script plots the fits images in fisheye and Hammer projections.
#
# Input:
#   (1) reading in the mask to get the x,y center and the fisheye view radius
#	(2) all the processed fisheye fit images
#
# Output:
#   (1) *fisheye.png
#   (2) *hammer.png
#
# History:
#	Li-Wei Hung -- Created
#
#------------------------------------------------------------------------------#
import copy
import matplotlib as mpl
from matplotlib.transforms import Transform
import numpy as n
import warnings

from astropy.io import fits
from astropy.time import Time, TimeDelta
from glob import glob
from matplotlib import pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from skimage.transform import rotate

# Local Source
import colormaps
import process_input as p
import upper_hammer

#------------------------------------------------------------------------------#


def main():
    """
    This script plots fits images in fisheye and Hammer projection. See the 
    script description for detail.
    """
    #--------------------------------------------------------------------------#
    #						  Generate Polar Coordinates				   	   #
    #--------------------------------------------------------------------------#
    # Mask - read in the fisheye mask center coordinates and radius
    mask = fits.open(p.mask, uint=False)[0].header
    xc, yc, r0 = int(mask['CENTERX']), int(
        mask['CENTERY']), int(mask['RADIUS'])
    X, Y = n.meshgrid(n.arange(-r0, r0), n.arange(-r0, r0))

    # Polar coordinates
    r = n.sqrt(X**2+Y**2) / r0
    theta = -n.arctan2(Y, X)

    # Fisheye takes r in degree
    r_deg = 90 * r
    theta_f = theta + n.pi/2

    # Hammer plot requires the values to be sorted
    r_str = n.pi/2 - r * n.pi/2
    inds = n.argsort(theta[:, 0])
    theta_s = theta[inds, :]
    r_s = r_str[inds, :]

    #--------------------------------------------------------------------------#
    #						  Define Plot settings						   	   #
    #--------------------------------------------------------------------------#
    # General plot settings
    plt.close('all')
    plt.style.use('dark_background')
    plt.rcParams['image.cmap'] = 'NPS_mag'
    cmap = copy.copy(mpl.cm.get_cmap("NPS_mag"))
    cmap.set_bad(color='black')

    colorbar = plt.imread(p.calibration+'colorbar.jpg')  # colorbar
    logo = plt.imread(p.calibration+'AH_logo1.jpg').copy()  # arrowhead
    logo[n.where(logo[:, :] == [8, 6, 7])] = 0  # set background to black

    # Fisheye plot setting
    fig0 = plt.figure('fisheye', figsize=(6, 7))
    ax0 = fig0.add_subplot(111, projection='polar')
    fig0.tight_layout(rect=(0.02, 0.05, 0.98, 0.9))
    ax0.set_rlim(0, 90)
    ax0.set_yticklabels([])
    ax0.tick_params(colors='darkgray')
    ax0.set_theta_zero_location('N')
    imagebox0 = OffsetImage(logo, zoom=0.25)
    imagebox0.image.axes = ax0
    ab0 = AnnotationBbox(imagebox0, (0.06, 0.06), xycoords='figure fraction',
                         frameon=False)
    ax0.add_artist(ab0)
    imagebox1 = OffsetImage(colorbar, zoom=0.2)
    imagebox1.image.axes = ax0
    ab1 = AnnotationBbox(imagebox1, (0.2, 0.892), xycoords='figure fraction',
                         frameon=False)
    ax0.add_artist(ab1)
    for mag in range(14, 25, 2):
        ax0.text(-0.11+0.043*(mag-14), 1.08, mag, color='gray', fontsize=9,
                 ha='left', va='center', transform=ax0.transAxes)
    ax0.text(0.12, 1.15, r'V mags arcsec$^{-2}$', color='gray', fontsize=9,
             ha='center', va='top', transform=ax0.transAxes)
    ax0.text(-0.01, -0.095, 'U.S. National Park Service', color='gray',
             fontsize=9, ha='left', va='center', transform=ax0.transAxes)
    ax0.text(-0.01, -0.13, 'Night Skies Program', color='gray',
             fontsize=9, ha='left', va='center', transform=ax0.transAxes)

    # Hammer plot setting
    fig1 = plt.figure('hammer', figsize=(15, 5.6))
    ax1 = fig1.add_subplot(111, projection="upper_hammer")
    fig1.tight_layout(rect=(0.03, -0.55, 0.98, 0.97))
    imagebox2 = OffsetImage(logo, zoom=0.26)
    imagebox2.image.axes = ax1
    ab2 = AnnotationBbox(imagebox2, (0.05, 0.08), xycoords='figure fraction',
                         frameon=False)
    ax1.add_artist(ab2)
    imagebox3 = OffsetImage(colorbar, zoom=0.47)
    imagebox3.image.axes = ax1
    ab3 = AnnotationBbox(imagebox3, (0.19, 0.92), xycoords='figure fraction',
                         frameon=False)
    ax1.add_artist(ab3)
    for mag in range(14, 25):
        ax1.text(-0.017+0.0348*(mag-14), 1.06, mag, color='gray', fontsize=10,
                 ha='left', va='top', transform=ax1.transAxes)
    ax1.text(0.14, 1.04, r'V mags arcsec$^{-2}$', color='gray', fontsize=10,
             ha='center', va='top', transform=ax1.transAxes)
    ax1.text(0.032, 0.41, 'U.S. National Park Service', color='gray',
             fontsize=10, ha='left', va='center', transform=ax1.transAxes)
    ax1.text(0.043, 0.385, 'Night Skies Program', color='gray',
             fontsize=10, ha='left', va='center', transform=ax1.transAxes)

    # Suppressing a MatPlotLib benign warning about pcolormesh shading
    warnings.filterwarnings("ignore", category=UserWarning)

    #--------------------------------------------------------------------------#
    #				Plot the image in fisheye and Hammer projections		   #
    #--------------------------------------------------------------------------#
    for f in glob(p.data_cal+'img-00*-sky*.fit'):

        print('projecting ' + f[len(p.data_cal):])
        imgf = fits.open(f, uint=False)[0]
        hdr = imgf.header
        img = imgf.data[yc-r0:yc+r0, xc-r0:xc+r0]
        img_hammer = rotate(img.astype('float32'), -90, cval=n.nan)[inds, :]
        t = Time(hdr['DATE-OBS'])+TimeDelta(p.UTCoffset*3600, format='sec')
        date = str(t.datetime.date())
        hour = str(round(t.datetime.hour+t.datetime.minute/60, 1))+' hours'

        # plot fisheye
        ax0.pcolormesh(theta_f, r_deg, img, shading='flat', vmin=14, vmax=24)
        ax0.grid(True, color='gray', linestyle='dotted', linewidth=.5)
        ax0.text(0.5, 1.2, hdr['LOCATION'], color='gray', fontsize=14,
                 ha='center', va='center', transform=ax0.transAxes)
        dtext = ax0.text(1.08, 1.12, date, color='gray', fontsize=12,
                         ha='right', va='center', transform=ax0.transAxes)
        htext = ax0.text(1.08, 1.08, hour, color='gray', fontsize=12,
                         ha='right', va='center', transform=ax0.transAxes)
        ax0.text(1.08, -0.095, 'Observer: '+hdr['OBSERVER'], color='gray', fontsize=9,
                 ha='right', va='center', transform=ax0.transAxes)
        ax0.text(1.08, -0.13, 'Processer: '+hdr['PROCESS'], color='gray', fontsize=9,
                 ha='right', va='center', transform=ax0.transAxes)
        fig0.savefig(f[:-4]+'_fisheye.png', dpi=200)
        # plt.show(block=False)
        dtext.remove()
        htext.remove()

        # plot hammer
        ax1.pcolormesh(theta_s, r_s, img_hammer, vmin=14, vmax=24)
        ax1.grid(True)
        ax1.text(0.5, 1.07, hdr['LOCATION'], color='gray', fontsize=16,
                 ha='center', va='center', transform=ax1.transAxes)
        obstime = date + '  ' + hour
        ta = ax1.text(1.0, 1.07, obstime, color='gray', fontsize=16,
                      ha='right', va='center', transform=ax1.transAxes)
        ax1.text(1, 0.41, 'Observer: '+hdr['OBSERVER'], color='gray', fontsize=10,
                 ha='right', va='center', transform=ax1.transAxes)
        ax1.text(1, 0.385, 'Processer: '+hdr['PROCESS'], color='gray', fontsize=10,
                 ha='right', va='center', transform=ax1.transAxes)
        fig1.savefig(f[:-4]+'_hammer.png')
        ta.remove()


if __name__ == '__main__':
    main()
