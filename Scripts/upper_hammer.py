#-----------------------------------------------------------------------------#
#upper_hammer.py
#
#NPS Night Skies Program
#
#Last updated: 2020/10/21
#
#This script 
#
#Input: 
#   (1) 
#
#Output:
#   (1) 
#
#History:
#	Li-Wei Hung -- Created 
#
#-----------------------------------------------------------------------------#

import matplotlib.projections as mprojections
import matplotlib.spines as mspines
import numpy as n

from matplotlib.axes import Axes
from matplotlib.patches import Wedge
from matplotlib.projections.geo import HammerAxes
from matplotlib.ticker import FixedLocator
from matplotlib.transforms import Affine2D


class UpperHammerAxes(HammerAxes):
	name = 'upper_hammer'
	def cla(self):
		HammerAxes.cla(self)
		Axes.set_xlim(self, -n.pi, 2*n.pi)
		Axes.set_ylim(self, 0, n.pi / 2.0)
		self.xaxis.set_ticks_position('top')
		self.tick_params(axis='x', length=0)
		self.tick_params(colors='grey', labelsize=11)
		self.grid(color='gray', linestyle='dotted', linewidth=.5)
		self.set_yticks(n.linspace(0, n.pi/2, 7)) #include latitude 90 label
		self.set_xticklabels(['N','30°','60°','E','120°','150°','S','210°','240°','W','300°','330°','N'])

	
	def _gen_axes_patch(self):
		return Wedge((0.5, 0.5), 0.5, 0, 180)
	
	def _gen_axes_spines(self):
		pass
		path = Wedge((0, 0), 1.0, 0, 180).get_path()
		spine = mspines.Spine(self, 'circle', path)
		spine.set_color('gray')
		spine.set_patch_circle((0.5, 0.5), 0.5)
		return {'wedge':spine}

	def _set_lim_and_transforms(self):
		super()._set_lim_and_transforms()
		
		self._xaxis_text2_transform = \
			Affine2D().scale(1.0, 0.0) + \
			self.transData + \
			Affine2D().translate(0.0, -12.0)
		
		yaxis_stretch = Affine2D().scale(n.pi*2, 1.02).translate(-n.pi, 0)
		yaxis_space = Affine2D().scale(1.0, 1.05)
		
		yaxis_text_base = \
			yaxis_stretch + \
			self.transProjection + \
			(yaxis_space +
			self.transAffine +
			self.transAxes)
		
		self._yaxis_text1_transform = \
			yaxis_text_base + \
			Affine2D().translate(-6, -4)
		
	def get_yaxis_text1_transform(self, pad):
		return self._yaxis_text1_transform, 'center', 'right'
		
	def set_longitude_grid(self, degrees):
		"""
		Set the number of degrees between each longitude grid.
	
		This is an example method that is specific to this projection
		class -- it provides a more convenient interface to set the
		ticking than set_xticks would.
		"""
		grid = n.linspace(-180, 180, int(360/degrees+1))
		self.xaxis.set_major_locator(FixedLocator(n.deg2rad(grid)))
		self.xaxis.set_major_formatter(self.ThetaFormatter(degrees))
	
		
mprojections.register_projection(UpperHammerAxes)

if __name__ == '__main__':
	import matplotlib.pyplot as plt
	#plt.close('all')
	fig = plt.figure(figsize=(15,5))
	ax = fig.add_subplot(111, projection='upper_hammer')
	ax.grid(True)
	plt.tight_layout(rect=(0.03,-0.6,0.98,0.97))
	plt.show()