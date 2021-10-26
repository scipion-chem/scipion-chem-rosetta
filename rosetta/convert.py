# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Alberto Manuel Parra Pérez (amparraperez@gmail.com)
# *
# * Biocomputing Unit, CNB-CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import os
from .objects import GridAGD

def adt2agdGrid(adtGrid, agdfile=None, outDir=None):
  e_map = adtGrid.getFileName()
  if agdfile == None:
    if outDir != None:
      agdName = e_map.split('/')[-1].replace('.e.map', '.agd')
      agdfile = os.path.join(outDir, agdName)
    else:
      agdfile = e_map.replace('.e.map', '.agd')

  x_center, y_center, z_center = adtGrid.getMassCenter()
  npts = (adtGrid.getRadius() * 2) / adtGrid.getSpacing()

  with open(agdfile, "w") as agd:
    #https://docs.eyesopen.com/toolkits/python/oechemtk/grids.html
    agd.write("Title:\n")
    agd.write("Mid: %12.6f %12.6f %12.6f\n" % (x_center, y_center, z_center))
    agd.write("Dim: %6d %6d %6d\n" % (npts, npts, npts))
    agd.write("Spacing: %12.6f\n" % adtGrid.getSpacing())
    agd.write("Values:\n")

    with open(e_map, "r") as emap:
      for line in emap.readlines():
        if not line.startswith(("GRID", "MACROMOLECULE", "SPACING", "CENTER", "NELEMENTS")):
          agd.write("%-12.6e\n" % float(line.strip()))

  return GridAGD(agdfile)
