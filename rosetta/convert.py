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

  x_center, y_center, z_center = adtGrid.massCenter.get()
  npts = (adtGrid.radius.get() * 2) / adtGrid.spacing.get()

  with open(agdfile, "w") as agd:
    agd.write("Title: \n")
    agd.write("Mid:   %s     %s    %s\n" % (round(x_center, 3), round(y_center, 3), round(z_center, 3)))
    agd.write("Dim:     %s     %s     %s\n" % (round(npts, None), round(npts, None), round(npts, None)))
    agd.write("Spacing:     %s\n" % (adtGrid.spacing.get()))
    agd.write("Values:\n")

    with open(e_map, "r") as emap:
      for line in emap.readlines():
        if not line.startswith(("GRID", "MACROMOLECULE", "SPACING", "CENTER", "NELEMENTS")):
          agd.write(line)

  return GridAGD(agdfile)
