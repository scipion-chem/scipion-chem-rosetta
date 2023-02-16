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

import math
try:
    from autodock.objects import GridADT
except:
    print('Autodock plugin cannot be imported, so ADT grid cannot be calculated')

class GridAGD(GridADT):
    """ Represent a grid file in agd (ASCIII) format """
    def __init__(self, filename=None, **kwargs):
        super().__init__(filename, **kwargs)
        if filename != None:
            self.parseFile()

    def parseFile(self):
        with open(self.getFileName()) as f:
            for line in f:
                if line.startswith('Mid:'):
                    self.setMassCenter(list(map(float, line.split()[1:])))
                elif line.startswith('Dim:'):
                    npts=float(line.split()[1])
                    self.setNumberOfPoints(npts)
                elif line.startswith('Spacing:'):
                    self.setSpacing(float(line.split()[1]))
        self.setRadius(math.sqrt(npts * self.getSpacing()))
