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
import pyworkflow.object as pwobj
import pwem.objects.data as data
try:
    from autodock.objects import GridADT
except:
    print('Autodock plugin cannot be imported, so ADT grid cannot be calculated')

class RosettaRaysProtein(data.EMFile):
    """ Represent a RAY file """
    def __init__(self, filename=None, **kwargs):
        data.EMFile.__init__(self, filename, **kwargs)


class RosettaRaysStruct(data.AtomStruct):
    """ Represent a RAY pdb file """
    def __init__(self, filename=None, pseudoatoms=False, **kwargs):
        data.AtomStruct.__init__(self, filename, pseudoatoms, **kwargs)


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


class DarcScore(data.EMObject):
    """Score given by DARC to each small molecule file (and its conformers) """
    def __init__(self, **kwargs):
        data.EMObject.__init__(self, **kwargs)
        self.ID = pwobj.String(kwargs.get('ID', None))
        self.scoreDarc = pwobj.Float(kwargs.get('scoreDarc', None))

        # Variable Columns depending on if it is used a electrostatic grid or not
        #  Total_Energy
        #  Interface_Energy
        #  Interface_Hb
        #  Total_Pack
        #  Interface_Unsat
        #  Thetalig

    def getID(self):
        return self.ID.get()

    def getScoreDarc(self):
        return self.scoreDarc.get()

    def getTotal_Energy(self):
        return self.Total_Energy.get()

    def getInterface_Energy(self):
        return self.Interface_Energy.get()

    def getInterface_Hb(self):
        return self.Interface_Hb.get()

    def getTotal_Pack(self):
        return self.Total_Pack.get()

    def getInterface_Unsat(self):
        return self.Interface_Unsat.get()

    def getThetalig(self):
        return self.Thetalig.get()


class SetScores(data.EMSet):
    """ Set of docking score to each small molecule"""
    ITEM_TYPE = DarcScore
    FILE_TEMPLATE_NAME = 'DarcScores%s.sqlite'

    def __init__(self, **kwargs):
        data.EMSet.__init__(self, **kwargs)
