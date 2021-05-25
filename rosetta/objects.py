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


import pyworkflow.object as pwobj
import pwem.objects.data as data
from pwem.objects.data import AtomStruct

class RaysProtein(data.EMFile):
    """ Represent a RAY file """
    def __init__(self, filename=None, **kwargs):
        data.EMFile.__init__(self, filename, **kwargs)


class RaysStruct(data.AtomStruct):
    """ Represent a RAY pdb file """
    def __init__(self, filename=None, pseudoatoms=False, **kwargs):
        data.AtomStruct.__init__(self, filename, pseudoatoms, **kwargs)


class GridAGD(data.EMFile):
    """ Represent a grid file in agd (ASCIII) format """
    def __init__(self, filename=None, **kwargs):
        data.EMFile.__init__(self, filename, **kwargs)


class DarcScore(data.EMObject):
    """Score given by DARC to each small molecule file (and its conformers) """
    def __init__(self, zincID, scoreDarc, **kwargs):
        data.EMObject.__init__(self)
        self.zincID = pwobj.String(zincID)
        self.scoreDarc = pwobj.Float(scoreDarc)
        # Variable Columns depending on if it is used a electrostatic grid or not
        #  Total_Energy
        #  Interface_Energy
        #  Interface_Hb
        #  Total_Pack
        #  Interface_Unsat
        #  Thetalig
        # self.smallMoleculeFile = pwobj.String(kwargs.get('smallMolFilename', None))

    def getzincID(self):
        return self.zincID.get()

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
