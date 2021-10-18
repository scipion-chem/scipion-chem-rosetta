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


from pathlib import Path

from pyworkflow.tests import *
from pyworkflow.protocol import *

from pwem.protocols import ProtImportPdb, ProtSetFilter
from rosetta.protocols import RosettaProteinPreparation, RosettaProtDARC
from pwchem.protocols import ProtChemImportSmallMolecules, ProtChemOBabelPrepareLigands



class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
      cls.ds = DataSet.getDataSet('model_building_tutorial')
      cls.dsLig = DataSet.getDataSet("smallMolecules")

      setupTestProject(cls)
      cls._runImportPDB()
      cls._runImportSmallMols()

      cls._runPrepareLigandsOBabel()
      cls._runPrepareLigandsADT()
      cls._runPrepareReceptorADT()

    @classmethod
    def _runImportSmallMols(cls):
      cls.protImportSmallMols = cls.newProtocol(
        ProtChemImportSmallMolecules,
        filesPath=cls.dsLig.getFile('mol2'))
      cls.launchProtocol(cls.protImportSmallMols)

    @classmethod
    def _runImportPDB(cls):
      cls.protImportPDB = cls.newProtocol(
        ProtImportPdb,
        inputPdbData=1,
        pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1_noHETATM.pdb'))
      cls.launchProtocol(cls.protImportPDB)

    @classmethod
    def _runPrepareLigandsOBabel(cls):
      cls.protOBabel = cls.newProtocol(
        ProtChemOBabelPrepareLigands,
        inputType=0, method_charges=0,
        inputSmallMols=cls.protImportSmallMols.outputSmallMolecules,
        doConformers=True, method_conf=0, number_conf=2,
        rmsd_cutoff=0.375)

      cls.launchProtocol(cls.protOBabel)

    @classmethod
    def _runPrepareReceptorADT(cls):
      cls.protPrepareReceptor = cls.newProtocol(
        RosettaProteinPreparation,
        inputStructure=cls.protImportPDB.outputPdb,
        repair=3)

      cls.launchProtocol(cls.protPrepareReceptor)

    def _runSetFilter(self, inProt, number, property):
      protFilter = self.newProtocol(
        ProtSetFilter,
        operation=ProtSetFilter.CHOICE_RANKED,
        threshold=number, rankingField=property)
      protFilter.inputSet.set(inProt)
      protFilter.inputSet.setExtended('outputPockets')

      self.launchProtocol(protFilter)
      return protFilter


    def _runDARC(self, pocketsProt=None):
        if pocketsProt == None:
            protAutoDock = self.newProtocol(
                RosettaProtDARC,
                wholeProt=True,
                inputAtomStruct=self.protPrepareReceptor.outputStructure,
                inputLibrary=self.protOBabel.outputSmallMolecules,
                radius=37, gaRun=2,
                numberOfThreads=8)
            self.launchProtocol(protAutoDock)
            smOut = getattr(protAutoDock, 'outputSmallMolecules', None)
            self.assertIsNotNone(smOut)

        else:
            protAutoDock = self.newProtocol(
                RosettaProtDARC,
                wholeProt=False,
                inputPockets=pocketsProt.outputPockets,
                inputLibrary=self.protPrepareLigand.outputSmallMolecules,
                pocketRadiusN=2, gaRun=2,
                numberOfThreads=8)
            self.launchProtocol(protAutoDock)
            smOut = getattr(protAutoDock, 'outputSmallMolecules_1', None)
            self.assertIsNotNone(smOut)

        return protAutoDock


class TestDARC(TestImportBase):

    def test_1(self):
        """ Complete Docking with GRID and search conformers on the fly
        """
        print("\n Complete Docking with GRID and search conformers on the fly \n")
        protDARC = self._runDARC()
        self.launchProtocol(protDARC)
