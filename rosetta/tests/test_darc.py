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


from pyworkflow.tests import *

from pwem.protocols import ProtImportPdb, ProtSetFilter
from rosetta.protocols import RosettaProteinPreparation, RosettaProtDARC
from pwchem.protocols import ProtChemImportSmallMolecules, ProtChemOBabelPrepareLigands, \
  ProtChemRDKitPrepareLigands, ProtDefineStructROIs

ADT = False

class TestImportBase(BaseTest):
    '''Initial protocols needed to run DARC. Paralelized and checking with plugins are available'''
    @classmethod
    def setUpClass(cls):
      cls.ds = DataSet.getDataSet('model_building_tutorial')
      cls.dsLig = DataSet.getDataSet("smallMolecules")

      setupTestProject(cls)
      cls._runImportPDB()
      cls._runImportSmallMols()
      cls._waitOutput(cls.protImportPDB, 'outputPdb', sleepTime=5)
      cls._waitOutput(cls.protImportSmallMols, 'outputSmallMolecules', sleepTime=5)

      cls._runPrepareLigandsOBabel()
      cls._runPrepareLigandsRDKit()
      cls._runPrepareReceptor()

      cls._waitOutput(cls.protOBabel, 'outputSmallMolecules', sleepTime=5)
      cls._waitOutput(cls.protPrepareLigandRDKit, 'outputSmallMolecules', sleepTime=5)
      cls._waitOutput(cls.protPrepareReceptor, 'outputStructure', sleepTime=5)

      cls.pocketProt = cls._runPocketFinder()
      cls._waitOutput(cls.pocketProt, 'outputStructROIs', sleepTime=5)

    @classmethod
    def _runImportSmallMols(cls):
      cls.protImportSmallMols = cls.newProtocol(
        ProtChemImportSmallMolecules,
        filesPath=cls.dsLig.getFile('mol2'))
      cls.proj.launchProtocol(cls.protImportSmallMols, wait=False)

    @classmethod
    def _runImportPDB(cls):
      protImportPDB = cls.newProtocol(
        ProtImportPdb,
        inputPdbData=1, pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1.pdb'))
      cls.launchProtocol(protImportPDB, wait=False)
      cls.protImportPDB = protImportPDB

    @classmethod
    def _runPrepareLigandsOBabel(cls):
      cls.protOBabel = cls.newProtocol(
        ProtChemOBabelPrepareLigands,
        inputType=0, method_charges=0,
        doConformers=True, method_conf=0, number_conf=2,
        rmsd_cutoff=0.375)
      cls.protOBabel.inputSmallMolecules.set(cls.protImportSmallMols)
      cls.protOBabel.inputSmallMolecules.setExtended('outputSmallMolecules')

      cls.proj.launchProtocol(cls.protOBabel, wait=False)

    @classmethod
    def _runPrepareLigandsRDKit(cls):
      cls.protPrepareLigandRDKit = cls.newProtocol(
        ProtChemRDKitPrepareLigands,
        doConformers=True, numConf=2)
      cls.protPrepareLigandRDKit.inputSmallMolecules.set(cls.protImportSmallMols)
      cls.protPrepareLigandRDKit.inputSmallMolecules.setExtended('outputSmallMolecules')

      cls.proj.launchProtocol(cls.protPrepareLigandRDKit, wait=False)

    @classmethod
    def _runPrepareReceptor(cls):
      cls.protPrepareReceptor = cls.newProtocol(
        RosettaProteinPreparation,
        inputAtomStruct=cls.protImportPDB.outputPdb,
        rchains=True,
        chain_name='{"model": 0, "chain": "C", "residues": 141}')

      cls.proj.launchProtocol(cls.protPrepareReceptor, wait=False)

    @classmethod
    def _runPocketFinder(cls):
        protPocketFinder = cls.newProtocol(
          ProtDefineStructROIs,
          inputAtomStruct=cls.protPrepareReceptor.outputStructure,
          inROIs='1) Residues: {"model": 0, "chain": "C", "index": "99-99", "residues": "I"}')

        cls.proj.launchProtocol(protPocketFinder, wait=False)
        return protPocketFinder

    def _runDARC(self, adtLigs=False, pocketsProt=None):
        if adtLigs:
            protLigs = self.protPrepareLigandRDKit
        else:
            protLigs = self.protOBabel

        if pocketsProt == None:
            protDARC = self.newProtocol(
                RosettaProtDARC,
                fromReceptor=0,
                target_residue='99:C',
                numberOfThreads=8)

            protDARC.inputAtomStruct.set(self.protPrepareReceptor)
            protDARC.inputAtomStruct.setExtended('outputStructure')
            protDARC.inputSmallMolecules.set(protLigs)
            protDARC.inputSmallMolecules.setExtended('outputSmallMolecules')

            protDARC.use_electro.set(False)

            self.launchProtocol(protDARC)
            self.assertIsNotNone(getattr(protDARC, 'outputSmallMolecules', None))

        else:
            protDARC = self.newProtocol(
                RosettaProtDARC,
                fromReceptor=1,
                numberOfThreads=8)

            protDARC.inputStructROIs.set(self.pocketProt)
            protDARC.inputStructROIs.setExtended('outputStructROIs')
            protDARC.inputSmallMolecules.set(protLigs)
            protDARC.inputSmallMolecules.setExtended('outputSmallMolecules')

            protDARC.use_electro.set(False)

            self.launchProtocol(protDARC)
            self.assertIsNotNone(getattr(protDARC, 'outputSmallMolecules', None))

        return protDARC


class TestDARC(TestImportBase):

    def test_1(self):
        """ Complete Docking from whole protein and shape only
        """
        print("\n Complete Docking from whole protein and shape only \n")
        protDARC = self._runDARC()

    def test_2(self):
        """ Complete Docking from protein pockets and shape only
        """
        print("\n Complete Docking from protein pockets and shape only \n")
        protDARC = self._runDARC(adtLigs=ADT, pocketsProt=self.pocketProt)

