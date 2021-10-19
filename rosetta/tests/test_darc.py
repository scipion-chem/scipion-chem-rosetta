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
from pwchem.protocols import ProtChemImportSmallMolecules, ProtChemOBabelPrepareLigands

pocketFinder=None
try:
    from autodock.protocols import ProtChemADTPrepareLigands, Autodock_GridGeneration, ProtChemAutoLigand
    ADT = True
    try:
        from p2rank.protocols import P2RankFindPockets
        pocketFinder = 'p2rank'
    except:
        pocketFinder='autoligand'
except:
    print('Autodock plugin cannot be imported, so ADT grid cannot be calculated')
    ADT = False
    try:
        from fpocket.protocols import FpocketFindPockets
        pocketFinder='fpocket'
    except:
        print('Cannot import any pocket finder (autoligand, p2rank, fpocket).'
              'Test on pocket will not be performed')

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
      cls._runPrepareReceptor()
      cls._waitOutput(cls.protPrepareReceptor, 'outputSmallMolecules', sleepTime=5)
      cls._waitOutput(cls.protPrepareReceptor, 'outputStructure', sleepTime=5)

      if ADT:
          cls._runPrepareLigandsADT()
          cls._runGridGeneration()

      if pocketFinder != None:
          if pocketFinder == 'autoligand':
              cls._waitOutput(cls.protGridADT, 'outputGrid', sleepTime=5)
          pocketProt = cls._runPocketFinder()
          cls._waitOutput(pocketProt, 'outputPockets', sleepTime=5)
          cls._runSetFilter(inProt=pocketProt, number=2, property='_score')

    @classmethod
    def _runImportSmallMols(cls):
      cls.protImportSmallMols = cls.newProtocol(
        ProtChemImportSmallMolecules,
        filesPath=cls.dsLig.getFile('mol2'))
      cls.proj.launchProtocol(cls.protImportSmallMols, wait=False)

    @classmethod
    def _runImportPDB(cls):
      cls.protImportPDB = cls.newProtocol(
        ProtImportPdb,
        inputPdbData=0,
        pdbId='4erf')
      cls.proj.launchProtocol(cls.protImportPDB, wait=False)

    @classmethod
    def _runPrepareLigandsOBabel(cls):
      cls.protOBabel = cls.newProtocol(
        ProtChemOBabelPrepareLigands,
        inputType=0, method_charges=0,
        doConformers=True, method_conf=0, number_conf=2,
        rmsd_cutoff=0.375)
      cls.protOBabel.inputSmallMols.set(cls.protImportSmallMols)
      cls.protOBabel.inputSmallMols.setExtended('outputSmallMolecules')

      cls.proj.launchProtocol(cls.protOBabel, wait=False)

    @classmethod
    def _runPrepareLigandsADT(cls):
      cls.protPrepareLigandADT = cls.newProtocol(
        ProtChemADTPrepareLigands,
        doConformers=True, method_conf=0, number_conf=2,
        rmsd_cutoff=0.375)
      cls.protPrepareLigandADT.inputSmallMols.set(cls.protImportSmallMols)
      cls.protPrepareLigandADT.inputSmallMols.setExtended('outputSmallMolecules')

      cls.proj.launchProtocol(cls.protPrepareLigandADT, wait=False)

    @classmethod
    def _runPrepareReceptor(cls):
      cls.protPrepareReceptor = cls.newProtocol(
        RosettaProteinPreparation,
        inputAtomStruct=cls.protImportPDB.outputPdb,
        rchains=True,
        chain_name='{"Chain": "C", "Number of residues": 93, "Number of chains": 3}')

      cls.proj.launchProtocol(cls.protPrepareReceptor, wait=False)

    @classmethod
    def _runPocketFinder(cls):
        if pocketFinder == 'autoligand':
          protPocketFinder = cls.newProtocol(
            ProtChemAutoLigand,
            inputAtomStruct=cls.protPrepareReceptor.outputStructure,
            prevGrid=True, inputGrid=cls.protGridADT.outputGrid,
            nFillPoints=10)
        elif pocketFinder == 'p2rank':
          protPocketFinder = cls.newProtocol(
            P2RankFindPockets,
            inputAtomStruct=cls.protPrepareReceptor.outputStructure)
        elif pocketFinder == 'fpocket':
          protPocketFinder = cls.newProtocol(
            FpocketFindPockets,
            inputAtomStruct=cls.protPrepareReceptor.outputStructure)

        cls.proj.launchProtocol(protPocketFinder, wait=False)
        return protPocketFinder

    @classmethod
    def _runSetFilter(cls, inProt, number, property):
      cls.protFilter = cls.newProtocol(
        ProtSetFilter,
        operation=ProtSetFilter.CHOICE_RANKED,
        threshold=number, rankingField=property)
      cls.protFilter.inputSet.set(inProt)
      cls.protFilter.inputSet.setExtended('outputPockets')

      cls.proj.launchProtocol(cls.protFilter, wait=False)
      return cls.protFilter

    @classmethod
    def _runGridGeneration(cls):
        cls.protGridADT = cls.newProtocol(
          Autodock_GridGeneration,
          inputAtomStruct=cls.protPrepareReceptor.outputStructure,
          radius=24.0,
          spacing=0.6)

        cls.proj.launchProtocol(cls.protGridADT, wait=False)
        return cls.protGridADT


    def _runDARC(self, ADTLigs=False, pocketsProt=None, gridProt=None):
        if ADTLigs:
            protLigs = self.protPrepareLigandADT
        else:
            protLigs = self.protOBabel

        if pocketsProt == None:
            protDARC = self.newProtocol(
                RosettaProtDARC,
                fromPockets=False,
                target_residue='99:C',
                numberOfThreads=8)

            protDARC.inputAtomStruct.set(self.protPrepareReceptor)
            protDARC.inputAtomStruct.setExtended('outputStructure')
            protDARC.inputLigands.set(protLigs)
            protDARC.inputLigands.setExtended('outputSmallMolecules')

            if gridProt == None:
                protDARC.shape_only.set(True)
            else:
                protDARC.shape_only.set(False)
                protDARC.grid.set(gridProt)
                protDARC.grid.setExtended('outputGrid')

            self.launchProtocol(protDARC)
            smOut = getattr(protDARC, 'outputSmallMolecules', None)
            self.assertIsNotNone(smOut)

        else:
            protDARC = self.newProtocol(
                RosettaProtDARC,
                fromPockets=True,
                numberOfThreads=8)

            protDARC.inputPockets.set(self.protFilter)
            protDARC.inputPockets.setExtended('outputPockets')
            protDARC.inputLigands.set(protLigs)
            protDARC.inputLigands.setExtended('outputSmallMolecules')

            if gridProt == None:
                protDARC.shape_only.set(True)
            else:
                protDARC.shape_only.set(False)
                protDARC.grid.set(gridProt)
                protDARC.grid.setExtended('outputGrid')

            self.launchProtocol(protDARC)
            smOut = getattr(protDARC, 'outputSmallMolecules', None)
            self.assertIsNotNone(smOut)

        return protDARC


class TestDARC(TestImportBase):

    def test_1(self):
        """ Complete Docking from whole protein and shape only
        """
        print("\n Complete Docking from whole protein and shape only \n")
        protDARC = self._runDARC()
        self.launchProtocol(protDARC)

    def test_2(self):
        """ Complete Docking from protein pockets and shape only
        """
        print("\n Complete Docking from protein pockets and shape only \n")
        if pocketFinder!=None:
            self._waitOutput(self.protFilter, 'outputPockets', sleepTime=5)
            if ADT:
                protDARC = self._runDARC(ADTLigs=True, pocketsProt=self.protFilter)
            else:
                protDARC = self._runDARC(pocketsProt=self.protFilter)
            self.launchProtocol(protDARC)
        else:
            print('Cannot import any pocket finder (autoligand, p2rank, fpocket).'
                  'Test on pocket will not be performed')

    def test_3(self):
        """ Complete Docking from whole protein and ADT electrostatics
        """
        print("\n Complete Docking from whole protein and ADT electrostatics \n")
        if ADT:
            self._waitOutput(self.protGridADT, 'outputGrid', sleepTime=10)
            protDARC = self._runDARC(ADTLigs=True, gridProt=self.protGridADT)
            self.launchProtocol(protDARC)
        else:
            print('Autodock cannot be imported, docking with electrostatics cannot be made')

    def test_4(self):
        """ Complete Docking from protein pockets and ADT electrostatics
        """
        print("\n Complete Docking from protein pockets and ADT electrostatics \n")
        if ADT:
          protDARC = self._runDARC(gridProt=self.protGridADT,
                                   pocketsProt=self.protFilter)
          self.launchProtocol(protDARC)
        else:
          print('Autodock cannot be imported, docking with electrostatics cannot be made')
