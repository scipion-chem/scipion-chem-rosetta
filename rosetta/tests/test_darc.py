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
from pathlib import Path

from pyworkflow.tests import *
from pyworkflow.protocol import *
from pwem.protocols.protocol_import import ProtImportPdb

from rosetta.protocols.protocol_target_preparation import RosettaProteinPreparation as PrepareProtein
from rosetta.protocols.protocol_ligand_preparation import Rosetta_ligand_preparation as PrepareLigand
from rosetta.protocols.protocol_generate_rays import Rosetta_make_rayFile as makeRay
from rosetta.protocols.protocol_generate_grid import Autodock_GridGeneration as makeGrid
from rosetta.protocols.protocol_darc import Rosetta_darc as darc

from bioinformatics.protocols.protocol_import_smallMolecules import ProtBioinformaticsImportSmallMolecules as importSM




class TestImportBase(BaseTest):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        path_test = Path(__file__).parent
        cls.path_data = os.path.join(path_test, "data")

    @classmethod
    def _importPDB(cls, path):
        inputPdbData = 1  # file
        args = {'inputPdbData': inputPdbData,
                'pdbFile': path
                }

        protocol = cls.newProtocol(ProtImportPdb, **args)
        cls.launchProtocol(protocol)
        pdb = protocol.outputPdb
        return pdb


    @classmethod
    def _prepareProtein(cls, target):

        args = {'inputpdb': target,
                'addH': True,
                'waters': True,
                'HETATM': True,
                "rchains": True,
                "chain_name": '{"Chain": "A", "Number of residues": 92, "Number of chains": 3}',
                "cseed": True
                }

        protocol = cls.newProtocol(PrepareProtein, **args)
        cls.launchProtocol(protocol)
        pdb = protocol.outputStructure
        return pdb


    @classmethod
    def _makegrid(cls, target):
        radius = 37
        spacing = 0.5
        args = {'inputpdb': target,
                'radius': radius,
                'spacing': spacing
                }

        protocol = cls.newProtocol(makeGrid, **args)
        cls.launchProtocol(protocol)
        gridAGD = protocol.outputGrid

        return gridAGD

    @classmethod
    def _makeray(self, prep_target):

        # Generate rays
        target_residue = 61

        args = {'inputpdb': prep_target,
                'gpuList': 0,  # Not GPU
                'target_residue': target_residue,
                'electrostatics': False,
                'cseed': True}

        protocol = self.newProtocol(makeRay, **args)
        self.launchProtocol(protocol)
        ray_structure = protocol.outputStructure
        ray_txt = protocol.outputRay_TXT

        return ray_structure, ray_txt


    @classmethod
    def _importSmallM(cls, path):
        inputPdbData = 1  # file
        args = {'multiple': inputPdbData,
                'filesPath': path,
                'filesPattern': '*'}

        protocol = cls.newProtocol(importSM, **args)
        cls.launchProtocol(protocol)
        setofSM = protocol.outputSmallMols
        return setofSM

    @classmethod
    def _prepareLigand(self, path_data):
        ligand_path = os.path.abspath(os.path.join(path_data, "smallMolecules", "mol2"))

        # Import SetOfSmallMolecules
        smallM = self._importSmallM(ligand_path)

        args = {'inputType': 0,  # SmallMolecules
                'inputSmallMols': smallM,
                'method_charges': 0,  # gasteiger
                'conformer': True,
                "method_conf": 0,
                "number_conf": 100,
                "rmsd_cutoff": 0.375,
                }

        protocol = self.newProtocol(PrepareLigand, **args)
        self.launchProtocol(protocol)
        small = protocol.outputSmallMols
        return small

    @classmethod
    def NonZero(cls, path):
        return os.path.isfile(path) and os.path.getsize(path) > 0



class TestDARC(TestImportBase):

    def test_1(self):
        """ Complete Docking with GRID and search conformers on the fly
        """
        print("\n Complete Docking with GRID and search conformers on the fly \n")

        prot_path = os.path.abspath(os.path.join(self.path_data, "4erf.pdb"))

        # Import PDB as Scipion object and prepare it
        target = self._importPDB(prot_path)
        prep_target = self._prepareProtein(target)

        # Generate grid and rays
        grid = self. _makegrid(prep_target)
        ray_structure, ray_txt= self._makeray(target)

        # Import 4 molecules as Scipion object and prepare it
        set_small = self._prepareLigand(self.path_data)

        args = {'protein': prep_target,
                'gpuList': 1,  # Not GPU
                'ligands': set_small,
                'ray_file': ray_txt,
                'shape_only': False,
                'grid': grid,
                'search_conformers': True,
                'cseed': True}


        protocol = self.newProtocol(darc, **args)
        self.launchProtocol(protocol)
        scores = protocol.DarcScores



        self.assertIsNotNone(scores,
                             "There was a problem with the entire docking)- Please check it")


        self.assertTrue(scores.getSize() == 4,
                        "There was a problem with Rosetta DARC. The number of complexes correctly evaluated "
                        "is less than the number given ")

        n_columns = len(list(scores.getFirstItem().getAttributes()))
        self.assertTrue(n_columns == 2,
                        "There is a incorrect number of columns. Therefore there are a lost of information")  # 4 fixed columns + 2 given

        list_ID = ["ZINC00001099", "ZINC00001453", "ZINC00000480", "ZINC00001019"]  # Molecule ID
        list_results = [239, 403, 496, 613]  # expected score
        i = 0

        for sc in scores:
            self.assertTrue(sc.getID() == list_ID[i],
                            "The molecule ID and the order does not correspond to the expected")
            self.assertTrue(int(sc.getScoreDarc())==list_results[i],
                            "The score and the order does not correspond to the expected ")

            i += 1




