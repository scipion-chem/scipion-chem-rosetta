# **************************************************************************
# *
# * Name:     test of protocol_generate_rays.py
# *
# * Authors:    Alberto M. Parra Pérez (amparraperez@gmail.com)
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

from rosetta.protocols.protocol_target_preparation import RosettaProteinPreparation as Prepare
from rosetta.protocols.protocol_generate_rays import Rosetta_make_rayFile as makeRay
from rosetta.protocols.protocol_generate_grid import Autodock_GridGeneration as makeGrid



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

        protocol = cls.newProtocol(Prepare, **args)
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
    def NonZero(cls, path):
        return os.path.isfile(path) and os.path.getsize(path) > 0



class TestRaysGeneration(TestImportBase):

    def test_1(self):
        """ Generate the rays with 1 target residue and without grid
        """
        print("\n Generate the rays with 1 target residue and without grid \n")

        prot_path = os.path.abspath(os.path.join(self.path_data, "4erf.pdb"))

        # Import PDB as Scipion object and prepare it
        target = self._importPDB(prot_path)
        prep_target = self._prepareProtein(target)

        # Generate rays
        target_residue = 61

        args = {'inputpdb': prep_target,
                'gpuList': 0,  # Not GPU
                'target_residue': target_residue,
                'cseed': True}

        protocol = self.newProtocol(makeRay, **args)
        self.launchProtocol(protocol)
        ray_structure = protocol.outputStructure
        ray_txt = protocol.outputRay_TXT



        self.assertIsNotNone(ray_structure,
                             "There was a problem with rays generation (PDB-RAYS)- Please check it")
        ray_structure_path = ray_structure.getFileName()
        self.assertTrue(ray_structure_path.endswith(".pdb"),
                        "Format ray structure is wrong. It should be pdb")
        self.assertTrue(self.NonZero(ray_structure_path),
                        "Structure ray pdb file is empty")

        self.assertIsNotNone(ray_txt,
                             "There was a problem with rays generation (TXT-RAYS) - Please check it")
        ray_txt_path = ray_txt.getFileName()
        self.assertTrue(ray_txt_path.endswith(".txt"),
                        "Format ray text file is wrong. It should be txt")
        self.assertTrue(self.NonZero(ray_txt_path),
                        "Structure ray text file is empty")



    def test_2(self):
        """ Generate the rays with 2 target residues and with grid
        """
        print("\n Generate the rays with 2 target residues and with electrostatic grid \n")

        prot_path = os.path.abspath(os.path.join(self.path_data, "4erf.pdb"))

        # Import PDB as Scipion object and prepare it
        target = self._importPDB(prot_path)
        prep_target = self._prepareProtein(target)

        grid = self._makegrid(prep_target)

        # Generate rays
        target_residue = 61
        additional_res = '57,56'

        args = {'inputpdb': prep_target,
                'gpuList': 0,  # Not GPU
                'target_residue': target_residue,
                'multiple_origin': True,
                'target_residues': additional_res,
                'electrostatics': True,
                'grid': grid,
                'cseed': True}

        protocol = self.newProtocol(makeRay, **args)
        self.launchProtocol(protocol)
        ray_structure = protocol.outputStructure
        ray_txt = protocol.outputRay_TXT
        agd_ray = protocol.outputGRID_AGD


        self.assertIsNotNone(ray_structure,
                             "There was a problem with rays generation (PDB-RAYS)- Please check it")
        ray_structure_path = ray_structure.getFileName()
        self.assertTrue(ray_structure_path.endswith(".pdb"),
                        "Format ray structure is wrong. It should be pdb")
        self.assertTrue(self.NonZero(ray_structure_path),
                        "Structure ray pdb file is empty")

        self.assertIsNotNone(ray_txt,
                             "There was a problem with rays generation (TXT-RAYS) - Please check it")
        ray_txt_path = ray_txt.getFileName()
        self.assertTrue(ray_txt_path.endswith(".txt"),
                        "Format ray text file is wrong. It should be txt")
        self.assertTrue(self.NonZero(ray_txt_path),
                        "Structure ray text file is empty")

        self.assertIsNotNone(agd_ray,
                             "There was a problem with rays generation (GRID post RAYS)- Please check it")
        agd_ray_path = agd_ray.getFileName()
        self.assertTrue(agd_ray_path.endswith(".agd"),
                        "Format grid file is wrong. It should be agd")
        self.assertTrue(self.NonZero(agd_ray_path),
                        "Grid file file is empty")


