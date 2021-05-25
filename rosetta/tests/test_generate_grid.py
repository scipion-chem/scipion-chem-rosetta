# **************************************************************************
# *
# * Name:     test of protocol_generate_grid.py
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
from pwem.protocols.protocol_import import ProtImportPdb

from rosetta.protocols.protocol_target_preparation import RosettaProteinPreparation as Prepare
from rosetta.protocols.protocol_generate_grid import Autodock_GridGeneration as Grid
from rosetta.objects import GridAGD





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



class TestGridGeneration(TestImportBase):

    def test_1(self):
        """ Generate the grid of a prepared protein by DARC
        """
        print("\n Generate the grid of a prepared protein by DARC \n")

        prot_path = os.path.abspath(os.path.join(self.path_data, "4erf.pdb"))

        # Import PDB as Scipion object and prepare it
        target = self._importPDB(prot_path)
        prep_target = self._prepareProtein(target)

        # Generate grid
        radius = 37
        spacing = 0.5
        args = {'inputpdb': prep_target,
                'radius': radius,
                'spacing': spacing
                }

        protocol = self.newProtocol(Grid, **args)
        self.launchProtocol(protocol)
        grid = protocol.outputGrid
        grid_path = os.path.abspath(grid.getFileName())

        self.assertIsNotNone(grid,
                             "There was a problem with grid generation - Please check it")

        self.assertTrue(grid_path.endswith(".agd"),
                             "Format grid is incorrect")

        """
        with open(grid_path) as f:
            for i, l in enumerate(f):
                pass
        number_lines_agd = i+1  # 857380
        """
        number_lines_agd = int(os.popen('wc -l %s' %grid_path).read().split()[0])  # Using bash

        expected_radius = ((radius + 10) / spacing)
        if expected_radius % 2 == 0:
            expected_lines = (expected_radius+1)**3 + 5
        else:
            expected_lines = (expected_radius)**3 + 5

        # Test number of lines according the selected radius
        self.assertTrue(number_lines_agd == int(expected_lines),
                             "Length grid is incorrect. It should be %s" % int(expected_lines))
