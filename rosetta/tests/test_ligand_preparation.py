# ***************************************import numpy as np***********************************
# *
# * Name:     test of protocol_ligand_preparation.py
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
from collections import Counter

from pyworkflow.tests import *
from pathlib import Path
from bioinformatics.protocols.protocol_import_smallMolecules import ProtBioinformaticsImportSmallMolecules as importSM


from rosetta.protocols.protocol_ligand_preparation import Rosetta_ligand_preparation as Prepare




class TestImportBase(BaseTest):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        path_test = Path(__file__).parent
        cls.path_data = os.path.join(path_test, "data")

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



class TestLigandPreparation(TestImportBase):

    def test_1(self):
        """ Prepare a set of 4 ligands with conformer generation (genetic algotithm)
        """
        print("\n Prepare a set of 4 ligands with conformer generation (genetic algotithm) \n")

        ligand_path = os.path.abspath(os.path.join(self.path_data, "smallMolecules", "mol2"))

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

        protocol = self.newProtocol(Prepare, **args)
        self.launchProtocol(protocol)
        small_1 = protocol.outputSmallMols


        self.assertIsNotNone(small_1,
                             "There was a problem with the import")
        self.assertTrue(small_1.getSize()==4,
                        "There was a problem with the import or conversion and the SetOfSmallMolecules is empty")

        n_columns = len(list(small_1.getFirstItem().getAttributes()))
        self.assertTrue(n_columns == 5,
                        "There is a incorrect number of columns. Therefore there are a lost of files")  # 4 fixed columns + 5 given

        list_conformers = [100, 100, 36, 48] # number of expected conformers
        i = 0

        for mol in small_1:
            self.assertTrue((mol.getFileName()).endswith("_withH.mol2"),
                            "The format of first molecule is wrong. It must be in mol2 format")
            try:
                self.assertTrue((mol.getConformersFileName()).endswith("_conformers.mol2"),
                            "The format of conformers molecules is wrong. It must be in mol2 format")


                with open(os.path.abspath(mol.getConformersFileName())) as f:
                    lines = dict(Counter([line for line in f if "@<TRIPOS>MOLECULE" in line]))
                    n_conformers = list(lines.values())[0]
                    self.assertTrue(n_conformers==list_conformers[i],
                                    "Number of conformers is wrong. It is %s and it should be %s" %(n_conformers, list_conformers[i]))
                    i += 1

            except:
                self.assertTrue((mol.getConformersFileName()).endswith("Not available"),
                                "Something was wrong in column of _ConformessFile")


            try:
                self.assertTrue((mol.getParamsFileName()).endswith(".params"),
                                "The format of params file is wrong. It must be in params format")
            except:
                self.assertTrue((mol.getParamsFileName()).endswith("Not available"),
                                "Something was wrong in column of _ParamsFile")

            try:
                self.assertTrue((mol.getPDBFileName()).endswith(".pdb"),
                                "The format of pdb file is wrong. It must be in pdb format")
            except:
                self.assertTrue((mol.getPDBFileName()).endswith("Not available"),
                                "Something was wrong in column of _PDBFile")




    def test_2(self):
        """ Prepare a set of 4 ligands with conformer generation (confab algotithm)
        """
        print("\n Prepare a set of 4 ligands with conformer generation (confab algotithm) \n")

        ligand_path = os.path.abspath(os.path.join(self.path_data, "smallMolecules", "pdb"))

        # Import SetOfSmallMolecules
        smallM = self._importSmallM(ligand_path)

        args = {'inputType': 0,  # SmallMolecules
                'inputSmallMols': smallM,
                'method_charges': 1,  # mmff94
                'conformer': True,
                "method_conf": 1,
                "number_conf": 1000,
                "rmsd_cutoff": 0.5,
                }

        protocol = self.newProtocol(Prepare, **args)
        self.launchProtocol(protocol)
        small_1 = protocol.outputSmallMols


        self.assertIsNotNone(small_1,
                             "There was a problem with the import")
        self.assertTrue(small_1.getSize() == 4,
                        "There was a problem with the import or conversion and the SetOfSmallMolecules is empty")

        n_columns = len(list(small_1.getFirstItem().getAttributes()))
        self.assertTrue(n_columns == 5,
                        "There is a incorrect number of columns. Therefore there are a lost of files")  # 4 fixed columns + 5 given

        for mol in small_1:
            self.assertTrue((mol.getFileName()).endswith("_withH.mol2"),
                            "The format of first molecule is wrong. It must be in mol2 format")
            try:
                self.assertTrue((mol.getConformersFileName()).endswith("_conformers.mol2"),
                                "The format of conformers molecules is wrong. It must be in mol2 format")

            except:
                self.assertTrue((mol.getConformersFileName()).endswith("Not available"),
                                "Something was wrong in column of _ConformersFile")

            try:
                self.assertTrue((mol.getParamsFileName()).endswith(".params"),
                                "The format of params file is wrong. It must be in params format")
            except:
                self.assertTrue((mol.getParamsFileName()).endswith("Not available"),
                                "Something was wrong in column of _ParamsFile")

            try:
                self.assertTrue((mol.getPDBFileName()).endswith(".pdb"),
                                "The format of pdb file is wrong. It must be in pdb format")
            except:
                self.assertTrue((mol.getPDBFileName()).endswith("Not available"),
                                "Something was wrong in column of _PDBFile")