# **************************************************************************
# *
# * Name:     test of protocol_target_preparation.py
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
from pyworkflow.tests import *
from pathlib import Path
from pwem.protocols.protocol_import import ProtImportPdb
from pwem.convert.atom_struct import AtomicStructHandler
from rosetta.protocols.protocol_target_preparation import RosettaProteinPreparation as Prepare

import numpy as np



class TestImportBase(BaseTest):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        path_test = Path(__file__).parent

    @classmethod
    def _importPDB(cls):
        args = {'inputPdbData': 0,
                'pdbId': '4erf'
                }

        protocol = cls.newProtocol(ProtImportPdb, **args)
        cls.launchProtocol(protocol)
        pdb = protocol.outputPdb
        return pdb

    @classmethod
    def _StructureInformation(cls, protein):
        structureHandler = AtomicStructHandler()
        structureHandler.read(protein)
        chains, residues = structureHandler.getModelsChains()

        chains = list(chains.items())[0][1] # Chains
        n_chains = len(chains)

        n_residues = dict(chains).values() # Residues, without including ligands, waters and others
        n_residues_total = sum(n_residues)


        n_waters = [] # Waters in each chain
        not_waters = []
        residues = dict(list(residues.items())[0][1])
        for chain in residues.values():
            n = 0
            no = 0
            for i in chain:
                if i[1] == 'HOH':
                    n += 1
                else:
                    no += 1
            n_waters.append(n)
            not_waters.append(no)

        n_waters_total = sum(n_waters)
        n_ligands = sum(not_waters) - n_residues_total

        atoms = structureHandler.getStructure().get_atoms()
        coords = [np.append(atom.get_coord(), 1) for atom in atoms]
        n_atoms = len(coords)

        return n_chains, n_residues_total, list(n_residues), n_waters_total, n_waters, n_ligands, n_atoms





class TestTargetPreparation(TestImportBase):


    def test_1(self):
        """ Remove from protein structure only waters
        """
        print("\n Removing water from a protein structure \n")

        # Import PDB as Scipion object
        target = self._importPDB()

        args = {'inputAtomStruct': target,
                'addH': False,
                'waters': True,
                'HETATM': False,
                "rchains": False,
                "cseed": True
                }

        protocol = self.newProtocol(Prepare, **args)
        self.launchProtocol(protocol)
        pdb_out = protocol.outputStructure
        pdb_out_path = os.path.abspath(pdb_out.getFileName())

        n_chains_0, n_residues_total_0, n_residues_chain_0, \
        n_waters_total_0, n_waters_0, n_ligands_0, n_atoms_0 = self._StructureInformation(target.getFileName())
        n_chains, n_residues_total, n_residues_chain, \
        n_waters_total, n_waters, n_ligands, n_atoms = self._StructureInformation(pdb_out_path)


        self.assertIsNotNone(pdb_out,
                             "There was a problem with output generation - Please check it")

        self.assertTrue(n_chains == n_chains_0 and n_chains==3,
                        "Some protein CHAIN was removed without having to")

        self.assertTrue(n_residues_total_0 == n_residues_total,
                        "Some protein RESIDUES was removed without having to")

        self.assertTrue(n_ligands_0 == n_ligands and n_ligands==3,
                        "Some LIGANDS was removed without having to")

        self.assertTrue(n_waters_total == 0,
                        "All the WATERS were not removed correctly")




    def test_2(self):
        """ Remove from protein structure only ligands and 2 chains.
            We use a cif format to check the conversion to pdb
        """
        print("\n Removing ligands and chains from a protein structure in cif format \n")

        # Import PDB as Scipion object
        target = self._importPDB()

        args = {'inputAtomStruct': target,
                'addH': False,
                'waters': False,
                'HETATM': True,
                "rchains": True,
                "chain_name": '{"model": 0, "chain": "A", "residues": 92}',
                "cseed": True
                }

        protocol = self.newProtocol(Prepare, **args)
        self.launchProtocol(protocol)
        pdb_out = protocol.outputStructure
        pdb_out_path = os.path.abspath(pdb_out.getFileName())


        n_chains_0, n_residues_total_0, n_residues_chain_0, \
        n_waters_total_0, n_waters_0, n_ligands_0, n_atoms_0 = self._StructureInformation(target.getFileName())
        n_chains, n_residues_total, n_residues_chain, \
        n_waters_total, n_waters, n_ligands, n_atoms = self._StructureInformation(pdb_out_path)


        self.assertIsNotNone(pdb_out,
                             "There was a problem with output generation - Please check it")

        end = (os.path.basename(pdb_out.getFileName()).endswith(".pdb"))
        self.assertTrue(end == True,
                             "The cif to pdb conversion failed. Please check it")

        self.assertTrue(n_chains == 1,
                        "The number of chains is wrong since it have to be one (chain A)")

        self.assertTrue(n_residues_chain_0[0] == n_residues_chain[0],
                        "Some protein RESIDUES was removed without having to")

        self.assertTrue(n_waters_0[0] == n_waters[0],
                        "Some WATERS was removed without having to")

        self.assertTrue(n_ligands == 0,
                        "All the LIGANDS were not removed correctly")



    def test_3(self):
        """ Adding Hydrogens.
            We use a cif format to check the conversion to pdb
        """
        print("\n Adding Hydrogens in protein in a cif format \n")

        # Import PDB as Scipion object
        target = self._importPDB()

        args = {'inputAtomStruct': target,
                'addH': True,
                'waters': False,
                'HETATM': False,
                "rchains": False,
                "cseed": True
                }

        protocol = self.newProtocol(Prepare, **args)
        self.launchProtocol(protocol)
        pdb_out = protocol.outputStructure
        pdb_out_path = os.path.abspath(pdb_out.getFileName())


        n_chains_0, n_residues_total_0, n_residues_chain_0, \
        n_waters_total_0, n_waters_0, n_ligands_0, n_atoms_0 = self._StructureInformation(target.getFileName())
        n_chains, n_residues_total, n_residues_chain, \
        n_waters_total, n_waters, n_ligands, n_atoms = self._StructureInformation(pdb_out_path)


        self.assertIsNotNone(pdb_out,
                             "There was a problem with output generation - Please check it")

        end = (os.path.basename(pdb_out.getFileName()).endswith(".pdb"))
        self.assertTrue(end == True,
                             "The cif to pdb conversion failed. Please check it")

        self.assertTrue(n_chains == n_chains_0,
                        "Some protein CHAIN was removed without having to")

        self.assertTrue(n_residues_total_0 == n_residues_total,
                        "Some protein RESIDUES was removed without having to")

        self.assertTrue(n_waters_total==0,
                        "Score programa removed all waters. If fail here, check if it works properly")

        self.assertTrue(n_ligands_0 == n_ligands,
                        "Some LIGANDS was removed without having to")

        self.assertTrue(n_atoms_0 < n_atoms,
                        "Hydrogens not added")




    def test_4(self):
        """ Clean and add hydrogens
        """
        print("\n Add hydrogens and completed clean of protein \n")

        # Import PDB as Scipion object
        target = self._importPDB()

        args = {'inputAtomStruct': target,
                'addH': True,
                'waters': True,
                'HETATM': True,
                "rchains": True,
                "chain_name": '{"model": 0, "chain": "A", "residues": 92}',
                "cseed": True
                }

        protocol = self.newProtocol(Prepare, **args)
        self.launchProtocol(protocol)
        pdb_out = protocol.outputStructure
        pdb_out_path = os.path.abspath(pdb_out.getFileName())


        n_chains_0, n_residues_total_0, n_residues_chain_0, \
        n_waters_total_0, n_waters_0, n_ligands_0, n_atoms_0 = self._StructureInformation(target.getFileName())
        n_chains, n_residues_total, n_residues_chain, \
        n_waters_total, n_waters, n_ligands, n_atoms = self._StructureInformation(pdb_out_path)


        self.assertIsNotNone(pdb_out,
                             "There was a problem with output generation - Please check it")

        self.assertTrue(n_chains == 1,
                        "The number of chains is wrong since it have to be one (chain A)")

        self.assertTrue(n_residues_chain_0[0] == n_residues_chain[0],
                        "Some protein RESIDUES was removed without having to")

        self.assertTrue(n_waters[0] == 0,
                        "All the WATERS were not removed correctly")

        self.assertTrue(n_ligands == 0,
                        "All the LIGANDS were not removed correctly")

        self.assertTrue((n_atoms_0/3) < n_atoms,
                        "Hydrogens not added")


