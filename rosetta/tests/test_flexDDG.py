# **************************************************************************
# *
# * Name:     test of protocol_flexDDG.py
# *
# * Authors:    Judith Maestro Ciria
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

from pyworkflow.tests import BaseTest, setupTestProject
from pwem.protocols.protocol_import import ProtImportPdb

from rosetta.protocols.protocol_flexDDG import ProtRosettaFlexDDG


class TestFlexDDG(BaseTest):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def test_single_mutation(self):
        args = {'inputPdbData': 1,
                'pdbFile': '/home/jmaestro/flexDDG_test/inputs/complex_7ZFC_Nb126.pdb'}
        protImport = self.newProtocol(ProtImportPdb, **args)
        self.launchProtocol(protImport)
        target = protImport.outputPdb

        args = {'inputAtomStruct': target,
                'chainsToMove': 'B',
                'toMutateList': 'VB32A',
                'nstruct': 1,
                'numberBackrubTrials': 10,
                'backrubTrajectoryStride': 10,
                'maxMinimizationIter': 5,
                'absScoreConvergenceThresh': 200.0}
        protFlexDDG = self.newProtocol(ProtRosettaFlexDDG, **args)
        self.launchProtocol(protFlexDDG)

        outputStats = getattr(protFlexDDG, 'outputStats', None)
        self.assertIsNotNone(outputStats, 'No outputStats was generated')
        self.assertEqual(len(outputStats), 1, 'Expected exactly 1 mutation result')

        item = list(outputStats)[0]
        self.assertEqual(item.mutation.get(), 'VB32A')
        print('ddG (GAM) = %.3f, ddG (raw) = %.3f' % (item.ddg.get(), item.ddgRaw.get()))

    def test_chainsToMove_from_wizard_json(self):
        """ The chain-selection wizard stores JSON (e.g. {"chain": "B", "model": 0, ...}) in
        chainsToMove instead of a plain chain letter. This reproduces that scenario. """
        args = {'inputPdbData': 1,
                'pdbFile': '/home/jmaestro/flexDDG_test/inputs/complex_7ZFC_Nb126.pdb'}
        protImport = self.newProtocol(ProtImportPdb, **args)
        self.launchProtocol(protImport)
        target = protImport.outputPdb

        args = {'inputAtomStruct': target,
                'chainsToMove': '{"model": 0, "chain": "B", "residues": 121}',
                'toMutateList': 'VB32A',
                'nstruct': 1,
                'numberBackrubTrials': 10,
                'backrubTrajectoryStride': 10,
                'maxMinimizationIter': 5,
                'absScoreConvergenceThresh': 200.0}
        protFlexDDG = self.newProtocol(ProtRosettaFlexDDG, **args)
        self.launchProtocol(protFlexDDG)

        outputStats = getattr(protFlexDDG, 'outputStats', None)
        self.assertIsNotNone(outputStats, 'No outputStats was generated')
        self.assertEqual(len(outputStats), 1, 'Expected exactly 1 mutation result')
