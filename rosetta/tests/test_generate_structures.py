# **************************************************************************
# *
# * Name:     test of protocol_target_preparation.py
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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
from pwem.protocols.protocol_import import ProtImportPdb, ProtImportVolumes
from rosetta.protocols import ProtRosettaGenerateStructures


class TestGenerateStructures(BaseTest):

    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        setupTestProject(cls)

        cls._runImportPDB()
        cls._runImportVolume()

    @classmethod
    def _runImportPDB(cls):
      cls.protImportPDB = cls.newProtocol(
        ProtImportPdb,
        inputPdbData=1,
        pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1.pdb'))
      cls.launchProtocol(cls.protImportPDB)

    @classmethod
    def _runImportVolume(cls):
      args = {'filesPath': cls.ds.getFile('volumes/emd_3488.map'),
        'samplingRate': 1.05, 'setOrigCoord': True,
        'x': 0.0, 'y': 0.0, 'z': 0.0
      }
      cls.protImportVolume = cls.newProtocol(ProtImportVolumes, **args)
      cls.launchProtocol(cls.protImportVolume)

    @classmethod
    def _runGenerateStructures(cls):
      cls.protGenStructures = cls.newProtocol(
        ProtRosettaGenerateStructures,
        inputStructure=cls.protImportPDB.outputPdb,
        inputVolume=cls.protImportVolume.outputVolume,
        resolution=1.05, numMods=2)

      cls.launchProtocol(cls.protGenStructures)

    def test(self):
        self._runGenerateStructures()

        self.assertIsNotNone(getattr(self.protGenStructures, 'outputAtomStructs', None))

