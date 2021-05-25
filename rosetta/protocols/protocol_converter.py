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

from pyworkflow.protocol.params import PointerParam, EnumParam
from pwem.objects.data import AtomStruct
from pwem.protocols import EMProtocol

from bioinformatics.objects import SetOfSmallMolecules, SmallMolecule

from rosetta import Plugin


def inputArg(fn):

    if fn.endswith('.pdb'):    # Protein Data Bank
        args = "-ipdb"
    elif fn.endswith('.cif'):  # cif (crystallography information)
        args = "-icif"

    elif fn.endswith('.sdf'):  # MDL MOL FORMAT
        args = "-isdf"
    elif fn.endswith('.sd'):   # MDL MOL FORMAT
        args = "-isd"

    elif fn.endswith('.mol2'): # Sybyl Mol2 format (3D)
        args = "-imol2"
    elif fn.endswith('.smi') or fn.endswith('.smiles'):  # Smiles format (2D)
        args = "-ismi"
    else:
        error = " Input format was not recognize (The allowed format are pdb, cif, mol2, sdf and smi)"
        raise Exception(error)
    return args + " %s" %fn



def outputArg(fnRoot, format, protocol):

    if format == 0:
        fnOut = protocol._getExtraPath(fnRoot + ".pdb")
        args = " -O %s" % (fnRoot + ".pdb")
    elif format == 1:
        fnOut = protocol._getExtraPath(fnRoot + ".cif")
        args = " -O  %s" % (fnRoot + ".cif")
    elif format == 2:
        fnOut = protocol._getExtraPath(fnRoot + ".mol2")
        args = " -O  %s" % (fnRoot + ".mol2")
    elif format == 3:
        fnOut = protocol._getExtraPath(fnRoot + ".sdf")
        args = " -O  %s" % (fnRoot + ".sdf")
    elif format == 4:
        fnOut = protocol._getExtraPath(fnRoot + ".smi")
        args = " -O  %s" % (fnRoot + ".smi")

    return fnOut, args


class ConvertStructures(EMProtocol):
    """Convert a set of input ligands or a receptor structure to a specific file format"""
    _label = 'Convert format'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputType', EnumParam, default=0,
                       choices=["Small molecules", 'Target structure'],
                       label='Input type')

        form.addParam('inputSmallMols', PointerParam, pointerClass="SetOfSmallMolecules", condition='inputType==0',
                      label='Set of small molecules:', allowsNull=False)

        form.addParam('outputFormatSmall', EnumParam, default=2, condition='inputType==0',
                       choices=['PDB', 'cif', 'Mol2', 'SDF', 'Smiles'],
                       label='Output format',
                       help = "If you try to convert a 2D format (ex. smi) to 3D format,"
                              "m you can do this but is wrong")


        form.addParam('inputStructure', PointerParam, pointerClass= "AtomStruct",
                      condition='inputType==1',
                      label='Input structure:', allowsNull=False)

        form.addParam('outputFormatTarget', EnumParam, default=2,
                      condition='inputType==1',
                       choices=['PDB', 'cif', 'Mol2'],
                       label='Output format')


        # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertStep')

    def convertStep(self):

        if self.inputType==0:  # Small molecules
            outputSmallMolecules = SetOfSmallMolecules().create(path=self._getPath(), suffix='SmallMols')

            for mol in self.inputSmallMols.get():

                fnSmall = os.path.abspath(mol.smallMoleculeFile.get())
                fnRoot = os.path.splitext(os.path.split(fnSmall)[1])[0]

                args = inputArg(fnSmall)
                fnOut, argout = outputArg(fnRoot, self.outputFormatSmall.get(), self)

                args += argout


                Plugin.runOPENBABEL(self, args=args, cwd=os.path.abspath(self._getExtraPath()))


                smallMolecule = SmallMolecule(smallMolFilename=fnOut)
                outputSmallMolecules.append(smallMolecule)

            if len(outputSmallMolecules) > 0:
                self._defineOutputs(outputSmallMols=outputSmallMolecules)
                self._defineSourceRelation(self.inputSmallMols, outputSmallMolecules)

        else:
            fnStructure = self.inputStructure.get().getFileName()
            args = inputArg(os.path.abspath(fnStructure))
            fnRoot = os.path.splitext(os.path.split(fnStructure)[1])[0]

            fnOut, argout = outputArg(fnRoot, self.outputFormatTarget.get(), self)

            args += argout

            Plugin.runOPENBABEL(self, args=args, cwd=os.path.abspath(self._getPath()))


            target = AtomStruct(filename=fnOut)
            self._defineOutputs(outputStructure=target)
            self._defineSourceRelation(self.inputStructure, target)



    def _summary(self):
        summary=[]
        if self.inputType.get()==0:
            if self.outputFormatSmall.get() == 0:
                summary.append('Converted to PDB')
            elif self.outputFormatTarget.get() == 1:
                summary.append('Converted to Cif')
            elif self.outputFormatSmall.get() == 2:
                summary.append('Converted to Mol2')
            elif self.outputFormatSmall.get() == 3:
                summary.append('Converted to SDF')
            elif self.outputFormatSmall.get() == 4:
                summary.append('Converted to Smiles')

        else:
            if self.outputFormatTarget.get() == 0:
                summary.append('Converted to PDB')
            elif self.outputFormatTarget.get() == 1:
                summary.append('Converted to Cif')
            elif self.outputFormatTarget.get() == 2:
                summary.append('Converted to Mol2')
        return summary