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

"""
This protocol will remove any water, ligand or HETATM that a pdb has and
it will remove redundant chains.

In addition, this protocol will use a Rosetta suite program, named score,
to build the missing atoms in the protein, including the hydrogens.
This is necessary to subsequently generate the electrostatic potential grid.

"""


from pyworkflow import Config
from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pyworkflow.protocol.params import LEVEL_ADVANCED

from pwem.protocols import EMProtocol
from pwem.objects.data import AtomStruct
from pwem.convert.atom_struct import cifToPdb

import os, shutil, json, glob

from rosetta import Plugin
from rosetta.constants import *

from pwchem.utils import clean_PDB


class RosettaProteinPreparation(EMProtocol):
    """
    This protocol will remove any water, ligand or HETATM that a pdb has and
    it will remove redundant chains.

    In addition, this protocol will use a Rosetta suite program (named score)
    to build the missing atoms in the protein, including the hydrogens.
    This is necessary to subsequently generate the electrostatic potential grid.
    """

    _label = 'DARC Protein Preparation'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """
        form.addSection(label=Message.LABEL_INPUT)

        form.addParam("inputAtomStruct", params.PointerParam, pointerClass="AtomStruct",
                      label="Atomic structure",
                      important=True, allowsNull=False,
                      help="Select the atomic structure of the docking target protein")

        form.addParam("addH", params.BooleanParam,
                      label="Add hydrogens",
                      default=True, important=True,
                      help="Build the missing atoms in the protein structure, including the hydrogens.\n"
                           "This is necessary to subsequently generate the electrostatic potential grid")

        clean = form.addGroup("Clean atomic structure")
        clean.addParam("waters", params.BooleanParam,
                       label='Remove waters',
                       default=True, important=True,
                       help='Remove all waters molecules from a pdb file')

        clean.addParam("HETATM", params.BooleanParam,
                       label='Remove ligands HETATM',
                       default=True, important=True,
                       help='Remove all ligands and HETATM contained in the protein')

        clean.addParam("rchains", params.BooleanParam,
                       label='Remove redundant chains',
                       default=False, important=True,
                       help='Remove redundant chains in the proteins')

        clean.addParam("chain_name", params.StringParam,
                       label="Conserved chain",
                       important=True,
                       condition="rchains==True",
                       help="Select the chain on which you want to carry out the "
                            "molecular docking. You must know the protein and structure "
                            "file that you loaded. \n\nFor example, the protein mdm2 "
                            "(4ERF) has a C1 symmetry, which indicates that its chains "
                            "are at least 95% equal, so you would write A, C or E "
                            "(names of the chains).")

        advanced = form.addGroup("Runs option",  expertLevel=LEVEL_ADVANCED)
        advanced.addParam("cseed", params.BooleanParam,
                       label='Use Constant Seed',
                       default=False,
                       help='Use this option to get reproducible results')

        advanced.addParam("seed", params.IntParam,
                       label='Set seed',
                       default=1111111,
                       condition="cseed",
                       help='Set a integer number as constant seed. The default one is 1111111 ')



    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('remove_WLHC')
        if self.addH.get():
            self._insertFunctionStep('score_optH')
        self._insertFunctionStep('createOutput')


    def remove_WLHC(self):
        """ Clean the pdb file from waters and ligands
        """
        # Get a PDB format file to the protein structure
        pdb_ini = self.inputAtomStruct.get().getFileName()
        filename = os.path.splitext(os.path.basename(pdb_ini))[0]
        fnPdb = self._getExtraPath('%s_clean.pdb' % filename)

        if self.rchains.get():
            chain = json.loads(self.chain_name.get())  # From wizard dictionary
            chain_id = chain["Chain"].upper().strip()
        else:
            chain_id = None
        cleanedPDB = clean_PDB(self.inputAtomStruct.get(), fnPdb,
                               self.waters.get(), self.HETATM.get(), chain_id)


    def score_optH(self):
        """ Rosetta will build any missing atoms including hydrogens in
            input protein PDB files.
            We will use the program score.linuxccrelease to do that
        """

        # Take the cleaned pdb file

        pdb_file = glob.glob(self._getExtraPath("*_clean.pdb"))[0]
        name_pdbfile = os.path.splitext(os.path.basename(pdb_file))[0]

        self.name_protein = name_pdbfile.split("_")[0]
        pdb_file = os. path.abspath(pdb_file)

        # Create the args of the program
        args = ""
        args += " -in:file:"
        args += "s %s" % pdb_file # PDB file to add missing atoms

        args += " -out:output -no_optH false"

        args += " -out:file:scorefile %s.sc" % name_pdbfile

        if self.cseed.get():
            args += " -run:constant_seed"
            args += " -run:jran %s" %self.seed.get()

        # Execute the program score to construct missing atoms
        # It creates different files by default in the protocol path:
        #   - <PDB_INPUT_NAME>_0001.pdb
        #   - <PDB_INPUT_NAME>_0001.sc (scorefile (default by program ->default.sc))
        Plugin.runRosettaProgram(Plugin.getProgram(SCORE), args, cwd=self._getPath())

        #Move and rename the files to Path from the ExtraPath
        scoresFile = self._getPath("%s_0001.pdb" % name_pdbfile)
        pdbOut = self._getPath("%s.pdb" % self.name_protein)
        self.cleanScores(scoresFile, pdbOut)
        shutil.move(self._getPath("%s.sc" % name_pdbfile), self._getExtraPath("%s.sc" % self.name_protein))

        self._store()


    def createOutput(self):
        """ Create a object atomStruct as a output from cleaned pdb with missing atoms
        """

        pdb_file = glob.glob(self._getExtraPath("*_clean.pdb"))[0]
        name_pdbfile = os.path.splitext(os.path.basename(pdb_file))[0]
        name_protein = name_pdbfile.split("_")[0]

        if self.addH.get():
            pdb_file_out = self._getPath('%s.pdb' % name_protein)
        else:
            pdb_file_out = pdb_file

        if os.path.exists(os.path.expanduser(pdb_file_out)):
            target = AtomStruct(filename=pdb_file_out)
            self._defineOutputs(outputStructure=target)
            self._defineSourceRelation(self.inputAtomStruct, target)



    # --------------------------- UTILS functions ------------------------------

    def _validate(self):
        """ Validate if the structure input fields are complete
        """
        errors = []

        if self.inputAtomStruct.get() is None:
            errors.append("A pdb file was not entered in the Atomic structure field. Please enter it.")

        # Check that the program exists
        program = Plugin.getProgram(SCORE)
        if not os.path.exists(os.path.expanduser(program)):
            errors.append("Cannot find " + program)

            # If there is any error at this point it is related to config variables
            errors.append("Check configuration file: " +
                          Config.SCIPION_LOCAL_CONFIG)
            errors.append("and set ROSETTA_HOME variables properly.")
            if program is not None:
                errors.append("Current values:")
                errors.append("ROSETTA_HOME = %s" % Plugin.getVar(ROSETTA_HOME))
                errors.append("SCORE = %s" % SCORE)


        return errors


    # --------------------------- INFO functions -----------------------------------

    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []

        if self.isFinished():

            pdb_file = glob.glob(self._getExtraPath("*_clean.pdb"))[0]
            name_pdbfile = os.path.splitext(os.path.basename(pdb_file))[0]
            name_protein = name_pdbfile.split("_")[0]

            filename = os.path.basename(self._getPath('%s.pdb' % name_protein))
            filename_original = os.path.basename(self._getExtraPath('%s.pdb' % name_protein))
            filename_clean = os.path.basename(self._getExtraPath('%s_clean.pdb' % name_protein))
            filename_sc = os.path.basename(self._getExtraPath('%s.sc' % name_protein))

            if self.rchains.get():
                chain = json.loads(self.chain_name.get())  # From wizard dictionary
                chain_total = chain["Number of chains"]
                summary.append("This protocol has removed waters, ligands and chains that have not been selected"
                               " (a total of %s chains was deleted). "
                               "Furthermore, it has added missing atoms to the structure, such as hydrogens. \n\n"
                               "The generated pdb file is in the protocol folder and its name is *% s* \n\n"
                               "Also, in the extra folder you can find:\n\n"
                               "- The original pdb from the target protein --> *%s* \n"
                               "- The cleaned pdb file without missing atoms --> *%s* \n"
                               "- The file with different scores about the protein and his structure --> *%s*  \n"
                               % (str(int(chain_total)-1), filename, filename_original,
                                  filename_clean, filename_sc))
            else:
                summary.append("This protocol has removed waters and ligands."
                               "Furthermore, it has added missing atoms to the structure, such as hydrogens. \n\n"
                               "The generated pdb file is in the protocol folder and its name is *%s* \n\n"
                               "Also, in the extra folder you can find:\n\n"
                               "- The original pdb from the target protein --> *%s* \n"
                               "- The cleaned pdb file without missing atoms --> *%s* \n"
                               "- The file with different scores about the protein and his structure --> *%s*  \n"
                               % ( filename, filename_original, filename_clean, filename_sc))
        else:
            summary.append("Protocol not yet completed")
        return summary


    def _methods(self):
        methods = []
        methods.append("In this protocol, the atomic structure file is parsed in pdb format (if a file is \n "
                       "provided in cif format, it is transformed into pdb with the *ciftopdb* function of \n"
                       "the pwem package) to eliminate waters (HOH), ligands or hetatm atoms and the \n"
                       "protein chains not selected \n\n")

        program = Plugin.getProgram(SCORE)

        methods.append("In addition, to complete the protein in terms of atoms lost in side chains such \n"
                       "as hydrogens, the Rosetta score program is used, which also makes a physical \n"
                       "assessment of the protein using the score function REF2015 (protein.sc file).\n"
                       "The program is located in : *%s*" % program)

        return methods


    def _citations(self):
        return ['LeaverFay2011']



    def cleanScores(self, scoresFile, pdbOut):
        with open(pdbOut, 'w') as f:
            with open(scoresFile) as fIn:
                for line in fIn:
                    if self.isPDBLine(line):
                        f.write(line)

    def isPDBLine(self, line):
        options = ['REMARK', 'ATOM', 'HETATM', 'HEADER', 'EXPDTA', 'TER', 'END']
        for opt in options:
          if line.startswith(opt):
              return True
        return False



