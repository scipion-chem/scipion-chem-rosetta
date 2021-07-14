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
Protocol Steps:
0. Input: Get the ligands through the ZINC protocol or import molecules (setofSmallMolecules)
1. Convert to mol2 and pdb format for images (openbabel)
2. Add hydrogens (obabel -h) and remove water
    2.1 Assign charges (Add gasteiger charges (ADT or computeGasteigerCharges de RDKit or babel --partialcharges mmff94 or gasteiger)
3. Generate low energy conformers (openbabel with Confab or RDKIT AllChem.EmbedMolecule)
4. Obtain the .params file from the conformers (rosetta bach_molfile_to_params)
"""

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params
from pyworkflow.protocol.params import LEVEL_ADVANCED
from pyworkflow.utils import Message
import pyworkflow.object as pwobj

import os
import re
import glob
import shutil

from pwchem import Plugin as Pbio
from pwchem.objects import SetOfSmallMolecules, SmallMolecule

from rosetta import Plugin
from rosetta.constants import *



class Rosetta_ligand_preparation(EMProtocol):
    """
    Prepare a set of molecules for use in a docking program (for example, Rosetta DARC).
    Sets the partial atomic charges, generates low-energy conformers, and generates the
    required parameter file for Rosetta DARC. This is done with OpenBabel.
    """

    _label = 'DARC Ligand preparation'
    _dic_method = {0: "gasteiger",
                  1: "mmff94",
                  2: "qeq",
                  3: "qtpie",
                  4: "eqeq",
                  5: "eem",
                  6: "none"}


    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSmallMols', params.PointerParam, pointerClass="SetOfSmallMolecules",
                      label='Set of small molecules:', allowsNull=False,
                      help='It must be in pdb or mol2 format, you may use the converter')

        H_C = form.addGroup("Charges assignation")
        H_C.addParam('method_charges', params.EnumParam,
                      choices=list(self._dic_method.values()),
                      default=0,
                      label='Assign charges method', allowsNull=False,
                      help='Choose the method to add partial charges to small molecule atoms. OpenBabel is used.')

        H_C.addParam('ph', params.BooleanParam,
                      default=False,
                      label='Do you want to set pH of the ligand?:',  # Next step --> PRopKA
                      help='Set the pH of the ligand environment using OpenBabel program. \n'
                           'In this way, the ligand will present more or less hydrogens depending on the pH. '
                           'Note that openbabel (program that add H) has not a general model of pH-dependent '
                           'protonation. It has a set of rules in a file called phmodel.txt,'
                           ' particularly for amino acids. \n\n For more accurate pKa determination and H addition'
                           ', you will probably want semiempirical quantum calculations or a more complete model '
                           '(It is not trivial).',
                      expertLevel=LEVEL_ADVANCED)

        H_C.addParam('phvalue', params.FloatParam, condition="ph",
                      default=7.4,
                      label='pH value:',
                      help='Set the pH of the ligand environment.',
                      expertLevel=LEVEL_ADVANCED)


        form.addParam('conformer', params.BooleanParam, default = False,
                     label='Do you want to generate conformers? ', allowsNull=False,
                     help='You can produce conformers of the ligand in order to do a better rigid docking')


        conformers = form.addGroup("Conformers generation", condition = "conformer")
        conformers.addParam('method_conf', params.EnumParam,
                            choices=["OpenBabel Genetic Algorithm", "OpenBabel Confab"],
                            default=0,
                            label='Method of conformers generation',
                            help='Method of conformers generation. If Confab fails due to the impossibility '
                                 'of assigning a force fields (there is a possibility that it may occur), you should'
                                 'use Genetic Algorithm generator ')


        conformers.addParam('number_conf', params.IntParam,
                            default=200,
                            label='Max. number of conformers:',
                            help='Set the number of conformers generated by OpenBabel from the same molecule.')

        conformers.addParam('rmsd_cutoff', params.FloatParam, condition = "method_conf != 0",
                            default=0.5,
                            label='RMSD cutoff:',
                            help='Set the number of conformers generated by OpenBabel from the same molecule.',
                            expertLevel=LEVEL_ADVANCED)



    # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('addcharges')
        if self.conformer.get():
            self._insertFunctionStep('conformer_generation')
        self._insertFunctionStep('get_params_file')
        self._insertFunctionStep('createOutput')



    def addcharges(self):
        """ Assign the charges using a method available
            in the open-access and free program openbabel
        """
        for mol in self.inputSmallMols.get():
            fnSmall = mol.smallMoleculeFile.get()  # File paths
            fnMol = os.path.split(fnSmall)[1]      # Name of complete file
            fnRoot = os.path.splitext(fnMol)[0]    # Molecule name: ID
            fnFormat = os.path.splitext(fnMol)[1]  # Format file


            # 0. Convert mol2 or pdb format to sdf to ensure the correct assignation of charges
            if fnFormat == ".mol2":
                args = " -imol2 %s --partialcharge none -O %s.sdf" % (os.path.abspath(fnSmall), fnRoot)
                Plugin.runOPENBABEL(self, args=args, cwd=os.path.abspath(self._getTmpPath()))

            elif fnFormat == ".pdb":
                args = " -ipdb %s --partialcharge none -O %s.sdf" % (os.path.abspath(fnSmall), fnRoot)
                Plugin.runOPENBABEL(self, args=args, cwd=os.path.abspath(self._getTmpPath()))
            else:
                raise Exception("Molecules must be in pdb or mol2 format")


        # Run over all sdf files generated
        for filesdf in glob.glob(self._getTmpPath("*")):
            fnRoot = re.split(".s", os.path.split(filesdf)[1])[0]

            # 1. Add all hydrogens or add hydrogens depending on the desirable pH with babel (-p)
            # 2. Add and calculate partial charges with different methods
            index_method = self.method_charges.get()
            cmethod = self._dic_method[index_method]

            # With a given pH
            if self.ph.get():
                args = " -isdf %s -p %s --partialcharge %s -O %s_withH.mol2" % (os.path.abspath(filesdf),
                                                                                str(self.phvalue.get()), cmethod,
                                                                                fnRoot)

                Plugin.runOPENBABEL(self, args=args, cwd=os.path.abspath(self._getExtraPath()))

            else:
                args = " -isdf %s -h --partialcharge %s -O %s_withH.mol2" % (os.path.abspath(filesdf), cmethod, fnRoot)
                Plugin.runOPENBABEL(self, args=args, cwd=os.path.abspath(self._getExtraPath()))



    def conformer_generation(self):
        """ Generate a number of conformers of the same small molecule in mol2 format with
            openbabel using two different algorithm
        """

        # 3. Generate mol2 conformers file for each molecule with OpenBabel

        for file in glob.glob(self._getExtraPath("*_withH.mol2")):
            fnRoot = re.split("_", os.path.split(file)[1])[0]  # ID or filename without _withH.mol2

            if self.method_conf.get() == 0:  # Genetic algorithm
                args =" %s --conformer --nconf %s --score rmsd --writeconformers " \
                      "-O %s_conformers.mol2" %(os.path.abspath(file), str(self.number_conf.get()), fnRoot)
            else:  # confab
                args = " %s -O %s_conformers.mol2 --confab --original --verbose " \
                       "--conf %s --rcutoff %s" % (os.path.abspath(file), fnRoot, str(self.number_conf.get()),
                                                  str(self.rmsd_cutoff.get()))

            Plugin.runOPENBABEL(self, args=args, cwd=os.path.abspath(self._getExtraPath()))



    def get_params_file(self):
        """ Generate params file that DarC will use to dock the ligand in the target protein
        """
        # 1 Create a list.txt with the path of conformers mol2 files
        with open(self._getExtraPath("molfile_list.txt"), "w+") as file:
            if self.conformer.get():
                for path in glob.glob(self._getExtraPath("*_conformers.mol2")):
                    file.write(os.path.abspath(path) + "\n")
            else:
                for path in glob.glob(self._getExtraPath("*_withH.mol2")):
                    file.write(os.path.abspath(path) + "\n")


        # 2. Launch batch_molfile_to_params.py for each file. It will generate a pdb file and params file
        args = ""
        database_path = os.path.join(Plugin.getHome(), ROSETTA_DATABASE_PATH)
        args += " -d %s " % database_path
        mol2params_path = os.path.join(Plugin.getHome(), ROSETTA_PARAMS_PATH, PARAMS_FILE)
        args += "--script_path %s " % mol2params_path
        mollist = os.path.abspath(self._getExtraPath("molfile_list.txt"))
        args += "%s" % mollist

        # Execute the program bach_molfile_to_params to create the params file that will be used by DARC programs
        # It creates a directory called params with several directories, each one called as the mol2 file.
        # Inside of each one, we can find:
        #   - 000.params
        #   - 000_conformers.pdb
        #   - log.txt
        Plugin.runRosettaProgram(Plugin.getProgram(BATCH_PARAMS_FILE, path=ROSETTA_PARAMS_PATH), args=args,
                                 cwd=os.path.abspath(self._getPath()))



    def createOutput(self):
        """Create a set of Small Molecules as output with the path to:
              - Path to small molecule with H (mol2 format)
              - Path to conformers file (mol2 format)
              - Path to params file
              - Image of the conformers pdb if it is possible
        """

        outputSmallMolecules = SetOfSmallMolecules().create(path=self._getPath(), suffix='_for_DARC')


        for fnSmall in glob.glob(self._getExtraPath("*_withH.mol2")):
            fnRoot = re.split("_", os.path.split(fnSmall)[1])[0] # ID or filename without _withH.mol2

            smallMolecule = SmallMolecule(smallMolFilename=fnSmall)

            if self.conformer.get(): # Conformers path if it is exist
                try:
                    if os.path.exists(os.path.join(self._getExtraPath(), "%s_conformers.mol2" % fnRoot)):
                        smallMolecule._ConformersFile = pwobj.String(self._getExtraPath("%s_conformers.mol2" % fnRoot))

                    # Path to getPath/params/<fnRoot>_conformers
                    path_params = os.path.join(self._getPath(), "params", "%s_conformers" %fnRoot)

                    path, dirs, files = next(os.walk(path_params))
                    file_count = len(files)

                    if os.path.exists(path_params) and file_count==3:
                        for file in glob.glob(os.path.join(path_params,"*.params")):
                            if os.path.isfile(file) and file.endswith('.params'):
                                fout = os.path.abspath(os.path.join(path_params, "%s.params" % fnRoot))
                                shutil.move(file, fout)
                                smallMolecule._ParamsFile = pwobj.String(fout)
                            else:
                                smallMolecule._ParamsFile = pwobj.String("Not available")

                        for file in glob.glob(os.path.join(path_params,"*.pdb")):
                            if os.path.isfile(file) and file.endswith('.pdb'):
                                fout = os.path.abspath(os.path.join(path_params, "%s.pdb" % fnRoot))
                                shutil.move(file, fout)
                                smallMolecule._PDBFile = pwobj.String(fout)
                            else:
                                smallMolecule._PDBFile = pwobj.String("Not available")

                    else:
                        smallMolecule._ParamsFile = pwobj.String("Not available")
                        smallMolecule._PDBFile = pwobj.String("Not available")

                    smallMolecule._PDBLigandImage = pwobj.String("Not available")

                    outputSmallMolecules.append(smallMolecule)
                except:
                    pass


            else:
                try:
                    smallMolecule._ConformersFile = pwobj.String("Not generated")

                    # Path to getPath/params/<fnRoot>_conformers
                    path_params = os.path.join(self._getPath(), "params", "%s_withH" % fnRoot)

                    path, dirs, files = next(os.walk(path_params))
                    file_count = len(files)

                    if os.path.exists(path_params) and file_count==3:
                        for file in glob.glob(os.path.join(path_params, "*.params")):
                            if os.path.isfile(file) and file.endswith('.params'):
                                smallMolecule._ParamsFile = pwobj.String(file)
                            else:
                                smallMolecule._ParamsFile = pwobj.String("Not available")

                        for file in glob.glob(os.path.join(path_params, "*.pdb")):
                            if os.path.isfile(file) and file.endswith('.pdb'):
                                smallMolecule._PDBFile = pwobj.String(file)
                                fnOut = self._getExtraPath("%s.png" % fnRoot)
                                args = Pbio.getPluginHome('utils/rdkitUtils.py') + " draw %s %s" % (file, fnOut)
                                try:
                                    Pbio.runRDKit(self, "python3", args)
                                    smallMolecule._PDBLigandImage = pwobj.String(fnOut)
                                except:
                                    smallMolecule._PDBLigandImage = pwobj.String("Not available")

                            else:
                                smallMolecule._PDBFile = pwobj.String("Not available")
                                smallMolecule._PDBLigandImage = pwobj.String("Not available")


                    else:
                        smallMolecule._ParamsFile = pwobj.String("Not available")
                        smallMolecule._PDBFile = pwobj.String("Not available")
                        smallMolecule._PDBLigandImage = pwobj.String("Not available")


                    outputSmallMolecules.append(smallMolecule)
                except:
                    pass


        # Get a list of all the file paths that ends with .pdb in protocol path to remove
        fileList = glob.glob(self._getPath("*.pdb"))
        for filePath in fileList:
            try:
                os.remove(filePath)
            except:
                print("Error while deleting file : ", filePath)


        if outputSmallMolecules is not None:
            self._defineOutputs(outputSmallMols=outputSmallMolecules)
            self._defineSourceRelation(self.inputSmallMols, outputSmallMolecules)



    # --------------------------- UTILS functions ------------------------------
    def _validate(self):
        """ Validate if the inputs are in mol2 or pdb format
        """
        errors = []

        if not self.inputSmallMols.get() is None:
            for mol in self.inputSmallMols.get():
                filename = mol.smallMoleculeFile.get()

                if (filename.endswith(".pdb") or filename.endswith(".mol2")):
                    pass
                else:
                    errors.append("The input %s format it is not correct. "
                                  "Please convert it in pdb or mol2 format with the converter." % filename)


        return errors
