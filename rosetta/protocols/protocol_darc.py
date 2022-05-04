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
This protocol uses a Rosetta suite program (named make_ray_files) to generate
a RAY file for the input protein. To generate this ray-file we need to input
the protein in PDB format and specify a target residue (or more than one)
at the interface.

The output will be a file named ray_<PDBname>_0001_<TargetResidue>.txt
"""

from pyworkflow.utils import Message, createLink
from pyworkflow.utils.path import makePath
from pyworkflow.protocol import params
import pyworkflow.object as pwobj
from pyworkflow.protocol.params import (LEVEL_ADVANCED, GPU_LIST)
from pwem.protocols import EMProtocol
from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import splitConformerFile, runOpenBabel, generate_gpf, calculate_centerMass
from pwchem import Plugin as pwchem_plugin

#Soft importing autodock for electrostatic grid generation
try:
    from autodock import Plugin as autodock_plugin
    from autodock.objects import GridADT
    ADTGrid = True
except:
    print('Autodock plugin cannot be imported, so ADT grid cannot be calculated')
    ADTGrid = False

import shutil
import os, re
import glob

from rosetta import Plugin
from rosetta.constants import *
from ..convert import adt2agdGrid
from rosetta.utils.batchParamsToMol_script import getBatchMolToParamsPath


class RosettaProtDARC(EMProtocol):
    """
    This protocol uses a Rosetta suite program (named make_ray_files) to generate
    a RAY file for the input protein. To generate this ray-file we need to input
    the protein in PDB format and specify a target residue (or more than one)
    at the interface.

    The output will be a file named ray_<PDBname>_0001_<TargetResidue>.txt
    """
    _label = 'DARC'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParamsRays(self, form):
        form.addParam("inputAtomStruct", params.PointerParam, pointerClass="AtomStruct",
                    label="Protein structure", condition='not fromPockets',
                    important=True,
                    allowsNull=False,
                    help="Select the atomic structure of the prepared protein file")

        form.addParam("target_residue", params.StringParam,
                      label="Target residue",
                      important=True,  condition='not fromPockets',
                      allowsNull=False,
                      help="Specify the number id of protein residue that you want to "
                           "use to generate the shape of pocket, around of that residue,"
                           " (preferably in the protein surface) using the rays method.")

        form.addParam("multiple_target", params.BooleanParam,
                      label='Multiple target residues',
                      default=False, condition='not fromPockets',
                      important=True,
                      help='Use to define the bound site for the small molecules')

        form.addParam('target_residues', params.StringParam, condition='multiple_target and not fromPockets',
                      label='Additional target residues',
                      help='Write the additional number of protein residues used to throw the rays and create'
                           'the docking surface. The numbers have to be separated by semicolon, comma or space.\n\n'
                           'Examples: \n\n'
                           '1- 99,61 \n\n'
                           '2- 99;61 \n\n'
                           '3- 99 61  \n\n')

        form.addParam("change_origin", params.BooleanParam,
                      label='Change origin for ray generation',
                      default=False, condition='not fromPockets',
                      important=True,
                      help='If you want change the default () origin of rays generation, select "yes"')

        form.addParam("origin_residue", params.IntParam,
                      label="Residue origin",
                      important=True,
                      condition='change_origin and not fromPockets',
                      allowsNull=True,
                      help="Specify the number id of protein residue that you want to"
                           " use as origin points for casting rays. This origin will be"
                           " the center of mass of that residue. ")

        return form


    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Ligands')
        group.addParam("inputLigands", params.PointerParam, pointerClass="SetOfSmallMolecules", #PDB of the ligands
                      label="Ligands to dock",
                      important=True,
                      allowsNull=False,
                      help="Select the ligand or ligand conformers for molecular docking")

        group = form.addGroup('Receptor')
        group.addParam('fromPockets', params.BooleanParam,
                      label="Use pockets to guide the docking",
                      default=True, important=True,
                      help="Whether to center the docking with pockets or with a residue  id")
        group.addParam('inputPockets', params.PointerParam, pointerClass="SetOfPockets",
                       label='Input pockets:', condition='fromPockets',
                       help="The protein pockets to dock in")
        group.addParam('mergeOutput', params.BooleanParam, default=True, expertLevel=LEVEL_ADVANCED,
                       label='Merge outputs from pockets:', condition='fromPockets',
                       help="Merge the outputs from the different pockets")
        group = self._defineParamsRays(group)

        form.addSection('Docking')
        if ADTGrid:
            form.addParam("shape_only", params.BooleanParam,
                          label="Use only shape of protein and ligand",
                          default=True,
                          important=True,
                          help="")
            form.addParam("grid", params.PointerParam, pointerClass="GridADT",
                          condition="not shape_only",
                          label="Grid file", important=True, allowsNull=False,
                          help="Select the grid file")

        else:
            self.shape_only = params.Boolean(True)
        
        form.addParam("search_conformers", params.BooleanParam,
                      label="Search conformers on the fly",
                      default=False,
                      important=True,
                      help="")

        form.addParam("minimize_output", params.BooleanParam,
                      label="Minimize output complex",
                      default=False,
                      important=True,
                      help="")

        # Advanced parameters =======================================
        advanced = form.addGroup("Advanced parameters", expertLevel=LEVEL_ADVANCED)
        advanced.addParam('num_runs', params.IntParam,
                          default=100,
                          label='Runs for PSO:',
                          help='Set the number of runs used during Particle Swarm Optimization (PSO). '
                               'Default value is 100 runs.',
                          expertLevel=LEVEL_ADVANCED)

        advanced.addParam('num_particles', params.IntParam,
                          default=100,
                          label='Particles for PSO:',
                          help='Set the number of particles used during Particle Swarm Optimization (PSO). '
                               'Default value is 100 particles.',
                          expertLevel=LEVEL_ADVANCED)

        advanced.addParam('missing_weight', params.FloatParam,
                          default=5.48,
                          label='Weight for missing points in the pockets:',
                          help='Set the weight of those rays that do not reach the protein pocket'
                               ' or the marked area around an atom or several in the ray generation step',
                          expertLevel=LEVEL_ADVANCED)

        advanced.addParam('steric_weight', params.FloatParam,
                          default=0.61,
                          label='Weight for missing points in the ligand:',
                          help='Set the weight of those rays that do not reach the small molecule',
                          expertLevel=LEVEL_ADVANCED)

        advanced.addParam('extra_weight', params.FloatParam,
                          default=5.47,
                          label='Weight for ligand out of pocket:',
                          help='Set the weight of the small molecule when it is moving away from the pocket',
                          expertLevel=LEVEL_ADVANCED)

        runs = form.addGroup("Runs option",  expertLevel=LEVEL_ADVANCED)
        runs.addParam("cseed", params.BooleanParam,
                       label='Use Constant Seed',
                       default=False,
                       help='Use this option to get reproducible results')

        runs.addParam("seed", params.IntParam,
                       label='Set seed',
                       default=1111111,
                       condition="cseed",
                       help='Set a integer number as constant seed. The default one is 1111111 ')

        form.addParallelSection(threads=4, mpi=1)
        form.addHidden(GPU_LIST, params.StringParam,
                       label='Choose GPU ID',
                       default='1',
                       help="GPU may have several cores. Set it one if "
                            "you don't know what we are talking about but you have a GPU."
                            "For DARC, first core index is 1, second 2, and so on. Write 0 if you do not want"
                            "to use GPU")

 # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        self.originalReceptorFile = self.getOriginalReceptorFile()
        # Insert processing steps
        cId = self._insertFunctionStep('convertInputStep', prerequisites=[])
        raysSteps = []
        if self.fromPockets:
            for pocket in self.inputPockets.get():
                gId = self._insertFunctionStep('generateRaysStep', pocket.clone(), prerequisites=[cId])
                raysSteps.append(gId)
        else:
          gId = self._insertFunctionStep('generateRaysStep', prerequisites=[cId])
          raysSteps.append(gId)

        darcSteps = []
        usedBases = []
        for mol in self.inputLigands.get():
            if not mol.getMolBase() in usedBases:
                if self.fromPockets:
                    for pocket in self.inputPockets.get():
                        dId = self._insertFunctionStep('darcStep', mol.clone(), pocket.clone(), prerequisites=raysSteps)
                        darcSteps.append(dId)
                else:
                    dId = self._insertFunctionStep('darcStep', mol.clone(), prerequisites=raysSteps)
                    darcSteps.append(dId)
                usedBases.append(mol.getMolBase())

        self._insertFunctionStep('createOutputStep', prerequisites=darcSteps)

    def convertInputStep(self):
        #Converting the ADT grid to the Rosetta agd format
        if not self.shape_only:
            adtGridName = self.grid.get().getFileName().split('/')[-1]
            self.agdGrid = adt2agdGrid(self.grid.get(), self._getExtraPath(adtGridName.replace('.e.map', '.agd')))


        # Generate params file that DARC will use to dock the ligand in the target protein
        writtenFiles = []
        with open(self._getExtraPath("molfile_list.txt"), "w+") as file:
            for mol in self.inputLigands.get():
                confFile = mol.getConformersFileName()

                if confFile != None and not confFile in writtenFiles:
                    writtenFiles.append(confFile)
                    if not 'mol2' in confFile and not 'sdf' in confFile:
                        confFile = self.convertConformersFile(confFile, outExt='mol2')
                    else:
                        confFile = os.path.abspath(confFile)
                    file.write(os.path.abspath(confFile) + "\n")

                elif confFile == None:
                    writtenFiles.append(confFile)
                    molFile = mol.getFileName()
                    if not 'mol2' in molFile and not 'sdf' in molFile:
                        confFile = self.convertFile(molFile)
                    else:
                        confFile = os.path.abspath(molFile)
                    file.write(os.path.abspath(confFile) + "\n")

        # 2. Launch batch_molfile_to_params.py for each file. It will generate a pdb file and params file
        database_path = os.path.join(Plugin.getHome(), ROSETTA_DATABASE_PATH)
        args = " -d %s" % database_path
        mol2params_path = os.path.join(Plugin.getHome(), ROSETTA_PARAMS_PATH, PARAMS_FILE)
        args += " --script_path %s" % mol2params_path
        mollist = os.path.abspath(self._getExtraPath("molfile_list.txt"))
        args += " %s" % mollist

        # Execute the program bach_molfile_to_params to create the params file that will be used by DARC programs
        # It creates a directory called params with several directories, each one called as the mol2 file.
        # Inside of each one, we can find:
        #   - 000.params
        #   - 000_conformers.pdb
        #   - log.txtç+
        batchParamsToMol_script = getBatchMolToParamsPath()
        self.runJob('chmod', ' 755 {}'.format(batchParamsToMol_script),
                         cwd=os.path.abspath(self._getExtraPath()))
        Plugin.runRosettaProgram(batchParamsToMol_script, args=args,
                                 cwd=os.path.abspath(self._getExtraPath()))

    def generateRaysStep(self, pocket=None):
        """Generate the txt and pdb file with the protein pocket mapping around a given residue
        """
        if pocket != None:
            rayDir = self._getExtraPath('pocket_{}'.format(pocket.getObjId()))
            makePath(rayDir)
            args = self.getRaysArgs(outDir=rayDir, pocket=pocket)
            # Generate 2 file with different formats (pdb (rays are hetatm) and txt).
            # Run Make Ray Files w/wo GPU
            if GPU_LIST == 0:
                Plugin.runRosettaProgram(Plugin.getProgram(MAKE_RAY_FILES), args, cwd=rayDir)
            else:
                args += " -gpu %s" % str(self.gpuList.get())
                Plugin.runRosettaProgram(Plugin.getProgram(MAKE_RAY_FILES_GPU), args, cwd=rayDir)
        else:
            rayDir = self._getExtraPath('pocket_1')
            makePath(rayDir)
            args = self.getRaysArgs(outDir=rayDir)
            # Generate 2 file with different formats (pdb (rays are hetatm) and txt).
            # Run Make Ray Files w/wo GPU
            if GPU_LIST == 0:
              Plugin.runRosettaProgram(Plugin.getProgram(MAKE_RAY_FILES), args, cwd=rayDir)
            else:
              args += " -gpu %s" % str(self.gpuList.get())
              Plugin.runRosettaProgram(Plugin.getProgram(MAKE_RAY_FILES_GPU), args, cwd=rayDir)


    def darcStep(self, ligand, pocket=None):
        """ Launch a docking process with Rosetta DARC for each ligand
        """

        # Add protein file where the program will generate the rays (REQUIRED)
        pdb_file = self.getOriginalReceptorFile()

        # Save compound with errors during docking
        compound_Error = []

        rayDir = self.getPocketDir(pocket)  # Run DARC for each ligand in the set of small molecules (and his conformers)

        # Create the args of the program and add protein file
        args = ""
        args += " -protein %s" % os.path.abspath(pdb_file)


        paramsDir = self.getParamsDir(ligand)
        # Add ligand file
        ligand_pdb = self.getRosettaConfFile(paramsDir, ligand)
        newLigandPDB = self.changeParamFileCode(ligand_pdb, ligand)
        args += " -ligand %s" % os.path.abspath(newLigandPDB)

        # Add params ligand file which path is in the set
        params_file = glob.glob(os.path.join(paramsDir, "*.params"))[0]
        newLigandParams = self.changeParamFileCode(params_file, ligand)
        args += " -extra_res_fa %s" % os.path.abspath(newLigandParams)

        # Add protein ray file
        ray_file = self.getRayFile(rayDir)
        args += " -ray_file %s" % os.path.abspath(ray_file)

        # Shape only or with electrostatics charges
        if self.shape_only:
            args += " -darc_shape_only"
        else:
            #args += " -add_electrostatics"

            args += " -espGrid_file %s" % os.path.abspath(self.getAGDFile())

        # Search conformers on the fly. DARC 2.0. Optimize conformer during docking
        if self.search_conformers.get():
            args += " -search_conformers True"

        # Minimize the best scoring DARC output model and give some metrics more and
        # Calculate ligand theta value during this minimization
        if self.minimize_output.get():
            args += " -minimize_output_complex True"
            args += " -calculate_thetaLig True"

        # Print Darc output ligand model with protein as PDB file
        args += " -print_output_complex True"

        # Append ligand file name to output files, instead of ligand code
        args += " -use_ligand_filename"

        # Use advanced options
        args += " -num_runs %s" % self.num_runs.get()
        args += " -num_particles %s" % self.num_particles.get()
        args += " -missing_point_weight %s" % self.missing_weight.get()
        args += " -steric_weight %s" % self.steric_weight.get()
        args += " -extra_point_weight %s" % self.extra_weight.get()

        if self.cseed.get():
            args += " -run:constant_seed"
            args += " -run:jran %s" % self.seed.get()

        try:
            # Run DARC w/wo GPU
            if GPU_LIST == 0:
                Plugin.runRosettaProgram(Plugin.getProgram(DARC), args, cwd=rayDir)
            else:
                args += " -gpu %s" % str(self.gpuList.get())
                Plugin.runRosettaProgram(Plugin.getProgram(DARC_GPU), args, cwd=rayDir)
        except:
            compound_Error.append(ligand_pdb)

    def createOutputStep(self):
        """Create a set of darc score for each small molecule and ID"""
        if self.checkSingleOutput():
            outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())

        outDirs = self.getAllPocketDirs()
        for outDir in outDirs:
            savedMols = []
            scoresDic = self.parseScores(outDir)
            gridId = self.getGridId(outDir)
            if not self.checkSingleOutput():
                outputSet = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix=gridId)

            pdbFiles = self.getLigandFiles(outDir)
            for mol in self.inputLigands.get():
                molName, molBase = mol.getUniqueName(), mol.getMolBase()
                for pFile in pdbFiles:
                    if molBase in pFile and not molBase in savedMols:
                        if not self.minimize_output or 'mini_' in pFile:
                          newMol = SmallMolecule()
                          newMol.copy(mol)
                          newMol.cleanObjId()
                          newMol.setGridId(gridId)
                          newMol.setMolClass('Rosetta')
                          newMol.setDockId(self.getObjId())
                          newMol._energy = pwobj.Float(scoresDic[molBase])

                          newPDBFile = self._getPath(newMol.getUniqueName() + '_1.pdb')
                          shutil.copy(os.path.join(outDir, pFile), newPDBFile)
                          newMol.poseFile.set(newPDBFile)
                          outputSet.append(newMol)
                          savedMols.append(molBase)

            if not self.checkSingleOutput():
              outputSet.proteinFile.set(self.getOriginalReceptorFile())
              outputSet.setDocked(True)
              self._defineOutputs(**{'outputSmallMolecules_{}'.format(gridId): outputSet})
              self._defineSourceRelation(self.inputLigands, outputSet)

        if self.checkSingleOutput():
          outputSet.setDocked(True)
          outputSet.proteinFile.set(self.getOriginalReceptorFile())
          self._defineOutputs(outputSmallMolecules=outputSet)
          self._defineSourceRelation(self.inputLigands, outputSet)


############################## UTILS ########################
    def getRaysArgs(self, outDir, pocket=None):
        # Add protein file where the program will generate the rays (REQUIRED)
        pdb_file = self.getOriginalReceptorFile()
        pdb_file_extra = os.path.join(outDir, os.path.basename(pdb_file))
        createLink(pdb_file, pdb_file_extra)

        args = "-protein %s" % os.path.basename(pdb_file)

        # Add Database path
        database_path = os.path.join(Plugin.getHome(), ROSETTA_DATABASE_PATH)
        args += " -database %s" % database_path

        # Add the specific residue that will be the center of ray generation (REQUIRED)
        # To use multiple origin points for casting rays or not
        if pocket == None:
            if self.multiple_target.get():
                residues_string = self.target_residues.get()
                res = list(set(list(filter(None, re.split(",|;| ", residues_string.upper()))) +
                          [self.target_residue.get()]))
                target_residues = ",".join(res)
                # res = list(map(int, res))
                args += " -central_relax_pdb_num %s" % (self.target_residue.get(), target_residues)
                args += " -multiple_origin"
            else:
                args += " -central_relax_pdb_num %s" % self.target_residue.get()
        else:
            residues = self.getPocketResidues(pocket)
            residueNum = re.findall(r'\d+', residues)[0]
            args += " -central_relax_pdb_num %s" % (residues)
            #args += " -multiple_origin"

        # By default the PocketGrid expands if pocket points are identified near the edge, this flag disables the autoexpanding feature.
        args += " -pocket_static_grid"

        # To use the mass center of a amino acid as origin of throwing tha rays
        if pocket == None and self.change_origin.get():
            args += " -set_origin "
            args += "-origin_res_num %i" % (self.origin_residue.get())  # number of residue in a chain

        # To include electrostatics calculations
        if not self.shape_only:
            grid_file = self.getAGDFile()
            args += " -espGrid_file %s" % os.path.abspath(grid_file)

        else:
            # Used to ignore electrostatic score and perform shape only calculation
            args += " -darc_shape_only"

        if self.cseed.get():
            args += " -run:constant_seed"
            args += " -run:jran %s" % self.seed.get()

        return args


    def switchResidueFormat(self, residue):
      '''From A_100 to 100:A'''
      return '{}:{}'.format(residue.split('.')[1], residue.split('.')[0])

    def getPocketResidues(self, pocket):
        resStr = ''
        for residue in pocket.getMostCentralResidues(n=1):
            resStr += self.switchResidueFormat(residue) + ','
        return resStr[:-1]

    def getOriginalReceptorFile(self):
        if not hasattr(self, 'originalReceptorFile'):
            if not self.fromPockets:
                pdb_file = self.inputAtomStruct.get().getFileName()
            else:
                pdb_file = self.inputPockets.get().getProteinFile()
        else:
            pdb_file = self.originalReceptorFile
        return pdb_file

    def getReceptorName(self):
        return self.getOriginalReceptorFile().split('/')[-1].split('.')[0]

    def getRayFile(self, rayDir=None):
        if rayDir==None:
            rayDir = self._getExtraPath()
        for file in os.listdir(rayDir):
            if file.startswith('ray_') and file.endswith('.txt'):
                return os.path.join(rayDir, file)
        return None

    def getPocketDir(self, pocket=None):
      if pocket is None:
          return self._getExtraPath('pocket_1')
      for lDir in os.listdir(self._getExtraPath()):
          if lDir == 'pocket_{}'.format(pocket.getObjId()):
              return self._getExtraPath(lDir)


    def getAllPocketDirs(self):
      dirs = []
      for lDir in os.listdir(self._getExtraPath()):
          if lDir.startswith('pocket_'):
              dirs.append(self._getExtraPath(lDir))
      return dirs

    def getGridId(self, outDir):
        return outDir.split('_')[-1]


    def getLigandFiles(self, outDir):
        lFiles = []
        for file in os.listdir(outDir):
            if not self.minimize_output.get() and file.startswith('LIGAND_'):
                lFiles.append(file)
            elif self.minimize_output.get() and file.startswith('mini_LIGAND'):
                lFiles.append(file)
        return lFiles

    def getParamsDir(self, mol):
        for dir in os.listdir(self._getExtraPath('params')):
            if mol.getMolBase() in dir:
                return self._getExtraPath('params/{}'.format(dir))

    def changeParamFileCode(self, file, ligand):
        if file.endswith('.params'):
            sep = '.'
        elif file.endswith('.pdb'):
            sep = '_'

        ligandDir, ligandFn = os.path.split(file)
        ligandCode = ligandFn.split(sep)[0]

        newLigandFile = os.path.join(ligandDir, ligandFn.replace(ligandCode, ligand.getMolBase()))
        if file != newLigandFile:
            shutil.copy(file, newLigandFile)
        return newLigandFile

    def checkSingleOutput(self):
        return self.mergeOutput.get() or len(self.getAllPocketDirs()) == 1

    def convertFile(self, inFile, outExt='mol2', outDir=None):
        if outDir == None:
            outDir = self._getTmpPath()
        confName, ext = os.path.splitext(inFile)
        oFile = inFile.replace(ext[1:], outExt)
        args = ' -i{} {} -o{} -O {}'.format(ext[1:], inFile, outExt, oFile)
        runOpenBabel(protocol=self, args=args, cwd=outDir)
        return os.path.join(outDir, oFile)

    def convertConformersFile(self, confFile, outExt='mol2'):
        confName, ext = os.path.splitext(confFile)
        outDir = self._getTmpPath(os.path.basename(confName))
        makePath(outDir)
        outDir = splitConformerFile(confFile, outDir)

        formattedFiles = []
        for inFile in os.listdir(outDir):
            oFile = self.convertFile(inFile, outExt=outExt, outDir=outDir)
            formattedFiles.append(oFile)

        outConfFile = os.path.join(outDir, os.path.basename(confFile.replace(ext[1:], outExt)))
        with open(outConfFile, 'w') as f:
            for oFile in formattedFiles:
                f.write(open(oFile).read() + '\n')

        return outConfFile

    def getPDBReceptor(self):
        strFile = self.getOriginalReceptorFile()
        inName, inExt = os.path.splitext(os.path.basename(strFile))
        if inExt != '.pdb':
            pdbFile = os.path.abspath(self._getTmpPath('pdbInput.pdb'))
            if not os.path.exists(pdbFile):
                args = ' -i{} {} -opdb -O {}'.format(inExt[1:], os.path.abspath(strFile), pdbFile)
                runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())
        else:
            pdbFile = strFile
        return pdbFile

    def getPDBQTReceptor(self):
      strFile = self.getOriginalReceptorFile()
      inName, inExt = os.path.splitext(os.path.basename(strFile))
      if inExt != '.pdbqt':
          pdbqtFile = os.path.abspath(self._getTmpPath('{}.pdbqt'.format(self.getReceptorName())))
          if not os.path.exists(pdbqtFile):
              args = ' -v -r %s -o %s' % (strFile, pdbqtFile)
              self.runJob(pwchem_plugin.getMGLPath('bin/pythonsh'),
                          autodock_plugin.getADTPath('Utilities24/prepare_receptor4.py') + args)
      else:
          pdbqtFile = strFile
      return pdbqtFile

    def getAGDFile(self):
        return self.agdGrid.getFileName()

    def parseScores(self, outDir):
        scoresDic = {}
        with open(os.path.join(outDir, 'darc_score.sc')) as fIn:
            for line in fIn:
                if not self.minimize_output:
                    code, score = line.split()[0], line.split()[1]
                else:
                    code, score = line.split()[1], line.split()[2]
                ligCode = code[len(self.getReceptorName())+1:]
                ligCode = '_'.join(ligCode.split('_')[:-1])
                scoresDic[ligCode] = score
        return scoresDic

    def getRosettaConfFile(self, dir, ligand):
        for file in os.listdir(dir):
            if '.pdb' in file and not ligand.getMolBase() in file:
                return os.path.join(dir, file)
