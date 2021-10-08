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
from pyworkflow.protocol import params
from pyworkflow.protocol.params import (LEVEL_ADVANCED, GPU_LIST)
import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from pwchem.objects import SetOfSmallMolecules, SmallMolecule

from operator import itemgetter
import shutil
import os, re
import csv
import glob

from rosetta import Plugin
from rosetta.constants import *
from rosetta.objects import (DarcScore, SetScores)
from ..convert import adt2agdGrid



class Rosetta_darc(EMProtocol):
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

        form.addParam("target_residue", params.IntParam,
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

        form.addHidden(GPU_LIST, params.StringParam,
                       label='Choose GPU ID',
                       default='1',
                       help="GPU may have several cores. Set it one if "
                            "you don't know what we are talking about but you have a GPU."
                            "For DARC, first core index is 1, second 2, and so on. Write 0 if you do not want"
                            "to use GPU")

        group = form.addGroup('Ligands')
        group.addParam("ligands", params.PointerParam, pointerClass="SetOfSmallMolecules", #PDB of the ligands
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
        group = self._defineParamsRays(group)

        form.addSection('Docking')
        form.addParam("shape_only", params.BooleanParam,
                      label="Use only shape of protein and ligand",
                      default=True,
                      important=True,
                      help="")

        form.addParam("grid", params.PointerParam, pointerClass="GridADT",
                      condition="shape_only==False",
                      label="Grid file",
                      important=True,
                      allowsNull=False,
                      help="Select the grid file")

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

 # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        # Insert processing steps
        cId = self._insertFunctionStep('convertInputStep', prerequisites=[])
        gridSteps = []
        if self.fromPockets:
            for pocket in self.inputPockets.get():
                rayDir = self._getExtraPath('pocket_{}'.format(pocket.clone().getObjId()))
                gId = self._insertFunctionStep('generateRaysStep', pocket.clone(), prerequisites=[cId])
                gridSteps.append(gId)
        else:
          gId = self._insertFunctionStep('generateRaysStep', prerequisites=[cId])
          gridSteps.append(gId)

        darcSteps = []
        for mol in self.ligands.get():
            dId = self._insertFunctionStep('darcStep', mol.clone(), prerequisites=gridSteps)
            darcSteps.append(dId)
        else:
            dId = self._insertFunctionStep('darcStep', mol.clone(), prerequisites=gridSteps)
            darcSteps.append(dId)


        self._insertFunctionStep('createOutputStep', prerequisites=darcSteps)

    def convertInputStep(self):
        if not self.shape_only.get():
          adtGridName = self.grid.get().getFileName().split('/')[-1]
          self.agdGrid = adt2agdGrid(self.grid.get(), self._getExtraPath(adtGridName.replace('.e.map', '.agd')))


    def generateRaysStep(self, pocket=None):
        """Generate the txt and pdb file with the protein pocket mapping around a given residue
        """
        if pocket != None:
            rayDir = self._getExtraPath('pocket_{}'.format(pocket.clone().getObjId()))
            os.mkdir(rayDir)
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
            os.mkdir(rayDir)
            args = self.getRaysArgs(outDir=rayDir)
            # Generate 2 file with different formats (pdb (rays are hetatm) and txt).
            # Run Make Ray Files w/wo GPU
            if GPU_LIST == 0:
              Plugin.runRosettaProgram(Plugin.getProgram(MAKE_RAY_FILES), args, cwd=rayDir)
            else:
              args += " -gpu %s" % str(self.gpuList.get())
              Plugin.runRosettaProgram(Plugin.getProgram(MAKE_RAY_FILES_GPU), args, cwd=rayDir)


    def darcStep(self, ligand):
        """ Launch a docking process with Rosetta DARC for each ligand
        """

        # Add protein file where the program will generate the rays (REQUIRED)
        pdb_file = self.getProteinFile()

        # Save compound with errors during docking
        compound_Error = []

        for rayDir in self.getPocketDirs():  # Run DARC for each ligand in the set of small molecules (and his conformers)
            # Create the args of the program and add protein file
            args = ""
            args += " -protein %s" % os.path.abspath(pdb_file)

            # Add ligand file
            ligand_pdb = ligand.getPDBFileName()
            args += " -ligand %s" % os.path.abspath(ligand_pdb)

            # Add params ligand file which path is in the set
            params_file = ligand.getParamsFileName()
            args += " -extra_res_fa %s" % os.path.abspath(params_file)

            # Check if params_file and ligand_pdb is available.
            # If it is not we skip the molecule and start again with another oner
            if ligand_pdb == "Not available" or params_file == "Not available":
                continue

            # Add protein ray file
            ray_file = self.getRayFile(rayDir)
            args += " -ray_file %s" % os.path.abspath(ray_file)

            # Shape only or with electrostatics charges
            if self.shape_only.get():
                args += " -darc_shape_only"
            else:
                #args += " -add_electrostatics"

                args += " -espGrid_file %s" % os.path.abspath(self.agdGrid.getFileName())

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

        outDirs = self.getPocketDirs()
        outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())
        for outDir in outDirs:
            gridId = self.getGridId(outDir)
            pdbFiles = self.getLigandFiles(outDir)
            for mol in self.ligands.get():
                molName = mol.getMolName()
                for pFile in pdbFiles:
                    if molName.split('_')[0] in pFile:
                        newMol = SmallMolecule()
                        newMol.copy(mol)
                        newMol.cleanObjId()
                        newMol.setGridId(gridId)

                        newPDBFile = self._getPath(newMol.getUniqueName() + '_1.pdb')
                        shutil.copy(os.path.join(outDir, pFile), newPDBFile)
                        newMol.poseFile.set(newPDBFile)
                        outputSet.append(newMol)
            self._defineOutputs(outputSmallMolecules=outputSet)


############################## UTILS ########################
    def getRaysArgs(self, outDir, pocket=None):
        # Add protein file where the program will generate the rays (REQUIRED)
        pdb_file = self.getProteinFile()
        pdb_file_extra = os.path.join(outDir, os.path.basename(pdb_file))
        createLink(pdb_file, pdb_file_extra)

        args = " -protein %s" % os.path.basename(pdb_file_extra)

        # Add Database path
        database_path = os.path.join(Plugin.getHome(), ROSETTA_DATABASE_PATH)
        args += " -database %s" % database_path

        # Add the specific residue that will be the center of ray generation (REQUIRED)
        # To use multiple origin points for casting rays or not
        if pocket == None:
            if self.multiple_target.get():
                residues_string = self.target_residues.get()
                res = list(filter(None, re.split(",|;| ", residues_string.upper())))
                target_residues = ",".join(res)
                # res = list(map(int, res))
                args += " -central_relax_pdb_num %s,%s" % (self.target_residue.get(), target_residues)
                args += " -multiple_origin"
            else:
                args += " -central_relax_pdb_num %s" % self.target_residue.get()
        else:
            residues = self.getPocketResidues(pocket)
            args += " -central_relax_pdb_num %s" % (residues)
            args += " -multiple_origin"

        # By default the PocketGrid expands if pocket points are identified near the edge, this flag disables the autoexpanding feature.
        args += " -pocket_static_grid"

        # To use the mass center of a amino acid as origin of throwing tha rays
        if pocket == None and self.change_origin.get():
            args += " -set_origin "
            args += "-origin_res_num %i" % (self.origin_residue.get())  # number of residue in a chain

        # To include electrostatics calculations
        if not self.shape_only.get():
            grid_file = self.agdGrid.getFileName()
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
        for residue in pocket.getMostCentralResidues():
            resStr += self.switchResidueFormat(residue) + ','
        return resStr[:-1]

    def getProteinFile(self):
        if not self.fromPockets:
          pdb_file = self.inputAtomStruct.get().getFileName()
        else:
          pdb_file = self.inputPockets.get().getProteinFile()
        return pdb_file

    def getRayFile(self, rayDir=None):
        if rayDir==None:
            rayDir = self._getExtraPath()
        for file in os.listdir(rayDir):
            if file.startswith('ray_') and file.endswith('.txt'):
                return os.path.join(rayDir, file)
        return None

    def getPocketDirs(self):
      dirs = []
      for lDir in os.listdir(self._getExtraPath()):
          if lDir.startswith('pocket_'):
              dirs.append(self._getExtraPath(lDir))
      return dirs

    def getScoreFiles(self):
        sFiles = []
        for file in os.listdir(self._getPath()):
            if file.startswith('darc_score'):
                sFiles.append(self._getPath(file))
        return sFiles

    def getGridId(self, outDir):
        if self.fromPockets:
            return outDir.split('_')[-1]
        else:
            return '0'

    def getLigandFiles(self, outDir):
        lFiles = []
        for file in os.listdir(outDir):
            if not self.minimize_output.get() and file.startswith('LIGAND_'):
                lFiles.append(file)
            elif self.minimize_output.get() and file.startswith('mini_LIGAND'):
                lFiles.append(file)
        return lFiles


