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

#IMPORT
from pyworkflow.protocol import params
from pyworkflow.protocol.params import (LEVEL_ADVANCED, GPU_LIST)
import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from pyworkflow.utils import Message

from operator import itemgetter
import shutil
import os
import csv
import glob

from rosetta import Plugin
from rosetta.constants import *
from rosetta.objects import (DarcScore, SetScores)



class Rosetta_darc(EMProtocol):
    """
    """
    _label = 'DARC'

    # -------------------------- DEFINE param functions ----------------------
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

        form.addParam("protein", params.PointerParam, pointerClass="AtomStruct",
                      label="Protein structure",
                      important=True,
                      allowsNull=False,
                      help="Select the atomic structure of the prepared protein file")

        form.addParam("ligands", params.PointerParam, pointerClass="SetOfSmallMolecules", #PDB of the ligands
                      label="Ligand to dock",
                      important=True,
                      allowsNull=False,
                      help="Select the ligand or ligand conformers for molecular docking")

        form.addParam("ray_file", params.PointerParam, pointerClass="RaysProtein", #txt
                      label="Ray file",
                      important=True,
                      allowsNull=False,
                      help="Select the ligand params file")

        form.addParam("shape_only", params.BooleanParam,
                      label="Use only shape of protein and ligand",
                      default=True,
                      important=True,
                      help="")

        form.addParam("grid", params.PointerParam, pointerClass="GridAGD",
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
                      default=True,
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

 # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('darc')
        self._insertFunctionStep('organize_files')
        self._insertFunctionStep('createOutput')

    def darc(self):
        """
        """

        # Add protein file where the program will generate the rays (REQUIRED)
        pdb_file = self.protein.get().getFileName()
        #pdb_file_extra = os.path.join(self._getExtraPath(), os.path.basename(pdb_file))
        #createLink(pdb_file, pdb_file_extra)

        # Add set of pdb of the small molecules (REQUIRED)
        ligand_set = self.ligands.get()


        for ligand in ligand_set:  # Run DARC for each ligand in the set of small molecules (and his conformers)
            # Create the args of the program and add protein file
            args = ""
            args += " -protein %s" % os.path.abspath(pdb_file)

            # Add ligand file
            #ligand_pdb = self.convertmol2_to_pdb(os.path.abspath(ligand.getFileName()))
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
            ray_file = self.ray_file.get().getFileName()
            args += " -ray_file %s" % os.path.abspath(ray_file)

            # Shape only or with electrostatics charges
            if self.shape_only.get():
                args += " -darc_shape_only"
            else:
                #args += " -add_electrostatics"
                grid_file = self.grid.get().getFileName()
                args += " -espGrid_file %s" % os.path.abspath(grid_file)


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


            # Run DARC w/wo GPU
            if GPU_LIST == 0:
                Plugin.runRosettaProgram(Plugin.getProgram(DARC), args, cwd=self._getExtraPath())
            else:
                args += " -gpu %s" % str(self.gpuList.get())
                Plugin.runRosettaProgram(Plugin.getProgram(DARC_GPU), args, cwd=self._getExtraPath())


    def organize_files(self):
        # DARC build multiple files and their organization is better

        # Move the darc_score.sc to main path of the Protocol
        score_extra = self._getExtraPath("darc_score.sc")
        scores = self._getPath("darc_score.sc")

        # Create a ligand directory and Ligand-protein complex (ligand selected from all conformers)
        os.mkdir(os.path.join(self._getExtraPath(), "ligands"))
        os.mkdir(os.path.join(self._getExtraPath(), "darc_complex"))
        if self.minimize_output.get():
            os.mkdir(os.path.join(self._getExtraPath(), "minimization"))
            os.mkdir(os.path.join(self._getExtraPath(), "minimization", "ligands"))
            os.mkdir(os.path.join(self._getExtraPath(), "minimization", "darc_complex"))

        # Sort by DARC score
        reader = csv.reader(open(score_extra), delimiter="\t")
        score_list =[]
        if self.minimize_output.get():
            header="TAG\tDARC_score\tTotal_Energy\tInterface_Energy\tInterface_HB\tTotal_packstats\tInterface_unsat\tthetaLig\n"
        else:
            header = "TAG\tDARC_score\n"

        with open(scores, "w+") as sc:
            sc.write(header)
            add = []
            if self.minimize_output.get():
                for row in reader:
                    temp = []
                    for x in row[1:]:
                        temp.extend(x.split())
                    add.append(temp)

                for row in sorted(add, key=itemgetter(1)):
                    sc.write('\t'.join(str(x) for x in row))
                    sc.write('\n')
                    score_list.append(row)
            else:
                for row in sorted(reader, key=itemgetter(1)):
                    sc.write('\t'.join(str(x) for x in row))
                    sc.write('\n')


        os.remove(score_extra)

        # Move files to specific directory
        for file in glob.glob(self._getExtraPath("*")):
            basename = os.path.basename(file)
            if basename.startswith("DARC") and os.path.isfile(file):
                shutil.move(file, os.path.join(self._getExtraPath(),"darc_complex",basename))
            elif basename.startswith("LIGAND") and os.path.isfile(file):
                shutil.move(file, os.path.join(self._getExtraPath(), "ligands", basename))

            if self.minimize_output.get():
                if basename.startswith("mini_LIGAND") and os.path.isfile(file):
                    shutil.move(file, os.path.join(self._getExtraPath(), "minimization", "ligands", basename))
                elif basename.startswith("mini_") and os.path.isfile(file):
                    shutil.move(file, os.path.join(self._getExtraPath(), "minimization", "darc_complex", basename))


    def createOutput(self):
        """Create a set of darc score for each small molecule and ZINC ID"""

        scores = self._getPath("darc_score.sc")
        outputDarcScore = SetScores().create(path=self._getPath(), suffix='')

        with open(scores) as sc:
            lines = sc.readlines()
            for line in lines[1:]:
                sc = line.split("\t")
                zincID = sc[0]
                scoreDarc = sc[1]

                dscore = DarcScore(zincID=zincID, scoreDarc=scoreDarc)

                if self.minimize_output.get():
                    dscore.Total_Energy = pwobj.Float(sc[2])
                    dscore.Interface_Energy = pwobj.Float(sc[3])
                    dscore.Interface_Hb = pwobj.Float(sc[4])
                    dscore.Total_Pack = pwobj.Float(sc[5])
                    dscore.Interface_Unsat = pwobj.Float(sc[6])
                    dscore.Thetalig = pwobj.Float(sc[7])

                outputDarcScore.append(dscore)

        if outputDarcScore is not None:
            print(outputDarcScore.getSize())
            self._defineOutputs(DarcScores=outputDarcScore)
            self._defineSourceRelation(self.protein, outputDarcScore)






    def convertmol2_to_pdb(self, filename):
        # Use babel to do this
        fnMol = os.path.split(filename)[1]  # Name of complete file
        fnRoot = os.path.splitext(fnMol)[0]  # Name without format
        fnFormat = os.path.splitext(fnMol)[1]  # Format file

        if fnFormat == ".mol2":
            # Convert mol2 format to pdb
            args = " -imol2 %s -O %s.pdb" % (os.path.abspath(filename), fnRoot)
            Plugin.runOPENBABEL(self, args=args, cwd=os.path.abspath(self._getTmpPath()))

        abs_path = os.path.abspath(self._getTmpPath("%s.pdb" % fnRoot))

        return abs_path

