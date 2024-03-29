# -*- coding: utf-8 -*-
# **************************************************************************
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

"""
This protocol uses a Rosetta suite program (rosetta_scripts) to generate a set of possible atomic structures which fit
an electronic density map.

Adaptation of shared code from Gabriel Lander

The output will be a file named ray_<PDBname>_0001_<TargetResidue>.txt
"""

import os, sys

from pyworkflow.utils import Message
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import toPdb
from pwem.objects import SetOfAtomStructs, AtomStruct

from pwem.emlib.image import ImageHandler
from pwem.convert import Ccp4Header

from rosetta import Plugin
from rosetta.constants import *


class ProtRosettaGenerateStructures(EMProtocol):
    """
    This protocol uses a Rosetta suite program (rosetta_scripts) to generate a set of possible atomic
    structures which fit an electronic density map.

    The output will be a file named ray_<PDBname>_0001_<TargetResidue>.txt
    """
    _label = 'Generate structures'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('General')
        group.addParam("inputStructure", params.PointerParam, pointerClass="AtomStruct",
                      label="Reference atomic structure: ",
                      important=True, allowsNull=False,
                      help="Select the reference atomic structure")
        group.addParam("inputVolume", params.PointerParam, pointerClass="Volume",
                      label="Volume related to the atomic structure: ",
                      help="Select an electron density map for the input structure."
                           "If not provided, the input structure associated volume will be used")
        group.addParam('resolution', params.FloatParam,
                      label='Resolution (A):', default=3.0,
                      help='Set the resolution of the input volume.')
        group.addParam('numMods', params.IntParam,
                       label='Number of output structures:', default=10,
                       help='Set the number of output Rosetta structures to be generated')
        group.addParam('skipIdealize', params.BooleanParam, default=True,
                       label='Skip Rosetta Idealize: ',
                       help='Skip Rosetta idealize to the input atomic structure')

        group = form.addGroup('Symmetry')
        group.addParam('sym', params.StringParam, label='Symmetry type: ', default='C1',
                       help='Type of symmetry to use in the input.\n'
                            'Supported symmetries: C, D, T')
        group.addParam('asu', params.StringParam, label='ASU chains: ',
                       default='A', condition='sym!="C1"',
                       help='Chains in the ASU (comma delimited, i.e. A,C,F)')

        group.addParam('symChains', params.StringParam, label='Symmetric chains: ',
                       condition='sym!="C1"',
                       help='Symmetry-related chains, 2 chains for C sym, 3 for D sym\n'
                            'Comma-separated instance of symmetric chains in the complex\n'
                            'i.e. C2 where ASU is chain A and B; symmetric are C and D respectively (hemoglobin):\n'
                            'ASU chains: A,B\nSymmetric chains: A,C (or B,D)\n\n'
                            'i.e. C3 where ASU is chain A and B; symmetric are C,D and E,F respectively:\n'
                            'ASU chains: A,B\nSymmetric chains: A,C,D (or B,E,F)\n\n')

        group = form.addGroup('Other parameters')
        group.addParam('hydrogen', params.BooleanParam, label='Include hydrogens: ',
                       default=False, help='Include hydrogens in structures')
        group.addParam('membrane', params.BooleanParam, label='Membrane protein: ',
                       default=False, help='Whether the input protein is placed into a membrane')

        form.addParallelSection(threads=4, mpi=1)
        form.addHidden(params.USE_GPU, params.BooleanParam, default=True,
                       label="Use GPU for execution: ",
                       help="This protocol has both CPU and GPU implementation.\
                                                         Select the one you want to use.")

        form.addHidden(params.GPU_LIST, params.StringParam, default='0', label="Choose GPU IDs",
                       help="Add a list of GPU devices that can be used")

 # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('prepareInputStep')
        self._insertFunctionStep('runRosettaScript')
        self._insertFunctionStep('createOutputStep')

    def prepareInputStep(self):
        #Fix volume
        inVol = self._getInputVolume()
        mrcFile = self._getExtraPath('inpVolume.mrc')
        ImageHandler().convert(inVol, mrcFile)
        Ccp4Header.fixFile(mrcFile, mrcFile, inVol.getOrigin(force=True).getShifts(),
                           inVol.getSamplingRate(), Ccp4Header.START)

        #Convert structure to pdb
        pdbFile = self.inputStructure.get().getFileName()
        name, ext = os.path.splitext(pdbFile)
        self.pdbfile = self._getExtraPath(os.path.basename(pdbFile))
        os.link(pdbFile, self.pdbfile)
        if ext != '.pdb':
          outFile = self.pdbfile.replace(ext, '.pdb')
          toPdb(self.pdbfile, outFile)
          self.pdbfile = outFile

        #Clean pdb
        self.pdbfile = self.cleanPDB(self.pdbfile)

        #Run rosetta Idealize
        if not self.skipIdealize.get():
            self.pdbfile = self.runRosettaIdealize()

        if self.isSymmetric():
            self.pdbfile, self.symfile = self.generateRosettaSym()

        # get number of residues in ASU
        self.residues = self.getResList(self.pdbfile)
        print("ASU contains %i residues" % (len(self.residues)))

    def runRosettaScript(self):
      rosettaDir = Plugin.getRosettaDir()
      program, programGPU = 'rosetta_scripts.static.linuxgccrelease', 'rosetta_scripts.opencl.linuxgccrelease'
      xmlRosettaFile = os.path.abspath(self._getExtraPath('multicycle.xml'))

      args = " -database {}/main/database".format(rosettaDir)
      args += " -in::file::s %s " % self.pdbfile
      args += " -parser::protocol {} ".format(xmlRosettaFile)
      if self.isSymmetric():
        args += " -parser::script_vars symmdef=%s " % self.symfile
        args += " -score_symm_complex false "
      if self.nucleic:
        args += " -relax::dna_move "
      args += " -ignore_unrecognized_res "
      args += " -edensity::mapreso %.3f " % self.resolution.get()
      args += " -edensity::cryoem_scatterers "
      args += " -in::file::centroid_input "
      args += " -out::suffix _rev2 "
      args += " -cryst::crystal_refine "
      args += " -restore_talaris_behavior "
      args += " -nstruct %i" % self.numMods.get()

      self.writeRosettaXML(xmlRosettaFile)

      print('Launching Rosetta scripts')
      if not getattr(self, params.USE_GPU):
        program = Plugin.getProgram(program)
        print(program, args)
        print('---------------------------\n')
        sys.stdout.flush()
        Plugin.runRosettaProgram(program, args, cwd=self._getPath())
      else:
        programGPU = Plugin.getProgram(programGPU)
        args += " -gpu %s" % str(getattr(self, params.GPU_LIST).get())
        print(programGPU, args)
        print('---------------------------\n')
        sys.stdout.flush()
        Plugin.runRosettaProgram(programGPU, args, cwd=self._getPath())


    def createOutputStep(self):
        outputSet = SetOfAtomStructs.create(self._getPath())
        for file in os.listdir(self._getPath()):
            if '_rev2_' in file:
                pdbFile = self._getPath(file)
                aStr = AtomStruct(filename=pdbFile)
                outVol = self._getInputVolume().clone()
                outVol.setLocation(self._getExtraPath('inpVolume.mrc'))
                aStr.setVolume(outVol)
                outputSet.append(aStr)

        self._defineOutputs(outputAtomStructs=outputSet)
        self._defineSourceRelation(self.inputStructure, outputSet)


    def _validate(self):
        errors = []

        # Check that the input volume exist
        if self._getInputVolume() is None:
            errors.append("Error: You should provide a volume.\n")
        return errors

###################################### UTILS ####################

    def cleanPDB(self, pdbfile):
        self.nucleic = False
        hydrogen = self.hydrogen.get()

        outf = self._getExtraPath(os.path.splitext(os.path.basename(pdbfile))[0] + "_clean.pdb")
        outpdb = open(outf, 'w')

        chainids = []
        # get rid of anything that is not ATOM, and get rid of hydrogens
        for line in open(pdbfile):
          if line.startswith("ATOM"):
            if hydrogen is not True and line[76:78].strip() == "H": continue
            if line[17:20].strip() in ['A', 'U', 'C', 'G', 'DA', 'DT', 'DG', 'DC']:
              self.nucleic = True
            if line[21] not in chainids: chainids.append(line[21])
            outpdb.write(line)

        outpdb.close()

        pdbfile = os.path.abspath(outf)
        self.chainIds = chainids

        # check symmetry
        sym = self.sym.get().upper()
        asu = self.asu.get().split(',')

        if self.isSymmetric():
          # number of asymmetric units
          try:
            # C symmetry
            num = int(float(sym[1:]))
            # if D sym, multiply by 2
            if sym[0] == 'D':
              num *= 2
          except:
            # Tetrahedral symmetry
            print("setting Tetrahedral symmetry")
            if sym[0] == 'T':
              num = 12
            else:
              exit("\nsymmetry=%s not handled by script\n" % sym)

          print("number of symmetry-related ASUs (sym=%s): %i" % (sym, num))
          print("number of chains: %i" % len(chainids))

          self.symchains = self.symChains.get().split(',')

          # create symmetry file for Rosetta
          if len(asu) == 1:
            outsymf = os.path.splitext(os.path.basename(outf))[0] + "_sym.pdb"
            self.keepChainSubset(outf, outsymf, self.symchains)
            pdbfile = os.path.abspath(outsymf)

        return pdbfile

    def runRosettaIdealize(self):
        program, programGPU = 'idealize_jd2.static.linuxgccrelease', 'idealize_jd2.opencl.linuxgccrelease'
        rosettaDir = Plugin.getRosettaDir()

        # cmd+= " -database $ROSETTA3_DB"
        args = " -database {}/main/database".format(rosettaDir)
        args += " -in::path ./"
        args += " -in::file::s %s" % self.pdbfile
        args += " -ignore_unrecognized_res"
        args += " -no_optH"
        args += " -out::path ./"
        args += " -out::path::pdb ./"
        args += " -chainbreaks "
        args += " -overwrite"
        args += " > idealize.out"

        print("\nRunning Rosetta idealize:")
        print('Idealize', args)
        print("----------------------------")

        if params.GPU_LIST == 0:
            Plugin.runRosettaProgram(program, args, cwd=self._getExtraPath())
        else:
            args += " -gpu %s" % str(self.gpuList.get())
            Plugin.runRosettaProgram(programGPU, args, cwd=self._getExtraPath())

        tmpfile = os.path.splitext(self.pdbfile)[0] + "_0001.pdb"
        outfile = os.path.splitext(self.pdbfile)[0] + "_ideal.pdb"
        if not os.path.isfile(tmpfile):
          exit("\nError: Rosetta idealize did not run\n")

        # rename idealize output with "_ideal.pdb"
        os.system("mv %s %s" % (tmpfile, outfile))

        return outfile

    def generateRosettaSym(self):
      # Generate and edit Rosetta symmetry file
      sym = self.sym.get().upper()
      asu = self.asu.get().split(',')

      symfile = os.path.abspath(self._getExtraPath(sym.lower() + ".symm"))

      print('SymChains: ', self.symchains)
      cmd = "{}/main/source/src/apps/public/symmetry/make_symmdef_file.pl".format(Plugin.getRosettaDir())
      args = " -m NCS"
      args += " -a %s" % (self.symchains[0][0])
      args += " -i %s" % (self.symchains[1][0])
      # if D sym also add C chain for symmetry
      if sym.upper()[0] == 'D':
        args += " %s" % (self.symchains[-1])
      args += " -p %s" % self.pdbfile
      args += " -r 1000"
      args += " > %s" % symfile

      print("\nRunning Rosetta symmetry:")
      print(args)
      print("----------------------------")
      self.runJob(cmd, args)

      pdbfile = os.path.splitext(self.pdbfile)[0] + "_INPUT.pdb"

      if not os.path.isfile(symfile) or not os.path.isfile(pdbfile):
        exit("\nError: Rosetta symmetry file not generated\n")

      symedit = os.path.abspath(self._getExtraPath(sym.lower() + "_edit.symm"))
      fout = open(symedit, 'w')

      for line in open(symfile):
        if line.startswith("set_dof JUMP0_to_com x"):
          fout.write("set_dof JUMP0_to_com x y z\n")
        elif line.startswith("set_dof JUMP0_to_subunit angle_x"):
          fout.write("set_dof JUMP0_to_subunit angle_x angle_y angle_z\n")
        elif line.startswith("set_dof JUMP0_0 x"):
          continue
        else:
          fout.write(line)
      fout.close()

      # if multiple subunits in ASU, don't use pdb from Rosetta
      if len(asu) > 1:
        pdbfile = os.path.splitext(self.pdbfile)[0] + "_asu_INPUT.pdb"
        self.keepChainSubset(self.pdbfile, pdbfile, asu)

      return pdbfile, symedit
    
    def writeRosettaXML(self, xmlRosettaFile):
        # write protocol xml
        xml = open(xmlRosettaFile, 'w')

        # determine electron density weight
        # according to Rosetta:
        # this value should be 50-60 for 3A, 20-25 for 3.5-4A, and 50-60 for 4.5-5A

        # estimating this with a linear fit to these values
        eweight = 70
        if self.resolution.get() > 2.5:
          eweight = 524 * 2.71828 ** (-.808 * self.resolution.get())

        # list of options that need sym flag added
        symflags = ["<cen weights", "<dens_soft", "<dens weights"]

        for line in generateStructuresXML.split('\n'):
          # for add symmetry to xml file
          if self.isSymmetric():
            if line.strip().startswith(tuple(symflags)):
              l = line.strip()[:-1]
              l += " symmetric='1'>"
              line = "\t\t%s\n" % l

            elif line.strip().startswith("<SetupForDensityScoring"):
              line = '\t\t<SetupForSymmetry name="setupsymm" definition=\"%s\"/>\n' % self.symfile

            elif line.strip().startswith("<SwitchResidueTypeSetMover"):
              continue

            elif line.strip().startswith("<MinMover"):
              l = "<SymMinMover name='cenmin' scorefxn='cen' type='lbfgs_armijo_nonmonotone' max_iter='200' tolerance='0.00001' bb='1' chi='1' jump='ALL'/>"
              line = "\t\t%s\n" % l

            elif line.strip().startswith('<Add mover="setupdens"'):
              line = "\t\t<Add mover='setupsymm'/>\n"

          # if more than 1000 residues in ASU, use localrelax instead of fastrelax
          if len(self.residues) > 1000 and line.strip().startswith("<FastRelax"):
            line = "\t\t<LocalRelax name='relaxcart' scorefxn='dens' max_iter='100' ncyc='1' ramp_cart='0' K='16' nexp='2'/>\n"
          if line.strip().endswith("(REPLACE WITH EWEIGHT)"):
            line = "\t\t\t<Reweight scoretype='elec_dens_fast' weight='%i'/>\n" % eweight
          elif line.strip().startswith("<LoadDensityMap"):
            mapFile = os.path.abspath(self._getExtraPath('inpVolume.mrc'))
            line = "\t\t<LoadDensityMap name='loaddens' mapfile=\"%s\"/>\n" % mapFile
          elif line.strip().endswith("MEMBRANE PROTEIN") and not self.membrane.get():
            continue

          xml.write(line.rstrip('\r\n') + "\n")
        xml.close()

    ################### SMALL UTILS ###################
    def _getInputVolume(self):
        if self.inputVolume.get() is None:
            fnVol = self.inputStructure.get().getVolume()
        else:
            fnVol = self.inputVolume.get()
        return fnVol

    def keepChainSubset(self, pdbin, pdbout, chains):
      outpdb = open(pdbout, 'w')
      for line in open(pdbin):
        if line.startswith("ATOM") or line.startswith("HETATM"):
          if line[21] in chains: outpdb.write(line)

      outpdb.close()
      if not os.path.isfile(pdbout):
        exit("\nError generating chain subset file: %s\n" % pdbout)

    def getResList(self, pdb):
      reslist = []
      f = open(pdb)
      for line in f:
        if line.startswith("ATOM") or line.startswith("HETATM"):
          chain = line[21]
          res = line[22:26].strip()
          resid = "%s.%s" % (res, chain)
          if resid not in reslist:
            reslist.append(resid)
      f.close()
      return reslist

    def isSymmetric(self):
        return self.sym.get() != 'C1'
