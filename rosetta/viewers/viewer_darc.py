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


import os, subprocess
import glob

from tkinter import messagebox

from pyworkflow.protocol.params import LabelParam, EnumParam, BooleanParam
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from pwem.viewers import TableView, Chimera, ChimeraView
import pyworkflow.utils as pwutils
from pwchem.viewers import PyMolViewer

from rosetta.protocols.protocol_darc import Rosetta_darc


def errorWindow(tkParent, msg):
    try:
        # if tkRoot is null the error message may be behind
        # other windows
        messagebox.showerror("Error",  # bar title
                              msg,  # message
                              parent=tkParent)
    except:
        print(("Error:", msg))


class DARCViewer(ProtocolViewer):
    """ Viewer for Rosetta program DARC
    """
    _label = 'DARC Viewer'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [Rosetta_darc]
    DARC_SCORE = 'darc_score.sc'


    def __init__(self,  **kwargs):
        ProtocolViewer.__init__(self,  **kwargs)
        # Get the complex files with the best ligand between all conformers given to darc
        self.complexDic, self.complexNames = self.getFilesDic(pattern="DARC*")
        self.ligansDic, self.ligandsNames = self.getFilesDic(pattern="LIGAND*")
        self.ligandsDic_mini, self.ligandsNames_mini = self.getFilesDic(pattern="mini_LIGAND*")
        self.complexDic_mini, self.complexNames_mini = self.getFilesDic(pattern="mini_*", subsPattern="mini_LIGAND*")

    def _defineParams(self, form):
        form.addSection(label='Visualize with PyMol')
        form.addParam('singleLigand', BooleanParam,
                      default=False,
                      label='Display single ligand: ',
                      help='Display the target with a single ligand docked')
        form.addParam('displayPymol', EnumParam, condition='not singleLigand',
                      choices=self.getChoices(pymol=True), default=0,
                      label='Display docking on pocket result: ',
                      help='Docking results are grouped by their pocket, choose the one to visualize')
        form.addParam('displayPymolSingle', EnumParam, condition='singleLigand',
                      choices=self.getChoicesSingle(), default=0,
                      label='Display single ligand: ',
                      help='Display this single ligand with the target')


        form.addSection(label='Visualize with ChimeraX')
        # group = form.addGroup('Overall results')
        view = form.addGroup("Visualization in Chimerax")
        view.addParam('show_complex', LabelParam,
                      label="Visualize complex protein-ligand")

        view.addParam('complex_v', EnumParam,
                      choices= self.complexNames,
                      default=0,
                      label="Choose a complex: ",
                      help="")

        view.addParam('show_ligand', LabelParam,
                      label="Visualize only ligand conformation")

        view.addParam('ligand_v', EnumParam,
                      choices= self.ligandsNames,
                      default=0,
                      label="Choose a ligand: ",
                      help="")

        if self.protocol.minimize_output.get():
            view_mini = form.addGroup("Visualization of complexes after minimization")
            view_mini.addParam('show_complex_mini', LabelParam,
                          label="Visualize complex protein-ligand")

            view_mini.addParam('complex_v_mini', EnumParam,
                          choices= self.complexNames_mini,
                          default=0,
                          label="Choose a complex: ",
                          help="")

            view_mini.addParam('show_ligand_mini', LabelParam,
                          label="Visualize only ligand conformation")

            view_mini.addParam('ligand_v_mini', EnumParam,
                          choices= self.ligandsNames_mini,
                          default=0,
                          label="Choose a ligand: ",
                          help="")

    def _getVisualizeDict(self):
        return {
            'complex_v': self._visualizeComplex,
            'ligand_v': self._visualizeLigand,
            'complex_v_mini': self._visualize_miniComplex,
            'ligand_v_mini': self._visualize_miniLigand,
            'displayPymol': self._viewPymol,
            'displayPymolSingle': self._viewSinglePymol
        }


    def _visualizeComplex(self, e=None):
        """Visualize a complex protein-small molecule in Chimera"""
        basename = self.complexNames[self.complex_v.get()]
        path = os.path.abspath(self.complexDic[basename])
        self.displayChimera(file_path= path)

    def _visualizeLigand(self, e=None):
        """Visualize a complex protein-small molecule in Chimera"""
        basename = self.ligandsNames[self.ligand_v.get()]
        path = os.path.abspath(self.ligansDic[basename])
        self.displayChimera(file_path= path)

    def _visualize_miniComplex(self, e=None):
        """Visualize a complex protein-small molecule in Chimera"""
        basename = self.complexNames_mini[self.complex_v_mini.get()]
        path = os.path.abspath(self.complexDic_mini[basename])
        self.displayChimera(file_path= path)

    def _visualize_miniLigand(self, e=None):
        """Visualize a complex protein-small molecule in Chimera"""
        basename = self.ligandsNames_mini[self.ligand_v_mini.get()]
        path = os.path.abspath(self.ligandsDic_mini[basename])
        self.displayChimera(file_path= path)



    # ------------------- UTILS functions --------------------------

    def displayChimera(self, file_path=None):

        if file_path == None:
            errorWindow(self.getTkRoot(), "Any structure available")
            return

        #pwutils.runJob(None, Chimera.getProgram(), file_path, env=Chimera.getEnviron())
        subprocess.Popen([Chimera.getProgram(), file_path], env=Chimera.getEnviron())

    def getFilesDic(self, pattern, subsPattern=None):
        filesDic, names = {}, []
        for pDir in self.protocol.getPocketDirs():
            gridId = self.protocol.getGridId(pDir)
            files = glob.glob(os.path.join(pDir, pattern))
            #Substract files with the subspattern
            if subsPattern != None:
                subFiles = glob.glob(os.path.join(pDir, subsPattern))
            else:
                subFiles = []
            for f in files:
                if not f in subFiles:
                    names.append('g{}_'.format(gridId) + os.path.basename(f))
                    filesDic[names[-1]] = f
        names.sort()
        return filesDic, names



    def getChoicesSingle(self):
      self.outputLigands = {}
      for oAttr in self.protocol.iterOutputAttributes():
        if 'outputSmallMolecules' in oAttr[0]:
          molSet = getattr(self.protocol, oAttr[0])
          for mol in molSet:
            curMol = mol.clone()
            pName = curMol.getUniqueName()
            self.outputLigands[pName] = curMol

      outputLabels = list(self.outputLigands.keys())
      outputLabels.sort()
      return outputLabels

    def getChoices(self, pymol=True):
      outputLabels = []
      for oAttr in self.protocol.iterOutputAttributes():
        outputLabels.append(oAttr[0])

      outputLabels.sort()
      if pymol and len(outputLabels) > 1:
        outputLabels = ['All'] + outputLabels
      return outputLabels

    def _viewSinglePymol(self, e=None):
      ligandLabel = self.getEnumText('displayPymolSingle')
      mol = self.outputLigands[ligandLabel].clone()
      pmlFile = self.protocol._getExtraPath('pocket_{}'.format(mol.getGridId())) + '/{}.pml'.format(ligandLabel)
      self.writePmlFile(pmlFile, self.buildPMLDockingSingleStr(mol, ligandLabel, addTarget=True))

      pymolV = PyMolViewer(project=self.getProject())
      pymolV.visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

    def _viewPymol(self, e=None):
      if self.getEnumText('displayPymol') != 'All':
        outName = self.getEnumText('displayPymol')
        pmlFile = self.protocol._getPath('{}.pml'.format(outName))
        pmlStr = self.buildPMLDockingsStr(outName)
      else:
        pmlFile = self.protocol._getPath('allOutputMols.pml')
        outName = self.getChoices()[1]
        pmlStr = self.buildPMLDockingsStr(outName)
        for outName in self.getChoices()[2:]:
          pmlStr += self.buildPMLDockingsStr(outName, addTarget=False)

      self.writePmlFile(pmlFile, pmlStr)
      pymolV = PyMolViewer(project=self.getProject())
      pymolV.visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

    def buildPMLDockingSingleStr(self, mol, molName, addTarget=True):
      pmlStr = ''
      if addTarget:
        pmlStr = 'load {}\n'.format(os.path.abspath(self.protocol.getProteinFile()))

      pdbFile = os.path.abspath(mol.getPoseFile())
      pmlStr += 'load {}, {}\nhide spheres, {}\nshow sticks, {}\n'.format(pdbFile, molName, molName, molName)
      return pmlStr

    def buildPMLDockingsStr(self, outName, addTarget=True):
      molSet = getattr(self.protocol, outName)
      pmlStr = ''
      if addTarget:
        pmlStr = 'load {}\n'.format(os.path.abspath(self.protocol.getProteinFile()))
      if not os.path.exists(self.getPDBsDir(outName)):
        molDic = {}
        for i, mol in enumerate(molSet):
          pFile = os.path.abspath(mol.getPoseFile())
          pName = mol.getUniqueName()
          molDic[pName] = pFile

        molKeys = list(molDic.keys())
        molKeys.sort()
        for molName in molKeys:
          pFile = molDic[molName]
          pmlStr += 'load {}, {}\ndisable {}\nhide spheres, {}\nshow sticks, {}\n'.format(
            pFile, molName, molName, molName, molName).format(pFile, molName, molName)

      else:
        molFiles = list(os.listdir(self.getPDBsDir(outName)))
        molFiles.sort()
        for pdbFile in molFiles:
          pFile = os.path.abspath(os.path.join(self.getPDBsDir(outName), pdbFile))
          pName = pdbFile.split('.')[0]
          pmlStr += 'load {}, {}\ndisable {}\n'.format(pFile, pName, pName)
      return pmlStr

    def writePmlFile(self, pmlFile, pmlStr):
      with open(pmlFile, 'w') as f:
        f.write(pmlStr)
        f.write('zoom')
      return pmlFile, pmlStr

    def getPDBsDir(self, outName):
      return self.protocol._getExtraPath(outName)


