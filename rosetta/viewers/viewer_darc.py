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


import glob

from tkinter import messagebox

from pyworkflow.protocol.params import LabelParam, EnumParam
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from pwem.viewers import ChimeraView
from pwchem.utils.utilsViewer import *
from pwchem.viewers import SmallMoleculesViewer

from rosetta.protocols.protocol_darc import RosettaProtDARC

def errorWindow(tkParent, msg):
    try:
        # if tkRoot is null the error message may be behind
        # other windows
        messagebox.showerror("Error",  # bar title
                              msg,  # message
                              parent=tkParent)
    except:
        print(("Error:", msg))


class DARCViewer(SmallMoleculesViewer):
    """ Viewer for Rosetta program DARC
    """
    _label = 'Viewer DARC docking'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [RosettaProtDARC]

    def __init__(self,  **kwargs):
        super().__init__(**kwargs)
        # Get the complex files with the best ligand between all conformers given to darc (CHIMERA)
        self.complexDic, self.complexNames = self.getFilesDic(pattern="DARC*")
        self.ligansDic, self.ligandsNames = self.getFilesDic(pattern="LIGAND*")
        self.ligandsDic_mini, self.ligandsNames_mini = self.getFilesDic(pattern="mini_LIGAND*")
        self.complexDic_mini, self.complexNames_mini = self.getFilesDic(pattern="mini_*", subsPattern="mini_LIGAND*")

    def _defineParams(self, form):
        super()._defineParams(form)
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
        viewDic = super()._getVisualizeDict()
        darcDic = {
            'complex_v': self._visualizeComplex,
            'ligand_v': self._visualizeLigand,
            'complex_v_mini': self._visualize_miniComplex,
            'ligand_v_mini': self._visualize_miniLigand,
        }
        viewDic.update(darcDic)
        return viewDic


    def _visualizeComplex(self, e=None):
        """Visualize a complex protein-small molecule in Chimera"""
        basename = self.complexNames[self.complex_v.get()]
        path = os.path.abspath(self.complexDic[basename])
        return self.displayChimera(file_path= path)

    def _visualizeLigand(self, e=None):
        """Visualize a complex protein-small molecule in Chimera"""
        basename = self.ligandsNames[self.ligand_v.get()]
        path = os.path.abspath(self.ligansDic[basename])
        return self.displayChimera(file_path= path)

    def _visualize_miniComplex(self, e=None):
        """Visualize a complex protein-small molecule in Chimera"""
        basename = self.complexNames_mini[self.complex_v_mini.get()]
        path = os.path.abspath(self.complexDic_mini[basename])
        return self.displayChimera(file_path= path)

    def _visualize_miniLigand(self, e=None):
        """Visualize a complex protein-small molecule in Chimera"""
        basename = self.ligandsNames_mini[self.ligand_v_mini.get()]
        path = os.path.abspath(self.ligandsDic_mini[basename])
        return self.displayChimera(file_path= path)



    # ------------------- UTILS functions --------------------------

    def displayChimera(self, file_path=None):

        if file_path == None:
            errorWindow(self.getTkRoot(), "Any structure available")
            return

        view = ChimeraView(file_path)
        return [view]

    def getFilesDic(self, pattern, subsPattern=None):
        filesDic, names = {}, []
        for pDir in self.protocol.getAllPocketDirs():
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
