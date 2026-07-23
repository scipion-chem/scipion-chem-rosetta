# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Judith Maestro Ciria
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

import json

from rosetta.protocols import ProtRosettaFlexDDG
from rosetta.constants import AA_THREE_TO_ONE

from pwem.wizards import EmWizard
import pwem.convert as emconv
from pwchem.wizards import SelectChainWizardQT, SelectMultiChainWizard

SelectChainWizardQT().addTarget(protocol=ProtRosettaFlexDDG,
                              targets=['mutChain'],
                              inputs=['inputAtomStruct'],
                              outputs=['mutChain'])

SelectChainWizardQT().addTarget(protocol=ProtRosettaFlexDDG,
                              targets=['ROIChain'],
                              inputs=['inputAtomStruct'],
                              outputs=['ROIChain'])

SelectMultiChainWizard().addTarget(protocol=ProtRosettaFlexDDG,
                                   targets=['chainsToMove'],
                                   inputs=['inputAtomStruct'],
                                   outputs=['chainsToMove'])


class AddMutationsFlexDDG(EmWizard):
    _targets = [(ProtRosettaFlexDDG, ['addMutation'])]

    def getPositions(self, form):
        protocol = form.protocol
        allRanPos = protocol.RangPositions.get().split(", ")
        return allRanPos

    def getaaTo(self, form):
        protocol = form.protocol
        return "X" if protocol.mutSaturation else str(protocol.mutResidue.get())

    def getchainResidues(self, form):
        protocol = form.protocol
        structureHandler = emconv.AtomicStructHandler()
        structureHandler.read(protocol.inputAtomStruct.get().getFileName())
        structureHandler.getStructure()
        modelsLength, modelsFirstResidue = structureHandler.getModelsChains()
        chainResidues = {}

        for modelID, chains in modelsFirstResidue.items():
            for chainID, residues in chains.items():
                if chainID not in chainResidues:
                    chainResidues[chainID] = {}
                for residue in residues:
                    res_id = residue[0]
                    res_type = residue[1]
                    chainResidues[chainID][res_id] = res_type

        return chainResidues

    def getchain(self, form):
        protocol = form.protocol
        chain = json.loads(protocol.mutChain.get())['chain']
        return chain

    def getROIOrigen(self, form):
        protocol = form.protocol
        ROIOrigin = protocol.ROIOrigin.get()
        return ROIOrigin

    def getSructROI(self, form):
        protocol = form.protocol
        inputStructROI = protocol.inputStructROI.get()
        return inputStructROI

    def getROIChain(self, form):
        protocol = form.protocol
        chainStr = protocol.ROIChain.get()
        if chainStr and chainStr.strip():
            return json.loads(chainStr)['chain']
        return None

    def getMutations(self, form):
        aaTo = self.getaaTo(form)
        chainResidues = self.getchainResidues(form)
        ROIOrigin = self.getROIOrigen(form)
        mutations = []

        if ROIOrigin == 0:
            allRanPos = self.getPositions(form)
            chain = self.getchain(form)
            for ranPos in allRanPos:
                ran = ranPos.split("-")
                for chain, residues_dict in chainResidues.items():
                    if chain == self.getchain(form):
                        for pos in range(int(ran[0]), int(ran[1]) + 1):
                            if pos in residues_dict:
                                aaFrom = AA_THREE_TO_ONE[residues_dict[pos]]
                                mutation = '{}{}{}{}'.format(aaFrom, chain, pos, aaTo)
                                mutations.append(mutation)
        else:
            structROI = self.getSructROI(form)
            roiChain = self.getROIChain(form)
            for item in structROI:
                chain_res = item.getDecodedCResidues()
                for roi in chain_res:
                    res = roi.split("_")
                    chain = res[0]
                    pos = int(res[1])
                    if roiChain and chain != roiChain:
                        continue
                    for ch, residues_dict in chainResidues.items():
                        if ch == chain:
                            if pos in residues_dict:
                                aaFrom = AA_THREE_TO_ONE[residues_dict[pos]]
                                mutation = '{}{}{}{}'.format(aaFrom, chain, pos, aaTo)
                                mutations.append(mutation)

        return mutations

    def show(self, form, *params):
        protocol = form.protocol
        mutations = self.getMutations(form)

        toMutateList = protocol.toMutateList.get()
        toMutateList += "\n" + "\n".join(mutations)
        form.setVar('toMutateList', toMutateList.strip())


class ClearMutationsFlexDDG(EmWizard):
  _targets = [(ProtRosettaFlexDDG, ['clearLabel'])]

  def show(self, form, *params):
    form.setVar('toMutateList', '')
