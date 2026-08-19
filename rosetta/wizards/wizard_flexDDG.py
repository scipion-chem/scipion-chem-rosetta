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

from rosetta.protocols import ProtRosettaFlexDDG

from pwchem.wizards import SelectChainWizardQT, SelectMultiChainWizard, AddMutationsWizard, ClearMutationsWizard

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


class AddMutationsFlexDDG(AddMutationsWizard):
    _targets = [(ProtRosettaFlexDDG, ['addMutation'])]


class ClearMutationsFlexDDG(ClearMutationsWizard):
  _targets = [(ProtRosettaFlexDDG, ['clearLabel'])]
