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

"""
Wrapper around the Flex ddG protocol (Barlow et al. 2018, J. Phys. Chem. B) from
https://github.com/Kortemme-Lab/flex_ddG_tutorial

Computes the change in binding free energy (interface ddG) between two protein
chains upon a point mutation, using Rosetta's backrub ensemble-based protocol.
"""

import os
import re
import json
import math
import sqlite3
import statistics

from pyworkflow.constants import BETA
from pyworkflow.object import Object, Float, String
import pyworkflow.protocol.params as params
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pyworkflow.utils import Message

from pwem.protocols import EMProtocol
from pwem.objects.data import SetOfStats
import pwem.convert as emconv

from pwchem.utils import cleanPDB

from rosetta import Plugin
from rosetta.constants import *


class ProtRosettaFlexDDG(EMProtocol):
    """
    This protocol computes the change in binding free energy (interface ddG) between two
    chains of a protein complex upon a point mutation, using the Flex ddG protocol
    (backrub ensemble sampling + talaris2014 GAM-reweighted scoring) as implemented in Rosetta.
    """
    _label = 'Flex ddG'
    _devStatus = BETA
    stepsExecutionMode = STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _addMutationForm(self, form):
        form.addParam('multiPosition', params.BooleanParam, default=False,
                      label='Use a set of ROIs.',
                      help='Mutate and calculate the change in binding free energy (ΔΔG) '
                           'over a set of Regions Of Interest (ROIs).\nTo calculate ΔΔG '
                           'at specific positions, select "No" and directly specify the '
                           'mutations in "List of mutations".')

        form.addParam('ROIOrigin', params.EnumParam, default=0, condition='multiPosition',
                       label='Source of ROIs: ', choices=['Manual', 'SetOfStructROIs'],
                       help='Select the source of the regions of interest.')

        form.addParam('mutChain', params.StringParam, allowsNull=False,
                      label='Chain to mutate', condition='ROIOrigin==0 and multiPosition',
                      help='Specify the protein chain to mutate.')

        form.addParam('RangPositions', params.StringParam, allowsNull=False,
                      label='Range of positions: ', condition='ROIOrigin==0 and multiPosition',
                      help='Specify the first and last index of each position range, separating '
                           'each range with a comma, i.e., "[FIRST_1]-[LAST_1], [FIRST_2]-[LAST_2]". '
                           'For example, "1-30, 50-70" will select for mutation all residues between '
                           'positions 1 and 30, and between positions 50 and 70 in the corresponding '
                           'chain.')

        form.addParam('inputStructROI', params.PointerParam, pointerClass="SetOfStructROIs",
                      label='Input structural ROI', condition='ROIOrigin==1 and multiPosition',
                      allowsNull=False, help='Select the source of the ROIs.')

        form.addParam('ROIChain', params.StringParam, default='', allowsNull=True,
                      label='Chain to filter (optional)', condition='ROIOrigin==1 and multiPosition',
                      help='Restrict the mutations generated from the ROIs to this chain only. '
                           'If left empty, mutations for every chain present in the ROIs will be added.')

        form.addParam('mutSaturation', params.BooleanParam, default=True,
                       label='Saturation mutagenesis', condition='multiPosition',
                       help='Perform saturation mutagenesis, that is, replace each position '
                            'with each of the 20 protein-forming aminoacids (ACDEFGHIKLMNPQRSTVWY).')

        form.addParam('mutResidue', params.StringParam, allowsNull=False,
                      label="Residue to introduce", condition='multiPosition and not mutSaturation',
                      help='Define the substitute residue which will be introduced with its '
                           'one-letter code.')

        form.addParam('addMutation', params.LabelParam,
                      label='Add defined mutations', condition='multiPosition',
                      help='Add the defined mutations to the list of mutations below.')

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Structure template')
        group.addParam('inputAtomStruct', params.PointerParam, pointerClass="AtomStruct",
                      label='Input atomic structure', allowsNull=False,
                      help='The atomic structure must contain the protein-protein (or protein-nanobody) '
                           'complex, with the two interacting partners in different chains.')

        group.addParam('chainsToMove', params.StringParam, allowsNull=False,
                      label='Chain(s) defining one side of the interface',
                      help='Chain(s) that are moved away from the rest of the complex to compute the '
                           'unbound state energy, i.e. one side of the protein-protein interface. '
                           'A single chain (e.g. "B") or several chains making up one side of the '
                           'interface (e.g. "L,H") can be specified.\nIt is recommended that this is the '
                           'chain where the mutation(s) are introduced.')

        group = form.addGroup('Define mutation')
        self._addMutationForm(group)
        group.addParam('toMutateList', params.TextParam, width=70,
                      default='', label='List of mutations:',
                      help='The syntax of a mutation is "[aaFrom][Chain][Position][aaTo]". For example, '
                           'the mutation "CA182Y", mutates position 182 of chain A that is a (C)ystein to '
                           'a t(Y)rosine.\nTo perform saturation mutagenesis (the amino acid is replaced '
                           'by each of the 20 protein-forming aminoacids), in the [aaTo] parameter specify '
                           '"X". For example, CA182X mutates Cys182 of the chain A to all protein-forming '
                           'aminoacids (ACDEFGHIKLMNPQRSTVWY).\nEach mutation is scored independently.')
        group.addParam('clearLabel', params.LabelParam,
                       label='Clear mutation list',
                       help='Clear mutations list')

        section = form.addSection(label='Flex ddG parameters')
        section.addParam('nstruct', params.IntParam, default=35,
                       label='Number of backrub trajectories (nstruct)',
                       help='Number of independent backrub trajectories generated per mutation. '
                            'ddG is averaged over these replicates. The Rosetta-benchmarked value is 35. '
                            'Lower this (e.g. 3) only for a quick functional test, as the ddG estimate '
                            'will be much noisier.')
        section.addParam('numberBackrubTrials', params.IntParam, default=35000,
                       label='Number of backrub trials',
                       help='Number of backrub sampling steps per trajectory. The benchmarked value is '
                            '35000. Lower this (e.g. 1000) only for a quick functional test.')
        section.addParam('backrubTrajectoryStride', params.IntParam, default=35000,
                       label='Backrub trajectory stride',
                       help='After every N backrub steps, the ddG calculation is performed. Leave this '
                            'equal to the number of backrub trials for the fastest run (a single ddG '
                            'checkpoint at the end of the trajectory), which is what this protocol assumes '
                            'when parsing results.')
        section.addParam('maxMinimizationIter', params.IntParam, default=5000,
                       label='Max. minimization iterations',
                       help='Maximum number of minimization gradient descent steps. The benchmarked value '
                            'is 5000.')
        section.addParam('absScoreConvergenceThresh', params.FloatParam, default=1.0,
                       label='Abs. score convergence threshold',
                       help='Maximum allowed change in total score after minimization. If exceeded, '
                            'another minimization cycle is run. The benchmarked value is 1.0.')

        form.addParallelSection(threads=16, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        convertId = self._insertFunctionStep(self.convertInputStep, prerequisites=[])

        # Flex ddG needs one independent rosetta_scripts call per (mutation, replicate) job, but
        # exposing one Scipion step per job does not scale (a saturation scan can easily reach
        # tens of thousands of jobs). Instead, the jobs are split into as many chunks as available
        # threads, and each chunk is processed sequentially by a single step, so the number of
        # Scipion steps stays around "number of threads" regardless of how many jobs there are.
        nChunks = max(1, self.numberOfThreads.get())
        chunkStepIds = []
        for it in range(nChunks):
            stepId = self._insertFunctionStep(self.runFlexDDGChunkStep, it, nChunks,
                                              prerequisites=[convertId])
            chunkStepIds.append(stepId)

        self._insertFunctionStep(self.createOutputStep, prerequisites=chunkStepIds)

    def convertInputStep(self):
        fnPDB = self._getExtraPath('complex.pdb')
        cleanPDB(self.inputAtomStruct.get().getFileName(), fnPDB)

        xmlFile = self._getExtraPath(FLEXDDG_XML_FILE)
        with open(xmlFile, 'w') as f:
            f.write(flexDDGXML)

    def runFlexDDGChunkStep(self, it, nChunks):
        mutations = self._getMutationsList()
        for mutIdx, structIdx in self._getJobChunk(it, nChunks):
            mutation = mutations[mutIdx]
            try:
                self._runSingleFlexDDG(mutation, structIdx)
            except Exception as e:
                aaFrom, chain, position, aaTo = mutation
                print('Flex ddG failed for mutation %s%s%d%s, replicate %d: %s'
                     % (aaFrom, chain, position, aaTo, structIdx, e))

    def _runSingleFlexDDG(self, mutation, structIdx):
        mutDir = self._getMutationDir(mutation, structIdx)
        os.makedirs(mutDir, exist_ok=True)

        resfile = os.path.join(mutDir, FLEXDDG_RESFILE)
        self._writeResfile(resfile, mutation)

        args = '-s %s' % os.path.abspath(self._getExtraPath('complex.pdb'))
        args += ' -parser:protocol %s' % os.path.abspath(self._getExtraPath(FLEXDDG_XML_FILE))
        args += ' -parser:script_vars chainstomove=%s' % self._getChainsToMove()
        args += ' mutate_resfile_relpath=%s' % os.path.abspath(resfile)
        args += ' number_backrub_trials=%d' % self.numberBackrubTrials.get()
        args += ' max_minimization_iter=%d' % self.maxMinimizationIter.get()
        args += ' abs_score_convergence_thresh=%.1f' % self.absScoreConvergenceThresh.get()
        args += ' backrub_trajectory_stride=%d' % self.backrubTrajectoryStride.get()
        args += ' -restore_talaris_behavior'
        args += ' -in:file:fullatom'
        args += ' -ignore_unrecognized_res'
        args += ' -ignore_zero_occupancy false'
        args += ' -ex1 -ex2'
        args += ' -out:path:all %s' % os.path.abspath(mutDir)
        args += ' -nstruct 1'

        program = Plugin.getProgram(ROSETTA_SCRIPTS)
        Plugin.runRosettaProgram(program, args, cwd=mutDir)

    def createOutputStep(self):
        mutations = self._getMutationsList()
        outputSet = SetOfStats.create(self.getPath())

        for mutation in mutations:
            aaFrom, chain, position, aaTo = mutation
            mutName = '%s%s%d%s' % (aaFrom, chain, position, aaTo)

            ddgValues, gamValues = [], []
            for structIdx in range(1, self.nstruct.get() + 1):
                mutDir = self._getMutationDir(mutation, structIdx)
                db3File = os.path.join(mutDir, FLEXDDG_DB_FILE)
                if not os.path.isfile(db3File):
                    continue
                try:
                    result = self._extractDdG(db3File)
                except Exception as e:
                    print('Could not extract ddG from %s: %s' % (db3File, e))
                    continue
                ddgValues.append(result['ddG_raw'])
                gamValues.append(result['ddG_gam'])

            if not ddgValues:
                print('No successful Flex ddG runs found for mutation %s' % mutName)
                continue

            item = Object()
            item.setObjLabel(label=mutName)
            item.mutation = String(mutName)
            item.ddg = Float(statistics.mean(gamValues))
            item.ddgRaw = Float(statistics.mean(ddgValues))
            item.ddgStd = Float(statistics.stdev(gamValues) if len(gamValues) > 1 else 0.0)
            item.nstruct = Float(len(ddgValues))
            outputSet.append(item)

        self._defineOutputs(outputStats=outputSet)
        self._defineTransformRelation(self.inputAtomStruct, outputSet)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = []

        structureHandler = emconv.AtomicStructHandler()
        structureHandler.read(self.inputAtomStruct.get().getFileName())
        structureHandler.getStructure()
        modelsLength, modelsFirstResidue = structureHandler.getModelsChains()

        validChains = set()
        chainResidues = {}
        for modelID, chains in modelsFirstResidue.items():
            for chainID, residues in chains.items():
                filteredResidues = [res for res in residues if res[1] != 'HOH']
                validChains.add(chainID)
                if chainID not in chainResidues:
                    chainResidues[chainID] = filteredResidues

        if not self.chainsToMove.get() or not self.chainsToMove.get().strip():
            errors.append('You must specify the chain(s) that define one side of the interface '
                          '("Chain(s) defining one side of the interface").')
        else:
            for ch in self._getChainsToMove().split(','):
                if ch.strip() not in validChains:
                    errors.append('Chain "%s" (from "Chain(s) defining one side of the interface") is '
                                  'not present in the PDB file. The PDB file contains the following '
                                  'chains: %s.' % (ch.strip(), ", ".join(sorted(validChains))))

        if not self.toMutateList.get().strip():
            errors.append('You have not added any mutation to the list. Do so using the "Add defined '
                          'mutations" wizard once you have defined it.')
        else:
            for line in self.toMutateList.get().strip().split('\n'):
                line = line.strip()
                if not line:
                    continue
                match = self._mutationPattern().match(line)
                if not match:
                    errors.append('The mutation "%s" does not have the 4 necessary parameters. '
                                  'Mutation format must be "[aaFrom][Chain][Position][aaTo]".' % line)
                    continue

                aaFrom, chain, position, aaTo = match.groups()
                if chain not in validChains:
                    errors.append('The chain "%s" of the mutation "%s" is not present in the PDB file. '
                                  'The PDB file contains the following chains: %s.'
                                  % (chain, line, ", ".join(validChains)))
                elif not position.isdigit():
                    errors.append('The position of the mutation "%s" must be an integer.' % line)
                elif aaFrom not in AA_THREE_TO_ONE.values():
                    errors.append('The wild-type aminoacid of the mutation "%s" does not exist or is '
                                  'not written with its one-letter code.' % line)
                elif aaTo != 'X' and aaTo not in AA_THREE_TO_ONE.values():
                    errors.append('The mutant aminoacid of the mutation "%s" does not exist or is not '
                                  'written with its one-letter code.' % line)
                else:
                    intPosition = int(position)
                    residuesDict = {res[0]: res[1] for res in chainResidues.get(chain, [])}
                    if intPosition not in residuesDict:
                        if residuesDict:
                            firstResidue = next(iter(residuesDict))
                            lastResidue = list(residuesDict)[-1]
                            errors.append('Position "%d" in chain "%s" for mutation "%s" is out of range. '
                                          'The chain "%s" has positions from %s to %s.'
                                          % (intPosition, chain, line, chain, firstResidue, lastResidue))
                    elif AA_THREE_TO_ONE[residuesDict[intPosition]] != aaFrom:
                        errors.append('The wild-type aminoacid "%s" at position "%d" in chain "%s" for '
                                      'mutation "%s" does not match the PDB file. The aminoacid at that '
                                      'position is %s (%s).'
                                      % (aaFrom, intPosition, chain, line, residuesDict[intPosition],
                                         AA_THREE_TO_ONE[residuesDict[intPosition]]))

        program = Plugin.getProgram(ROSETTA_SCRIPTS)
        if not os.path.exists(os.path.expanduser(program)):
            errors.append('Cannot find Rosetta rosetta_scripts binary: %s' % program)

        return errors

    def _summary(self):
        summary = []
        if self.isFinished():
            summary.append('Flex ddG was run for %d mutation(s), %d backrub trajector%s each.'
                           % (len(self._getMutationsList()), self.nstruct.get(),
                              'y' if self.nstruct.get() == 1 else 'ies'))
        return summary

    def _methods(self):
        methods = []
        methods.append('Prediction of the binding free energy change (ΔΔG) for protein-protein '
                       'interactions due to point mutation(s), using the Flex ddG protocol (backrub '
                       'ensemble sampling) implemented in Rosetta.\nThe reported ddg is the GAM-reweighted '
                       'talaris2014 score, as recommended by Barlow et al. 2018; ddgRaw is the unweighted '
                       'talaris2014 total_score ddG.')
        return methods

    def _citations(self):
        return ['Barlow2018', 'LeaverFay2011']

    # --------------------------- UTILS functions ------------------------------
    def _getJobChunk(self, it, nChunks):
        """ Deterministically returns the (mutIdx, structIdx) jobs assigned to chunk "it" out of
        "nChunks" total chunks, striping the full job list so chunks stay evenly sized. """
        mutations = self._getMutationsList()
        jobs = [(mutIdx, structIdx) for mutIdx in range(len(mutations))
               for structIdx in range(1, self.nstruct.get() + 1)]
        return jobs[it::nChunks]

    def _getChainsToMove(self):
        """ Returns a comma-separated string of chain IDs (e.g. "B" or "L,H") suitable for
        Rosetta's InterfaceDdGMover chain_name argument.
        The "Chain(s) defining one side of the interface" field can be filled either by typing
        chain IDs directly (e.g. "B" or "L,H") or through the chain-selection wizard, which stores
        a JSON dictionary such as {"chain": "B", ...} or {"model-chain": "0-B, 0-H"}. """
        value = self.chainsToMove.get().strip()
        try:
            chainJson = json.loads(value)
        except (ValueError, TypeError):
            return value

        if 'chain' in chainJson:
            return chainJson['chain'].upper().strip()
        elif 'model-chain' in chainJson:
            modelChains = chainJson['model-chain'].upper().strip()
            chains = [x.split('-')[1] for x in modelChains.split(',')]
            return ','.join(c.strip() for c in chains)
        return value

    def _mutationPattern(self):
        return re.compile(r'([A-Za-z]+)([A-Za-z]+)([^a-zA-Z]+)([A-Za-z]+)')

    def _getMutationsList(self):
        """ Parses self.toMutateList into a deduplicated list of (aaFrom, chain, position, aaTo)
        tuples, expanding saturation ('X') mutations into the 20 canonical aminoacids. """
        mutations = []
        seen = set()
        for line in self.toMutateList.get().strip().split('\n'):
            line = line.strip()
            if not line:
                continue
            match = self._mutationPattern().match(line)
            if not match:
                continue
            aaFrom, chain, position, aaTo = match.groups()
            if not position.isdigit():
                continue
            position = int(position)

            aaTos = CANONICAL_AAS if aaTo == 'X' else [aaTo]
            for a in aaTos:
                key = (aaFrom, chain, position, a)
                if key not in seen:
                    seen.add(key)
                    mutations.append(key)
        return mutations

    def _getMutationDir(self, mutation, structIdx):
        aaFrom, chain, position, aaTo = mutation
        mutName = '%s%s%d%s' % (aaFrom, chain, position, aaTo)
        return self._getExtraPath(mutName, '%02d' % structIdx)

    def _writeResfile(self, resfile, mutation):
        aaFrom, chain, position, aaTo = mutation
        with open(resfile, 'w') as f:
            f.write('NATAA\nstart\n%d %s PIKAA %s\n' % (position, chain, aaTo))

    def _extractDdG(self, db3File):
        """ Extracts the interface ddG from a Flex ddG ddG.db3 sqlite database, assuming a single
        ddG checkpoint (backrub_trajectory_stride == number_backrub_trials). Returns the raw
        talaris2014 ddG and the GAM-reweighted ddG (Barlow et al. 2018). """
        conn = sqlite3.connect(db3File)
        c = conn.cursor()
        rows = c.execute('''
            SELECT batches.name, structure_scores.struct_id, score_types.score_type_name,
                   structure_scores.score_value
            FROM structure_scores
            INNER JOIN batches ON batches.batch_id = structure_scores.batch_id
            INNER JOIN score_types ON score_types.batch_id = structure_scores.batch_id
                                  AND score_types.score_type_id = structure_scores.score_type_id
        ''').fetchall()
        conn.close()

        maxStruct = {}
        for batchName, structId, term, value in rows:
            state = batchName[:-9] if batchName.endswith('_dbreport') else batchName
            maxStruct[state] = max(maxStruct.get(state, 0), structId)

        stateScores = {}
        for batchName, structId, term, value in rows:
            state = batchName[:-9] if batchName.endswith('_dbreport') else batchName
            if structId == maxStruct[state]:
                stateScores.setdefault(state, {})[term] = value

        required = ['bound_wt', 'unbound_wt', 'bound_mut', 'unbound_mut']
        for state in required:
            if state not in stateScores:
                raise RuntimeError('Missing state "%s" in %s' % (state, db3File))

        commonTerms = set(stateScores['bound_wt']) & set(stateScores['unbound_wt']) & \
                     set(stateScores['bound_mut']) & set(stateScores['unbound_mut'])

        ddgRaw = ((stateScores['bound_mut']['total_score'] - stateScores['unbound_mut']['total_score']) -
                 (stateScores['bound_wt']['total_score'] - stateScores['unbound_wt']['total_score']))

        gamTotal = 0.0
        for term in FLEXDDG_GAM_PARAMS:
            if term not in commonTerms:
                continue
            termDdg = ((stateScores['bound_mut'][term] - stateScores['unbound_mut'][term]) -
                      (stateScores['bound_wt'][term] - stateScores['unbound_wt'][term]))
            gamTotal += self._gamTransform(termDdg, term)

        return {'ddG_raw': ddgRaw, 'ddG_gam': gamTotal}

    def _gamTransform(self, x, term):
        a, b = FLEXDDG_GAM_PARAMS[term]
        return -1.0 * math.exp(a) + 2.0 * math.exp(a) / (1.0 + math.exp(-1.0 * x * math.exp(b)))
