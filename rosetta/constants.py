# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Alberto Manuel Parra Pérez (amparraperez@gmail.com)
# *           Judith Maestro Ciria 
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

# Path to the binaries inside the ROSETTA folder
ROSETTA_BINARIES_PATH = "main/source/bin"
ROSETTA_DATABASE_PATH = "main/database"
ROSETTA_PARAMS_PATH = "main/source/scripts/python/public"


# Name of programs for linux
SCORE = 'score.static.linuxgccrelease'  # rescores PDBs and silent files, extracts, PDBs from silent files,
                                        # assembles PDBs into silent files.

PARAMS_FILE = 'molfile_to_params.py'
BATCH_PARAMS_FILE = 'batch_molfile_to_params.py'

MAKE_RAY_FILES = 'make_ray_files.static.linuxgccrelease'  # create a ray file to map the pocket or interface
MAKE_RAY_FILES_GPU = 'make_ray_files.opencl.linuxgccrelease'  # create a ray file to map the pocket or interface
                                                              # using GPU

ROSETTA_SCRIPTS = 'rosetta_scripts.static.linuxgccrelease'
ROSETTA_SCRIPTS_GPU = 'rosetta_scripts.opencl.linuxgccrelease'


# ------------------------------------ Flex ddG ------------------------------------
AA_THREE_TO_ONE = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
                   'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
                   'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
                   'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}

CANONICAL_AAS = 'ACDEFGHIKLMNPQRSTVWY'

FLEXDDG_RESFILE = 'nataa_mutations.resfile'
FLEXDDG_XML_FILE = 'ddG-backrub.xml'
FLEXDDG_DB_FILE = 'ddG.db3'

# GAM reweighting parameters fitted on the ZEMu benchmark set (Barlow et al. 2018, J. Phys. Chem. B)
# used to reweight the raw talaris2014 ddG score terms into the final, benchmarked ddG prediction.
FLEXDDG_GAM_PARAMS = {
    'fa_sol':      (6.940, -6.722),
    'hbond_sc':    (1.902, -1.999),
    'hbond_bb_sc': (0.063,  0.452),
    'fa_rep':      (1.659, -0.836),
    'fa_elec':     (0.697, -0.122),
    'hbond_lr_bb': (2.738, -1.179),
    'fa_atr':      (2.313, -1.649),
}

# RosettaScripts protocol capture from https://github.com/Kortemme-Lab/flex_ddG_tutorial (ddG-backrub.xml)
flexDDGXML = '''<ROSETTASCRIPTS>
  <SCOREFXNS>
    <ScoreFunction name="fa_talaris2014" weights="talaris2014"/>
    <ScoreFunction name="fa_talaris2014_cst" weights="talaris2014">
      <Reweight scoretype="atom_pair_constraint" weight="1.0"/>
      <Set fa_max_dis="9.0"/>
    </ScoreFunction>
  </SCOREFXNS>

  <!-- ### Only required input file (other than PDB) - mutation resfile ### -->
  <!-- #### All residues must be set to be NATAA packable at top of resfile ### -->
  <TASKOPERATIONS>
    <ReadResfile name="res_mutate" filename="%%mutate_resfile_relpath%%"/>
  </TASKOPERATIONS>

  <RESIDUE_SELECTORS>
    <Task name="resselector" fixed="0" packable="0" designable="1" task_operations="res_mutate"/>
    <Neighborhood name="bubble" selector="resselector" distance="8.0"/>
    <PrimarySequenceNeighborhood name="bubble_adjacent" selector="bubble" lower="1" upper="1"/>
    <StoredResidueSubset name="restore_neighbor_shell" subset_name="neighbor_shell"/>
    <Not name="everythingelse" selector="restore_neighbor_shell"/>
  </RESIDUE_SELECTORS>
  <TASKOPERATIONS>
    <OperateOnResidueSubset name="repackonly" selector="restore_neighbor_shell">
      <RestrictToRepackingRLT/>
    </OperateOnResidueSubset>
    <OperateOnResidueSubset name="norepack" selector="everythingelse">
      <PreventRepackingRLT/>
    </OperateOnResidueSubset>
    <UseMultiCoolAnnealer name="multicool" states="6"/>
    <ExtraChiCutoff name="extrachizero" extrachi_cutoff="0"/>
    <InitializeFromCommandline name="commandline_init"/>
    <RestrictToRepacking name="restrict_to_repacking"/>
  </TASKOPERATIONS>

  <FILTERS>
  </FILTERS>

  <MOVERS>
    <StoreResidueSubset name="neighbor_shell_storer" subset_name="neighbor_shell" residue_selector="bubble_adjacent" />

    <AddConstraintsToCurrentConformationMover name="addcst" use_distance_cst="1" coord_dev="0.5" min_seq_sep="0" max_distance="9" CA_only="1" bound_width="0.0" cst_weight="0.0"/>
    <ClearConstraintsMover name="clearcst"/>
    <MinMover name="minimize" scorefxn="fa_talaris2014_cst" chi="1" bb="1" type="lbfgs_armijo_nonmonotone" tolerance="0.000001" max_iter="%%max_minimization_iter%%" abs_score_convergence_threshold="%%abs_score_convergence_thresh%%"/>

    <PackRotamersMover name="repack" scorefxn="fa_talaris2014" task_operations="commandline_init,repackonly,norepack,multicool"/>
    <PackRotamersMover name="mutate" scorefxn="fa_talaris2014" task_operations="commandline_init,res_mutate,norepack,multicool"/>

    <ReportToDB name="dbreport" batch_description="interface_ddG" database_name="ddG.db3">
      <ScoreTypeFeatures/>
      <ScoreFunctionFeatures scorefxn="fa_talaris2014"/>
      <StructureScoresFeatures scorefxn="fa_talaris2014"/>
    </ReportToDB>

    <ReportToDB name="structreport" batch_description="interface_ddG_struct" database_name="struct.db3">
      <PoseConformationFeatures/>
      <PdbDataFeatures/>
      <JobDataFeatures/>
      <ResidueFeatures/>
      <PoseCommentsFeatures/>
      <ProteinResidueConformationFeatures/>
      <ResidueConformationFeatures/>
    </ReportToDB>

    <SavePoseMover name="save_wt_bound_pose" restore_pose="0" reference_name="wt_bound_pose"/>
    <SavePoseMover name="save_backrub_pose" restore_pose="0" reference_name="backrubpdb"/>
    <SavePoseMover name="restore_backrub_pose" restore_pose="1" reference_name="backrubpdb"/>

    <InterfaceDdGMover name="int_ddG_mover" wt_ref_savepose_mover="save_wt_bound_pose" chain_name="%%chainstomove%%" db_reporter="dbreport" scorefxn="fa_talaris2014"/>

    <ScoreMover name="apply_score" scorefxn="fa_talaris2014_cst" verbose="0"/>

    <!-- This ParsedProtocol allows the ddG calculation to take place multiple times along the backrub trajectory, if desired -->
    <ParsedProtocol name="finish_ddg_post_backrub">
      <Add mover_name="save_backrub_pose"/>
      <Add mover_name="structreport"/>

      <Add mover_name="repack"/>

      <Add mover_name="addcst"/>
      <Add mover_name="minimize"/>
      <Add mover_name="clearcst"/>

      <Add mover_name="save_wt_bound_pose"/>
      <Add mover_name="structreport"/>
      <Add mover_name="restore_backrub_pose"/>

      <Add mover_name="mutate"/>

      <Add mover_name="addcst"/>
      <Add mover_name="minimize"/>
      <Add mover_name="clearcst"/>
      <Add mover_name="structreport"/>

      <Add mover_name="int_ddG_mover"/>
    </ParsedProtocol>

    <BackrubProtocol name="backrub" mc_kt="1.2" ntrials="%%number_backrub_trials%%" pivot_residue_selector="restore_neighbor_shell" task_operations="restrict_to_repacking,commandline_init,extrachizero" recover_low="0" trajectory_stride="%%backrub_trajectory_stride%%" trajectory_apply_mover="finish_ddg_post_backrub"/>

  </MOVERS>
  <APPLY_TO_POSE>
  </APPLY_TO_POSE>
  <PROTOCOLS>
    <Add mover_name="addcst"/>
    <Add mover_name="apply_score"/> <!-- Necessary to initialize neighbor graph -->
    <Add mover_name="neighbor_shell_storer"/>

    <Add mover_name="minimize"/>
    <Add mover_name="clearcst"/>

    <Add mover_name="backrub"/>
  </PROTOCOLS>
  <OUTPUT />
</ROSETTASCRIPTS>'''

generateStructuresXML = '''<ROSETTASCRIPTS>
	<SCOREFXNS>
		<ScoreFunction name="cen" weights="score4_smooth_cart">
			<Reweight scoretype="elec_dens_fast" weight="15"/>
		</ScoreFunction>

		<ScoreFunction name="dens_soft" weights="soft_rep">
			<Reweight scoretype="cart_bonded" weight="0.5"/>
			<Reweight scoretype="pro_close" weight="0.0"/>
			<Reweight scoretype="fa_sol" weight="0.0"/> #REMOVE THIS LINE IF NOT A MEMBRANE PROTEIN
			<Reweight scoretype="elec_dens_fast" weight="25"/>
		</ScoreFunction>
		
		<ScoreFunction name="dens" weights="talaris2014_cart">
			<Reweight scoretype="elec_dens_fast" weight="25"/>
			<Reweight scoretype="fa_sol" weight="0.0"/> #REMOVE THIS LINE IF NOT A MEMBRANE PROTEIN
		<Set scale_sc_dens_byres="R:0.76,K:0.76,E:0.76,D:0.76,M:0.76, C:0.81,Q:0.81,H:0.81,N:0.81,T:0.81,S:0.81,Y:0.88,W:0.88, A:0.88,F:0.88,P:0.88,I:0.88,L:0.88,V:0.88"/> #These values were empirically determined by the ROSETTA group and SHOULD NOT BE CHANGED
		</ScoreFunction>
	</SCOREFXNS>
	<MOVERS>
		<SetupForDensityScoring name="setupdens"/>
		<LoadDensityMap name="loaddens" mapfile="../noMask0723.mrc"/>
		<SwitchResidueTypeSetMover name="tocen" set="centroid"/>
		<MinMover name="cenmin" scorefxn="cen" type="lbfgs_armijo_nonmonotone" max_iter="200" tolerance="0.00001" bb="1" chi="1" jump="ALL"/>
		#This section will use Z-scores to determing poorly fit areas for rebuilding with increasingly more strict cutoffs
		<CartesianSampler name="cen5_50" automode_scorecut="-0.5" scorefxn="cen" mcscorefxn="cen" fascorefxn="dens_soft" strategy="auto" fragbias="density" rms="2.0" ncycles="200" fullatom="0" bbmove="1" nminsteps="25" temp="4"/>
		<CartesianSampler name="cen5_60" automode_scorecut="-0.3" scorefxn="cen" mcscorefxn="cen" fascorefxn="dens_soft" strategy="auto" fragbias="density" rms="1.5" ncycles="200" fullatom="0" bbmove="1" nminsteps="25" temp="4"/>
		<CartesianSampler name="cen5_70" automode_scorecut="-0.1" scorefxn="cen" mcscorefxn="cen" fascorefxn="dens_soft" strategy="auto" fragbias="density" rms="1.5" ncycles="200" fullatom="0" bbmove="1" nminsteps="25" temp="4"/>
		<CartesianSampler name="cen5_80" automode_scorecut="0.0" scorefxn="cen" mcscorefxn="cen" fascorefxn="dens_soft" strategy="auto" fragbias="density" rms="1.0" ncycles="200" fullatom="0" bbmove="1" nminsteps="25" temp="4"/>
		<FastRelax name="relaxcart" scorefxn="dens" repeats="1" cartesian="1"/>
	</MOVERS>
	<PROTOCOLS>
		<Add mover="setupdens"/>
		<Add mover="loaddens"/>
		<Add mover="cenmin"/>
		<Add mover="relaxcart"/>
		<Add mover="cen5_50"/>
		<Add mover="relaxcart"/>
		<Add mover="cen5_60"/>
		<Add mover="relaxcart"/>
		<Add mover="cen5_70"/>
		<Add mover="relaxcart"/>
		<Add mover="cen5_80"/>
		<Add mover="relaxcart"/>
		<Add mover="relaxcart"/>
	</PROTOCOLS>
	<OUTPUT scorefxn="dens"/>
</ROSETTASCRIPTS>'''
