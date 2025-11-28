# PyMOL Visualization Script - mol_001 Protein-Ligand Interactions
# =====================================================================
# Molecule: 4-((6,7-dimethoxyquinazolin-4-yl)amino)cyclohexane-1-carboxylic acid
# SMILES: COc1cc2ncnc(NC3CCC(C(=O)O)CC3)c2cc1OC
# Vina Affinity: -8.78 kcal/mol
# Target: EGFR Kinase (PDB: 1M17)
# =====================================================================

# Reinitialize PyMOL
reinitialize

# =============================================================================
# LOAD STRUCTURES
# =============================================================================

# Load EGFR receptor (use PDB for better residue info)
fetch 1M17, EGFR
# Alternative: load local file
# load C:/De Novo Drug Discovery AI/docking/vina_results/1M17_receptor.pdbqt, EGFR

# Load mol_001 docked pose
load C:/De Novo Drug Discovery AI/docking/vina_results/pdbqt/mol_001_docked.pdbqt, mol_001

# Alternative if PDB pose exists:
# load C:/De Novo Drug Discovery AI/docking/vina_results/poses_pdb/mol_001_pose1.pdb, mol_001

# =============================================================================
# CLEAN UP RECEPTOR
# =============================================================================

# Remove waters and other small molecules
remove resn HOH
remove resn AQ4  # Remove native erlotinib ligand for clarity

# =============================================================================
# PROTEIN DISPLAY SETTINGS
# =============================================================================

# Hide everything first
hide everything

# Show protein as cartoon
show cartoon, EGFR
color gray80, EGFR
set cartoon_transparency, 0.7, EGFR

# Show surface for context (optional - can be slow)
# show surface, EGFR
# set surface_color, white, EGFR
# set transparency, 0.85, EGFR

# =============================================================================
# LIGAND DISPLAY - mol_001
# =============================================================================

# Show mol_001 as sticks
show sticks, mol_001
set stick_radius, 0.25, mol_001

# Color by element with orange carbons
color orange, mol_001 and elem C
color red, mol_001 and elem O
color blue, mol_001 and elem N
color white, mol_001 and elem H

# Show mesh around ligand
# create mol_001_mesh, mol_001
# show mesh, mol_001_mesh
# color orange, mol_001_mesh
# set mesh_width, 0.5

# =============================================================================
# BINDING SITE DEFINITION
# =============================================================================

# Select residues within 5A of ligand
select binding_pocket, byres (EGFR within 5 of mol_001)

# Key EGFR ATP-binding site residues
select key_residues, EGFR and resi 718+719+720+721+726+743+744+745+766+767+768+769+770+771+772+773+774+775+776+777+778+779+788+790+791+792+793+794+795+796+797+830+831+832+833+834+835+836+837+838+839+853+854+855+856+857

# Gatekeeper and hinge region
select hinge_region, EGFR and resi 790+791+792+793+794+795
select gatekeeper, EGFR and resi 790  # T790 - important for resistance

# Catalytic residues
select catalytic, EGFR and resi 793+797+835+836+837

# =============================================================================
# DISPLAY BINDING SITE
# =============================================================================

# Show binding pocket as lines
show lines, binding_pocket

# Show key residues as sticks
show sticks, key_residues
color lightpink, key_residues and elem C
util.cnc key_residues

# Highlight specific interaction partners
select THR830, EGFR and resi 830
select THR790, EGFR and resi 790
select MET793, EGFR and resi 793
select LYS745, EGFR and resi 745
select ASP855, EGFR and resi 855
select GLU762, EGFR and resi 762

# Color specific residues
color tv_orange, THR830 and elem C
color tv_blue, MET793 and elem C
color tv_yellow, LYS745 and elem C
color tv_green, THR790 and elem C

# Show interaction partners as sticks
show sticks, THR830 or MET793 or LYS745 or THR790 or ASP855

# =============================================================================
# HYDROGEN BOND ANALYSIS
# =============================================================================

# Detect H-bonds using geometry-based method (mode=2)
# Cutoff: 3.5 Angstroms

# H-bonds between mol_001 and protein
dist hbonds_all, mol_001, EGFR, mode=2, cutoff=3.5

# Specific H-bond to THR830 (documented interaction)
dist hbond_THR830, mol_001, THR830, mode=2, cutoff=3.5
color orange, hbond_THR830
set dash_width, 3, hbond_THR830
set dash_gap, 0.2, hbond_THR830
set dash_color, orange, hbond_THR830

# H-bonds to hinge region (MET793 - key for quinazoline binding)
dist hbond_hinge, mol_001, MET793, mode=2, cutoff=3.5
color cyan, hbond_hinge
set dash_width, 3, hbond_hinge

# General H-bond styling
set dash_round_ends, 1
hide labels, hbonds_all
color yellow, hbonds_all

# =============================================================================
# HYDROPHOBIC CONTACTS
# =============================================================================

# Select hydrophobic residues near ligand
select hydrophobic_contacts, EGFR and resi 718+726+751+766+775+790+797+835+838+855 and (byres (EGFR within 4 of mol_001))

# Show as thin sticks
show sticks, hydrophobic_contacts
color palegreen, hydrophobic_contacts and elem C

# =============================================================================
# POLAR CONTACTS ANALYSIS
# =============================================================================

# Find all polar contacts (H-bonds + salt bridges)
dist polar_contacts, mol_001, binding_pocket, mode=2, cutoff=3.5
color magenta, polar_contacts
set dash_width, 2, polar_contacts

# =============================================================================
# ELECTROSTATIC SURFACE (Optional)
# =============================================================================

# Uncomment to show electrostatic potential surface
# set surface_quality, 1
# set surface_type, 0
# show surface, binding_pocket
# set_color pos_color, [0.0, 0.0, 1.0]
# set_color neg_color, [1.0, 0.0, 0.0]
# ramp_new e_pot, EGFR, [-5, 0, 5], [neg_color, white, pos_color]
# set surface_color, e_pot, EGFR

# =============================================================================
# LABELING
# =============================================================================

# Label key residues
label key_residues and name CA, "%s%s" % (resn, resi)
set label_color, black
set label_size, 14
set label_position, (0, 0, 3)

# Label specific interaction partners more prominently
label THR830 and name CA, "THR830\n(H-bond)"
label MET793 and name CA, "MET793\n(Hinge)"
label THR790 and name CA, "THR790\n(Gatekeeper)"
label LYS745 and name CA, "LYS745"

# =============================================================================
# LIGAND PHARMACOPHORE HIGHLIGHTS
# =============================================================================

# Quinazoline nitrogen atoms (H-bond acceptors)
select quin_N, mol_001 and elem N
show spheres, quin_N
set sphere_scale, 0.3, quin_N
color blue, quin_N

# Carboxylic acid (H-bond donor/acceptor)
select COOH, mol_001 and (name O* or name OXT)
show spheres, COOH
set sphere_scale, 0.25, COOH
color red, COOH

# Methoxy groups
select OMe, mol_001 and name O and neighbor (name C and neighbor (name C))
show spheres, OMe
set sphere_scale, 0.2, OMe
color firebrick, OMe

# =============================================================================
# DISTANCE MEASUREMENTS
# =============================================================================

# Measure key distances manually
# NH to MET793 backbone carbonyl (typical quinazoline H-bond)
# dist d1, mol_001 and elem N, MET793 and name O, mode=0

# COOH to THR830 hydroxyl
# dist d2, mol_001 and name O, THR830 and name OG1, mode=0

# =============================================================================
# VIEW SETTINGS
# =============================================================================

# Center on mol_001
center mol_001
zoom mol_001, 12

# Set clipping planes
clip slab, 40

# Background
bg_color white

# Depth cue and fog
set depth_cue, 0
set fog, 0

# Specular lighting
set spec_reflect, 1.5
set spec_power, 200
set specular, 1

# =============================================================================
# RAY TRACING SETTINGS
# =============================================================================

set ray_shadows, 1
set ray_trace_mode, 1
set antialias, 2
set ray_opaque_background, 1

# =============================================================================
# CREATE GROUPS FOR EASY TOGGLING
# =============================================================================

group Ligand, mol_001
group Protein, EGFR
group Binding_Site, binding_pocket key_residues
group Interactions, hbonds_all hbond_THR830 hbond_hinge polar_contacts
group Key_Residues, THR830 THR790 MET793 LYS745 ASP855

# =============================================================================
# SAVE VIEWS
# =============================================================================

# View 1: Overview
set_view (\
     0.980226040,    0.113399386,   -0.162182897,\
    -0.060825769,    0.944879711,    0.321765423,\
     0.189751863,   -0.305610687,    0.932927132,\
     0.000000000,    0.000000000,  -65.000000000,\
    22.000000000,    0.300000000,   52.800003052,\
    50.000000000,   80.000000000,  -20.000000000 )

# =============================================================================
# SCENE CREATION (Optional)
# =============================================================================

# Scene 1: Full complex
scene overview, store

# Scene 2: Focus on H-bonds
hide cartoon
zoom mol_001, 8
scene hbond_focus, store

# Scene 3: Binding pocket surface
show cartoon
show surface, binding_pocket
set transparency, 0.7, binding_pocket
scene pocket_surface, store

# Return to overview
scene overview, recall

# =============================================================================
# OUTPUT INSTRUCTIONS
# =============================================================================

print ""
print "============================================================"
print "mol_001 - EGFR INTERACTION ANALYSIS"
print "============================================================"
print ""
print "Structure: 4-((6,7-dimethoxyquinazolin-4-yl)amino)cyclohexane-1-carboxylic acid"
print "Vina Affinity: -8.78 kcal/mol"
print "QED Score: 0.87"
print ""
print "KEY INTERACTIONS:"
print "  - H-bond: NH -> THR830-OG1 (~3.0 A)"
print "  - H-bond: Quinazoline N1 -> MET793 (hinge region)"
print "  - Hydrophobic: Cyclohexyl in hydrophobic pocket"
print "  - Methoxy groups: van der Waals with Leu788, Val726"
print ""
print "GROUPS AVAILABLE:"
print "  Ligand        - mol_001 molecule"
print "  Protein       - EGFR receptor"
print "  Binding_Site  - Pocket residues"
print "  Interactions  - H-bonds and contacts"
print "  Key_Residues  - Important binding partners"
print ""
print "SCENES AVAILABLE:"
print "  scene overview        - Full view"
print "  scene hbond_focus     - H-bond closeup"
print "  scene pocket_surface  - Surface representation"
print ""
print "TO RENDER IMAGE:"
print "  ray 2400, 2400"
print "  png mol_001_EGFR_interaction.png, dpi=300"
print ""
print "============================================================"
