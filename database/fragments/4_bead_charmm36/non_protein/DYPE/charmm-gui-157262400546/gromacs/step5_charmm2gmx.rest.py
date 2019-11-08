"""
Generated by CHARMM-GUI (http://www.charmm-gui.org)

step5_charmm2gmx_rest.py

This program is for generating GROMACS restraints file.

Correspondance: jul316@lehigh.edu or wonpil@lehigh.edu
Last update: August 20, 2015
"""

import sys, os

infile  = sys.argv[1]

funcForPosres = 1
funcForDihres = 1

posres = []
dihres = []

prot = False
memb = False
mice = False
carb = False
hmmm = False
hchl = False

for line in open(infile.lower(), 'r'):
    line = line.strip()
    if len(line) > 0:
        segments = line.split()
        moltype = segments[0]
        restype = segments[1]
        if moltype == 'PROT' and restype == 'POS':
            atom1 = int(segments[2])
            stat1 = segments[3]
            posres.append([moltype, atom1, stat1])
            prot  = True
        if moltype == 'MEMB' and restype == 'POS':
            atom1 = int(segments[2])
            posres.append([moltype, atom1])
            memb  = True
        if moltype == 'MEMB' and restype == 'DIH':
            atom1 = int(segments[2])
            atom2 = int(segments[3])
            atom3 = int(segments[4])
            atom4 = int(segments[5])
            phi0  = float(segments[6])
            dphi  = float(segments[7])
            dihres.append([moltype, atom1, atom2, atom3, atom4, phi0, dphi])
            memb  = True
        if moltype == 'MICE' and restype == 'POS':
            atom1 = int(segments[2])
            posres.append([moltype, atom1])
            mice  = True
        if moltype == 'CARB' and restype == 'DIH':
            atom1 = int(segments[2])
            atom2 = int(segments[3])
            atom3 = int(segments[4])
            atom4 = int(segments[5])
            phi0  = float(segments[6])
            dphi  = float(segments[7])
            dihres.append([moltype, atom1, atom2, atom3, atom4, phi0, dphi])
            carb  = True
        if ( moltype == 'HMMM' or moltype == 'HCHL' ) and restype == 'POS':
            atom1 = int(segments[2])
            posres.append([moltype, atom1])
            hmmm  = True
            if moltype == 'HCHL': hchl = True

posres.sort()
dihres.sort()

if len(posres) == 0 and len(dihres) == 0: exit()

if not os.path.exists('gromacs/restraints'):
    os.mkdir('gromacs/restraints')
outfile = open('gromacs/restraints/'+infile.upper().split('_')[1]+'_rest.itp', 'w')

fc_memb  = {}

fc_memb['bb']    = [4000.0, 4000.0, 2000.0, 1000.0,  500.0,  200.0,   50.0,    0.0]
fc_memb['sc']    = [2000.0, 2000.0, 1000.0,  500.0,  200.0,   50.0,    0.0,    0.0]
fc_memb['lpos']  = [1000.0, 1000.0, 1000.0,  400.0,  200.0,   40.0,    0.0,    0.0]
fc_memb['ldih']  = [1000.0, 1000.0,  400.0,  200.0,  200.0,  100.0,    0.0,    0.0]
fc_memb['carb']  = [1000.0, 1000.0,  400.0,  200.0,  200.0,  100.0,    0.0,    0.0]
fc_memb['hmmm']  = [  20.0,   20.0,   20.0,   20.0,   20.0,   20.0,   20.0,   20.0]
fc_memb['mpos']  = [1000.0, 1000.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0]

for step in range(8):
    if step == 7:
        outfile.write('#ifdef STEP7\n')
    else:
        outfile.write('#ifdef STEP6_%d\n' % step)
    if prot:
        outfile.write('#define fc_bb %.1f\n' % fc_memb['bb'][step])
        outfile.write('#define fc_sc %.1f\n' % fc_memb['sc'][step])
    if memb:
        outfile.write('#define fc_lpos %.1f\n' % fc_memb['lpos'][step])
        outfile.write('#define fc_ldih %.1f\n' % fc_memb['ldih'][step])
    if mice: outfile.write('#define fc_mpos %.1f\n' % fc_memb['mpos'][step])
    if carb: outfile.write('#define fc_cdih %.1f\n' % fc_memb['carb'][step])
    if hmmm:
        if hchl: outfile.write('#define fc_hmmm %.1f\n' % (fc_memb['hmmm'][step]*10))
        else:    outfile.write('#define fc_hmmm %.1f\n' % fc_memb['hmmm'][step])
    outfile.write('#endif\n')

if len(posres) > 0:
    outfile.write('\n[ position_restraints ]\n')
    for iposres in posres:
        moltype = iposres[0]
        atom1 = iposres[1]
        if moltype == 'PROT':
            stat1 = iposres[2]
            if stat1 == 'BB': fc = 'fc_bb'
            if stat1 == 'SC': fc = 'fc_sc'
            outfile.write('%5d %5d %8s %8s %8s\n' % (atom1, funcForPosres, fc, fc, fc))
        elif moltype == 'MICE':
            fc = 'fc_mpos'
            outfile.write('%5d %5d %8s %8s %8s\n' % (atom1, funcForPosres, fc, fc, fc))
        else:
            if moltype == 'MEMB': fc = 'fc_lpos'
            if moltype == 'HMMM' or moltype == 'HCHL': fc = 'fc_hmmm'
            outfile.write('%5d %5d %8.1f %8.1f %8s\n' % (atom1, funcForPosres, 0, 0, fc))

if len(dihres) > 0:
    outfile.write('\n[ dihedral_restraints ]\n')
    for idihres in dihres:
        moltype = idihres[0]
        atom1   = idihres[1]
        atom2   = idihres[2]
        atom3   = idihres[3]
        atom4   = idihres[4]
        phi0    = idihres[5]
        dphi    = idihres[6]
        if moltype == 'MEMB': fc = 'fc_ldih'
        if moltype == 'CARB': fc = 'fc_cdih'
        outfile.write('%5d %5d %5d %5d %5d %8.1f %8.1f %8s\n' % (atom1, atom2, atom3, atom4, funcForDihres, phi0, dphi, fc))

