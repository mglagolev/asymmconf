#!/usr/bin/python3
from simple_lammps import angle_string, dihedral_string
from simple_lammps import atom_string, vel_string, bond_string
import math
import sys
import argparse

parser = argparse.ArgumentParser(
    description = 'Create an initial configuration for multiblock copolymer'
                + ' from alternating flexible and helical blocks')

parser.add_argument(
    '--box', metavar = 'SIZE', type = float, nargs = 1,
    help = "Size of the cubic box")

parser.add_argument(
    '--nchains', metavar = 'CHAINS', type = int, nargs = 1,
    default = [1], help = "Number of polymer chains")

parser.add_argument(
    '--nblocks', metavar = 'BLOCKS', type = int, nargs = 1,
     help = "Number of blocks in a macromolecule")

parser.add_argument(
    '--n1', metavar = 'LENGTH', type = int, nargs = 1,
     help = "Length of the odd blocks")

parser.add_argument(
    '--n2', metavar = 'LENGTH', type = int, nargs = 1,
     help = "Length of the even blocks")

parser.add_argument(
    '--pendants2', metavar = 'PATTERN', type = int, nargs = '*',
    default = [0], help = "Repeating pattern containing the number of"
    + " pendants for the monomer units of the even block")

args = parser.parse_args()

box = args.box[0]
nchains = args.nchains[0]
nblocks = args.nblocks[0]
nblock1 = args.n1[0]
nblock2 = args.n2[0]
pendants_pattern = args.pendants2

ncols = int(math.sqrt(nchains))
nrows = math.ceil(nchains / ncols)

n_odd = nblocks % 2 + int(nblocks / 2)
n_even = int(nblocks / 2)
chain_length = nblock1 * n_odd + nblock2 * n_even 

xspacing = box / ncols
yspacing = box / nrows
lbond = min(1., box / float(chain_length))

lbond1, lbond2, lbond12 = lbond, lbond, lbond
xoffset, yoffset, zoffset = xspacing / 2., yspacing / 2., lbond / 2.


#Create lammps data

atoms = []
bonds = []
angles = []
dihedrals = []
out = ""


atomid = 0
for nchain in range(nchains):
    nrows = int(math.sqrt(nchains))
    ixmol = nchain % nrows
    iymol = int(nchain / nrows)
    x = xoffset + ixmol * xspacing
    y = yoffset + iymol * yspacing
    z = zoffset

    nblock = 0
    while True:
        this_block1_list = []
        this_block2_list = []
        nblock += 1
        if nblock > nblocks:
            break
        for nmono in range(nblock1):
            atomid += 1
            if nmono == 0:
                if nblock > 1:
                    z += lbond12
                    bonds.append(["3", atomid-1, atomid])
            else:
                z += lbond1
                bonds.append(["1", atomid-1, atomid])
            atoms.append([atomid, nchain + 1, "1", x, y, z, 0, 0, 0])
            this_block1_list.append(atomid)

        nblock += 1
        if nblock > nblocks:
            break
        for nmono in range(nblock2):
            atomid += 1
            this_block2_list.append(atomid)
            if nmono == 0:
                if nblock > 1:
                    z += lbond12
                    bonds.append(["3", atomid-1, atomid])
            else:
                z += lbond2
                bonds.append(["2", this_block2_list[-2],
                                   this_block2_list[-1]])
            atoms.append([atomid, nchain + 1, "2", x, y, z, 0, 0, 0])
            if len(this_block2_list) >= 3:
                angles.append(["1", this_block2_list[-3],
                                    this_block2_list[-2],
                                    this_block2_list[-1]])
            if len(this_block2_list) >= 4:
                dihedrals.append(["1", this_block2_list[-4],
                                       this_block2_list[-3],
                                       this_block2_list[-2],
                                       this_block2_list[-1]])
            pendant_index = nmono % len(pendants_pattern)
            for npendant in range(pendants_pattern[pendant_index]):
                atomid += 1
                atoms.append([atomid, nchain + 1, "3", x, y+lbond1*(npendant+1), z,
                                   0, 0, 0 ])
                bonds.append(["4", atomid-1, atomid])

natomtypes = len(set([atom[2] for atom in atoms]))
nbondtypes = len(set([bond[0] for bond in bonds]))
nangletypes = len(set([angle[0] for angle in angles]))
ndihedraltypes = len(set([dihedral[0] for dihedral in dihedrals]))

#Headers
out = "LAMMPS Description\n\n"\
    + str(atomid).rjust(8) + " atoms\n"\
    + str(len(bonds)).rjust(8) + " bonds\n"\
    + str(len(angles)).rjust(8) + " angles\n"\
    + str(len(dihedrals)).rjust(8) + " dihedrals\n"\
    + "       " + str(natomtypes) + " atom types\n"\
    + "       " + str(nbondtypes) + " bond types\n"\
    + "       " + str(nangletypes) + " angle types\n"\
    + "       " + str(ndihedraltypes) + " dihedral types\n"\
    + "\n"\
    + "0 " + str(box) + " xlo xhi"\
    + "0 " + str(box) + " ylo yhi"\
    + "0 " + str(box) + " zlo zhi"\

#Atoms
out += "\nAtoms\n\n"
# atom id, mol id, atom type, x, y, z, ix, iy, iz
for atom in atoms:
    out += atom_string(*atom) + "\n"

#Velocities
out += "\nVelocities\n\n"
# atom id, vx, vy, vz
for iatom in range(1, atomid + 1):
    out += vel_string(iatom, 0., 0., 0.) + "\n"

#Bonds
out += "\nBonds\n\n"
# bond id, bond type, atom1 id, atom2 id
for ibond in range(len(bonds)):
    bond = bonds[ibond]
    out += bond_string(ibond + 1, bond[0], bond[1], bond[2]) + "\n"

#Angles
if len(angles) > 0:
    out += "\nAngles\n\n"
# angle id, angle type, atom1 id, atom2 id, atom3 id
for iangle in range(len(angles)):
    angle = angles[iangle]
    out += angle_string(iangle + 1, angle[0], angle[1],
                        angle[2], angle[3]) + "\n"

#Dihedrals
if len(dihedrals) > 0:
    out += "\nDihedrals\n\n"
# dihedral id, dihedral type, atom1, atom2, atom3, atom4
for idihedral in range(len(dihedrals)):
    dihedral = dihedrals[idihedral]
    out += dihedral_string(idihedral + 1, dihedral[0], dihedral[1],
                           dihedral[2], dihedral[3], dihedral[4]) + "\n"


sys.stdout.write(out)