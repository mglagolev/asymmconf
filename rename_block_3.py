#!/usr/bin/python3

import MDAnalysis as mda
import argparse

parser = argparse.ArgumentParser(
    description = 'Rename the second A block in an ABA triblock copolymer')

parser.add_argument(
    'input', metavar = 'LAMMPS_DATA', action = "store", help = "input file")

parser.add_argument(
    'output', metavar = 'LAMMPS_DATA', action = "store", help = "output file")

parser.add_argument(
    '--newtype', metavar = 'TYPE', type = str, nargs = 1,
    default = ["4"], help = "New type for block #3")

args = parser.parse_args()

u = mda.Universe(args.input, format = 'DATA')

for i in range(u.atoms.n_atoms):
    if i > 0:
        if u.atoms.resids[i] != u.atoms.resids[i-1]:
            block_1_type = u.atoms.types[i]
            block_1_ended = False
    if u.atoms.types[i] != block_1_type:
        block_1_ended = True
    if block_1_ended == True and u.atoms.types[i] == block_1_type:
        u.atoms.types[i] = args.newtype[0]

u.trajectory.ts.data['molecule_tag'] = u.atoms.resids
u.atoms.write(args.output)