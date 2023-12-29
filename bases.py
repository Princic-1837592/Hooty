"""
    Hooty: fast DNA distance matrix
    Copyright (C) 2024  Andrea Princic & Giacomo Chiappa

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

A = 0x0
C = 0x1
G = 0x2
T = 0x3
W = 0x4
S = 0x5
M = 0x6
K = 0x7
R = 0x8
Y = 0x9
B = 0xa
D = 0xb
H = 0xc
V = 0xd
N = 0xe
GAP = 0x10

BYTES_MAP = {
    "A": A,
    "C": C,
    "G": G,
    "T": T,
    "W": W,
    "S": S,
    "M": M,
    "K": K,
    "R": R,
    "Y": Y,
    "B": B,
    "D": D,
    "H": H,
    "V": V,
    "N": N,
    "-": GAP,
}

REPLACEMENTS = [
    (A, T),
    (G, C),
    (A, C),
    (G, T),
    (A, G),
    (C, T),
    (C, G, T),
    (A, G, T),
    (A, C, T),
    (A, C, G),
]
