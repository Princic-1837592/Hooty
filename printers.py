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

from math import inf, nan

SEPARATORS = {
    "comma": ",",
    "semicolon": ";",
    "tab": "\t",
}
DEFAULT_SEPARATOR = "semicolon"


def to_sv(matrix, species, groups, separator, precision):
    separator = SEPARATORS.get(separator, SEPARATORS[DEFAULT_SEPARATOR])
    header_groups = [[] for _ in range(groups[-1] + 2)]
    for specie, group in zip(species, groups):
        header_groups[group + 1].append(specie)
    header = list(map(lambda g: "-".join(g), header_groups))
    lines = [separator.join(header)]
    for row, group in zip(matrix, header[1:]):
        numbers = separator.join(map(lambda c: format_cell(c, precision), row))
        lines.append(f"{group}{separator}{numbers}")
    return "\n".join(lines)


def format_cell(cell, precision):
    vals = []
    for val in cell if type(cell) == tuple else (cell,):
        vals.append(format_value(val, precision))
    return " – ".join(vals)


def format_value(val, precision):
    if abs(val) == 0.0:
        val = f"0"
    elif val in (nan, inf, -inf):
        val = "/"
    else:
        val = f"{val * 100:.{precision}f}"
    return val
