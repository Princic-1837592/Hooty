from math import inf, nan


def to_csv(matrix, species, groups):
    return to_sv(matrix, ",", species, groups)


def to_tsv(matrix, species, groups):
    return to_sv(matrix, "\t", species, groups)


def to_sv(matrix, separator, species, groups):
    header_groups = [[] for _ in range(groups[-1] + 2)]
    for specie, group in zip(species, groups):
        header_groups[group + 1].append(specie)
    header = list(map(lambda g: "-".join(g), header_groups))
    lines = [separator.join(header)]
    for row, group in zip(matrix, header[1:]):
        numbers = separator.join(map(format_cell, row))
        lines.append(f"{group}{separator}{numbers}")
    return "\n".join(lines)


def format_cell(cell):
    vals = []
    for val in cell:
        vals.append(format_value(val))
    return " - ".join(vals)


def format_value(val):
    if abs(val) == 0.0:
        val = " 0.00"
    elif val in (nan, inf, -inf):
        val = "  /  "
    else:
        val = f"{val * 100:>5.2f}"
    return val
