from multiprocessing import cpu_count

from gooey import GooeyParser, Gooey

import printers
from main import main as internal_main


def main() -> int:
    parser = GooeyParser()
    parser.add_argument(
        "fasta_file",
        metavar="Fasta file",
        help="path to the fasta file",
        widget="FileChooser",
    )
    parser.add_argument(
        "species_file",
        metavar="Species file",
        help="path to the species file",
        widget="FileChooser",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="Output file",
        default=None,
        dest="output_file",
    )
    parser.add_argument(
        "-l",
        "--log",
        metavar="Log file",
        default=None,
        dest="log_file",
    )
    parser.add_argument(
        "-p",
        "--processes",
        metavar="Processes",
        default=cpu_count(),
        dest="processes",
        type=int,
        widget="IntegerField",
    )
    parser.add_argument(
        "-u",
        "--unambiguous",
        metavar="Unambiguous",
        default=False,
        action="store_true",
        help="ignore ambiguous bases",
        dest="unambiguous",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        metavar="Threshold",
        default=0.05,
        help="max percentage of ambiguous bases in a group",
        dest="threshold",
        type=float,
    )
    parser.add_argument(
        "--full-matrix",
        metavar="Full matrix",
        default=None,
        help="compute distances for all sequences",
        dest="full_matrix",
    )
    parser.add_argument(
        "-s",
        "--sep",
        metavar="Separator",
        choices=list(printers.SEPARATORS.keys()),
        default=printers.DEFAULT_SEPARATOR,
        help="separator for the output file",
        dest="separator",
    )
    return Gooey(
        program_name="Hooty",
        program_description="Compute distances between groups of sequences",
        use_cmd_args=True,
    )(internal_main)(parser)


if __name__ == "__main__":
    error_code = main()
    if error_code != 0:
        exit(error_code)
