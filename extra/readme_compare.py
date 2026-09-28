"""
Command-line script for generating the demo table for README files.

Requires the `tabulate` package (`pip install seqsim[extra]`).
"""

import seqsim
import tabulate


def main():
    demo1 = ["kitten", "sitting"]
    demo2 = [(1, 2, 3, 4), (3, 4, 2, 1)]

    # Collect results for all methods
    ret = []
    for method in sorted(seqsim.METHODS):
        func = seqsim.METHODS[method]
        row = [method, f"`{func.__module__.split('.')[-1]}.{func.__name__}`"]
        for demo in (demo1, demo2):
            row.append(seqsim.distance(demo, method=method))
            row.append(seqsim.distance(demo, method=method, normal=True))
        ret.append(row)

    print(
        tabulate.tabulate(
            ret,
            headers=[
                "Method",
                "Function",
                '"kitten" / "sitting"',
                "normalized",
                "(1, 2, 3, 4) / (3, 4, 2, 1)",
                "normalized",
            ],
            tablefmt="github",
            floatfmt=".4f",
        )
    )


if __name__ == "__main__":
    main()
