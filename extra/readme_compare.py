"""
Command-line script for generating the demo table for README files.
"""

import seqsim
import tabulate


def main():
    demo1 = ["kitten", "sitting"]
    demo2 = [(1, 2, 3, 4), (3, 4, 2, 1)]

    # Collect results for all methods
    ret = []
    for method in sorted(seqsim.METHODS):
        dist1 = seqsim.distance(demo1, method=method)
        dist2 = seqsim.distance(demo2, method=method)
        ret.append([method, dist1, dist2])

        dist1_norm = seqsim.distance(demo1, method=method, normal=True)
        dist2_norm = seqsim.distance(demo2, method=method, normal=True)
        ret.append([method + "_normal", dist1_norm, dist2_norm])

    print(
        tabulate.tabulate(
            ret,
            headers=["Method", '"kitten" / "sitting"', "(1, 2, 3, 4) / (3, 4, 2, 1)"],
            tablefmt="github",
        )
    )


if __name__ == "__main__":
    main()
