"""
Setup file for `seqsim`.
"""

# Import Python standard libraries
import pathlib

from setuptools import setup, find_packages

# The directory containing this file
LOCAL_PATH = pathlib.Path(__file__).parent

# The text of the README file
README_FILE = (LOCAL_PATH / "README.md").read_text(encoding="utf-8")

# Load requirements, so they are listed in a single place
with open("requirements.txt", encoding="utf-8") as fp:
    install_requires = [dep.strip() for dep in fp.readlines()]

# This call to setup() does all the work
setup(
    author="Tiago Tresoldi",
    author_email="tiago.tresoldi@lingfil.uu.se",
    classifiers=["License :: OSI Approved :: MIT License",
                 "Operating System :: OS Independent",
                 "Programming Language :: Python :: 3",
                 "Topic :: Software Development :: Libraries", ],
    description="Library for computing measures of similarity for sequences of hashable data types",
    # entry_points={"console_scripts": ["seqsim=seqsim.__main__:main"]},
    include_package_data=True,
    install_requires=install_requires,
    keywords=['sequence similarity', "sequence distance", 'string similarity', "string distance"],
    license="MIT",
    long_description=README_FILE,
    long_description_content_type="text/markdown",
    name="seqsim",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.6",
    test_suite="tests",
    tests_require=[],
    url="https://github.com/tresoldi/seqsim",
    version="0.2",  # remember to sync with __init__.py
    zip_safe=False,
)
