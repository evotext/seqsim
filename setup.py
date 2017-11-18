#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

with open('README.md') as readme_file:
    readme = readme_file.read()

requirements = []

test_requirements = []

setup(
    author="Tiago Tresoldi",
    author_email="tresoldi@gmail.com",
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Software Development :: Libraries',
        ],
    description="Package for calculating sequence similarity (string similarity).",
    include_package_data=True,
    install_requirements=requirements,
    keywords='sequence string similarity',
    license='MIT',
    long_description=readme,
    name='seqsim',
    package_dir={'seqsim' : 'seqsim'},
    packages=['seqsim'],
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/tresoldi/seqsim',
    version='0.0.1',
)
