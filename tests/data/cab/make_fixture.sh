#!/bin/sh
# Rebuild tests/data/cab_affinity.json from CAB's own functions.
#
#   sh tests/data/cab/make_fixture.sh PATH/TO/common.js
#
# CAB's code is not distributed with seqsim (it carries no licence). The
# file used was common.js version 30, loaded by the page of version 34 of
# https://ada.geschkult.fu-berlin.de/manuscript_analysis/proximity.html,
# with SHA-256 b70c72091346770a904df9adec89817c0785dd1f91290099a3aaf7274ff98cc8.
# The script copies, by line number, CAB's functions baseOf, isBisWitness,
# isGapReading, groupReadings, typeWeight, locusRelevant, witnessReadings,
# intersects, bayesPosterior and appPairMatrix (and buildMsProfiles and
# pairMsStats for the second generator) into a node script, followed by the
# synthetic data generators of this folder, runs them with node, and
# converts their output with build_fixture.py.
set -e
COMMON=$1
HERE=$(cd "$(dirname "$0")" && pwd)
echo "b70c72091346770a904df9adec89817c0785dd1f91290099a3aaf7274ff98cc8  $COMMON" | sha256sum -c -
WORK=$(mktemp -d)
LINES='6,7p;12p;67p;101,104p;106p;112,147p'
{ echo "const A={};"; sed -n "$LINES" "$COMMON"; cat "$HERE/generator1.js"; } > "$WORK/h1.js"
{ echo "const A={};"; sed -n "$LINES;151,152p" "$COMMON"; cat "$HERE/generator2.js"; } > "$WORK/h2.js"
{ echo "const A={};"; sed -n "$LINES" "$COMMON"; cat "$HERE/generator3.js"; } > "$WORK/h3.js"
(cd "$WORK" && node h1.js && node h2.js && node h3.js)
python3 "$HERE/build_fixture.py" "$WORK" "$HERE/../cab_affinity.json"
