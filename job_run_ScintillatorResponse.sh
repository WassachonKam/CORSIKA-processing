#!/bin/bash
# Set up the IceCube CVMFS toolset (same idea as line 20 of the example)
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/setup.sh`

# Your IceTray build (the folder that contains env-shell.sh)
SOFT=/data/user/sverpoest/software/icetray/build/        # <-- change to your build
ENVSHELL=$SOFT/env-shell.sh

echo "Host: $(hostname)"
echo "Arguments: $@"

# Run the Python script inside the IceTray environment.
# "$@" forwards --particle/--energy/--zenith from the .sub file.
$ENVSHELL python3 /data/user/wkammeem/CORSIKA/submit_ScintillatorResponse.py "$@"