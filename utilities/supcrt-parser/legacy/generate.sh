#!/usr/bin/env bash

# This file generates the databases in all supported formats.
#
# WARNING: This file needs to be executed from the root directory of Supkrit!!!

for f in databases/linearized/*.pre.db; do ./databases/scripts/convert.py $f databases/generated; done

echo '** All files in this directory have been automatically moved to supkrit/Supkrit/databases **' > databases/generated/README.md
