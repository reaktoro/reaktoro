#!/bin/bash

mkdir -p output

python supcrt-convert-xml.py slop98.reaktoro.dat output/supcrt98.xml --exclude=organic-species.txt
python supcrt-convert-xml.py slop98.reaktoro.dat output/supcrt98-organics.xml
python supcrt-convert-xml.py slop07.reaktoro.dat output/supcrt07.xml --exclude=organic-species.txt
python supcrt-convert-xml.py slop07.reaktoro.dat output/supcrt07-organics.xml

echo "All supcrt databases generated successfully."
