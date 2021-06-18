#!/bin/bash

# Create an output dir if not existent
mkdir -p output

# Make sure output dir is empty
rm -rf output/*

# Convert the original SUPCRT databases to Reaktoro's xml format
python supcrt-convert-xml.py slop98.reaktoro.dat output/supcrt98.xml --exclude=organic-species.txt
python supcrt-convert-xml.py slop98.reaktoro.dat output/supcrt98-organics.xml
python supcrt-convert-xml.py slop07.reaktoro.dat output/supcrt07.xml --exclude=organic-species.txt
python supcrt-convert-xml.py slop07.reaktoro.dat output/supcrt07-organics.xml

# Go inside the output dir
cd output

# For each supcrt database file in xml format, zip it
for f in *.xml; do zip ${f%.*} $f; done

# For each zipped database file, create a C++ header containing
# a static  array of unsigned char that represents the contents
# of the zipped file as hexadecimals
for f in *.zip; do xxd -i $f >${f%.*}.hpp; done

# Move xml, zip, and hpp files to dedicated directories
mkdir {zip,hpp,xml}
mv *.zip zip
mv *.hpp hpp
mv *.xml xml

# Step out from the output dir
cd ..

echo "All supcrt databases generated successfully."
