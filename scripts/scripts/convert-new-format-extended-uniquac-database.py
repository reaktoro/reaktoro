# This script converts the Extended UNIQUAC database file to the new format. The
# new format is a dictionary with keys as the formula and values as the
# parameters. The old format is a list of lists where each list has the formula
# and the parameter. The script reads the input yaml file, converts the database
# to the new format, and saves the new yaml file.
#
# The script is called as:
#
# python convert-new-format-extended-uniquac-database.py input.yaml The script
#
# reads the input.yaml file, converts the database to the new format, and saves
# the input-new.yaml file.
#
# Author: Allan Leal
# Date: 07 May 2024

import oyaml as yaml
import sys

# Get the input file name from argument
input_file = sys.argv[1]

output_file = input_file.replace(".yaml", "-new.yaml")

# Load the input file as yaml
with open(input_file, 'r') as file:
    data = yaml.safe_load(file)

if "ActivityModelParams" not in data:
    print("ActivityModelParams not found in the input file. No conversion needed for input file: ", input_file)
    sys.exit(1)

if "ExtendedUNIQUAC" not in data["ActivityModelParams"]:
    print("ExtendedUNIQUAC not found in the input file. No conversion needed for input file: ", input_file)
    sys.exit(1)

params = data["ActivityModelParams"]["ExtendedUNIQUAC"]

# The new r, q, and u blocks
rnew = {}
qnew = {}
unew = {}

if "r" in params:
    r = params["r"]
    assert isinstance(r, list), "Expecting the database file with r block as a list in input file: " + input_file
    for formula, value in r:
        rnew[formula] = value

if "q" in params:
    q = params["q"]
    assert isinstance(q, list), "Expecting the database file with q block as a list in input file: " + input_file
    for formula, value in q:
        qnew[formula] = value

if "u" in params:
    u = params["u"]
    assert isinstance(u, list), "Expecting the database file with u block as a list in input file: " + input_file
    for formula1, formula2, values in u:
        if formula1 not in unew:
            unew[formula1] = {}
        values = "@[" + ", ".join([f"{x:.8g}".replace("1e", "1.0e") for x in values]) + "]@"
        unew[formula1][formula2] = values

params["r"] = rnew
params["q"] = qnew
params["u"] = unew

# Save the new yaml file
with open(output_file, 'w') as file:
    # yaml.dump(data, file)

    # Dump the modified data to a string
    new_yaml = yaml.dump(data)
    new_yaml = new_yaml.replace("@'", "").replace("'@", "").replace('@"', "").replace('"@', "")
    new_yaml = new_yaml.replace("[", "[ ").replace("]", " ]")

    # Print the new yaml string
    file.write(new_yaml)
