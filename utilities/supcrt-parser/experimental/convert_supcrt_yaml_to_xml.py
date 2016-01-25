import argparse, sys
from lxml import etree as xml

###############################################################################
# The following is needed to ensure that PyYAML uses OrderedDict instead of
# regular dict. This is needed to preserve the order of the YAML elements.
# This workaround is given at:
# http://stackoverflow.com/questions/5121931/in-python-how-can-you-load-yaml-
# mappings-as-ordereddicts/21048064#21048064
###############################################################################
import yaml
from yaml.representer import Representer
from yaml.constructor import Constructor, MappingNode, ConstructorError
import collections

_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG

def dict_representer(dumper, data):
    return dumper.represent_dict(data.iteritems())

def dict_constructor(loader, node):
    return collections.OrderedDict(loader.construct_pairs(node))

yaml.add_representer(collections.OrderedDict, dict_representer)
yaml.add_constructor(_mapping_tag, dict_constructor)
###############################################################################
from collections import OrderedDict


def xmlizer(key, value, root):
    if type(value) is list:
        children = []
        for entry in value:
            xmlizer(key, entry, root)
    elif type(value) in [OrderedDict, dict]:
        child = xml.Element(key)
        if value.has_key('value') and value.has_key('units'):
            child.text = str(value['value'])
            child.attrib['units'] = str(value['units'])
        else:
            for subkey, subvalue in value.iteritems():
                xmlizer(subkey, subvalue, child)
        root.append(child)
    else:
        child = xml.Element(key)
        child.text = str(value)
        root.append(child)


def main():
    # Create a command-line argument parser
    parser = argparse.ArgumentParser(prog='ReaktoroConverter')

    # Add the input argument
    parser.add_argument('input', type=str, \
        help='the relative path of the YAML database file, including its name')

    # Add the output argument (optional)
    parser.add_argument('output', type=str, nargs='?', \
        help='the relative path of the XML output file, including its name')

    # Parse the command-line arguments (remove the first argument, which is the name of this file
    args = parser.parse_args(sys.argv[1:])

    # Provide a default output file name if none was provided
    if args.output is None:
        args.output = os.path.splitext(args.input)[0] + '.xml'

    yaml_database = file(args.input)

    dictionary = yaml.load(yaml_database)

    doc = xml.Element('Database')

    for key, value in dictionary.iteritems():
        xmlizer(key, value, doc)

    root = doc.getroottree()
    root.write(args.output, xml_declaration=True, pretty_print=True)


if __name__ == '__main__':
    main()
