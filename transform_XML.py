#!/usr/bin/env python


"""
Author: Melanie van den Bosch
Master internship phylogenomics and dating of
salmonid trees
This script is to transform the BEAUTI.xml files
for each OG alignment.
"""

from sys import argv

def get_xml_filename(argv):
    xml_filename = argv[1]
    return xml_filename

def get_ali_filename(argv):
    ali_filename = argv[2]
    return ali_filename 

def parse_into_xml(xml_filename):
    with open(xml_filename) as xml:
        for line in xml:
            print line

def read_alignment(fasta_file):
    print "jemo"



if __name__ == "__main__":
    # Get input file names from cmd line arguments
    xml_filename = get_xml_filename(argv)
    #ali_filename = get_ali_filename(argv)

    # 
    #xml = "dummy_xml/secondary_constr_aa.xml"
    parse_into_xml(xml_filename)
