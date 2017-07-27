#!/usr/bin/env python


"""
Author: Melanie van den Bosch
Master internship phylogenomics and dating of
salmonid trees
This script is to transform the BEAUTI.xml files
for each OG alignment.
"""

from sys import argv
import re

def get_xml_filename(argv):
    """ Returns the xml filename from the command line

    Keyword arguments:
        argv -- The list with arguments for this script
    Returns:
        xml_filename -- string, name of the file provided as argument
    """
    xml_filename = argv[1]
    return xml_filename

def get_ali_filename(argv):
    """ Returns the alignment filename from the command line

    Keyword arguments:
        argv -- The list with arguments for this script
    Returns:
        xml_filename -- string, name of the file provided as argument
    """
    ali_filename = argv[2]
    return ali_filename 

def parse_into_xml(xml_filename, ali_filename):
    full_seq_new = ''
    ID_pattern = re.compile(r'id="(OG\d+)"')
    new_ID = get_ali_ID(ali_filename)
    seq_pattern = re.compile(r'(\s+<sequence id="seq_)(.+)("\s+taxon=")(.+)' \
                             '("\s* totalcount="\d+)(" value=")(.*)("\/>)')
    with open(xml_filename, "rt") as xml_in:
        with open("xml_out.xml", "wt") as xml_out:
            for line in xml_in:
                seq_match = seq_pattern.match(line)
                ID_match = ID_pattern.match(line)
                if seq_match:
                    s_beg = seq_match.group(1)
                    s_name = seq_match.group(2)
                    s_tax = seq_match.group(3)
                    s_taxname = seq_match.group(4)
                    s_count = seq_match.group(5)
                    s_val = seq_match.group(6)
                    s_ali = seq_match.group(7)
                    s_end = seq_match.group(8)
                    for ali_name, ali_seq in parse_alignment_info(ali_filename).iteritems():
                        # add the alignment lines into a variable
                        full_seq_new += s_beg + ali_name + s_tax + ali_name + s_count + s_val + \
                                        ali_seq + s_end + '\n'
                        
                    xml_out.write(s_beg.replace(s_beg, full_seq_new))
                    full_seq_new = ''
                
                elif ID_match:
                    old_ID = ID_match.group(1)
                    full_ID_new = 'id="' + new_ID + '"\n'
                    xml_out.write(old_ID.replace(old_ID, full_ID_new))
                
                else:
                    xml_out.write(line)
                

def parse_alignment_info(ali_filename):
    """ Opens the alignment file

    Keyword arguments:
        fasta_file -- string, the filename of the .fa alignment
    Returns:
        ali_content -- string, name of the file provided as argument
    """
    ali_dict = {}
    name_pattern = re.compile(r">.+_(\D+\|\w+.\w)")

    with open(ali_filename) as ali:
        for line in ali:
            if line.startswith(">"):
                name_match = name_pattern.match(line)
                if name_match:
                    name = name_match.group(1)
                    ali_dict[name] = ''
            else:
                ali_dict[name] += line.strip()
            for key in ali_dict:
                ali_dict[key] = "".join(ali_dict[key])
                
        return ali_dict

def get_ali_ID(ali_filename):
    ID_pattern = re.compile(r".+(OG\d+).fa")
    ID_match = ID_pattern.match(ali_filename)
    if ID_match:
        new_ID = ID_match.group(1)
    return new_ID
                                   
                



if __name__ == "__main__":
    # Get input file names from cmd line arguments
    #xml_filename = get_xml_filename(argv)
    #ali_filename = get_ali_filename(argv)
    ali_filename = "C:/Users/meeldurb/Google Drive/Master internship" \
                    " phylogenetics salmonids/Salmonid_genomics_resources/" \
                    "Orthologs_homeologs/orthogroups.03.06.2017/Alignments/OG0008390.fa"
    xml = "dummy_xml/secondary_constr_aa.xml"
    parse_into_xml(xml, ali_filename)
    ali = parse_alignment_info(ali_filename)
    #print ali
    
    
