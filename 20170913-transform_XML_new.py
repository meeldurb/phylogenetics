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
import os

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

def get_duplicates_filename(argv):
    """ Returns the filename of the duplicates from the command line

    Keyword arguments:
        argv -- The list with arguments for this script
    Returns:
        dup_filename -- string, name of the file provided as argument
    """
    dup_filename = argv[3]
    return dup_filename

def parse_into_xml(xml_filename, ali_filename, dup_filename):
    full_seq_new = ''
    MG3_1_lines = ''
    MG3_2_lines = ''
    ID_pattern = re.compile(r'(OG\d+_\d+\.)')
    new_ID = get_ali_ID(ali_filename)
    #print new_ID
    seq_pattern = re.compile(r'(\s+<sequence id="seq_)(.+)("\s+taxon=")(.+)' \
                             '("\s* totalcount="\d+)(" value=")(.*)("\/>)')
    MG1_1_pattern = re.compile(r'(\s+<taxon id=")(Omyk1)(" spec="Taxon"\/>)')
    MG1_2_pattern = re.compile(r'(\s+<taxon id=")(Ssal1)(" spec="Taxon"\/>)') 
    MG2_1_pattern = re.compile(r'(\s+<taxon id=")(Omyk2)(" spec="Taxon"\/>)') 
    MG2_2_pattern = re.compile(r'(\s+<taxon id=")(Ssal2)(" spec="Taxon"\/>)')
    MG3_1_pattern = re.compile(r'(\s+<taxon id=")(Eluc1)(" spec="Taxon"\/>)')
    MG3_2_pattern = re.compile(r'(\s+<taxon idref=")(Ssal3)("\/>)')
    
    dup_table = parse_dup_table(dup_filename)
    dup_ID = dup_table[str(new_ID)]
    #print dup_ID
    newdir = 'xml_outfiles_new/'
    if not os.path.exists(newdir):
        os.makedirs(newdir)
    xml_out_filename = newdir + new_ID + "xml"
    with open(xml_filename, "rt") as xml_in:
        with open(xml_out_filename, "wt") as xml_out:
            for line in xml_in:
                seq_match = seq_pattern.match(line)
                ID_match = ID_pattern.search(line)
                MG1_1_match = MG1_1_pattern.match(line)
                MG1_2_match = MG1_2_pattern.match(line)
                MG2_1_match = MG2_1_pattern.match(line)
                MG2_2_match = MG2_2_pattern.match(line)
                MG3_1_match = MG3_1_pattern.match(line)
                MG3_2_match = MG3_2_pattern.match(line)
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
                elif MG1_1_match:
                    p_name = MG1_1_match.group(2)
                    xml_out.write(line.replace(p_name, dup_ID[2]))
                elif MG1_2_match:
                    p_name = MG1_2_match.group(2)
                    xml_out.write(line.replace(p_name, dup_ID[1]))
                elif MG2_1_match:
                    p_name = MG2_1_match.group(2)
                    xml_out.write(line.replace(p_name, dup_ID[4]))
                elif MG2_2_match:
                    p_name = MG2_2_match.group(2)
                    xml_out.write(line.replace(p_name, dup_ID[3]))
                elif MG3_1_match:
                    p_beg = MG3_1_match.group(1)
                    p_name = MG3_1_match.group(2)
                    p_end = MG3_1_match.group(3)
                    for ali_name in parse_alignment_info(ali_filename):
                        spec_name = ali_name.split("|")[0]
                        if spec_name == "Eluc" or spec_name == "Okis" or spec_name == "Tthy":
                            MG3_1_lines += p_beg + ali_name + p_end + '\n'
                    xml_out.write(p_beg.replace(p_beg, MG3_1_lines))
                    MG3_lines = ''
                elif MG3_2_match:
                    p_beg = MG3_2_match.group(1)
                    p_name = MG3_2_match.group(2)
                    p_end = MG3_2_match.group(3)
                    for ali_name in parse_alignment_info(ali_filename):
                        spec_name = ali_name.split("|")[0]
                        if spec_name == "Omyk" or spec_name == "Ssal":
                            MG3_2_lines += p_beg + ali_name + p_end + '\n'
                    xml_out.write(p_beg.replace(p_beg, MG3_2_lines))
                    MG3_lines = ''
                elif ID_match:
                    old_ID = ID_match.group(1)
                    #print old_ID
                    #print new_ID
                    xml_out.write(line.replace(old_ID, new_ID))

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
    name_pattern = re.compile(r">(\D+\|\w+.\w+)")
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
    ID_pattern = re.compile(r".+(OG\d+_\d+\.)aln")
    ID_match = ID_pattern.match(ali_filename)
    if ID_match:
        new_ID = ID_match.group(1)
        return new_ID
                                   
def parse_dup_table(dup_filename):
    dup_dict = {}
    with open(dup_filename) as dup_table:
        for line in dup_table:
            clan_id, Eluc, Ssal_1, Omyk_1, Ssal_2, Omyk_2 = line.strip().split(';')
            #print clan_id
            dup_dict[clan_id] = [Eluc, Ssal_1, Omyk_1, Ssal_2, Omyk_2]
        return dup_dict



##def group_by_heading(xml_filename):
##    buffer = []
##    for line in xml_filename:
##        if line.startswith('                <taxonset id="monophyletic_group_1"'):
##            if buffer:
##                #yield buffer
##                buffer = [line]
##        else:
##            buffer.append(line)
##        yield buffer
##
##def try_out(xml_filename):
##    with open(xml_filename, 'r') as source:
##        for heading_and_lines in group_by_heading(source):
##            #print heading_and_lines
##            heading = heading_and_lines[0]
##            lines = heading_and_lines[1:3]
##        print heading
##        print lines
##        
            

    

if __name__ == "__main__":
    # Get input file names from cmd line arguments
    xml_filename = get_xml_filename(argv)
    ali_filename = get_ali_filename(argv)
    dup_filename = get_duplicates_filename(argv)

##    ali_filename = "C:/Users/meeldurb/Google Drive/" \
##                   "Master internship phylogenetics salmonids/" \
##                   "Salmonid_genomics_resources/" \
##                   "Orthologs_homeologs/orthogroups_2017/" \
##                   "cds_macse_nt_align/OG0000372_3.aln"
##    xml_filename = "dummy_xml/sec_constr_nt.xml"
##    dup_filename = "C:/Users/meeldurb/Dropbox/Melanie/" \
##                   "Master_internship_phylogenetics/" \
##                   "phylogenetics/RData/" \
##                   "20170913_duplicate_clans_filtered_final_nt.csv"
    #dupdict = parse_dup_table(dup_filename)
    parse_into_xml(xml_filename, ali_filename, dup_filename)


    
