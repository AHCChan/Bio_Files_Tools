HELP_DOC = """
BED POSTMERGE UNCOLLAPSE
(version 1.0)
by Angelo Chan

This is a program for unpacking a BED file after:
    1)  A number of BED files had their first 3 columns duplicated, and then a
        column inserted at the very end of the file with a different number for
        each file, starting from 1.
    2)  The BED files were concatenated into a single BED file.
    3)  The concatenated BED file was sorted.
    4)  The sorted BED file then underwent BEDTools Merge, and all columns
        except the first 3 were collapsed.

After unpacking, the new BED file will have:
    1)  The first 3 columns will be the genomic coordinates of the merged
        regions. (1 region per row)
    2)  For each file which went into this, either:
            A)  A number of blank columns, equal to the number of columns in the
                original file, if no locus from that file became part of the
                merged region.
            B)  All the original data for a locus, if only one locus from that
                file became part of the merged region.
            C)  The original data for multiple loci, if multiple loci from that
                file became part of the merged region. The number of columns
                will be the same as in the original file, with the contents of
                different rows separated commas.

This program will also produce a summary file, which will have:
    1)  The first 3 columns will be the genomic coordinates of the merged
        regions. (1 region per row)
    2)  For each file which went into this, the number of loci from that file
        which went into merged region.



USAGE:
    
    python27 BED__Postmerge_Uncollapse.py <input_path>
            [-a <output_path_all>] [-s <output_path_summary>]



MANDATORY:
    
    input_file
        
        The filepath of the input BED file. No headers allowed.

OPTIONAL:
    
    output_path_complete
        
        (DEFAULT path generation available)
        
        The filepath of the output file which contains all of the contents of
        the original files.
    
    output_path_summary
        
        (DEFAULT path generation available)
        
        The filepath of the output file which contains a summary of how many
        loci from each file went into each merged region.



EXAMPLE:
    
    python27 BED__Postmerge_Uncollapse.py data\merged_sorted_concat_added.bed -c
            data\uncollapsed__COMPLETE.bed -s data\uncollapsed__SUMMARY.bed

USAGE:
    
    python27 BED__Postmerge_Uncollapse.py <input_path>
            [-a <output_path_complete>] [-s <output_path_summary>]
"""

NAME = "BED__Postmerge_Uncollapse.py"



# Configurations ###############################################################

AUTORUN = True

WRITE_PREVENT = False # Completely prevent overwritting existing files
WRITE_CONFIRM = True # Check to confirm overwritting existing files

PRINT_ERRORS = True
PRINT_PROGRESS = True
PRINT_METRICS = True



# Minor Configurations #########################################################

FILEMOD__ALL =   "__UNPACKED_ALL"
FILEMOD__SUMMARY =   "__UNPACKED_SUMMARY"



# Imported Modules #############################################################

import _Controlled_Print as PRINT
from _Command_Line_Parser import *

from Table_File_Reader import *



# Strings ######################################################################

STR__use_help = "\nUse the -h option for help:\n\t python "\
"BED__Postmerge_Uncollapse.py -h"



STR__invalid_ID = """
ERROR: Invalid file ID: {S}
Please ensure that before concatenation and BEDTools Merge, all files have ID
numbers inserted into the last row. The first file should have an ID number of 1
and the following files should increment by 1."""



STR__metrics = "Rows in file: {N}"



STR__uncollapse_begin = "\nRunning Uncollapse_Postmerge_BED..."

STR__uncollapse_complete = "\nUncollapse_Postmerge_BED successfully finished."



# Apply Globals ################################################################

PRINT.PRINT_ERRORS = PRINT_ERRORS
PRINT.PRINT_PROGRESS = PRINT_PROGRESS
PRINT.PRINT_METRICS = PRINT_METRICS



# Functions ####################################################################

def Uncollapse_BED(path_in, path_out_all, path_out_summary):
    """
    Generate a series of FASTA files each containing a synthetic chromosome.
    
    @path_in
            (str - filepath)
            The filepath of the input file.
    @path_out_all
            (str - dirpath)
            The filepath of the output file containing all the original data.
    @path_out_summary
            (str - dirpath)
            The filepath of the output file containing the summary.
    
    Return a value of 0 if the function runs successfully.
    Return a value of 1 if there is a problem.
    
    Uncollapse_BED(str, str, str) -> int
    """
    # Setup reporting
    row_count = 0
    
    # I/O setup
    f = Table_Reader()
    f.Set_New_Path(path_in)
    f.Set_Delimiter("\t")
    f.Open()
    f.Close()
    o1 = open(path_out_all, "w")
    o2 = open(path_out_summary, "w")
    
    PRINT.printP(STR__uncollapse_begin)
    
    # Other setup
    f.Open()
    f.Read()
    first = f.Get()
    column_count = len(first[3:-1])
    cc_range = range(column_count)
    placeholder = column_count * "\t"
    indexes = range(1, column_count+1)
    empty_dict = {}
    for i in indexes: empty_dict[i] = None
    count_dict = {}
    for i in indexes: count_dict[i] = 0
    f.Close()
    
    # Main loop
    f.Open()
    while not f.EOF:
        f.Read()
        row_count += 1
        values = f.Get()
        # Setup
        error = False
        dict_ = dict(empty_dict)
        counts = dict(count_dict)
        # Coordinates of the merged region
        merge_str = "\t".join(values[:3])
        # Unpack
        data = values[3:]
        data_split = [i.split(",") for i in data]
        length = len(data_split[0])
        range_ = range(length)
        for i in range_:
            try:
                ID = int(data_split[-1][i])
            except:
                f.Close()
                o1.close()
                o2.close()
                PRINT.printE(STR__invalid_ID.format(S=data_split[-1][i]))
                return 1
            if dict_[ID]:
                for j in cc_range:
                    s = "," + data_split[j][i]
                    dict_[ID][j] = dict_[ID][j] + s
            else:
                temp = []
                for j in cc_range:
                    temp.append(data_split[j][i])
                dict_[ID] = temp
            counts[ID] += 1
        # Write
        o1.write(merge_str)
        o2.write(merge_str)
        for i in indexes:
            if dict_[i]:
                s = "\t" + "\t".join(dict_[i])
                o1.write(s)
            else:
                o1.write(placeholder)
            s = "\t" + str(counts[i])
            o2.write(s)
        o1.write("\n")
        o2.write("\n")
    
    # Finish
    f.Close()
    o1.close()
    o2.close()
    PRINT.printP(STR__uncollapse_complete)
    
    # Reporting
    PRINT.printM(STR__metrics.format(N=row_count))
    
    # Wrap up
    return 0



# Command Line Parsing #########################################################

def Parse_Command_Line_Input__Uncollapse_BED(raw_command_line_input):
    """
    Parse the command line input and call the Generate_Synthetic_Genome function
    with appropriate arguments if the command line input is valid.
    """
    PRINT.printP(STR__parsing_args)
    # Remove the runtime environment variable and program name from the inputs
    inputs = Strip_Non_Inputs(raw_command_line_input, NAME)
    
    # No inputs
    if not inputs:
        PRINT.printE(STR__no_inputs)
        PRINT.printE(STR__use_help)
        return 1
    
    # Help option
    if inputs[0] in LIST__help:
        print(HELP_DOC)
        return 0
    
    # Initial validation (Redundant in current version)
    if len(inputs) < 1:
        PRINT.printE(STR__insufficient_inputs)
        PRINT.printE(STR__use_help)
        return 1
    
    # Validate mandatory inputs
    path_in = inputs.pop(0)
    valid = Validate_Read_Path(path_in)
    if valid == 1:
        PRINT.printE(STR__IO_error_read.format(f = path_in))
        PRINT.printE(STR__use_help)
        return 1
    
    # Set up rest of the parsing
    path_out_a = Generate_Default_Output_File_Path_From_File(path_in,
            FILEMOD__ALL, True)
    path_out_s = Generate_Default_Output_File_Path_From_File(path_in,
            FILEMOD__SUMMARY, True)
    
    # Validate optional inputs (except output path)
    while inputs:
        arg = inputs.pop(0)
        flag = 0
        try: # Following arguments
            if arg in ["-a", "-s"]:
                arg2 = inputs.pop(0)
            else: # Invalid
                arg = Strip_X(arg)
                PRINT.printE(STR__invalid_argument.format(s = arg))
                PRINT.printE(STR__use_help)
                return 1
        except:
            PRINT.printE(STR__insufficient_inputs)
            PRINT.printE(STR__use_help)
            return 1
        if arg == "-a":
            path_out_a = arg2
        else: # arg == "-s"
            path_out_s = arg2
    
    # Validate output paths
    valid_out = Validate_Write_Path(path_out_a)
    if valid_out == 2: return 0
    if valid_out == 3:
        printE(STR__IO_error_write_forbid)
        return 1
    if valid_out == 4:
        printE(STR__In_error_write_unable)
        return 1
    valid_out = Validate_Write_Path(path_out_s)
    if valid_out == 2: return 0
    if valid_out == 3:
        printE(STR__IO_error_write_forbid)
        return 1
    if valid_out == 4:
        printE(STR__In_error_write_unable)
        return 1
    
    # Run program
    exit_state = Uncollapse_BED(path_in, path_out_a, path_out_s)
    
    # Exit
    if exit_state == 0: return 0
    else: return 1


    
def Validate_Write_Path(filepath):
    """
    Validates the filepath of the input file.
    Return 0 if the filepath is writtable.
    Return 1 if the user decides to overwrite an existing file.
    Return 2 if the user declines to overwrite an existing file.
    Return 3 if the file exists and the program is set to forbid overwriting.
    Return 4 if the program is unable to write to the filepath specified.
    
    Validate_Write_Path(str) -> int
    """
    try:
        f = open(filepath, "U")
        f.close()
    except: # File does not exist. 
        try:
            f = open(filepath, "w")
            f.close()
            return 0 # File does not exist and it is possible to write
        except:
            return 4 # File does not exist but it is not possible to write
    # File exists
    if WRITE_PREVENT: return 3
    if WRITE_CONFIRM:
        confirm = raw_input(STR__overwrite_confirm.format(f=filepath))
        if confirm not in LIST__yes: return 2
    # User is not prevented from overwritting and may have chosen to overwrite
    try:
        f = open(filepath, "w")
        f.close()
        if WRITE_CONFIRM: return 1 # User has chosen to overwrite existing file
        return 0 # Overwriting existing file is possible
    except:
        return 4 # Unable to write to specified filepath



# Main Loop ####################################################################

if AUTORUN and (__name__ == "__main__"):
    exit_code = Parse_Command_Line_Input__Uncollapse_BED(sys.argv)
