HELP_DOC = """
BED TO LINEAR COORDINATES
(version 2.0)
by Angelo Chan

This is a program for converting genomic coordinates (chromosome name, start
position, end position) into "linear" coordinates which can be graphed easily,
by adding the size of all "previous" chromosomes.

For example, with the following chromosome sizes:
    chr1    1000
    chr2    2000
    chr3    3000

The following coordinates:
    chr1    4000    4002
    chr1    5000    5004
    chr2    6000    6006
    chr2    7000    7008
    chr3    8000    8002
    chr3    9000    9004

... will produce:
    1   4000    4002    4001
    1   5000    5004    5002
    2   7000    7006    7003
    2   8000    8008    8004
    1   11000   11002   11001
    1   12000   12004   12002

The original coordinates have been increased by the size of the "previous"
chromosomes, and an adjusted midpoint coordinate will also be produced.

A new column is also created which alternates between 1 and 2 in accordance with
the chromosome, so that no two "adjacent" chromosomes will have the same
alternator number.

The newly created columns will be inserted to the right of all existing columns.

The chromsomes should be listed from top to bottom, in the order in which they
should eventually be graphed in subsequent analyses, from left to right.



USAGE:
    
    python27 BED_to_Linear_Coordinates <input_BED> <chr_sizes_file>
            [-o <output_path>] [-a <alternating_numbers>]



MANDATORY:
    
    input_BED
        
        The filepath of the input BED file. No headers allowed.
    
    chr_sizes_file
        
        The filepath of the input chromosome sizes file. This is a TSV file
        with the names of the chromosomes in the first column, and their
        corresponding sizes in the second column.

OPTIONAL:
    
    output_path
        
        (DEFAULT path generation available)
        
        The filepath of the output file.
    
    alternating_numbers
        
        (DEFAULT: 0,1)
        
        2 or more comma-separated numbers which will be assigned, alternatingly
        to "adjacent" chromosomes such that no 2 adjacent chromsomes will have
        the same number.
        
        The intended use of these numbers is for designating different GEOM
        shapes when the data is graphed using GGPLOT2 in R. Therefore, they
        should be integers between 0 and 25, inclusive.



EXAMPLES:
    
    python27 BED_to_Linear_Coordinates path\data.bed path\sizes.tsv

    python27 BED_to_Linear_Coordinates path\data.bed path\sizes.tsv
            -o path\data_graphable.bed -a 15,16,17,18

USAGE:
    
    python27 BED_to_Linear_Coordinates <input_BED> <chr_sizes_file>
            [-o <output_path>] [-a <alternating_numbers>]
"""

NAME = "BED_to_Linear_Coordinates.py"



# Configurations ###############################################################

AUTORUN = True

WRITE_PREVENT = False # Completely prevent overwritting existing files
WRITE_CONFIRM = True # Check to confirm overwritting existing files

PRINT_ERRORS = True
PRINT_PROGRESS = True
PRINT_METRICS = True



# Minor Configurations #########################################################

FILEMOD = "__LINEARIZED"



# Defaults #####################################################################
"NOTE: altering these will not alter the values displayed in the HELP DOC"

DEFAULT__alt_nums = [0,1]



# Imported Modules #############################################################

import _Controlled_Print as PRINT
from _Command_Line_Parser import *

from Table_File_Reader import *



# Strings ######################################################################

STR__use_help = "\nUse the -h option for help:\n\t python "\
"BED_to_Linear_Coordinates.py -h"



STR__no_size_str = """
ERROR: No size given for the chromosome: {s}
Please check your chromosome sizes file.
"""

STR__invalid_size_str = """
ERROR: Invalid chromosome size: {s}
Please check your chromosome sizes file.
"""

STR__invalid_data = """
ERROR: Invalid data in your BED file:
    {s}
"""

STR__invalid_alt_nums = """
ERROR: Invalid alternating numbers:
    {s}
Please specify at least 2integers, separated by commas.
"""


STR__metrics = "Rows in file: {N}"



STR__linearize_begin = "\nRunning BED_to_Linear_Coordinates..."

STR__linearize_complete = "\nBED_to_Linear_Coordinates successfully finished."



# Apply Globals ################################################################

PRINT.PRINT_ERRORS = PRINT_ERRORS
PRINT.PRINT_PROGRESS = PRINT_PROGRESS
PRINT.PRINT_METRICS = PRINT_METRICS



# Functions ####################################################################

def BED_to_Linear(path_BED, path_sizes, path_out, alt_numbers):
    """
    Create a new set of values which would allow genomic coordinate data to be
    plotted linearly.
    
    @path_in
            (str - filepath)
            The filepath of the input BED file.
    @path_sizes
            (str - filepath)
            The filepath of the chromosome sizes file. This is a TSV file with
            the names of the chromosomes in the first column, and their
            corresponding sizes in the second column.
    @path_out
            (str - filepath)
            The filepath of the output file.
    @alt_numbers
            (list<int>)
            2 or more comma-separated numbers which will be assigned,
            alternatingly to "adjacent" chromosomes such that no 2 adjacent
            chromsomes will have the same number.
    
    Return a value of 0 if the function runs successfully.
    Return a value of 1 if there is a problem.
    
    BED_to_Linear(str, str, str, list<int>) -> int
    """
    # Setup reporting
    row_count = 0
    
    # Get displacements
    processed = Process_Path_Sizes(path_sizes, alt_numbers)
    if not processed: return 1
    displacements, chr_nums = processed
    
    # I/O setup
    f = Table_Reader()
    f.Set_New_Path(path_BED)
    f.Set_Delimiter("\t")
    f.Open()
    f.Close()
    o = open(path_out, "w")
    
    PRINT.printP(STR__linearize_begin)
    
    # Main loop
    f.Open()
    while not f.EOF:
        f.Read()
        row_count += 1
        values = f.Get()
        raw = f.Get_Raw()
        raw = raw[:-1]
        # Get values
        chr_ = values[0]
        start = values[1]
        end = values[2]
        try:
            start = int(start)
            end = int(end)
            if start > end:
                raise Exception
        except:
            PRINT.printE(STR__invalid_data.format(s = raw))
            return 1
        # Process
        difference = end - start
        mid = start + difference
        offset = displacements[chr_]
        os_start = start + offset
        os_end = end + offset
        os_mid = mid + offset
        alt_chr = chr_nums[chr_]
        # String
        os_start = str(os_start)
        os_end = str(os_end)
        os_mid = str(os_mid)
        alt_chr = str(alt_chr)
        sb = (raw + "\t" + alt_chr + "\t" + os_start + "\t" + os_end + "\t" +
                os_mid + "\n")
        # Write
        o.write(sb)
    
    # Finish
    f.Close()
    o.close()
    PRINT.printP(STR__linearize_complete)
    
    # Reporting
    PRINT.printM(STR__metrics.format(N=row_count))
    
    # Wrap up
    return 0

def Process_Path_Sizes(path_sizes, alt_numbers):
    """
    Return a dictionary with the "displacement" values for each chromosome, and
    a dictionary with the "alternating number" for each chromosome.
    
    Return an empty list if there are any problems.
    
    Get_Displacements(str, list<int>) -> []
    Get_Displacements(str, list<int>) -> [dict<str:int>, dict<str:int>]
    """
    # Setup
    displacement = {}
    assigned_nums = {}
    numbers = list(alt_numbers)
    current_displacement = 0
    # Process
    f = open(path_sizes, "U")
    for line in f:
        values = line.split("\t")
        if values:
            if values[-1][-1] == "\n":
                values[-1] = values[-1][:-1]
            chr_ = values[0]
            try:
                size_str = values[1]
            except:
                PRINT.printE(STR__no_size_str.format(s = chr_))
                return []
            try:
                size = int(size_str)
            except:
                PRINT.printE(STR__invalid_size_str.format(s = size_str))
                return []
            # Update
            current_no = numbers.pop(0)
            assigned_nums[chr_] = current_no
            numbers.append(current_no) # Essentially rotate to the end
            displacement[chr_] = current_displacement
            current_displacement += size
    return [displacement, assigned_nums]



# Command Line Parsing #########################################################

def Parse_Command_Line_Input__BED_to_Linear(raw_command_line_input):
    """
    Parse the command line input and call the BED_to_Linear function with
    appropriate arguments if the command line input is valid.
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
    
    # Initial validation
    if len(inputs) < 2:
        PRINT.printE(STR__insufficient_inputs)
        PRINT.printE(STR__use_help)
        return 1
    
    # Validate mandatory inputs
    path_BED = inputs.pop(0)
    valid = Validate_Read_Path(path_BED)
    if valid == 1:
        PRINT.printE(STR__IO_error_read.format(f = path_in))
        PRINT.printE(STR__use_help)
        return 1
    path_sizes = inputs.pop(0)
    valid = Validate_Read_Path(path_sizes)
    if valid == 1:
        PRINT.printE(STR__IO_error_read.format(f = path_in))
        PRINT.printE(STR__use_help)
        return 1
    
    # Set up rest of the parsing
    path_out = Generate_Default_Output_File_Path_From_File(path_BED, FILEMOD,
            True)
    alt_nums = DEFAULT__alt_nums
    
    # Validate optional inputs (except output path)
    while inputs:
        arg = inputs.pop(0)
        flag = 0
        try: # Following arguments
            if arg in ["-o", "-a"]:
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
        if arg == "-o":
            path_out = arg2
        else: # arg == "-a"
            alt_nums = Validate_List_Of_Ints_NonNeg(arg2, ",")
            if len(alt_nums) < 2:
                PRINT.printE(STR__invalid_alt_nums.format(s = arg2))
                PRINT.printE(STR__use_help)
                return 1
    
    # Validate output paths
    valid_out = Validate_Write_Path(path_out)
    if valid_out == 2: return 0
    if valid_out == 3:
        printE(STR__IO_error_write_forbid)
        return 1
    if valid_out == 4:
        printE(STR__In_error_write_unable)
        return 1
    
    # Run program
    exit_state = BED_to_Linear(path_BED, path_sizes, path_out, alt_nums)
    
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
    exit_code = Parse_Command_Line_Input__BED_to_Linear(sys.argv)


