HELP_DOC = """
BED FILE CHROMOSOME NORMALIZER
(version 2.0)
by Angelo Chan

This is a program for normalizing genomic coordinates (chromosome name, start
position, end position) in a way which scales all chromosomes to the same size,
and the relative positions of all coordinates within that chromosome are
preserved.

The size (in normalized terms) is also calculated.

For example, with the following chromosome sizes:
    chr1    1000
    chr2    2000
    chr3    4000

The following coordinates:
    chr1    100     200
    chr1    500     700
    chr2    100     200
    chr2    500     700
    chr2    1000    1500
    chr3    2000    2400

..., scaling to 100, will produce the following coordinates:
    chr1    10  20
    chr1    50  70
    chr2    5   10
    chr2    25  35
    chr2    50  75
    chr3    50  60

Including the midpoint and size, the new columns will be:
    chr1    10  20  15  10
    chr1    50  70  60  20
    chr2    5   10  7   5
    chr2    25  35  30  10
    chr2    50  75  67  25
    chr3    50  60  55  10

The goal is to be able to graph the new coordinates in a straightforward manner
which also allows a sort of percentile comparison.

The three new columns will be added to the end of the existing data. You can use
Table_to_Table.py to rearrange the columns if you want to replace the original
coordinate columns with the new ones.



USAGE:
    
    python27 BED_Normalize_Chromosomes <input_BED> <chr_sizes_file>
            [-o <output_path>] [-t I|D <target>]



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
    
    I|D
        
        (DEFAULT: D)
        
        Whether to use integers or decimals for the normalized coordinates.
            
            I   Integers
            D   Decimals
    
    target
        
        (DEFAULT: 1000000)
        
        What number to normalize the chromosome sizes to.
        
        Note that using a low number may distort the results slightly due to
        decimal places. More so if using the [Integers] setting.



EXAMPLES:
    
    python27 BED_Normalize_Chromosomes path\data.bed path\sizes.tsv

    python27 BED_to_Linear_Coordinates path\data.bed path\sizes.tsv
            -o path\data_graphable.bed -t D 1

USAGE:
    
    python27 BED_Normalize_Chromosomes <input_BED> <chr_sizes_file>
            [-o <output_path>] [-t I|D <target>]
"""

NAME = "BED_Normalize_Chromosomes.py"



# Configurations ###############################################################

AUTORUN = True

WRITE_PREVENT = False # Completely prevent overwritting existing files
WRITE_CONFIRM = True # Check to confirm overwritting existing files

PRINT_ERRORS = True
PRINT_PROGRESS = True
PRINT_METRICS = True



# Minor Configurations #########################################################

FILEMOD = "__NORMALIZED"



# Defaults #####################################################################
"NOTE: altering these will not alter the values displayed in the HELP DOC"

DEFAULT__target_type = "DECIMALS"
DEFAULT__target = 1000000



# Imported Modules #############################################################

import _Controlled_Print as PRINT
from _Command_Line_Parser import * # 1.8

from Table_File_Reader import *



# Strings ######################################################################

STR__use_help = "\nUse the -h option for help:\n\t python "\
"BED_Normalize_Chromosomes.py -h"



STR__invalid_target_int = """
ERROR: Invalid target: {s}
Please specify a positive integer.
"""

STR__invalid_target_float = """
ERROR: Invalid target: {s}
Please specify a positive number.
"""



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

STR__missing_chr = """
ERROR: Chromosome not found in sizes file:
    {s}
"""



STR__metrics = "Rows in file: {N}"



STR__linearize_begin = "\nRunning BED_Normalize_Chr..."

STR__linearize_complete = "\nBED_Normalize_Chr successfully finished."



# Apply Globals ################################################################

PRINT.PRINT_ERRORS = PRINT_ERRORS
PRINT.PRINT_PROGRESS = PRINT_PROGRESS
PRINT.PRINT_METRICS = PRINT_METRICS



# Functions ####################################################################

def BED_Normalize_Chr(path_BED, path_sizes, path_out, target_type, target):
    """
    Generate a series of FASTA files each containing a synthetic chromosome.
    
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
    @target_type
            (str)
            What kind of number (integer/decimal) the normalized coordinates
            should be.
    @target
            (int)
            The size all chromosomes should be normalized to.
    
    Return a value of 0 if the function runs successfully.
    Return a value of 1 if there is a problem.
    
    BED_Normalize_Chr(str, str, str, str, int) -> int
    """
    # Setup reporting
    row_count = 0
    
    # Get displacements
    chr_sizes = Get_Chr_Sizes(path_sizes, target_type)
    if not chr_sizes: return 1
    
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
        if chr_ not in chr_sizes:
            PRINT.printE(STR__missing_chr.format(s = chr_))
            return 2
        try:
            start = int(start)
            end = int(end)
            if start > end:
                raise Exception
        except:
            PRINT.printE(STR__invalid_data.format(s = raw))
            return 3
        # Process
        sc_start, sc_end, sc_mid, size = Normalize_Coords(chr_, start, end,
                chr_sizes, target)
        # String
        sc_start = str(sc_start)
        sc_end = str(sc_end)
        sc_mid = str(sc_mid)
        size = str(size)
        sb = (raw + "\t" + chr_ + "\t" + sc_start + "\t" + sc_end + "\t" +
                sc_mid + "\t" + size + "\n")
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



def Get_Chr_Sizes(path_sizes, target_type):
    """
    Return a dictionary of the chromosome sizes.
    Return an empty dictionary if there is a problem.
    
    @path_sizes
            (str - filepath)
            The filepath of the chromosome sizes file. This is a TSV file with
            the names of the chromosomes in the first column, and their
            corresponding sizes in the second column. 
    @target_type
            (str)
            What kind of number (integer/decimal) the normalized coordinates
            should be.
    
    Get_Chr_Sizes(str) -> dict<str:int>
    """
    # Setup
    results = {}
    # Process
    f = open(path_sizes, "U")
    for line in f:
        values = line.split("\t")
        if values:
            if values[-1][-1] == "\n":
                values[-1] = values[-1][:-1]
            chr_ = values[0]
            # Get Size
            try:
                size_str = values[1]
            except:
                PRINT.printE(STR__no_size_str.format(s = chr_))
                return {}
            try:
                if target_type in LIST__integer:
                    size = int(size_str)
                elif target_type in LIST__decimal:
                    size = float(size_str)
                else:
                    raise Exception
            except:
                PRINT.printE(STR__invalid_size_str.format(s = size_str))
                return []
            # Append
            results[chr_] = size
    return results

def Normalize_Coords(chr_, start, end, chr_sizes, target):
    """
    Return a set of coordinates, normalized. Include the midpoint and size of
    the element in normalized terms.
    
    @chr_
            (str)
            The name of the chromosome.
    @start
            (int)
            The nucleotide number of the start (5') of the genetic element.
    @end
            (int)
            The nucleotide number of the end (3') of the genetic element.
    @chr_sizes
            (dict<str:int>)
            The dictionary of how big each chromosome is.
    @target
            (int)
            The size all chromosomes should be normalized to.
    
    Normalize_Coords(str,int,int,dict<str:int>,str,int/float) ->
            [str,int/float,int/float,int/float]
    """
    chr_size = chr_sizes[chr_]
    scaled_start = (start*target)/chr_size
    scaled_end = (end*target)/chr_size
    mid = (scaled_start + scaled_end) / 2
    size = scaled_end - scaled_start
    return [scaled_start, scaled_end, mid, size]



# Command Line Parsing #########################################################

def Parse_Command_Line_Input__BED_Normalize_Chr(raw_command_line_input):
    """
    Parse the command line input and call the BED_Normalize_Chr function with
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
        PRINT.printE(STR__IO_error_read.format(f = path_BED))
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
    target_type = DEFAULT__target_type
    target = DEFAULT__target
    
    # Validate optional inputs (except output path)
    while inputs:
        arg = inputs.pop(0)
        flag = 0
        try: # Following arguments
            if arg in ["-o"]:
                arg2 = inputs.pop(0)
            elif arg in ["-t"]:
                arg2 = inputs.pop(0)
                arg3 = inputs.pop(0)
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
        else: # arg == "-t"
            target_type = arg2
            if arg2 in LIST__integer:
                target = Validate_Int_Positive(arg3)
                if target == -1:
                    PRINT.printE(STR__invalid_target_int.format(s = arg3))
                    PRINT.printE(STR__use_help)
                    return 1
            elif arg2 in LIST__decimal: # Decimal
                target = Validate_Float_Positive(arg3)
                if target == -1:
                    PRINT.printE(STR__invalid_target_float.format(s = arg3))
                    PRINT.printE(STR__use_help)
                    return 1
            else:
                PRINT.printE(STR__invalid_int_dec.format(s = arg2))
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
    exit_state = BED_Normalize_Chr(path_BED, path_sizes, path_out, target_type,
            target)
    
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
    exit_code = Parse_Command_Line_Input__BED_Normalize_Chr(sys.argv)


