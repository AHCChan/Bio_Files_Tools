HELP_DOC = """
CHROMSOME SIZES TO BINS
(version 1.0)
by Angelo Chan

This is a program for generating a series of bins from a chromosome sizes file.

The size of the bin can be specified, as can what to do with the final segment
of the chromsome not covered by a full-sized bin.



USAGE:
    
    python27 Chr_Sizes_to_Bins.py <chr_sizes_file> [-o <output_path>] [-b
            <bin_size>] [-t F|A|S|H]



MANDATORY:
    
    chr_sizes_file
        
        The filepath of the input chromosome sizes file. This is a TSV file
        with the names of the chromosomes in the first column, and their
        corresponding sizes in the second column.

OPTIONAL:
    
    output_path
        
        (DEFAULT path generation available)
        
        The filepath of the output file.
    
    bin_size
        
        (DEFAULT: 1000)
        
        The size of each bin.
    
    (-t)

        (DEFAULT: A)
        
        What to do with the tail-end of the chromsome which cannot be covered by
        a full-sized bin.
            
            {F}ull        - Full-sized bin. Insert a full-sized bin which more
                            than covers the last nucleotides specified by the
                            chromosome size.
            {A}ccurate    - Insert a smaller-than-standard bin which covers
                            exactly the last nucleotides specified by the
                            chromosome size.
            {S}kip        - Skip the nucleotides at the end. Don't have a bin
                            for them.
            {H}alf        - If the remaining nucleotides span at least half the
                            length of a bin, insert a full-sized bin.



EXAMPLES:
    
    python27 Chr_Sizes_to_Bins.py path\sizes.tsv
    
    python27 Chr_Sizes_to_Bins.py path\sizes.tsv -o path\\bins.bed
    
    python27 Chr_Sizes_to_Bins.py path\sizes.tsv -o path\\bins.bed -b 1000000
    
    python27 Chr_Sizes_to_Bins.py path\sizes.tsv -o path\\bins.bed -t S

USAGE:
    
    python27 Chr_Sizes_to_Bins.py <chr_sizes_file> [-o <output_path>] [-b
            <bin_size>] [-t F|A|S|H]
"""

NAME = "Chr_Sizes_to_Bins.py"



# Configurations ###############################################################

AUTORUN = True

WRITE_PREVENT = False # Completely prevent overwritting existing files
WRITE_CONFIRM = True # Check to confirm overwritting existing files

PRINT_ERRORS = True
PRINT_PROGRESS = True
PRINT_METRICS = True



# Minor Configurations #########################################################

FILEMOD = "__BINNED"



# Defaults #####################################################################
"NOTE: altering these will not alter the values displayed in the HELP DOC"

DEFAULT__bin_size = 1000
DEFAULT__tail_mode = 2 # Accurate



# Imported Modules #############################################################

import _Controlled_Print as PRINT
from _Command_Line_Parser import *

from Table_File_Reader import *



# Enums ########################################################################

class MODE:
        FULL=1
        ACCURATE=2
        SKIP=3
        HALF=4



# Strings ######################################################################

STR__use_help = "\nUse the -h option for help:\n\t python "\
"Chr_Sizes_to_Bins.py -h"



STR__invalid_bin_size = """
ERROR: Invalid bin size:
    {s}
Please specify a positive integer.
"""

STR__invalid_mode = """
ERROR: Invalid mode specified:
    {s}
Please specify one of the following:
    FULL
    ACCURATE
    SKIP
    HALF
"""

STR__metrics = """
    Size of Genome: {A}
    Number of Bins: {B}
"""



STR__binning_begin = "\nRunning Bin_Chromosome_Sizes..."

STR__binning_complete = "\nBin_Chromosome_Sizes successfully finished."



# Lists ########################################################################

LIST__full = ["FULL", "Full", "full", "F", "f"]
LIST__accurate = ["ACCURATE", "Accurate", "accurate", "ACC", "Acc", "acc", "A",
        "a"]
LIST__skip = ["SKIP", "Skip", "skip", "S", "s"]
LIST__half = ["HALF", "Half", "half", "H", "h"]



# Apply Globals ################################################################

PRINT.PRINT_ERRORS = PRINT_ERRORS
PRINT.PRINT_PROGRESS = PRINT_PROGRESS
PRINT.PRINT_METRICS = PRINT_METRICS



# Functions ####################################################################

def Bin_Chromosome_Sizes(path_sizes, path_out, bin_size, tail_mode):
    """
    Create a new set of values which would allow genomic coordinate data to be
    plotted linearly.
    
    @path_sizes
            (str - filepath)
            The filepath of the chromosome sizes file. This is a TSV file with
            the names of the chromosomes in the first column, and their
            corresponding sizes in the second column.
    @path_out
            (str - filepath)
            The filepath of the output file.
    @bin_size
            (int)
            The size of the bins to be produced.
    @tail_mode
            (int) - Pseudo ENUM
            How to treat the tail end of the chromsome which doesn't span the
            length of a full-sized bin. Acceptable options are:
                1:  Insert a full-sized bin which more than covers the last
                    nucleotides specified by the chromosome size.
                2:  Insert a smaller-than-standard bin which covers exactly the
                    last nucleotides specified by the chromosome size.
                3:  Skip the nucleotides at the end. Don't have a bin for them.
                4:  If the remaining nucleotides span at least half the length
                    of a bin, insert a full-sized bin.
                
    
    Return a value of 0 if the function runs successfully.
    Return a value of 1 if there is a problem.
    
    Bin_Chromosome_Sizes(str, str, int, int) -> int
    """
    # Setup reporting
    genome_size = 0
    bin_count = 0
    
    # I/O setup
    f = Table_Reader()
    f.Set_New_Path(path_sizes)
    f.Set_Delimiter("\t")
    f.Open()
    f.Close()
    o = open(path_out, "w")
    
    PRINT.printP(STR__binning_begin)
    
    # Main loop
    f.Open()
    while not f.EOF:
        f.Read()
        chr_name = f[0]
        chr_size = f[1]
        chr_size = int(chr_size)
        genome_size += chr_size
        current = 0
        difference = chr_size - current
        write = True
        while difference > 0:
            if difference >= bin_size:
                start = current + 1
                end = current + bin_size
                current = end
                difference = chr_size - current
            else: # End of chromosome
                if tail_mode == MODE.FULL:
                    start = current + 1
                    end = current + bin_size
                elif tail_mode == MODE.ACCURATE:
                    start = current + 1
                    end = chr_size
                elif tail_mode == MODE.SKIP:
                    write = False
                elif tail_mode == MODE.HALF:
                    if difference > (bin_size/2):
                        start = current + 1
                        end = current + bin_size
                    else:
                        write = False
                else: # Bad Tail_Mode
                    f.Close()
                    o.close()
                    return 1
                difference = 0
            if write:
                bin_count += 1
                start = str(start)
                end = str(end)
                sb = chr_name + "\t" + start + "\t" + end + "\n"
                o.write(sb)
    
    # Finish
    f.Close()
    o.close()
    PRINT.printP(STR__binning_complete)
    
    # Reporting
    PRINT.printM(STR__metrics.format(A=genome_size, B=bin_count))
    
    # Wrap up
    return 0



# Command Line Parsing #########################################################

def Parse_Command_Line_Input__Chr_Sizes_to_Bins(raw_command_line_input):
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
    if len(inputs) < 1:
        PRINT.printE(STR__insufficient_inputs)
        PRINT.printE(STR__use_help)
        return 1
    
    # Validate mandatory inputs
    path_sizes = inputs.pop(0)
    valid = Validate_Read_Path(path_sizes)
    if valid == 1:
        PRINT.printE(STR__IO_error_read.format(f = path_in))
        PRINT.printE(STR__use_help)
        return 1
    
    # Set up rest of the parsing
    path_out = Generate_Default_Output_File_Path_From_File(path_sizes, FILEMOD,
            True)
    bin_size = DEFAULT__bin_size
    tail_mode = DEFAULT__tail_mode
    
    # Validate optional inputs (except output path)
    while inputs:
        arg = inputs.pop(0)
        flag = 0
        try: # Following arguments
            if arg in ["-o", "-b", "-t"]:
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
        elif arg == "-b":
            bin_size = Validate_Int_Positive(arg2)
            if bin_size == -1:
                PRINT.printE(STR__invalid_bin_size.format(s = arg2))
                PRINT.printE(STR__use_help)
                return 1
        else: # arg == "-t"
            if arg2 in LIST__full:
                tail_mode = MODE.FULL
            elif arg2 in LIST__accurate:
                tail_mode = MODE.ACCURATE
            elif arg2 in LIST__skip:
                tail_mode = MODE.SKIP
            elif arg2 in LIST__half:
                tail_mode = MODE.HALF
            else:
                PRINT.printE(STR__invalid_mode.format(s = arg2))
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
    exit_state = Bin_Chromosome_Sizes(path_sizes, path_out, bin_size, tail_mode)
    
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
    exit_code = Parse_Command_Line_Input__Chr_Sizes_to_Bins(sys.argv)


