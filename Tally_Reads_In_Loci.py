HELP_DOC = """
TALLY READS IN LOCI
(version 1.0)
by Angelo Chan

This is a program which takes a file containing loci of interest (BED file), a
data file containing alignment data (BED file), and tallying the number of reads
in each locus of interest.

This program also requires that both the loci file and the data file are sorted.
This program also requires a chromosome order file (TSV file) in which the
chromosomes are listed in the same order as they are sorted. A typical
chromosome sizes file will suffice.

(NOTE: SAM files can be converted into BED files using BED tools.)



USAGE:
    
    python27 Tally_Reads_In_Loci.py <loci_file> <chr_order_file> <data_file>
            [-o <output_path>] [-m <method>] [-t <threshold>]



MANDATORY:
    
    loci_file
        
        (BED file - BED3)
        The filepath of the input loci file. No headers allowed.
    
    chr_order_file
        
        (TSV file)
        The filepath of the chromosome order file. This can be a chromosome
        sizes file. Only the first column matters. No headers allowed.
    
    data_file
        
        (BED file - BED3)
        The filepath of the input data file. No headers allowed.

OPTIONAL:
    
    output_path
        
        (DEFAULT path generation available)
        
        The filepath of the output file.
    
    method
        
        (DEFAULT: M)
        
        How the score for each locus is calculated. The method also covers what
        to do when a read (data file) overlaps with more than one loci. Acceptable
        options are:
            
            {C}omplete    - The read contributes a score of 1 if it overlaps
                            completely with the locus.
            {F}raction    - The score contributed by a read to a loci is the
                            percentage of the read which overlaps.
                            (Ex. If the read is 100bp, and 75bp overlaps with
                            the loci of interest, it increases the loci score by
                            0.75.)
            {M}ajority    - Each read contributes a score of 1 if more than half
                            of the read overlaps with the locus.
            {N}ucleotides - Each read contributes a score to a locus equal to
                            the number of nucleotides which overlap with said
                            locus.
            {P}resent     - Each read contributes a score of 1 if there is any
                            overlap at all.
            {T}hreshold   - Each read contributes a score of 1 if it exceeds a
                            given threshold. If the threshold given is an
                            integer, the threshold required is the number of
                            overlapping nucleotides. If the threshold given is
                            a decimal number, the threshold required is the
                            percentage of the read which overlaps with the
                            locus.
        
        NOTE: Reads may contribute more than 1 point per read or more points
        than there are nucleotides in the read, if loci overlap. If this will be
        an issue, please ensure the integrity of your loci data. (Ex., by using
        BEDTools merge)
    
    threshold

        (DEFAULT: 0.5)
        
        Only relevant if the [threshold] option is selected for <method>.
        
        If the threshold given is an integer, the threshold required is the
        number of overlapping nucleotides. If the threshold given is a decimal
        number, the threshold required is the percentage of the read which
        overlaps with the locus.



EXAMPLES SCENARIO EXPLANATION:
    
    1:
    Basic use-case.
    
    2:
    Basic use-case with an output file specified.
    
    3:
    A count of the number of nucleotides which overlap with each locus.
    
    4:
    A count of the number of reads, for each locus, which overlaps with at least
    20bp.

EXAMPLES:
    
    python27 BED_to_Linear_Coordinates.py path\loci.bed path\chr_sizes.tsv
            path\aligned_reads.bed
    
    python27 BED_to_Linear_Coordinates.py path\loci.bed path\chr_sizes.tsv
            path\aligned_reads.bed -o path\loci_counts.bed
    
    python27 BED_to_Linear_Coordinates.py path\loci.bed path\chr_sizes.tsv
            path\aligned_reads.bed -o path\loci_nuc_numbers.bed -m N
    
    python27 BED_to_Linear_Coordinates.py path\loci.bed path\chr_sizes.tsv
            path\aligned_reads.bed -o path\loci_kmer20.bed -m T -t 20

USAGE:
    
    python27 Tally_Reads_In_Loci.py <loci_file> <chr_order_file> <data_file>
            [-o <output_path>] [-m <method>] [-t <threshold>]
"""

NAME = "Tally_Reads_In_Loci.py"



# Configurations ###############################################################

AUTORUN = True

WRITE_PREVENT = False # Completely prevent overwritting existing files
WRITE_CONFIRM = True # Check to confirm overwritting existing files

PRINT_ERRORS = True
PRINT_PROGRESS = True
PRINT_METRICS = True



# Minor Configurations #########################################################

FILEMOD = "__LOCI_TALLIES"



# Defaults #####################################################################
"NOTE: altering these will not alter the values displayed in the HELP DOC"

DEFAULT__method = 3 # METHOD.MAJORITY
DEFAULT__threshold = 0.5



# Imported Modules #############################################################

import _Controlled_Print as PRINT
from _Command_Line_Parser import * # 2.0

from Multitrack_BED_Reading_Coordinator import * #1.0



# Enums ########################################################################

class METHOD:
    COMPLETE=1
    FRACTION=2
    MAJORITY=3
    NUCLEOTIDES=4
    PRESENT=5
    THRESHOLD=6



# Strings ######################################################################

STR__use_help = "\nUse the -h option for help:\n\t python "\
"Tally_Reads_In_Loci.py -h"



STR__invalid_method = """
ERROR: Invalid method specified:
    {s}
Please specify one of:
    COMPLETE
    FRACTION
    MAJORITY
    NUCLEOTIDES
    PRESENT
    THRESHOLD
"""

STR__invalid_threshold = """
ERROR: Invalid threshold specified:
    {s}
Please specify either a positive integer, or a decimal number between 1 and 0.
"""

STR__metrics = """
                Loci: {A}
       Data Entriess: {B}
    
    Relevant Entries: {C}
         Total Score: {D}
       Average Score: {E}
"""



STR__tally_begin = "\nRunning Tally_The_Reads_In_Each_Loci..."

STR__tally_complete = "\nTally_The_Reads_In_Each_Loci successfully finished."



# Lists ########################################################################

LIST__complete = ["C", "c", "COMPLETE", "Complete", "complete", "ALL", "All",
        "all"]
LIST__fraction = ["F", "f", "FRACTION", "Fraction", "fraction"]
LIST__majority = ["M", "m", "MAJORITY", "Majority", "majority", "MOST", "Most",
        "most"]
LIST__nucleotides = ["N", "n", "NUCLEOTIDES", "Nucleotides", "nucleotides",
        "NUC", "Nuc", "nuc"]
LIST__present = ["P", "p", "PRESENT", "Present", "present"]
LIST__threshold = ["T", "t", "THRESHOLD", "Threshold", "threshold"]



# Apply Globals ################################################################

PRINT.PRINT_ERRORS = PRINT_ERRORS
PRINT.PRINT_PROGRESS = PRINT_PROGRESS
PRINT.PRINT_METRICS = PRINT_METRICS



# Functions ####################################################################

def Tally_The_Reads_In_Each_Loci(path_loci, path_chrs, path_data, path_output,
            method, threshold):
    """
    Create a new set of values which would allow genomic coordinate data to be
    plotted linearly.
    
    @path_loci
            (str - filepath)
            The filepath of the input loci file. No headers allowed.
    @path_chrs
            (str - filepath)
            The filepath of the chromosome order file. This can be a chromosome
            sizes file. Only the first column matters. No headers allowed.
    @path_data
            (str - filepath)
            The filepath of the input data file. No headers allowed.
    @path_output
            (str - filepath)
            The filepath of the output file.
    @method
            (int) - Pseudo ENUM
            An integer which indicates the method by which the score for each
            locus is calculated. The method also covers what to do when a read
            (data file) overlaps with more than one loci. Acceptable options
            are:
                1:  COMPLETE - The read contributes a score of 1 if it overlaps
                    completely with the locus.
                2:  FRACTION - The score contributed by a read to a loci is the
                    percentage of the read which overlaps.
                    (Ex. If the read is 100bp, and 75bp overlaps with the loci
                    of interest, it increases the loci score by 0.75.)
                3:  MAJORITY - Each read contributes a score of 1 if more than
                    half of the read overlaps with the locus.
                4:  NUCLEOTIDES - Each read contributes a score to a locus equal
                    to the number of nucleotides which overlap with said locus.
                5:  PRESENT - Each read contributes a score of 1 if there is an
                    overlap at all.
                6:  THRESHOLD - Each read contributes a score of 1 if it exceeds
                    a given threshold. If the threshold given is an integer, the
                    threshold required is the number of overlapping nucleotides.
                    If the threshold given is a decimal number, the threshold
                    required is the percentage of the read which overlaps with
                    the locus.
    @threshold
            (int/float)
            Only relevant if the [THRESHOLD] option is selected for @method.
            If the threshold given is an integer, the threshold required is the
            number of overlapping nucleotides. If the threshold given is a
            decimal number, the threshold required is the percentage of the read
            which overlaps with the locus.
    
    Return a value of 0 if the function runs successfully.
    Return a positive integer if there is a problem. The integer functions as an
    error code.
    
    Tally_Reads_In_Loci(str, str, str, str, int, int/float) -> int
    """
    # Setup reporting
    count_loci = 0
    count_entries = 0
    count_relevant = 0
    total_score = 0
    
    # I/O setup (1)
    MBC = Multitrack_BED_Coordinator()
    MBC.Set_Path_Loci(path_loci)
    MBC.Set_Path_Chrs(path_chrs)
    MBC.Add_Path_Data(path_data)
    
    if method == METHOD.COMPLETE:
        MBC.Set_Retain_Partial_Overlaps(False)
    else:
        MBC.Set_Retain_Partial_Overlaps(True)
    MBC.Set_Retain_Prior(False)
    MBC.Set_Retain_Skipped(False)
    MBC.Set_Retain_Remainder(False)
    MBC.Set_Retain_Remaining_Chrs(False)
    
    MBC.Toggle_Printing_E(True)
    
    # Setup
    if method == METHOD.MAJORITY:
        method = METHOD.THRESHOLD
        threshold = 0.5000000000000001 # Lowest detectable difference
    
    # Test
    rt = MBC.Open()
    if rt == 1:
        return 1 # MBC should print Error messages
    elif rt == 2:
        return 2 # MBC should print Error messages

    # I/O setup (2)
    o = open(path_output, "w")
    
    # Main loop
    PRINT.printP(STR__tally_begin)
    ### COMPLETE / PRESENT
    if method == METHOD.COMPLETE or method == METHOD.PRESENT:
        while not MBC.EOF:
            # Read
            rt = MBC.Read()
            if rt:
                MBC.Close()
                o.close()
                return rt + 2
            count_loci += 1
            reads = MBC.current[0]
            # Calculate
            score = 0
            for read in reads:
                count_relevant += 1
                score += 1
            total_score += score
            # Write
            sb = (MBC.cur_chr + "\t" + str(MBC.cur_start) + "\t" +
                    str(MBC.cur_end) + "\t" +
                    "\t".join(MBC.current_locus_data[3:]) + "\t" + str(score)
                    + "\n")
            o.write(sb)
    ### NUCLEOTIDE COUNT
    elif method == METHOD.NUCLEOTIDES:
        while not MBC.EOF:
            # Read
            rt = MBC.Read()
            if rt:
                MBC.Close()
                o.close()
                return rt + 2
            count_loci += 1
            reads = MBC.current[0]
            # Calculate
            score = 0
            for read in reads:
                count_relevant += 1
                # Left
                if MBC.cur_start < read[1]:
                    left = read[1]
                else:
                    left = MBC.cur_start
                # Right
                if MBC.cur_end > read[2]:
                    right = read[2]
                else:
                    right = MBC.cur_end
                # Overlap
                overlap = (right - left) + 1
                score += overlap
            total_score += score
            # Write
            sb = (MBC.cur_chr + "\t" + str(MBC.cur_start) + "\t" +
                    str(MBC.cur_end) + "\t" +
                    "\t".join(MBC.current_locus_data[3:]) + "\t" + str(score)
                    + "\n")
            o.write(sb)
    ### NUCLEOTIDE THRESHOLD
    elif (method == METHOD.THRESHOLD) and (type(threshold) == int):
        while not MBC.EOF:
            # Read
            rt = MBC.Read()
            if rt:
                MBC.Close()
                o.close()
                return rt + 2
            count_loci += 1
            reads = MBC.current[0]
            # Calculate
            score = 0
            for read in reads:
                # Left
                if MBC.cur_start < read[1]:
                    left = read[1]
                else:
                    left = MBC.cur_start
                # Right
                if MBC.cur_end > read[2]:
                    right = read[2]
                else:
                    right = MBC.cur_end
                # Overlap
                overlap = (right - left) + 1
                if overlap >= threshold:
                    count_relevant += 1
                    score += 1
            total_score += score
            # Write
            sb = (MBC.cur_chr + "\t" + str(MBC.cur_start) + "\t" +
                    str(MBC.cur_end) + "\t" +
                    "\t".join(MBC.current_locus_data[3:]) + "\t" + str(score)
                    + "\n")
            o.write(sb)
    ### FRACTION
    elif method == METHOD.FRACTION:
        while not MBC.EOF:
            # Read
            rt = MBC.Read()
            if rt:
                MBC.Close()
                o.close()
                return rt + 2
            count_loci += 1
            reads = MBC.current[0]
            locus_size = (MBC.cur_end - MBC.cur_start) + 1.0
            # Calculate
            score = 0
            for read in reads:
                # Left
                if MBC.cur_start < read[1]:
                    left = read[1]
                else:
                    left = MBC.cur_start
                # Right
                if MBC.cur_end > read[2]:
                    right = read[2]
                else:
                    right = MBC.cur_end
                # Overlap
                overlap = (right - left) + 1
                percentage = overlap / locus_size
                count_relevant += 1
                score += percentage
            total_score += score
            # Write
            sb = (MBC.cur_chr + "\t" + str(MBC.cur_start) + "\t" +
                    str(MBC.cur_end) + "\t" +
                    "\t".join(MBC.current_locus_data[3:]) + "\t" + str(score)
                    + "\n")
            o.write(sb)
    ### PERCENTAGE THRESHOLD (includes MAJORITY)
    else:
        while not MBC.EOF:
            # Read
            rt = MBC.Read()
            if rt:
                MBC.Close()
                o.close()
                return rt + 2
            count_loci += 1
            reads = MBC.current[0]
            # Calculate
            score = 0
            for read in reads:
                # Left
                if MBC.cur_start < read[1]:
                    left = read[1]
                else:
                    left = MBC.cur_start
                # Right
                if MBC.cur_end > read[2]:
                    right = read[2]
                else:
                    right = MBC.cur_end
                # Overlap
                overlap = (right - left) + 1
                read_size = (read[2] - read[1] ) + 1.0
                fraction = overlap / read_size
                if fraction >= threshold:
                    count_relevant += 1
                    score += 1
            total_score += score
            # Write
            sb = (MBC.cur_chr + "\t" + str(MBC.cur_start) + "\t" +
                    str(MBC.cur_end) + "\t" +
                    "\t".join(MBC.current_locus_data[3:]) + "\t" + str(score)
                    + "\n")
            o.write(sb)
    
    # Finish
    MBC.Read_Final()
    MBC.Close()
    o.close()
    PRINT.printP(STR__tally_complete)
    
    # Reporting
    metrics = [count_loci, MBC.files_data[0].line_no, count_relevant,
            total_score]
    Report_Metrics(metrics)
    
    # Wrap up
    return 0

def Report_Metrics(summary_metrics):
    """
    Print a report into the command line interface of the metrics of the
    operation.
    
    @summary_metrics
            (list<int>)
            A list of summary metrics for the data, including:
                * The number of loci in the loci file
                * The number of entries in the data file
                * The number of entries which were relevant
                * The total score of the output
    
    Report_Metrics([int, int, int, int]) -> None
    """
    # Unpacking
    loci, entries, relevant, score = summary_metrics
    # Average
    if relevant == 0:
        average = 0.0
    else:
        average = (float(score))/relevant
    # Strings
    loci = str(loci) + "   "
    entries = str(entries) + "   "
    relevant = str(relevant) + "   "
    if type(score) == int:
        score = str(score) + "   "
    else:
        score = str(average)
        score = Trim_Percentage_Str(score, 2)
    average = str(average)
    average = Trim_Percentage_Str(average, 2)
    # Pad Column
    col = [loci, entries, relevant, score, average]
    col = Pad_Column(col, 0, 0, " ", 0)
    loci, entries, relevant, score, average = col
    # Print
    PRINT.printM(STR__metrics.format(A = loci, B = entries, C = relevant,
            D = score, E = average))



# Command Line Parsing #########################################################

def Parse_Command_Line_Input__Tally_Reads_In_Loci(raw_command_line_input):
    """
    Parse the command line input and call the Tally_The_Reads_In_Each_Loci
    function with appropriate arguments if the command line input is valid.
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
    if len(inputs) < 3:
        PRINT.printE(STR__insufficient_inputs)
        PRINT.printE(STR__use_help)
        return 1
    
    # Validate mandatory inputs
    path_loci = inputs.pop(0)
    valid = Validate_Read_Path(path_loci)
    if valid == 1:
        PRINT.printE(STR__IO_error_read.format(f = path_loci))
        PRINT.printE(STR__use_help)
        return 1
    path_chrs = inputs.pop(0)
    valid = Validate_Read_Path(path_chrs)
    if valid == 1:
        PRINT.printE(STR__IO_error_read.format(f = path_chrs))
        PRINT.printE(STR__use_help)
        return 1
    path_data = inputs.pop(0)
    valid = Validate_Read_Path(path_data)
    if valid == 1:
        PRINT.printE(STR__IO_error_read.format(f = path_data))
        PRINT.printE(STR__use_help)
        return 1
    
    # Set up rest of the parsing
    path_out = Generate_Default_Output_File_Path_From_File(path_data, FILEMOD,
            True)
    method = DEFAULT__method
    threshold = DEFAULT__threshold
    
    # Validate optional inputs (except output path)
    while inputs:
        arg = inputs.pop(0)
        flag = 0
        try: # Following arguments
            if arg in ["-o", "-m", "-t"]:
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
        elif arg == "-m":
            if arg2 in LIST__complete:
                method = METHOD.COMPLETE
            elif arg2 in LIST__fraction:
                method = METHOD.FRACTION
            elif arg2 in LIST__majority:
                method = METHOD.MAJORITY
            elif arg2 in LIST__nucleotides:
                method = METHOD.NUCLEOTIDES
            elif arg2 in LIST__present:
                method = METHOD.PRESENT
            elif arg2 in LIST__threshold:
                method = METHOD.THRESHOLD
            else:
                PRINT.printE(STR__invalid_method.format(s = arg2))
                return 1
        else: # arg == "-t"
            try:
                threshold = int(arg2)
                if threshold < 1:
                    PRINT.printE(STR__invalid_threshold.format(s = arg2))
                    return 1
            except:
                try:
                    threshold = float(arg2)
                    if (threshold < 0) or (threshold > 1):
                        PRINT.printE(STR__invalid_threshold.format(s = arg2))
                        return 1
                except:
                    PRINT.printE(STR__invalid_threshold.format(s = arg2))
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
    exit_state = Tally_The_Reads_In_Each_Loci(path_loci, path_chrs, path_data,
            path_out, method, threshold)
    
    # Exit
    if exit_state == 0: return 0
    else: return exit_state



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
    exit_code = Parse_Command_Line_Input__Tally_Reads_In_Loci(sys.argv)


