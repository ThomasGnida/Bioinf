file = "Assignment_2_fasta"
#Alternatively for multiple fasta files or different formats more akin to original fasta files the first lines have
# to be changed to adapt to the fasta format and multiple files


IUPAC_set = { #Define IUPAC base nomenclature as set required for ambiguous bases
    'A': {'A'}, 'T': {'T'}, 'C': {'C'}, 'G': {'G'},
    'N': {'A', 'T', 'C', 'G'},
    'R': {'A', 'G'},
    'Y': {'C', 'T'},
    'S': {'G', 'C'},
    'W': {'A', 'T'},
    'K': {'G', 'T'},
    'M': {'A', 'C'},
    'B': {'C', 'G', 'T'},
    'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'},
    'V': {'A', 'C', 'G'},
}

def bases_are_compatible(b1, b2):
    '''
    :param b1: Base 1
    :param b2: Base 2
    :return: True if matches are compatible according to IUPAC
    '''
    if b1.upper() == b2.upper():
        return True
    combined = IUPAC_set.get(b1.upper(), {b1}) & IUPAC_set.get(b2.upper(), {b2})
    if combined:
        for base in [b1, b2]:
            if IUPAC_set.get(base.upper()) == combined:
                return True
    else:
        return False
    return False #Fallback

def compare_base_strings(s1, s2, threshold = 1, min_length = 10):
    '''
    :param s1: String of base sequence 1
    :param s2: String of base sequence 2
    :param threshold: Threshold for incorrect bases, default is 1
    :param min_length: Minimum length of bases, default is 10
    :return: True if strings match according to IUPAC or if there are more mismatches than defined by the threshold
    '''
    if len(s1) != len(s2):
        return False
    if len(s1) < min_length:
        return False
    incorrect_matches = 0
    for i in range(len(s1)):
        if not bases_are_compatible(s1[i], s2[i]):
            incorrect_matches += 1
    if incorrect_matches > threshold:
        return False

    return True

with open(file) as Input_file:
    next_line_is_sequence = False
    sequences = []
    final_sequence = ""
    expected_length = 0
    '''
    Read the file, after every >primer the next line contains a relevant sequence for task 1
    '''
    for line in Input_file.readlines():
        if next_line_is_sequence:
            sequences.append(line.strip())
            next_line_is_sequence = False
        if line.startswith(">primer"):
            next_line_is_sequence = True


    '''
    The following block uses splicing and the predefined functions compare_base_strings to check if the first j+1
    entries match the last j+1 entries of the current and next sequence. As the matching part of the last sequence is
    already removed from the previous string the last sequence can be appended to the string.
    '''
    overlapping_sequence = ""
    final_sequence = sequences[0]
    expected_length += len(sequences[0])
    for i in range(len(sequences)):
        if i == len(sequences) - 1:
            break
        for j in range(len(sequences[i])):
            if compare_base_strings(sequences[i+1][:j+1], sequences[i][-(j+1):], 1):
                overlapping_sequence = sequences[i+1][:j+1]

        final_sequence = final_sequence + sequences[i + 1][len(overlapping_sequence):]
        expected_length +=  len(sequences[i + 1][len(overlapping_sequence):])

    if len(final_sequence) != expected_length:
        print("Incorrect sequence length")
        print(expected_length, len(final_sequence))
    else:
        print("Task 1\n The Final sequence is:")
        print(final_sequence + "\n")



with open(file) as Input_file:
    #There are multiple adapter sequences from each distributor for different sets, e.g. CTGTCTCTTATACACATCT for
    # Illumina read1 and read2 for Illumina and Nextera prep kits. For simplifications 3 adapter sequences have
    # been provided that are used in the code.
    adapters_set = {
        "Illumina":"AGATCGGAAGAGC",
        "Nanopore":"TTTCTGTTGGTGCTG",
        "PacBio":"ATCTCTCTCAACA"
    }

    next_line_is_sequence = False
    sequences = []
    adapters = []
    final_sequence = ""

    '''
    For task 2 similar to task 1 the file is read line by line, if the line starts with >read the next line contains the
    relevant sequence. The line containing the >read also contains the device, this is assigned to an array of adapters
    that are read by splitting the line at the equal sign.
    Afterwards each sequence in the sequences array is compared for the the adapter that fits to the corresponding 
    device. Only if the adapter sequence matches exactly it is removed through splicing.
    '''
    for line in Input_file.readlines():
        if next_line_is_sequence:
            sequences.append(line.strip())
            next_line_is_sequence = False
        if line.startswith(">read"):
            adapters.append(line.split("=")[1].strip())
            next_line_is_sequence = True

    for i, seq in enumerate(sequences):
        adapter_seq = adapters_set.get(adapters[i])
        if compare_base_strings(adapter_seq, seq[:len(adapter_seq)], 0, 1):
            sequences[i] = seq[len(adapter_seq):]
        elif compare_base_strings(adapter_seq, seq[-len(adapter_seq):], 0, 1):
            sequences[i] = seq[:-len(adapter_seq)]

    print("Task 2 \n The trimmed sequences are:")
    for seq in sequences:
        print(seq)
