import re
import argparse
from pathlib import Path
import collections
import gzip


#log_file = "jobid_754307_100.out"


def parse_log(log_file, adapters):
    if log_file[-2:] == "gz": file = gzip.open(log_file, 'rt')
    else: file = open(log_file, 'r')
    data = file.read()


    # Split log file into entries based on double newlines
    entries = data.strip().split('\n\n')


    # Prepare output lines
    #output_lines = ["read_id\tstart_adapter\tend_adapter\tclassification"]
    #plus_strand_reads = []
    #minus_strand_reads = []
    #plus_strand_reads_5prime = []
    #minus_strand_reads_5prime = []


    #adapters choices ['5prime', '3prime', 'both']
    adapter_choice = [adapters]
    if adapters == 'both': adapter_choice = ['3prime','5prime']


    # initiate dictionary for reads
    read_dict = dict()
    total_counter = 0
    filtered_counter = 0
    for entry in entries:
        total_counter = total_counter + 1
        # Extract read_id
        read_id_match = re.search(r'^\S+', entry)
        read_id = read_id_match.group(0) if read_id_match else "unknown"


        # Extract all start adapters
        start_alignments_block = re.search(r"start alignments:(.*?)end:", entry, re.S)
        if start_alignments_block:
            start_adapters = re.findall(r"^\s+(.*?),", start_alignments_block.group(1), re.M)
        else:
            start_adapters = "no_adapter"
        start_adapter = ";".join(start_adapters) if start_adapters else "no_adapter"


        # Extract all end adapters
        end_alignments_block = re.search(r"end alignments:(.*)$", entry, re.S)
        if end_alignments_block:
            end_adapters = re.findall(r"^\s+(.*?),", end_alignments_block.group(1), re.M)
        else:
            end_adapters = "no_adapter"
        end_adapter = ";".join(end_adapters) if end_adapters else "no_adapter"


        # Split adapters into lists
        #start_adapters = start_adapter.split(';')
        #end_adapters = end_adapter.split(';')


        a_minion = 'SQK-NSK007'
        a5_sense = 'SMARTer_IIA_Oligo(+_strand_being_read)'
        a5_senseminion = 'SQK-NSK007;SMARTer_IIA_Oligo(+_strand_being_read)'
        a5_anti = 'SMARTer_IIA_rev_comp(-_strand_being_read)'
        a5_antiminion = 'SQK-NSK007;SMARTer_IIA_rev_comp(-_strand_being_read)'
        a3_sense = 'FLAM_rev_comp(+_strand_being_read)'
        a3_senseminion = 'SQK-NSK007;FLAM_rev_comp(+_strand_being_read)'
        a3_anti = 'FLAM_Fw(-_strand_being_read)'
        a3_antiminion = 'SQK-NSK007;FLAM_Fw(-_strand_being_read)'


        type = 'unexpected'
        strand = 'na'


        ## BOTH ADAPTERS IN PROPER DIRECTION
        if( (start_adapter == a5_sense or start_adapter == a5_senseminion) and 
        	(end_adapter == a3_sense or end_adapter == a3_senseminion) ): 
        	type = 'both'
        	strand = '+'
        if( (start_adapter == a3_anti or start_adapter == a3_antiminion) and 
        	(end_adapter == a5_anti or end_adapter == a5_antiminion) ): 
        	type = 'both'
        	strand = '-'
        ## 5' ADAPTER IN PROPER DIRECTION
        if( (start_adapter == a5_sense or start_adapter == a5_senseminion) and 
        	(end_adapter == a_minion or end_adapter == "no_adapter" ) ): 
        	type = '5prime'
        	strand = '+'
        if( (end_adapter == a5_anti or end_adapter == a5_antiminion) and 
        	(start_adapter == a_minion or start_adapter == "no_adapter" ) ): 
        	type = '5prime'
        	strand = '-'
        ## 3' ADAPTER IN PROPER DIRECTION
        if( (end_adapter == a3_sense or end_adapter == a3_senseminion) and 
        	(start_adapter == a_minion or start_adapter == "no_adapter" ) ): 
        	type = '3prime'
        	strand = '+'
        if( (start_adapter == a3_anti or start_adapter == a3_antiminion) and 
        	(end_adapter == a_minion or end_adapter == "no_adapter" ) ): 
        	type = '3prime'
        	strand = '-'


        ## condition on adapter parameter
        if type == adapters:
        	read_dict[read_id] = strand
        	filtered_counter = filtered_counter + 1


    print("*** PARSING PORECHOP ***", flush=True)
    print("Total number of sequenced reads = {}".format(total_counter), flush=True)
    print("Sequenced reads with {} adapters = {}".format(adapters, filtered_counter), flush=True)
    print(" ", flush=True)
    # return read dictionary where key is read name and value is strand
    return(read_dict)


if __name__ == "__main__":
    parser = arg.parse.ArgumentParser()
    # parameters
    group_param = parser.add_argument('Parameters')
    group_param.add_argument('--stranded', action = 'store_true', default = False, help = 'condition on stranded reads')
    group_param.add_argument('--adapters', type=str, default = '5prime', choices = ['5prime', '3prime', 'both'], help = 'position of adapter to condition upon')
    # inputs
    group_inputs = parser.add_argument('Inputs')
    group_inputs.add_argument('--porechop', type = str, help = 'porechop output file', required = True)


    args = parser.parse_args()

