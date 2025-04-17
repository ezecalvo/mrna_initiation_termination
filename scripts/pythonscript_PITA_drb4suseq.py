import re
import argparse
from pathlib import Path
import collections
import gzip
import pysam
import pybedtools

# custom imports
from porechop_parser import parse_log

def readBam(bamfile):
    #bamhere = pysam.Samfile(bamfile, 'rb')
    bamhere = pysam.AlignmentFile(bamfile, 'rb', require_index=False)
    return(bamhere)

def filterBamRNA(bamhere):
    # loop through reads in bam file and return directionary of read names, start/end coordinates
    startbed = []
    endbed = []
    lengthdict = dict()
    total_counter = 0
    for read in bamhere:
        total_counter = total_counter + 1
        #if total_counter % 100 == 0: 
        #    print("... {}".format(total_counter), flush=True)
        strand = '+'
        startpos = int(read.pos)
        endpos = int(read.aend)
        if read.is_reverse:
            startpos = int(read.aend)
            endpos = int(read.pos)
            strand = '-'
        if startpos < 1 or endpos < 1: continue
        # save relevant values
        readlength = read.query_length
        mappedlength = abs(endpos - startpos)
        lengthdict[str(read.qname)] = str(readlength) +':'+ str(mappedlength)
        startbed.append(str(bamhere.getrname(read.rname)) +'\t'+ str(startpos - 1) + '\t'+ str(startpos + 1) + '\t'+ str(read.qname) +'\t.\t'+ str(strand))
        endbed.append(str(bamhere.getrname(read.rname)) +'\t'+ str(endpos - 1) + '\t'+ str(endpos + 1) + '\t'+ str(read.qname) +'\t.\t'+ str(strand))
    print("*** READING READS ***", flush=True)
    print("Total number of mapped reads = {}".format(total_counter), flush=True)
    print("{} start and {} end coordinates".format(len(startbed), len(endbed)), flush=True)
    print(" ", flush=True)
    # return bed files of start and end coords
    return(startbed, endbed, lengthdict)


def filterBamDNA(bamhere, porechopreads):
    # porechopreads: read dictionary where key is read name and value is strand
    # strand = '+' means start coordinate read start, end coordinate is at read end
    # strand = '-' means start coordinate is at read end, end coordinate is at read start
    # loop through reads in bam file and return directionary of read names, start/end coordinates
    porechoplist = sorted([*porechopreads])
    #porechoplist = sorted([*porechopreads], key= lambda s: int(''.join(filter(str.isdigit, s))))
    startbed = []
    endbed = []
    lengthdict = dict()
    total_counter = 0
    filtered_counter = 0
    for read in bamhere:
        total_counter = total_counter + 1
        if total_counter % 100 == 0: 
            print("... {}, {} retained, {} in list".format(total_counter, filtered_counter, len(porechoplist)), flush=True)
            print("")
        if read.qname in porechoplist:
            filtered_counter = filtered_counter + 1
            porechoplist.remove(str(read.qname))
            #readname = str(read.qname)
            strand = porechopreads[str(read.qname)]
            startpos = int(read.pos)
            endpos = int(read.aend)
            if strand == '-':
                startpos = int(read.aend)
                endpos = int(read.pos)
            if startpos < 1 or endpos < 1: continue
            # save relevant values
            readlength = read.query_length
            mappedlength = abs(endpos - startpos)
            lengthdict[str(read.qname)] = str(readlength) +':'+ str(mappedlength)
            startbed.append(str(bamhere.getrname(read.rname)) +'\t'+ str(startpos - 1) + '\t'+ str(startpos + 1) + '\t'+ str(read.qname) +'\t.\t'+ str(strand))
            endbed.append(str(bamhere.getrname(read.rname)) +'\t'+ str(endpos - 1) + '\t'+ str(endpos + 1) + '\t'+ str(read.qname) +'\t.\t'+ str(strand))
    print("*** FILTERING READS ***", flush=True)
    print("Total number of mapped reads = {}".format(total_counter), flush=True)
    print("Primary alignments after adapter filter = {}".format(filtered_counter), flush=True)
    print(" ", flush=True)
    # return bed files of start and end coords
    return(startbed, endbed)


def intersectGene(endbed, genesbed, stranded):
    #https://daler.github.io/pybedtools/intersections.html
    # https://bedtools.readthedocs.io/en/latest/content/tools/map.html
    # make BedTool object from bedfile of read ends (from filterBam function)
    endbedscratch = pybedtools.BedTool('\n'.join(endbed), from_string=True)
    endbedsort = endbedscratch.sort()
    # make BedTool object from bedfile of genes (from file)
    genes = pybedtools.BedTool(genesbed)
    # perform bedtools intersect
    if stranded:
        bmap = endbedsort.map(genes, c=4, o='distinct', null='none',s=True)
    else:
        bmap = endbedsort.map(genes, c=4, o='distinct', null='none')
    # go through reads and get reads overlapping genes
    geneReads = dict()
    readlist = str(bmap).split('\n')[:-1]
    counter = 0
    for read in readlist:
        readhere = read.split('\t')
        gene = str(readhere[6]).split(',')
        readname = str(readhere[3])
        if 'none' not in gene and len(gene) == 1:
            counter = counter + 1
            geneReads[readname] = gene[0]
    print("*** INTERSECT 3' END OF READ WITH GENE ***", flush=True)
    print("Total number of filtered reads = {}".format(len(readlist)), flush=True)
    print("Filtered reads with 3' end in gene = {}".format(counter), flush=True)
    print(" ", flush=True)
    # return read dictionary where key is read name and value is strand
    return(geneReads)


def intersectAFE(startbed, AFEbed, stranded):
    #https://daler.github.io/pybedtools/intersections.html
    # https://bedtools.readthedocs.io/en/latest/content/tools/map.html
    # make BedTool object from bedfile of read ends (from filterBam function)
    startbedscratch = pybedtools.BedTool('\n'.join(startbed), from_string=True)
    startbedsort = startbedscratch.sort()
    # make BedTool object from bedfile of AFEs (from file)
    afes = pybedtools.BedTool(AFEbed)
    afesort = afes.sort()
    # add AFE# to gene column
    afes_string = str(afes).split('\n')
    afes_string_edited = []
    for entry in afes_string:
        if entry == '': break
        entryhere = entry.split('\t')
        entryhere[3] = entryhere[3] +':'+ entryhere[4]
        collapsedentry = "\t".join(entryhere)
        afes_string_edited.append(collapsedentry)
    afes_edited_bed = pybedtools.BedTool('\n'.join(afes_string_edited), from_string=True)
    afes_edited_sort = afes_edited_bed.sort()
    # perform bedtools intersect
    if stranded:
        bmap = startbedsort.map(afes_edited_sort, c=4, o='distinct', null='none',s=True)
    else:
        bmap = startbedsort.map(afes_edited_sort, c=4, o='distinct', null='none')
    # go through reads and get reads overlapping AFEs
    afeReads = dict()
    readlist = str(bmap).split('\n')[:-1]
    counter = 0
    for read in readlist:
        readhere = read.split('\t')
        geneafe = str(readhere[6]).split(',')
        readname = str(readhere[3])
        if 'none' not in geneafe and len(geneafe) == 1:
            counter = counter + 1
            afeReads[readname] = geneafe
    print("*** INTERSECT 5' END OF READ WITH AFEs ***", flush=True)
    print("Total number of filtered reads = {}".format(len(readlist)), flush=True)
    print("Filtered reads with 5' end in AFE = {}".format(counter), flush=True)
    print(" ", flush=True)
    # return read dictionary where key is read name and value is strand
    return(afeReads)


def intersectPAS(endbed, PASbed, stranded):
        #https://daler.github.io/pybedtools/intersections.html
    # https://bedtools.readthedocs.io/en/latest/content/tools/map.html
    # make BedTool object from bedfile of read ends (from filterBam function)
    endbedscratch = pybedtools.BedTool('\n'.join(endbed), from_string=True)
    endbedsort = endbedscratch.sort()
    # make BedTool object from bedfile of PASs (from file)
    ales = pybedtools.BedTool(PASbed)
    alesort = ales.sort()
    # perform bedtools intersect
    if stranded:
        bmap = endbedsort.map(alesort, c=4, o='distinct', null='none',s=True)
    else:
        bmap = endbedsort.map(alesort, c=4, o='distinct', null='none')
    # go through reads and get reads NOT overlapping PASs
    pasReads = []
    readlist = str(bmap).split('\n')[:-1]
    counter = 0
    for read in readlist:
        readhere = read.split('\t')
        pas = str(readhere[6])
        readname = str(readhere[3])
        if pas == 'none':
            counter = counter + 1
            pasReads.append(readname)
    print("*** INTERSECT 3' END OF READ WITH PASs ***", flush=True)
    print("Total number of filtered reads = {}".format(len(readlist)), flush=True)
    print("Filtered reads with 3' end not in PAS = {}".format(len(pasReads)), flush=True)
    print(" ", flush=True)
    # return read dictionary where key is read name and value is strand
    return(pasReads)


def conditionReads(geneOverlap, afeOverlap, pasOverlap, lengthdict, outname):
    finalReads = []
    outfile = gzip.open(outname +'_parsedReads.info.gz', 'wt')
    outfile.write('readname' +'\t'+ 'gene' +'\t'+ 'afe' +'\t'+ 'mapped_length' +'\t'+ 'read_length' +'\n')
    commonreads = set([*geneOverlap]) & set([*afeOverlap]) & set([*pasOverlap])
    for read in commonreads:
        gene = geneOverlap[read]
        afegene = afeOverlap[read][0].split(':')[0]
        afe = afeOverlap[read][0].split(':')[1]
        if gene == afegene:
            readlength = lengthdict[read].split(':')[0]
            mappedlength = lengthdict[read].split(':')[1]
            # output file with readname, gene, AFE name, mapped length, read length
            outfile.write(str(read) +'\t'+ str(gene) +'\t'+ str(afe) +'\t'+ str(mappedlength) +'\t'+ str(readlength) +'\n')
            finalReads.append(read)
    print("*** CONDITION ON READ OVERLAP CRITERIA ***", flush=True)
    print("****** (1) read start is in AFE", flush=True)
    print("****** (2) read end is NOT in PAS", flush=True)
    print("****** (3) read end is within a gene", flush=True)
    print("****** (4) read start AFE is in the same gene as read end", flush=True)
    #print("Total number of filtered reads = {}".format(len(readlist)))
    print("Reads meeting all these criteria = {}".format(len(finalReads)), flush=True)
    print(" ", flush=True)
    outfile.close()
    # return dictionary of final reads
    return(finalReads)

def filterBamFinal(bamhere, finalReads, outname):
    outbam = pysam.AlignmentFile(outname +"_parsedReads.bam", "wb", template=bamhere)
    for read in bamhere:
        readname = str(read.qname)
        if readname in [*finalReads]:
            # output bam file
            outbam.write(read)
    bamhere.close()
    outbam.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'NOTE: must have samtools and bedtools accessible on command line')
    # parameters
    group_param = parser.add_argument_group('Library type')
    group_param.add_argument('--librarytype', type=str, default='directRNA', choices = ['directRNA', 'cDNA'], help='type of library, determines whether to look for adapters or not (DEFAULT = directRNA)', required = False)
    group_param.add_argument('--stranded', action = 'store_true', default = False, help = 'consider read strand information in gene and exon level analyses (DEFAULT = False)', required = False)
    group_param.add_argument('--outbam', action = 'store_true', default = False, help = 'write out a bam with filtered reads', required = False)
    # shared inputs
    group_inputs = parser.add_argument_group('Inputs', 'shared inputs for all library types')
    group_inputs.add_argument('--bam', type = str, help = 'bam file with mapped LRS reads (including path), with index file in the same folder', required = True)
    group_inputs.add_argument('--genes', type = str, help = 'bed file with gene coordinates (including path)', required = True)
    group_inputs.add_argument('--AFE', type = str, help = 'bed file with AFE coordinates (including path)', required = True)
    group_inputs.add_argument('--PAS', type = str, help = 'bed file with PAS coordinates (including path)', required = True)
    # shared outputs
    group_outputs = parser.add_argument_group('Outputs', 'shared outputs for all library types')
    group_outputs.add_argument('--outname', type = str, help = 'name for output files (no suffix)', required = True)
    # cDNA library parameters
    group_cdna = parser.add_argument_group('cDNA libraries', 'extra arguments required for cDNA libraries')    
    group_cdna.add_argument('--porechop', type = str, help = 'porechop output file (including path)', required = False, default = "None")
    group_cdna.add_argument('--adapters', type=str, default = '5prime', choices = ['5prime', '3prime', 'both'], help = 'position of adapter to condition upon (DEFAULT = 5prime)', required = False)

    args = parser.parse_args()

    outname=args.outname

    ## read in bam files
    bamhere = readBam(args.bam)

    ## library type
    if args.librarytype == 'directRNA':
        startbed, endbed, lengthdict = filterBamRNA(bamhere)

    if args.librarytype == 'cDNA':
        if args.porechop == "None":
            sys.exit("ERROR! Need to include --porechop when --librarytype = cDNA")        
        ## parse the porechop output
        # returns read dictionary where key is read name and value is strand
        # strand = '+' means start coordinate read start, end coordinate is at read end
        # strand = '-' means start coordinate is at read end, end coordinate is at read start
        porechopreads = parse_log(args.porechop, args.adapters)
        ## filter bam file by porechop output, get dictionary of read starts/ends
        startbed, endbed, lengthdict = filterBamDNA(bamhere, porechopreads)

    ## intersect reads, get list of read names
    geneOverlap = intersectGene(endbed, args.genes, args.stranded)
    afeOverlap = intersectAFE(startbed, args.AFE, args.stranded)
    pasOverlap = intersectPAS(endbed, args.PAS, args.stranded)

    ## condition reads by intersects
    finalreads = conditionReads(geneOverlap, afeOverlap, pasOverlap, lengthdict, outname)

    ## write out filtered bam file
    if args.outbam:
        bamhere = readBam(args.bam)
        filterBamFinal(bamhere, finalreads, outname)

#adapters = '5prime'
#porechop='/pi/athma.pai-umw/analyses/valeria/drb_4su_PITA/K562_rep1_10m/porechop/logs/all_out.out.gz'
#bam='/pi/athma.pai-umw/analyses/valeria/drb_4su_PITA/K562_rep1_10m/minimap2/bams/K562_rep1_10m.merged.mapped.bam'
#bam='/pi/athma.pai-umw/analyses/athma/elongation/test.sorted.bam'
#genes='/pi/athma.pai-umw/genomes/hg38/hg38_ensembl_v95/Homo_sapiens.GRCh38.95.uniquegene.bed.bed'
#AFE='/pi/athma.pai-umw/analyses/ezequiel/TSS-TES/first_exonstmp.bed.gz'
#PAS='/pi/athma.pai-umw/analyses/ezequiel/TSS-TES/polyAtmp.bed.gz'


# load python3, miniconda3, samtools, bedtools modules

# directRNA, without bam file of filtered reads (fastest)
# python3 pythonscript_PITA_drb4suseq.py --librarytype directRNA --stranded --bam BAM --genes GENES --AFE AFE --PAS PAS --outname OUTNAME
# directRNA, with bam file of filtered reads (2nd fastest)
# python3 pythonscript_PITA_drb4suseq.py --librarytype directRNA --stranded --outbam --bam BAM --genes GENES --AFE AFE --PAS PAS --outname OUTNAME

# cDNA, without bam file of filtered reads (3rd fastest)
# python3 pythonscript_PITA_drb4suseq.py --librarytype cRNA --bam BAM --genes GENES --AFE AFE --PAS PAS --outname OUTNAME --porechop PORECHOP --adapters both
# cDNA, with bam file of filtered reads (slowest)
# python3 pythonscript_PITA_drb4suseq.py --librarytype cRNA --outbam --bam BAM --genes GENES --AFE AFE --PAS PAS --outname OUTNAME --porechop PORECHOP --adapters both

