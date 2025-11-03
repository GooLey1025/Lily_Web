import re
import argparse

parser = argparse.ArgumentParser(description="Convertrelative BLAST positions to absolute positions based on LG segmentation ranges.")
parser.add_argument("input", metavar="INPUT_FILE", help="Path to the input BLAST output file.")
parser.add_argument("output", metavar="OUTPUT_FILE", help="Path to the output file")
args=parser.parse_args()

# Input and output file paths
blast_result_file = args.input# Input BLAST output file
output_file = args.output   # Output file with absolute positions

with open(blast_result_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        fields = line.strip().split("\t")
        if len(fields) < 10:
            outfile.write("\t".join(fields) + "\n")
            continue

        sseqid = fields[1]  # Second column: target sequence ID
        # Try to extract segmentation range from sseqid: support name:start-end or name_start-end, optional trailing ':'
        match = re.match(r"^(.+?)[_:](\d+)-(\d+):?$", sseqid)

        if match:
            segment_start = int(match.group(2))
            # Convert relative positions to absolute positions
            try:
                sstart = int(fields[8])  # Ninth column: sstart
                send = int(fields[9])    # Tenth column: send
            except ValueError:
                outfile.write("\t".join(fields) + "\n")
                continue

            # seqkit subseq headers typically report 1-based inclusive coordinates in the suffix
            # Convert relative (1-based) BLAST positions by adding (segment_start - 1)
            offset = segment_start - 1
            absolute_sstart = sstart + offset
            absolute_send = send + offset

            fields[8] = str(absolute_sstart)
            fields[9] = str(absolute_send)

        # Write updated line to the output file (updated when segment suffix present)
        outfile.write("\t".join(fields) + "\n")
