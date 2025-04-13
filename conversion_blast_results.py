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
        sseqid = fields[1]  # Second column: target sequence ID
        sstart = int(fields[8])  # Ninth column: sstart
        send = int(fields[9])  # Tenth column: send

        # Only process sseqid starting with "LG"
        if sseqid.startswith("LG"):
            # Extract segmentation range from sseqid
            match = re.match(r"(\w+)_(\d+)-(\d+)", sseqid)
            if match:
                segment_start = int(match.group(2))  # Extract the left boundary of the segment
                segment_end = int(match.group(3))    # Extract the right boundary of the segment
            else:
                segment_start = 0
                segment_end = 0

            # Convert relative positions to absolute positions
            absolute_sstart = sstart + segment_start
            absolute_send = send + segment_start

            # Replace relative positions with absolute positions
            fields[8] = str(absolute_sstart)
            fields[9] = str(absolute_send)

        # Write updated line to the output file (LG updated, others unchanged)
        outfile.write("\t".join(fields) + "\n")
