#!/usr/bin/env python3

import os


def parse_name_and_length(line: str):
    parts = line.rstrip("\n").split("\t")
    if not parts:
        return None, None
    name = parts[0]
    length_str = parts[-1]
    try:
        length = int(length_str)
    except ValueError:
        return None, None
    return name, length


def main():
    lengths_path = "seq_lengths_per_sequence.txt"
    bed_path = "split_sequences.bed"

    if not os.path.exists(lengths_path):
        raise SystemExit(f"Missing {lengths_path}")

    with open(lengths_path, "r", encoding="utf-8") as f_in, open(
        bed_path, "w", encoding="utf-8"
    ) as f_out:
        for line in f_in:
            line = line.strip()
            if not line:
                continue
            name, length = parse_name_and_length(line)
            if name is None or length is None:
                continue

            # Split into contiguous chunks strictly less than 1GB (allow 1GB-1)
            max_segment_len = 1_000_000_000 - 1  # 999,999,999

            seg_start = 0
            while seg_start < length:
                seg_end = min(seg_start + max_segment_len, length)
                if seg_end > seg_start:
                    f_out.write(f"{name}\t{seg_start}\t{seg_end}\n")
                seg_start = seg_end


if __name__ == "__main__":
    main()


