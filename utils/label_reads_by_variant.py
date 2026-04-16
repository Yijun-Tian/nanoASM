import argparse
import pysam


def parse_variants(variant_str):
    """Parse variants into a list of tuples:
    (variant_id, ref_allele, alt_alleles, chromosome, start_coordinate)
    """
    variants = []

    for variant in variant_str.split(","):
        variant_data = variant.split(":")
        if len(variant_data) < 5:
            raise ValueError(
                "Each variant must include at least: id, ref, one or more alt alleles, chromosome, and start coordinate."
            )

        variant_id = variant_data[0]
        ref_allele = variant_data[1]
        alt_alleles = variant_data[2:-2]
        chromosome = variant_data[-2]
        start_coordinate = int(variant_data[-1])

        variants.append(
            (variant_id, ref_allele, alt_alleles, chromosome, start_coordinate)
        )

    return variants


def main(input_bam, output_bam, variants):
    bam_in = pysam.AlignmentFile(input_bam, "rb")
    bam_out = pysam.AlignmentFile(output_bam, "wb", header=bam_in.header)

    for variant_id, ref_allele, alt_alleles, chromosome, start in variants:
        stop = start + 1

        for read in bam_in.fetch(chromosome, start, stop):
            if read.query_sequence is None:
                continue

            read_positions = read.get_reference_positions(full_length=True)

            if start in read_positions:
                read_index = read_positions.index(start)
                read_base = read.query_sequence[read_index]
                tag_info = f"{variant_id}@{read_base}"

                matched = False
                if read_base == ref_allele:
                    matched = True
                else:
                    for alt in alt_alleles:
                        if read_base == alt:
                            matched = True
                            break

                if not matched:
                    tag_info = f"{variant_id}@{read_base}"

                if read.has_tag("VA"):
                    existing_tag = read.get_tag("VA")
                    new_tag = f"{existing_tag};{tag_info}"
                    read.set_tag("VA", new_tag)
                else:
                    read.set_tag("VA", tag_info)

            if (
                not read.is_secondary
                and not read.is_supplementary
                and read.has_tag("VA")
            ):
                bam_out.write(read)

    bam_in.close()
    bam_out.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Annotate BAM reads with variant allele tags."
    )
    parser.add_argument("-i", "--input_bam", required=True, help="Input BAM file")
    parser.add_argument("-o", "--output_bam", required=True, help="Output BAM file")
    parser.add_argument(
        "-v",
        "--variants",
        required=True,
        help="Variants in the format: id:ref:alt1[:alt2...]:chromosome:start[,more_variants]",
    )

    args = parser.parse_args()
    variant_list = parse_variants(args.variants)
    main(args.input_bam, args.output_bam, variant_list)
