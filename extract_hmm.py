import argparse
from pathlib import Path
import logging
from io import StringIO
from itertools import chain

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)


def argparser():
    parser = argparse.ArgumentParser(
        description="Extract listed HMM profiles from Pfam-A.hmm"
    )
    parser.add_argument(
        "--pfam", type=Path, required=True, help="Path to Pfam-A.hmm"
    )
    filter_method_args = parser.add_argument_group(
        title="How to filter accessions",
        argument_default="",
    )
    filter_method_args.add_argument(
        "--accessions",
        type=str,
        help="Path of accession file, or a string of accessions separated by ','",
    )
    filter_method_args.add_argument(
        "--filter_by_name",
        type=str,
        help="NAME string of the HMM profile, separated by ','",
    )
    filter_method_args.add_argument(
        "--filter_by_description",
        type=str,
        help="DESC string of the HMM profile, separated by ','",
    )
    filter_method_args.add_argument(
        "--filter_by_length",
        type=str,
        help="Length of the HMM profile, separated by ',', e.g. '100,200'",
    )
    filter_method_args.add_argument(
        "--combine_logic",
        type=str,
        choices=["OR", "AND"],
        default="OR",
        help="Logic to combine filters, default to 'OR'",
    )
    parser.add_argument(
        "--output", type=Path, required=True, help="Path to output file"
    )

    return parser


def get_hmm_profile_block(fileio: StringIO):
    block = []
    profile = {}
    for line in fileio:
        block.append(line)
        if line.startswith("NAME"):
            profile["NAME"] = line.split()[1]
        if line.startswith("ACC"):
            profile["ACC"] = line.split()[1]
        if line.startswith("DESC"):
            profile["DESC"] = line.split()[1]
        if line.startswith("LENG"):
            profile["LENG"] = int(line.split()[1])
        if line.startswith("NSEQ"):
            profile["NSEQ"] = int(line.split()[1])
        if line.startswith("EFFN"):
            profile["EFFN"] = float(line.split()[1])
        if line.startswith("GA"):
            profile["GA"] = float(line.split()[1])
        if line.startswith("TC"):
            profile["TC"] = float(line.split()[1])
        if line.startswith("NC"):
            profile["NC"] = float(line.split()[1])
        if line.startswith("//"):
            yield block, profile
            block = []


def extract_hmm_profiles_by(
    by: str, target_set: set, input_file: Path
) -> dict[str, dict[str, list]]:
    extracted = {}
    # extracted = {target: {acc: block}}
    # target from target_set
    # acc from profile["ACC"]
    # block is list of lines of one profile
    with input_file.open("r") as infile:
        for block, profile in get_hmm_profile_block(infile):
            t_string = profile[by]
            is_target = False
            match by:
                case "ACC":
                    for t in target_set:
                        if t_string.startswith(t):
                            target = t
                            is_target = True
                            break
                case "NAME" | "DESC":
                    for t in target_set:
                        if t in t_string:
                            target = t
                            is_target = True
                            break
                case "LENG" | "NSEQ" | "EFFN":
                    target_set = sorted(list(target_set))
                    target = f"{by}_{target_set[0]}_{target_set[1]}"
                    assert (
                        len(target_set) == 2
                    ), f"A range (a, b) is required when filter by {by}"
                    assert all(
                        [isinstance(t, int | float) for t in target_set]
                    ), f"A range (a, b) is required when filter by {by}"
                    if target_set[0] <= profile[by] < target_set[1]:
                        is_target = True
                case "GA":
                    raise NotImplementedError
                case "TC":
                    raise NotImplementedError
                case "NC":
                    raise NotImplementedError
                case _:
                    raise ValueError(f"Invalid filter method: {by}")

            if is_target:
                extracted.setdefault(target, {})[profile["ACC"]] = block
                is_target = False

    return extracted


def write_extracted_hmm_profiles(profiles: dict[str, list], output_file: Path):
    accessions = set()
    with output_file.open("w") as outfile:
        for acc, block in profiles.items():
            outfile.writelines(block)
            accessions.add(acc)
    return accessions


def main():
    parser = argparser()
    args = parser.parse_args()
    assert args.pfam.is_file(), "Pfam-A.hmm file not found"
    assert args.output.parent.is_dir(), "Output folder not found"
    pfam_file = args.pfam

    # accession
    extracted_acc = {}
    if args.accessions != "":
        if Path(args.accessions).is_file():
            logging.info("Parse accessions from file")
            accessions = set()
            with open(args.accessions, "r") as infile:
                for line in infile:
                    if line.startswith("Accession"):
                        continue
                    accessions.add(line.split("\t")[0].strip())
        else:
            logging.info("Accessions not found as file, use as string")
            accessions = set(args.accessions.split(","))
        assert (
            len(accessions) > 0
        ), f"No accessions found using {args.accessions}"
        logging.info("Filter by accessions...")
        extracted_acc.update(
            extract_hmm_profiles_by("ACC", accessions, pfam_file)
        )
        if len(accessions) > len(extracted_acc):
            not_found = set(accessions) - set(extracted_acc)
            logging.warning(
                "Some accessions were not found in the HMM file.\n"
                f"Accessions not found: {', '.join(not_found)}"
            )
        logging.info(f"Accessions found {len(extracted_acc)} profiles.")

    # name
    extracted_name = {}
    if args.filter_by_name != "":
        logging.info("Filter by names...")
        names = set(args.filter_by_name.split(","))
        extracted_name.update(extract_hmm_profiles_by("NAME", names, pfam_file))
        for name in names:
            if name not in extracted_name:
                logging.warning(f"Name {name} not found in the HMM file.")
            else:
                logging.info(
                    f"Name {name} found {len(extracted_name[name])} profiles."
                )

    # description
    extracted_desc = {}
    if args.filter_by_description != "":
        logging.info("Filter by descriptions...")
        descriptions = set(
            a.strip("\"'") for a in args.filter_by_description.split(",")
        )
        extracted_desc.update(
            extract_hmm_profiles_by("DESC", descriptions, pfam_file)
        )
        for desc in descriptions:
            if desc not in extracted_desc:
                logging.warning(
                    f"Description {desc} not found in the HMM file."
                )
            else:
                logging.info(
                    f"Description {desc} found {len(extracted_desc[desc])} profiles."
                )

    # length
    extracted_len = {}
    if args.filter_by_length != "":
        logging.info("Filter by lengths...")
        lengths = set([int(a) for a in args.filter_by_length.split(",")])
        extracted_len.update(
            extract_hmm_profiles_by("LENG", lengths, pfam_file)
        )
        if len(extracted_len) == 0:
            logging.error(
                f"Length {args.filter_by_length} not found in the HMM file."
            )
        else:
            for n in extracted_len:
                profiles = extracted_len[n]
                logging.info(f"{n} found {len(profiles)} profiles.")

    all_extracted = [
        extracted_acc,
        extracted_name,
        extracted_desc,
        extracted_len,
    ]
    # write to output file
    target_profiles = {}  # {acc: block}
    logging.info(
        f"Gather all filtered profiles together (logic = {args.combine_logic}), "
        "writing extracted profiles to output file..."
    )
    if args.combine_logic == "OR":
        for extracted in all_extracted:
            for target, profiles in extracted.items():
                target_profiles.update(profiles)
    else:
        accessions = set()
        for extracted in all_extracted:
            if len(extracted) == 0:
                continue  # skip empty extracted
            extracted_keys = set(
                chain.from_iterable(set(v.keys()) for v in extracted.values())
            )
            if not accessions:
                accessions = extracted_keys
            else:
                accessions &= extracted_keys

        for extracted in all_extracted:
            if len(extracted) == 0:
                continue  # skip empty extracted
            for acc in accessions:
                for target, profiles in extracted.items():
                    if acc in profiles:
                        target_profiles[acc] = profiles[acc]

    accessions = write_extracted_hmm_profiles(target_profiles, args.output)
    output_rename = args.output.parent / (
        (args.output.stem if args.output.suffix == ".hmm" else args.output.name)
        + f"_n{len(accessions)}.hmm"
    )
    args.output.rename(output_rename)

    logging.info(f"Total {len(accessions)} profiles extracted.")


if __name__ == "__main__":
    main()
