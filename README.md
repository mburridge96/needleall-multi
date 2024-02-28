# Overview

Simple cli tool for running EMBOSS needleall across multiple cores

# Usage
```bash
Cli tool for multi-process Emboss NeedleAll

Usage: needleall-multi [OPTIONS] --fasta <FASTA>

Options:
  -f, --fasta <FASTA>
          Path to target fasta file
  -w, --working-dir <WORKING_DIR>
          Working directory. Default value creates a wd in %Y%m%d%H%M%S format [default: date/]
  -o, --outfile <OUTFILE>
          Identities output file name [default: identities.tsv]
  -e, --errorfile <ERRORFILE>
          Needle all error file name [default: needle_error.error]
  -g, --gap-open-penalty <GAP_OPEN_PENALTY>
          Gap open penalty [default: 10.0]
  -g, --gap-extend-penalty <GAP_EXTEND_PENALTY>
          Gap extend penalty [default: 0.5]
  -t, --threshold <THRESHOLD>
          Threshold for result [default: -1.0]
  -c, --cpu-count <CPU_COUNT>
          Set cpu number for multi-processing [default: 0]
  -d, --debug-time
          Time debug
  -h, --help
          Print help
  -V, --version
          Print version
```
