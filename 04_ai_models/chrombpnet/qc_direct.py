# chromBPNet v1.0.1's `chrombpnet qc` creates <out>/auxiliary with exist_ok=False and then reads
# <out>/auxiliary/filtered.peaks.bed, so it can never run. Call the qc pipeline directly on the training directory.
import sys
import chrombpnet.parsers as parsers
import chrombpnet.pipelines as pipelines
sys.argv = ["chrombpnet", "qc"] + sys.argv[1:]
args = parsers.read_parser()
pipelines.chrombpnet_qc(args)
