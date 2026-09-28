#!/bin/bash
# Lean scPrinter / seq2PRINT environment in project space (the README recipe's RAPIDS stack is only for GPU chromVAR/UMAP).
set -euo pipefail
module load miniforge/25.11.0-0
export PIP_NO_CACHE_DIR=1      # pip's download cache would go to the small home directory
DATA=$(realpath -m ../results/scprinter_data)     # scPrinter's own downloads (genome, Tn5 bias, motifs; about 16 GB)
# the env goes to the first folder in envs_dirs of ~/.condarc: put a project-space folder there (see the main README)
conda create -y -n scprinter -c conda-forge -c bioconda python=3.11 macs2 macs3 bedtools htslib samtools pybedtools pybigwig \
    'numpy<2' pip ipykernel git weasyprint
conda activate scprinter
python -c "import sys; print(sys.version)"
pip install torch --index-url https://download.pytorch.org/whl/cu121
pip install "git+https://github.com/buenrostrolab/scPrinter@v1.2.0" modisco-lite
pip install "git+https://github.com/austintwang/finemo_gpu.git"
# scPrinter 1.2.0 imports tangermeme.tools.tomtom, which tangermeme >= 0.5 removed; conda's pyBigWig is built for numpy 1
pip install "tangermeme==0.4.4"
pip install --force-reinstall --no-deps "pyBigWig>=0.3.23"
# scPrinter 1.2.0 calls snapatac2.pp.import_data, renamed in snapATAC2 2.9; 2.8.0 needs numpy 1, so shap/zarr too
pip install "snapatac2==2.8.0" "shap<0.47" "zarr<3"
# imported by the training script but not declared as a dependency
pip install ema_pytorch
python -c "import scprinter, torch; print('scprinter', scprinter.__version__ if hasattr(scprinter,'__version__') else 'ok', '| torch', torch.__version__, '| CUDA build', torch.version.cuda)"
which seq2print_train macs3 bgzip tabix
conda env config vars set -n scprinter SCPRINTER_DATA=$DATA   # else scPrinter downloads into ~/.cache
conda env config vars set -n scprinter LD_LIBRARY_PATH=$CONDA_PREFIX/lib   # else the system libstdc++ (too old for scipy) loads first
