# Analysis Code for Cux2 and Atf4 paper
"DNA damage burden causes selective CUX2 neuron loss in neuroinflammation"
"Expansion of outer cortical layer Cux2+ neurons required selective adaptations for DNA repair"

This repository contains analysis scripts used to generate the results and figures reported in the associated manuscript.  
The code is provided for transparency and reproducibility and is not intended as a standalone software tool.

## System requirements

### Operating systems
- Linux with NVIDIA GPU support
- macOS (Apple silicon / ARM)

### Software dependencies
- Python (≥ 3.9)
- GPU-accelerated deep learning framework (e.g. PyTorch)
- NVIDIA CUDA Toolkit (Linux only)

### Tested platforms
- Ubuntu 20.04 with NVIDIA GPU and CUDA support
- macOS 13 on Apple silicon (ARM)

### Hardware requirements
GPU acceleration is required to run the full analysis efficiently.  
CPU-only execution is not supported or may be impractically slow.

## Installation
No formal installation is required beyond installing the dependencies listed above.  
The repository can be cloned using:
```bash
git clone https://github.com/USERNAME/REPOSITORY_NAME.git
cd REPOSITORY_NAME
```

## Usage

The scripts correspond to the analyses and figures reported in the manuscript.
Users may adapt the code to their own data by modifying input files and parameters as indicated in the scripts.

## Reproducibility

The provided scripts allow reproduction of the main quantitative results in the manuscript using the specified datasets.
Public datasets are referenced by accession numbers in the manuscript.

## Notes

No standalone demo dataset is provided, as the code is tailored to the analyses in this study.

## License

MIT License
