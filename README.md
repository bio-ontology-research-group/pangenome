# JASAPAGE: Phased genome assemblies and pangenome graphs of human populations of Japan and Saudi Arabia

# Saudi and Japanese Pangenome Analysis Pipeline

## Overview
This repository contains workflows and data for constructing pangenomes from Saudi and Japanese population genomic data using the minigraph-cactus pipeline. The workflows are implemented using Common Workflow Language (CWL), enabling reproducible and portable execution across different computing environments.

## Background

### Minigraph-Cactus Pipeline
Minigraph-cactus is a state-of-the-art pipeline that combines the efficiency of minigraph with the accuracy of Progressive Cactus for constructing pangenomes. This hybrid approach offers several advantages:

- Scalable construction of genome graphs from multiple assemblies
- Preservation of structural variants and complex genomic regions
- Efficient handling of large-scale genomic data
- Integration of both reference-based and reference-free approaches

The pipeline first uses minigraph to create initial graphs, followed by Progressive Cactus's sophisticated algorithms for multiple genome alignment, resulting in high-quality pangenome representations.

### Common Workflow Language (CWL)
CWL is a specification for describing analysis workflows and tools. Its key features include:

- Platform independence
- Explicit declaration of dependencies
- Built-in support for Docker/container technologies
- Scalability from single workstations to HPC environments
- Clear separation of tools from the workflow logic


### Workflows

## Installation

### Prerequisites
You need to install Docker or Singularity.
```bash
# Install required software
pip install cwltool
```

### Repository Setup
```bash
git clone https://github.com/JaSaPaGe/pangenome-cwl
```

## Usage

### DeepVariant workflow

```bash
cwl-runner variant-calling/main-vg.cwl inputs.yml
```

### Configuration
Create a YAML file with your inputs:
```yaml
graph:
  class: File
  location: reference/ksa_jpn_hg38.gbz
reads1:
  class: File
  location: input/sample_R1.fastq.gz
reads2:
  class: File
  location: input/sample_R2.fastq.gz
hapl:
  class: File
  location: reference/ksa_jpn_hg38.hapl
ref_prefix: GRCh38
model_type: WGS
ref:
  class: File
  location: reference/hg38.fa
  secondaryFiles:
    - class: File
      location: reference/hg38.fa.fai
output_cram: aligned_reads.cram
output_vcf: sample.vcf
output_gvcf: sample.gvcf
```

## Output Files
- `aligned_reads.cram`: Alignment against GRCh38
- `sample.gvcf`: Genomic VCF
- `sample.vcf`: Variant calls


## Citation
If you use this pipeline in your research, please cite:

## License
The data is available under CC0

## Contact
Please create GitHub issues if you encounter any problems. 
