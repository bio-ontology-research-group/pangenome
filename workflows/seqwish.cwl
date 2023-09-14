cwlVersion: v1.1
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/seqwish:0.7.9--h43eeafb_2"
  ResourceRequirement:
    coresMin: 4
    ramMin: $(8 * 1024)
    outdirMin: $(Math.ceil(inputs.seqsFA.size/(1024*1024*1024) + 20))
baseCommand: seqwish
arguments: [-t, $(runtime.cores),
            -k, $(inputs.kmerSize),
            -s, $(inputs.seqsFA),
            -p, $(inputs.seqsPAF),
            -g, $(inputs.seqsPAF.nameroot).gfa]

inputs:
  seqsFA: File
  seqsPAF: File
  kmerSize:
    type: int
    default: 16
outputs:
  seqwishGFA:
    type: File
    outputBinding:
      glob: $(inputs.seqsPAF.nameroot).gfa
