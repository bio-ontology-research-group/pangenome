cwlVersion: v1.1
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/seqwish:0.7.9--h43eeafb_2"
  ResourceRequirement:
    coresMin: 4
    ramMin: $(7 * 1024)
    outdirMin: $(Math.ceil(inputs.readsFA.size/(1024*1024*1024) + 20))
stdout: $(inputs.readsFA.nameroot).paf
baseCommand: seqwish
arguments: [-t, $(runtime.cores),
            -k, $(inputs.kmerSize),
            -s, $(inputs.readsFA),
            -p, $(inputs.readsPAF),
            -g, $(inputs.readsPAF.nameroot).gfa]

inputs:
  readsFA: File
  readsPAF: File
  kmerSize:
    type: int
    default: 16
outputs:
  seqwishGFA:
    type: File
    outputBinding:
      glob: $(inputs.readsPAF.nameroot).gfa
