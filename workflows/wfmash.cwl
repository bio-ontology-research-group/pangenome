cwlVersion: v1.1
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/wfmash:0.10.5--h94f6cfe_0"
  ResourceRequirement:
    coresMin: 8
    ramMin: $(4 * 1024)
stdout: $(inputs.seqsFA.nameroot).paf
baseCommand: wfmash
arguments: [-t, $(runtime.cores),
            -k, $(inputs.kmerSize),
            -s, $(inputs.readsFA),
            -p, $(inputs.readsPAF),
            -g, $(inputs.readsPAF.nameroot).gfa]

inputs:
  seqsFA: File
  stindex: File
  bgindex: File
outputs:
  wfmashPAF:
    type: File
    outputBinding:
      glob: $(inputs.seqsFA.nameroot).paf
