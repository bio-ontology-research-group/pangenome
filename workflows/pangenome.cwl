cwlVersion: v1.1
class: Workflow
requirements:
  ResourceRequirement:
    ramMin: $(8 * 1024)
    coresMin: 2
  MultipleInputFeatureRequirement: {}
  
inputs:
  seqsFA: File

outputs:
  graphGFA:
    type: File
    outputSource: seqwish/seqwishGFA
  alignment:
    type: File
    outputSource: wfmash/wfmashPAF

steps:
  wfmash:
    in:
      seqsFA: seqsFA
    out: [wfmashPAF]
    run: wfmash.cwl
  seqwish:
    in:
      seqsFA: seqsFA
      seqsPAF: wfmash/wfmashPAF
    out: [seqwishGFA]
    run: seqwish.cwl