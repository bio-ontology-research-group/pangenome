cwlVersion: v1.1
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/smoothxg:0.7.1--h40c17d1_0"
  ResourceRequirement:
    coresMin: 4
    ramMin: $(8 * 1024)
baseCommand: smoothxg
arguments: [-t, $(runtime.cores),
            -g, $(inputs.GFA),
            -r, $(inputs.nHaps),
            -o, $(inputs.GFA.nameroot).smooth.gfa,
            ]
            
inputs:
  GFA: File
  nHaps: int
outputs:
  smoothGFA:
    type: File
    outputBinding:
      glob: $(inputs.GFA.nameroot).smooth.gfa
