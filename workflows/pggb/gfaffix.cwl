cwlVersion: v1.1
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/gfaffix:0.1.5--h031d066_0"
  ResourceRequirement:
    coresMin: 4
    ramMin: $(8 * 1024)
baseCommand: gfaffix
arguments: [-o, $(inputs.GFA.nameroot).gfaffix.gfa,
            -t, $(inputs.GFA.nameroot).trans.gfa,
            $(inputs.GFA),
           ]
            
inputs:
  GFA: File
outputs:
  gfaffixGFA:
    type: File
    outputBinding:
      glob: $(inputs.GFA.nameroot).gfaffix.gfa
  gfatransGFA:
    type: File
    outputBinding:
      glob: $(inputs.GFA.nameroot).trans.gfa
