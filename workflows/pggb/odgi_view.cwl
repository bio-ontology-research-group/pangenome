cwlVersion: v1.1
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/odgi:0.8.3--py310h6cc9453_0"
  ResourceRequirement:
    coresMin: 4
    ramMin: $(8 * 1024)
baseCommand: odgi
arguments: [view, --idx, $(inputs.GFA), --to-gfa, -o,
            $(inputs.GFA.nameroot).gfa]
inputs:
  GFA: File
outputs:
  odgiView:
    type: File
    outputBinding:
      glob: $(inputs.GFA.nameroot).gfa