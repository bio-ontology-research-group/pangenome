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
arguments: [build, -g, $(inputs.GFA), -o, -,
            {shellQuote: false, valueFrom: "|"},
            odgi, sort, -i, -, -p, s, -o,
            $(inputs.GFA.nameroot).odgi]
inputs:
  GFA: File
outputs:
  odgiGraph:
    type: File
    outputBinding:
      glob: $(inputs.GFA.nameroot).odgi