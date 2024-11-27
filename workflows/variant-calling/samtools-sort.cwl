cwlVersion: v1.1
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
  InitialWorkDirRequirement:
    listing: [ $(inputs.aligned_reads) ]

  DockerRequirement:
    dockerPull: "coolmaksat/minimap2-samtools"
  ResourceRequirement:
    coresMin: 2
    ramMin: $(4 * 1024)
baseCommand: samtools
arguments: [view, -h, $(inputs.aligned_reads),
           {shellQuote: false, valueFrom: '|'},
            sed, -e, "s/GRCh38#0#//g",
           {shellQuote: false, valueFrom: '|'},
           samtools, sort, -@10, -O, CRAM, -T, /tmp/readsPrefix, -,
           ]
inputs:
  aligned_reads: File

outputs:
  aligned_reads_sorted:
    type: stdout
  
stdout: $(inputs.aligned_reads.nameroot + '.sorted.cram')
