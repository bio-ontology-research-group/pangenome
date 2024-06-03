cwlVersion: v1.1
class: Workflow
requirements:
  ResourceRequirement:
    ramMin: $(8 * 1024)
    coresMin: 2
  MultipleInputFeatureRequirement: {}
  
inputs:
  seqsFA: File
  nHaps: int

outputs:
  graphGFA:
    type: File
    outputSource: seqwish/seqwishGFA
  alignment:
    type: File
    outputSource: wfmash/wfmashPAF
  odgi_graph:
    type: File
    outputSource: odgi_view/odgiView

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
  smoothxg:
    in:
      GFA: seqwish/seqwishGFA
      nHaps: nHaps
    out: [smoothGFA]
    run: smoothxg.cwl
  gfaffix:
    in:
      GFA: smoothxg/smoothGFA
      nHaps: nHaps
    out: [gfaffixGFA,gfatransGFA]
    run: gfaffix.cwl
  odgi_build:
    in:
      GFA: gfaffix/gfaffixGFA
    out: [odgiGraph]
    run: odgi_build.cwl
  odgi_unchop:
    in:
      GFA: odgi_build/odgiGraph
    out: [odgiUnchop]
    run: odgi_unchop.cwl
  odgi_sort:
    in:
      GFA: odgi_unchop/odgiUnchop
    out: [odgiSorted]
    run: odgi_sort.cwl
  odgi_view:
    in:
      GFA: odgi_sort/odgiSorted
    out: [odgiView]
    run: odgi_view.cwl