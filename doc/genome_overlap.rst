Methods
=======

In order to assess the completeness of our de Novo transcriptome assembly (lamp10), we have compared the alignment of the generated transcripts against the existing genome annotations released with Pmarinus v7.0.75. First, we use blastn to align transcripts to the genome, using parameters `-evalue 0.000001`. Then, we use the coordinates from the annotation and the corresponding coordinates from the alignments to calculate the proportion of annotated sequence overlapped, proportion of transcripts overlapped, and the respective proportions of non-overlapped sequence and transcripts. We consider an annotated region to be overlapped by a trancript if it is at least 90% covered, with at least 98% identity [TODO: get better justification for these cutoffs other than "things Camille remembers reading"].

We give particular attention to alignments which entirely contain annotated regions, as these suggest extensions to existing annotations. When these alignments are from transcripts with homology evidence from other species, we consider them to represent putative extensions. Further, alignments which are entirely contained within an annotation suggest either an overly aggressive prediction in the genome, or an incompletely assembled transcript. 

Results
=======

Using our 90% overlap / 98% identity heuristic, We find [x]% of annotations to be covered by a transcript from lamp10, while [y]% of transcripts cover an annotation. Futher, []% of the genome is coveredby annotations, while [%] is covered by alignments from lamp10; []% of transcripts have any alignment to the genome.

We also find that []% of transcript alignments entirely contain an annotation, increasing the annotation size by [blah]%. []% of extensions are supported by homology to a known protein. []% of transcipt alignments are entirely contained by an annotation.

[figure: bar chart of proportions; histogram of extension percentages or raw sequence lengths]

Discussion
===========


