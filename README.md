# DRUG-seq toolbox

This is a toolbox for processing raw sequencing output from multiplexed RNA-seq experiments into a digital gene expression matrix that will contain integer counts of the number of transcripts per barcode. It provides a bash script that serves as a wrapper for multiple analysis steps, including demultiplexing of the raw data by molecular (UMI) and sample / well barcode, filtering barcodes by a list of expected barcodes, alignment of reads to a reference genome, collecting basic QC metrics and counting UMIs per barcode.

The pipeline uses some custom python scripts, and many tools from the [Drop-seq toolbox](https://github.com/broadinstitute/Drop-seq) (Mc Caroll lab, Harvard Medical school) as well as [Picard](https://broadinstitute.github.io/picard/) (Broad institute), which are all included in this toolbox.

Please refer to the [user manual](./man/Drugseq_toolbox_manual.pdf) for a description of how to use it.

Contact
This toolbox is provided by the [Snijder lab](https://www.snijderlab.org/) at ETH Zurich.

If you have questions or find a bug, son't hesitate to contact Rebekka Wegmann: wegmann@imsb.biol.ethz.ch

Notes and caution
This tool comes with no warranty and accurate function is not guaranteed.
