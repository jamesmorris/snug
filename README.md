snug
====

Interactive user interface written in Java built to streamline the quality control process of visually inspecting lists of variant sites called from short read sequencing data.

Getting Started
---------------

### Running the jar file
`java -jar snug.jar`

### Loading data
`File->Load files`

snug takes as input 3 files:

1. An indexed bam or a file containing the names of multiple indexed bams
2. A list of variants
3. An indexed reference sequence

#### Local/Remote


#### BAM files
+ All BAM files should be sorted and indexed.
+ In snug you can choose to either load a single BAM file or multiple BAMs.
+ To load multiple BAM files you need to create a BAM list file, which is just a text file containing the file name of each BAM you wish to load with each BAM on a new line.
+ In local mode BAM list files are assumed to contain the BAM file name and not the full path. As the bam list is assumed to be in the same directory as the bams.
+ In remote mode BAM list files should be local files that contain the full path to each BAM to be loaded on the remote server.
Bam list flies can also contain an optional second value for each bam which will be used as the bams name when it is displayed in snug (for example mother, father, child).

#### Variant files

The variants can be in any of the following common formats:

+ chr:position
+ chr-position
+ chr position
+ chr,position

#### Reference files

Results
-------
A results file containing the scores for each variant viewed will automaticaly be created in the same directory as the variants file with the same name plus a `.scores` extension.
The scores in the file correspond

+ -1 - No
+ 0  - Maybe
+ 1  - Yes


Display options
---------------
All these options are available under the `View` menu.

+ Collapsed view - reads are represented using lines rather than letters to reduce the space taken by each read.
+ Reference on top - display the reference sequence on the top of the display rather than the bottom.
+ Shade base by quality - 
+ Shade background by quality - 
+ Highlight low depth - highlight BAM files which have less than the minimum number of reads covering a variant.


FAQ
---
+ how do I index my BAM file?

`samtools index my.bam`
+ how do I index my reference file?

`samtools faidx my.fasta`
+ installing samtools link
+ 

Contact
-------
Please use the snug project issue tracker to report bugs or suggest features and inprovements.
