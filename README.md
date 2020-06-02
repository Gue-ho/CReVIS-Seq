# CReVIS-Seq

Lentiviruses have been widely used as a means of transferring exogenous DNAs into human cells to treat various genetic diseases. Lentiviral vectors are fundamentally integrated into the host genome but their integration sites are generally unpredictable, which may increase the uncertainty about their use in therapeutics. To determine viral integration sites in the host genome, several PCR-based methods have been developed. However, such methods require a well-designed primer set and are less accurate especially for cases in which a high copy number of viral DNAs are inserted in a clone or bulk cell populations. Here we describe CReVIS-seq, a highly efficient, genome-wide method for detecting viral insertion sites using high-throughput sequencing, which is based on in vitro circularization and subsequent cleavage of genomic DNAs in a CRISPR guide RNA-specific manner. Because CReVIS-seq does not require PCR amplification for the target enrichment step, it is free from PCR-derived biases and can be easily used to identify multiple target sites simultaneously, even in bulk cell populations. Furthermore, because of its versatile nature, CReVIS-seq can also be used to detect other structural variations that occur in various genomes. 

citation : Kim HS et al. CReVIS-seq: a highly accurate and multiplexable method for genome-wide mapping of lentiviral integration sites

# Usage

CReVIS-Seq requires Python3 and reference file which is indexed to bwa.

CRISPR-Sub can run with:

    python3 CReVIS_seq.py {LTR sequence} {target sequence} {NGS file} {output file name} {reference genome directory}

Example Code:
 
 You need to download and index the reference file to run bwa.
 
    python3 CReVIS_seq.py GGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA ACACTGACTAAAAGGGTCTG 10_1t.fastq -input2 10_2t.fastq result /directory/of/your/reference/genome/

# License
-------
CReVIS-Seq is licensed under the new BSD licence.

Copyright (c) 2019, Gue-ho Hwang and Sangsu Bae
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the name of the Hanyang University nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNERS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
