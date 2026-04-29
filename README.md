# Single-molecule variation in telomeric sequence and structure across humans

### Authors
Danilo Dubocanin**¹**, Mitchell R. Vollger**²**, Shane J Neph**³**, Maria Sara Del Rio Pisula**¹ ⁴**, Julian K Lucas**⁵**, Adriana E Sedeño-Cortés**³**, Ben J Mallory**³ ⁶**, Taylor D Real**³ ⁶**, Human Pangenome Reference Consortium, Nicolas Altemose**¹ ⁷** <sup>§</sup>, Andrew B. Stergachis**³ ⁶ ⁸** <sup>§</sup>

---

### Affiliations
1. **Department of Genetics**, Stanford University, Palo Alto, CA, USA
2. **Department of Human Genetics**, University of Utah, Salt Lake City, UT, USA
3. **Division of Medical Genetics**, Department of Medicine, University of Washington, Seattle, WA, USA
4. **Department of Biology**, Stanford University, Stanford, CA, USA
5. **UC Santa Cruz Genomics Institute**, University of California, Santa Cruz, Santa Cruz, CA, USA
6. **Department of Genome Sciences**, University of Washington, Seattle, WA, USA
7. **Biohub – San Francisco**, San Francisco, CA, USA
8. **Brotman Baty Institute for Precision Medicine**, Seattle, WA, USA

---

### Abstract
 The repetitive architectures of telomeric and subtelomeric regions have obscured studies of their genetic variation and chromatin organization across the human population. Here, we integrate near-complete diploid genome assemblies from 212 individuals with matched long-read sequencing data to construct an atlas of 316,146 telomere-spanning molecules across 12,080 chromosome-end-resolved telomere arrays. This atlas reveals that nearly every chromosome end harbors a structured and unique pattern of telomere variant repeats (TVR), or TVR code, with subtelomere-proximal TVR codes being heritable, somatically stable, and influenced by subtelomeric TAR1 regulatory elements. Despite ongoing cycles of telomere shortening and elongation in the germline, proximal TVR codes are maintained across the human population. These TVR codes expose rare telomerase-independent events that lengthen telomeres in the germline, including interchromosomal telomere exchange and recurrent internal duplications within telomere arrays. Furthermore, single-molecule chromatin fiber sequencing across 26,972 molecules spanning the telomere-subtelomere boundary confirms that TVR-rich regions adopt telomeric chromatin but introduce discrete discontinuities into otherwise compact telomeric chromatin fibers. Together, our results link chromosome-end sequence variation to telomere cap formation and telomerase-independent telomere extension mechanisms in the human germline.

---

### Data Availability

* **Processed data files for analysis:** Available via our [interactive web directory](https://s3.kopah.uw.edu/dubocd/index.html).
* **HPRC assemblies and standard reads:** ONT (R10) and PacBio (HiFi) reads are available at the [Human Pangenome raw sequencing data repository](https://data.humanpangenome.org/raw-sequencing-data).
* **Fiber-seq for HPRC samples:** Available via the [Human Pangenomics S3 directory](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/5ECA1D3E-1C37-44B7-BA20-576417F786F0--UCSC_HPRC_ONT_YEAR1_FIBERSEQ/).
* **Fiber-seq for CHM13 and HG002:** Available from GEO accessions [GSM7074431](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7074431) and [GSM7074433](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7074433), respectively.

---

### Dependencies

This analysis is performed using relatively standard bioinformatics and data science packages.

python (>=3.11) dependencies:
* numpy           1.24.4 
* pandas          1.5.3 
* tqdm            4.64.1 
* pybedtools      0.9.0
* matplotlib      3.7.1 
* seaborn         0.12.2  
* datasketch      1.6.5 
* edlib           1.3.9.post1 
* scipy           1.11.1 

other tools used in this analysis:

* **seqtk** available at https://github.com/lh3/seqtk. We used seqtk to call telomere-subtelomere junctions.
* **fibertools** available at https://github.com/fiberseq/fibertools-rs. A CLI tool for creating and interacting with Fiber-seq BAM files. We use various fibertools outputs to extract sequence and chromatin features from single molecules.

### Running the code 

To run all analyses in this paper, you can download all processed files from our [interactive web directory](https://s3.kopah.uw.edu/dubocd/index.html). Then, you can download the jupyter notebook corresponding to your analysis of interest and replace hard-coded file names with the location of the downloaded file on your system.
Alternatively, you can also download bam files from the [Human Pangenome raw sequencing data repository](https://data.humanpangenome.org/raw-sequencing-data), and process them as described in our methods section.
