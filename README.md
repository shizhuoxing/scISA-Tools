# scIsoSeq-Analysis-Tools

High-throughput single-cell RNA sequencing is a powerful technique advance gene expression regulation research to a higher resolution level for single cells. Recent technological advancements in high-throughput scRNA-Seq methods such Droplet-based single-cell sequencing allow thousands of cells to be captured and sequenced in a relatively short time frame and at a fraction of the cost. Such methods rely on capture of polyadenylated mRNA transcripts followed by cDNA synthesis, pooling, amplification, library construction, and cDNA sequencing.

But the most widely used single cell RNA sequencing is only generates short reads from one end of a cDNA template, limiting the reconstruction of highly diverse isofromsuch as alternative splicing and structural variation. Recent advances in long-read sequencing technologies present a potential solution to the shortcomings of short-read sequencing. Full-length cDNA reads can encompass the entire sequence of transcripts, but typically suffer from higher error rates and lower sequencing depth than short-read technologies.

Based on the foundation of our previous works, combining single-cell full-length RNA library construction and MTZL library construction technology, we were able to significantly improve the throughput of single-cell RNA sequencing on the PacBio sequencing platform, now we can obtain sufficient high-accuracy reads to accurately parse cell barcode and full-length transcript information. It's really amazing.

See our wiki for more information on MTZL:
https://github.com/shizhuoxing/BGI-Full-Length-RNA-Analysis-Pipeline/wiki

In this repo, I will provide tools for analyzing MTZL scIsoSeq data, including basic data processing, cellBC extraction, expression matrix generation, and so on.
