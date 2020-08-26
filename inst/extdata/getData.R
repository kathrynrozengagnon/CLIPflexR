require(SRAdb)
mySRA <- getSRAdbFile()
sra_con <- dbConnect( dbDriver("SQLite"), mySRA )
Fox3_Std <-  getFASTQfile("SRR1107535",sra_con = sra_con,destDir = getwd())
fl <- "/Users/thomascarroll/Desktop/Projects/Hatten/counts/SRR1107535.fastq.gz"
require(ShortRead)
f <- FastqSampler(fl, 10^6)
big <- yield(f) 
f <- FastqSampler(fl, 10^5)
med <- yield(f) 
f <- FastqSampler(fl, 10^4)
small <- yield(f) 
writeFastq(big,"Fox3_Std_big.fq.gz")
writeFastq(med,"Fox3_Std_med.fq.gz")
writeFastq(small,"Fox3_Std_small.fq.gz")

require(BSgenome.Hsapiens.UCSC.hg19)
require(magrittr)
DNAStringSet(list(chr2=BSgenome.Hsapiens.UCSC.hg19[["chr2"]][(10^7):((10^7)+10^5)],
             chr6=BSgenome.Hsapiens.UCSC.hg19[["chr6"]][(10^7):((10^7)+10^5)])) %>% 
writeXStringSet(filepath = "hg19Small.fa")
