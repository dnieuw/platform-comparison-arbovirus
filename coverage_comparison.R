library(ggridges)
library(data.table)
library(ggplot2)
library(stringr)

f2si<-function (number, rounding=F, digits=ifelse(rounding, NA, 6)) {
	lut <- c(1e-24, 1e-21, 1e-18, 1e-15, 1e-12, 1e-09, 1e-06, 
					 0.001, 1, 1000, 1e+06, 1e+09, 1e+12, 1e+15, 1e+18, 1e+21, 
					 1e+24, 1e+27)
	pre <- c("y", "z", "a", "f", "p", "n", "u", "m", "", "k", 
					 "M", "G", "T", "P", "E", "Z", "Y", NA)
	ix <- findInterval(number, lut)
	if (ix>0 && ix<length(lut) && lut[ix]!=1) {
		if (rounding==T && !is.numeric(digits)) {
			sistring <- paste0(round(number/lut[ix]), pre[ix])
		} else if (rounding == T || is.numeric(digits)) {
			sistring <- paste0(signif(number/lut[ix], digits), pre[ix])
		} else {
			sistring <- paste0(number/lut[ix], pre[ix])
		} 
	} else {
		sistring <- as.character(number)
	}
	return(sistring)
}

read_depthfile <- function(x){
	data <- fread(x)
	
	colnames(data) <- c("ref","position","coverage")
	
	sample_id <- str_remove(x, ".*/") %>% str_remove("\\..*")
	sample_info <- str_split(sample_id, "_", simplify=T)
	
	data[,c("sample_id","barcode", "approach", "platform") := 
			 	list(sample_id, sample_info[1], sample_info[2], sample_info[3])]
	return(data)
}

plot_ridges <- function(data){
	ggplot(data) + 
		geom_density_ridges_gradient(aes(x=position,
																		 y=CT,height=scaled_coverage,
																		 fill=coverage_comparison, colour=""),
																 stat="identity", 
																 scale=0.9) +
		scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(n = 11, "Spectral")),
												 trans="log10",
												 breaks=c(1,10,100,1000,10000,100000,1000000),
												 na.value = "gray"
		) +
		scale_colour_manual(values = alpha("black",0)) +              
		guides(colour=guide_legend("Coverage below\nthreshold")) +
		theme_minimal() +
		theme(panel.border = element_rect(color = "black", fill=NA),
					text = element_text(size=16)) + 
		labs(x="Genome position",
				 y="Scaled coverage depth",
				 fill="Coverage depth")
}

setwd("##PATH TO SNAKEMAKE OUTPUT FOLDER###")

ct_values <- fread("metadata_ct_value.csv")
setkey(ct_values, "barcode")

filelist <- list.files(path="depth", pattern = "*.tsv", full.names = T)
data_all <- rbindlist(lapply(filelist, read_depthfile))

#uncomment and comment below to get results before cleaning
# data <- data_all[grep("dedup|cleaned",sample_id, invert=T)]
data <- data_all[grep("cleaned",sample_id)]

data <- merge.data.table(data, ct_values, by = "barcode")

data[,scaled_coverage:=(coverage)/(max(coverage)),.(sample_id)]
data[is.infinite(scaled_coverage), scaled_coverage:=0]
data[,coverage_comparison:=coverage]
data[coverage<100 & platform=="Nanopore",coverage_comparison:=NA]
data[coverage<5 & platform!="Nanopore",coverage_comparison:=NA]
data[,scaled_position:=position/max(position),.(sample_id)]
data[,CT:=factor(CT,levels = c("CT33","CT29","CT25"), ordered = T)]

data[approach=="Meta", approach:="Metagenomic"]
data[approach=="Virocap", approach:="Capture"]
data[,approach:=factor(approach,levels = c("Amplicon","Metagenomic","Capture"), ordered = T)]
data[,platform:=factor(platform,levels = c("Illumina","Nanopore"), ordered = T)]

p <- plot_ridges(data)

p <- p + facet_grid(virus ~ approach + platform) +
	theme(strip.text = element_text(size=12),
				axis.text = element_text(size=10),
				axis.title = element_text(size=12),
				legend.text = element_text(size=10))

#uncomment and comment below to save results before cleaning
# ggsave("coverage_overview.pdf",p,width = 18, height=10)
# ggsave("coverage_overview.png",p,width = 18, height=10)

ggsave("coverage_overview_cleaned.pdf",p,width = 18, height=10)
ggsave("coverage_overview_cleaned.png",p,width = 18, height=10)

#Filter out positions outside of amplicon region
filtereddata <- data[(approach=="Amplicon"&
						 	((virus=="ZIKV" & position>=13 & position<=10653)|
						 	(virus=="YFV" & position>=123 & position<=10338)|
						 	(virus=="WNV" & position>=9 & position<=11011)|
						 	(virus=="USUV" & position<=11955)))|
						 		approach!="Amplicon"]

#Divide the positions with no 'NA' by the genome size to obtain the percentage of coverage above the threshold
plotdata <- filtereddata[,.(perc_coverage=round(100*sum(!is.na(coverage_comparison))/max(position),1)), .(approach, platform, virus, CT)]

#Some relabeling
plotdata[,CT:=factor(CT,levels = c("CT33","CT29","CT25"), ordered = T)]
plotdata[,approach:=factor(approach,levels = c("Amplicon","Metagenomic","Capture"), ordered = T)]
plotdata[,platform:=factor(platform,levels = c("Illumina","Nanopore"), ordered = T)]

#Generate the plot
p <- ggplot(plotdata) + 
	geom_bar(aes(x=CT, y=perc_coverage), position="dodge", stat="identity", fill = "#0ea1e6") +
	geom_text(aes(x=CT, y=perc_coverage/3, label=paste0(perc_coverage,"%")), position = position_dodge(0.9), hjust = 0) +
	coord_flip() +
	facet_grid(virus~approach+platform, scales = "free") +
	theme_minimal() +
	labs(x="CT value", y="Reliable genome coverage") +
	theme(panel.border = element_rect(color = "black", fill = NA))

#uncomment and comment below to save results before cleaning
# ggsave("perc_coverage_overview.pdf",p,width = 18, height=4)
# ggsave("perc_coverage_overview.png",p,width = 18, height=4)

ggsave("perc_coverage_overview_cleaned.pdf",p,width = 18, height=4)
ggsave("perc_coverage_overview_cleaned.png",p,width = 18, height=4)

####Read count figure
read_readstats <- function(x){
	data <- fread(x)
	data[,sample_id := str_remove(x, ".*/") %>% str_remove("\\..*")]
}

#Fill in path containing the readstats files resulting from running the Snakemake workflow
filelist <- list.files(path="readstats/", pattern = "*.tsv", full.names = T)

filelist <- filelist[grep("dedup",filelist, invert = T)]

#Read in readstats files and combine, parse the file names to extract which platform etc. 
read_stats <- rbindlist(lapply(filelist, read_readstats))
read_stats[,num_seqs := str_remove_all(num_seqs,",") %>% as.numeric()]
read_stats <- read_stats[,.("num_reads"=sum(num_seqs)),sample_id]
read_stats[,sample_id:=str_remove(sample_id, "_ZIKV|_USUV|_WNV|_YFV")]
read_stats[,sample_id:=str_replace(sample_id, "_mapped_dedup","_dedup")]
read_stats[,sample_id:=str_replace(sample_id, "_mapped_cleaned","_cleaned")]
read_stats[,c("barcode","approach","platform","stage"):=tstrsplit(sample_id, "_", fixed=T)]

#Add CT values
read_stats <- merge.data.table(read_stats, ct_values, by = "barcode")

#Some relabeling
read_stats[approach=="Meta", approach:="Metagenomic"]
read_stats[approach=="Virocap", approach:="Capture"]

read_stats[,CT:=factor(CT,levels = c("CT33","CT29","CT25"), ordered = T)]
read_stats[,approach:=factor(approach,levels = c("Amplicon","Metagenomic","Capture"), ordered = T)]
read_stats[,platform:=factor(platform,levels = c("Illumina","Nanopore"), ordered = T)]

read_stats[,percentage_reads:=round(.SD[,num_reads]/.SD[stage=="raw",num_reads]*100, digits = 1),.(barcode,approach,platform)]

#Generate the plot
p <- ggplot(read_stats) + 
	geom_bar(aes(x=CT, y=num_reads, fill=stage, group=stage), position="dodge", stat="identity") +
	geom_text(aes(x=CT, y=num_reads/3, label=paste(sapply(num_reads, f2si, rounding = T, digits = 3), paste0("(",percentage_reads,"%)")), group=stage), size = 3, position = position_dodge(0.9), hjust = 0) +
	coord_flip() +
	scale_fill_manual(values = c("#b50024","#e30057","#ff672b","#ffd52b"), name="Analysis process", labels=c("cleaned mapped","Mapped reads", "QC reads", "Raw reads"))+
	facet_grid(virus~approach+platform, scales = "free") +
	theme_minimal() +
	labs(x="CT value", y="Read count") +
	theme(panel.border = element_rect(color = "black", fill = NA),
				strip.text = element_text(size=12),
				axis.text = element_text(size=10),
				axis.title = element_text(size=12),
				legend.title = element_blank(),
				legend.text = element_text(size=10),
				legend.position="top")

ggsave("reads_overview.pdf",p,width = 18, height=10)
ggsave("reads_overview.png",p,width = 18, height=10)

read_vcffile <- function(x){
  data <- fread(x)
  
  sample_id <- str_remove(x, ".*/") %>% str_remove("\\..*")
  sample_info <- str_split(sample_id, "_", simplify=T)
  
  data[,c("sample_id","barcode", "approach", "platform") := 
         list(sample_id, sample_info[1], sample_info[2], sample_info[3])]
  return(data)
}

#Fill in path containing the depth files resulting from running the Snakemake workflow
filelist <- list.files(path="vcf", pattern = "*.tsv", full.names = T)

#read all files and combine
vcfdata_all <- rbindlist(lapply(filelist, read_vcffile))

vcfdata_all[approach=="Meta", approach:="Metagenomic"]
vcfdata_all[approach=="Virocap", approach:="Capture"]

#uncomment and comment below to get results before cleaning
# vcfdata <- vcfdata_all[grep("dedup|cleaned",sample_id, invert=T)]

vcfdata <- vcfdata_all[grep("cleaned",sample_id)]

vcfdata <- merge.data.table(data, vcfdata, all.x=T, by.x = c("sample_id","barcode","approach","platform","position"), by.y = c("sample_id","barcode","approach","platform","POS"))

major_variants <- vcfdata[MAJOR==1&!is.na(coverage_comparison),.N,.(virus,position,REF,ALT)][,.(virus,position,REF,ALT)]

all_combinations <- merge.data.table(major_variants, unique(vcfdata[,.(virus,platform,approach,CT)]), by=c("virus"), allow.cartesian = T)

plotdata <- merge.data.table(all_combinations, data[,.(virus,position,platform,approach,CT,coverage_comparison)], by=c("virus","position","platform","approach","CT"), all.x=T)
plotdata <- merge.data.table(plotdata, vcfdata[MAJOR==1,.(virus,position,REF,ALT,platform,approach,CT,coverage_comparison,AF)], by=c("virus","position","REF","ALT","platform","approach","CT","coverage_comparison"), all.x=T)

plotdata[is.na(coverage_comparison), AF:=0]
plotdata <- plotdata[!is.na(AF)]
plotdata[AF==0, AF:=NA]

plotdata[,CT:=factor(CT,levels = c("CT33","CT29","CT25"), ordered = T)]
plotdata[,approach:=factor(approach,levels = c("Amplicon","Metagenomic","Capture"), ordered = T)]
plotdata[,platform:=factor(platform,levels = c("Illumina","Nanopore"), ordered = T)]
plotdata[,mutation_annot:=paste0(position," ",REF,">",ALT)]
plotdata[str_length(REF)>2,mutation_annot:=paste0("del \u0394", position+1, "-", position+str_length(REF)-2)]
plotdata[str_length(REF)==2,mutation_annot:=paste0("del \u0394", position+1)]
plotdata[,mutation_annot:=factor(mutation_annot,levels = plotdata[,.N,.(position,mutation_annot)][order(position)][,mutation_annot], ordered = T)]

p <- ggplot(plotdata[!mutation_annot%in%plotdata[,.N,.(mutation_annot)][N==15,mutation_annot]]) + geom_tile(aes(y=paste(platform,approach),x=mutation_annot, fill=AF)) +
  scale_fill_gradientn(colors=c("blue","#ebebeb","red"), na.value = "gray", limits=c(0,1)) +
	theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill=NA),
        text = element_text(size=16),
        axis.text.x = element_text(angle = 60, hjust=1, size=8)) +
  facet_grid(CT~virus, scales = "free", space = "free") +
	labs(x="Variant", y="Platform / Approach", fill="Variant fraction")

p

#uncomment and comment below to save results before cleaning
# ggsave("variant_overview.pdf",p,width = 18, height=4)
# ggsave("variant_overview.png",p,width = 18, height=4)

ggsave("variant_overview_cleaned.pdf", p, width = 14, height=4)
ggsave("variant_overview_cleaned.png", p, width = 14, height=4)

p <- ggplot(plotdata[!mutation_annot%in%plotdata[,.N,.(mutation_annot)][N!=15,mutation_annot]]) + geom_tile(aes(y=paste(platform,approach),x=mutation_annot, fill=AF)) +
  scale_fill_gradientn(colors=c("blue","#ebebeb","red"), na.value = "gray", limits=c(0,1)) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill=NA),
        text = element_text(size=16),
        axis.text.x = element_text(angle = 60, hjust=1, size=8)) +
  facet_grid(CT~virus, scales = "free", space = "free") +
  labs(x="Variant", y="Platform / Approach", fill="Variant fraction")

p

ggsave("true_variants.pdf",p,width = 14, height=4)
ggsave("true_variants.png",p,width = 14, height=4)

