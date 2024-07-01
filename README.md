# epigenomics
Analysis code for epigenomics data

## Data processing

### 1.QC
```bash
fastp -w 8 -b 50 -B 50 -i sample1/sample1_r1.fq.gz -I sample1/sample1_r2.fq.gz -o QC/sample1/sample1_R1.clean.fq.gz -O QC/sample1/sample1_R2.clean.fq.gz --html QC/sample1/sample1.html --json QC/sample1/sample1.json
```

### 2.mapping
```bash
##step01
chromap --preset chip/atac -t 8 -x genome.fa.fai -r genome.fa -1  sample1/sample1_R1.clean.fq.gz -2  sample1/sample1_R2.clean.fq.gz --TagAlign -o mapping/sample1/sample1.tagAlign > mapping/sample1/sample1.stat 2>&1
        
##step02
awk '{{ORS=NR%2?FS:RS}}1'  mapping/sample1/sample1.tagAlign | awk  '{{OFS="\t"}}{{print $1, $2, $3, $7, $8, $9, $4, $5"@"$11, $6, $12}}' > mapping/sample1/sample1.bedpe
        
##step03
pairToBed -type neither -a mapping/sample1/sample1.bedpe -b genome.blacklist > mapping/sample1/sample1.blacklist.bedpe
        
##step04
pairToBed -type both -a mapping/sample1/sample1.blacklist.bedpe -b genome.whitelist |awk '{{OFS="\t"}}{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}}' | sort -k1,1 -k2,2n | uniq > mapping/sample1/sample1.whitelist.bedpe
        
#step05
get_random_lines_from_file(f'mapping/sample1/sample1.whitelist.bedpe', f'mapping/sample1/sample1.30M.tagAlign', '30M')
bedtotag = bedpetotagalign(f'mapping/sample1/sample1.30M.tagAlign', f'mapping/sample1/sample1.new.tagAlign')
```

### 3.spp
```bash
python3 encode_task_xcor.py --nth 4 --mito-chr-name MT/null --out-dir callpeak mapping/sample1/sample1.new.tagAlign

extsize = 150
files = glob.glob(f"callpeak/*.cc.fraglen.txt")
extsize=f1.readline().strip()


```

### 4.callpeak
```bash
###chip
#target_type=="narrow":
macs3 callpeak -n sample1 -t callpeak/sample1.new.tagAlign.no_chrM.tagAlign.gz -c chip.input.file -g genome.fa -q 0.01 --nomodel --shift 0 --extsize extsize --keep-dup all -B --SPMR --outdir callpeak/sample1/target_type

#target_type=="broad":
macs3 callpeak -n sample1 -t callpeak/sample1.new.tagAlign.no_chrM.tagAlign.gz -c chip.input.file -g genome.fa  --broad --broad-cutoff 0.1 --nomodel --shift 0  --keep-dup all -B --SPMR --outdir callpeak/sample1/target_type

#target_type=="unknown":    
macs3 callpeak -n sample1  -t callpeak/sample1.new.tagAlign.no_chrM.tagAlign.gz -c chip.input.file -g genome.fa -q 0.01 --nomodel --shift 0 --extsize extsize --keep-dup all -B --SPMR --outdir callpeak/sample1/target_type

macs3 callpeak -n sample1 -t callpeak/sample1.new.tagAlign.no_chrM.tagAlign.gz -c chip.input.file -g genome.fa  --broad --broad-cutoff 0.1 --nomodel --shift 0  --keep-dup all -B --SPMR --outdir callpeak/sample1/target_type


###atac
macs3 callpeak -n sample1 -t callpeak/sample1.new.tagAlign.no_chrM.tagAlign.gz -g genome.fa -q 0.01 --nomodel --shift 75 --extsize extsize --keep-dup all -B --SPMR --outdir callpeak/sample1/target_type
              
```

### 5.for -log10(p value) signal track
```bash
reads = os.popen(f"gzip -dc callpeak/sample1.new.tagAlign.no_chrM.tagAlign.gz | wc -l").read()
sval = str(int(reads)/1000000)
macs3 bdgcmp -t callpeak/sample1/*/sample1_treat_pileup.bdg -c callpeak/sample1/*/sample1_control_lambda.bdg --outdir callpeak -o sample1_ppois.bdg -m ppois -S sval

###fixed bdg
fixbdg = fix_bdg(f'callpeak/sample1/', genomesize)
```

### 6.peak matrix
```bash
#获取MACS peak list(将所有样本的peak合并到一个文件)
python3 get_sample_info.py -a scientific_name -s scientific_name
python3 get_peak_list.py -i callpeak -s sample_info_peak.txt -o peak_list

#将基因组划分为bin
python3 create_bins.py -a scientific_name -f genome.fa.fai -o scientific_name.bin

##提取每个样本中的bw信号(加到callpeak更合适？)
python3  extract_signal.py -p callpeak/sample -b scientific_name.bin

##这里需要把ChIP-Seq和ATAC-Seq的样本信息合并，然后在生成矩阵个
#peak matrix
python3 merge_peak.py -i callpeak -s sample_info_peak.txt -b scientific_name.bin -o peak_matrix
#signal matrix
python3 merge_signal.py -i callpeak -s sample_info_peak.txt -o signal_matrix
```

### 7.格式化
```bash
#处理peak list
##读取基因信息
python3 get_gene_info.py -i 01gene_details.txt -o genes.bed
bedtools intersect -loj -a peak_list -b genes.bed -wa -wb >peak_list_anno.bed
python3 format_peak_list.py -i peak_list_anno.bed -o peak_list_fmt.bed

#peak_matrix转化为peak_matrix_list
python3 format_peak_matrix.py -s sample_info_peak.txt -i peak_matrix -o peak_matrix.fmt

#signal_matrix转化为signal_matrix_list
python3 format_signal_matrix.py -s sample_info_peak.txt -i signal_matrix -o signal_matrix.fmt
```