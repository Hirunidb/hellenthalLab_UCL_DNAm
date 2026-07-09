# filtering out SNPs with Minimum Allele Frequency (MAF) >= 0.01
# SNP database downloaded from https://ftp.ensembl.org/pub/release-115/variation/gvf/homo_sapiens/1000GENOMES-phase_3.gvf.gz

cd /path/to/folder/

zcat 1000GENOMES-phase_3.gvf.gz   | awk -F'\t' '
    BEGIN { OFS="\t" }
    $0 !~ /^#/ && $3=="SNV" {
      maf_ok = 0
      n = split($9, a, ";")
      for (i = 1; i <= n; i++) {
        if (a[i] ~ /^(AFR|AMR|EAS|EUR|SAS)=/) {
          split(a[i], b, "=")
          split(b[2], vals, ",")
          for (j in vals) {
            if (vals[j] + 0 >= 0.01) maf_ok = 1
          }
        }
      }
      if (maf_ok) print "chr" $1 "_" $4
    }' > 0.1_snps_maf_0.01_chr_pos.txt