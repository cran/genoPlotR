# download fasta & gbk
FTP=ftp://ftp.ncbi.nih.gov/genomes/Bacteria/
for EXT in ptt rnt fna; do
    wget $FTP/Bartonella_bacilliformis_KC583/NC_008783.$EXT -O BB.$EXT
    wget $FTP/Bartonella_grahamii_as4aup/NC_012846.$EXT -O BG.$EXT
    wget $FTP/Bartonella_henselae_Houston-1/NC_005956.$EXT -O BH.$EXT
    wget $FTP/Bartonella_quintana_Toulouse/NC_005955.$EXT -O BQ.$EXT
done

# make databases; run blasts with BLAST+ 2.2.22
for ORG in BB BG BH BQ; do
    makeblastdb -in $ORG.fna -dbtype nucl -parse_seqids
done
blastn -query BB.fna -db BG.fna -evalue 1e-3 -outfmt 6 > BB_vs_BG.blastn.tab
blastn -query BG.fna -db BH.fna -evalue 1e-3 -outfmt 6 > BG_vs_BH.blastn.tab
blastn -query BH.fna -db BQ.fna -evalue 1e-3 -outfmt 6 > BH_vs_BQ.blastn.tab
