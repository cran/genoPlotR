# download gbk files for 4 Bartonella genomes
FTP=ftp://ftp.ncbi.nih.gov/genomes/Bacteria/
wget $FTP/Bartonella_bacilliformis_KC583/NC_008783.gbk -O BB.gbk
wget $FTP/Bartonella_grahamii_as4aup/NC_012846.gbk -O BG.gbk
wget $FTP/Bartonella_henselae_Houston-1/NC_005956.gbk -O BH.gbk
wget $FTP/Bartonella_quintana_Toulouse/NC_005955.gbk -O BQ.gbk

# run mauve: this supposes that you have progressiveMauve available in 
# your path. Mauve version used was 2.3.1
progressiveMauve --output=barto.xmfa --output-guide-tree=barto.guide_tree --backbone-output=barto.backbone --hmm-identity=0.8 --island-gap-size=10000 BB.gbk BG.gbk BH.gbk BQ.gbk
