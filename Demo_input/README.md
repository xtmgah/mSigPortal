############### Usage Example ###############<br><br>
#001 VCF<br>
python mSigPortal_Profiler_Extraction.py -f vcf -i Demo_input/demo_input_single.vcf.tar.gz -p Project -o Test_Output-sigle-VCF -g GRCh37 -t WGS    <br>
python mSigPortal_Profiler_Extraction.py -f vcf -i Demo_input/demo_input_single.vcf.gz -p Project -o Test_Output-sigle-VCF -g GRCh37 -t WGS -s True

#002 CSV<br>
python mSigPortal_Profiler_Extraction.py -f csv -i Demo_input/demo_input_multi.csv -p Project -o Test_Output-sigle-VCF -g GRCh37 -t WGS -F alt_allele_in_normal

#003 TSV<br>
python mSigPortal_Profiler_Extraction.py -f tsv -i Demo_input/demo_input_multi.tsv -p Project -o Test_Output-sigle-VCF -g GRCh37 -t WGS -F alt_allele_in_normal




