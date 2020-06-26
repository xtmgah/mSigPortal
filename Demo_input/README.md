############### Usage Example ###############<br><br>
#001 VCF<br>
python mSigPortal_Profiler_Extraction.py -f vcf -i Demo_input/demo_input_single.vcf.tar.gz -p Project -o Test_Output-sigle-VCF -g GRCh37 -t WGS    <br>
python mSigPortal_Profiler_Extraction.py -f vcf -i Demo_input/demo_input_single.vcf.gz -p Project -o Test_Output-sigle-VCF -g GRCh37 -t WGS -s True

#002 CSV<br>
python mSigPortal_Profiler_Extraction.py -f csv -i Demo_input/demo_input_multi.csv -p Project -o Test_Output-sigle-VCF -g GRCh37 -t WGS -F alt_allele_in_normal
python mSigPortal_Profiler_Extraction.py -f tsv -i Demo_input/demo_input_multi.tsv -p Project -o Test_Output_CSV_Split -g GRCh37 -t WGS -F germline_risk
python mSigPortal_Profiler_Extraction.py -f tsv -i Demo_input/demo_input_multi.tsv -p Project -o Test_Output_CSV_Split -g GRCh37 -t WGS -s True

#003 TSV<br>
python mSigPortal_Profiler_Extraction.py -f tsv -i Demo_input/demo_input_multi.tsv -p Project -o Test_Output-sigle-VCF -g GRCh37 -t WGS -F alt_allele_in_normal <br>
python mSigPortal_Profiler_Extraction.py -f tsv -F germline_risk -i Demo_input/demo_input_multi.tsv -p Project -o Test_Output -g GRCh37 -t WGS

#004 catalog_tsv<br>
python mSigPortal_Profiler_Extraction.py -f catalog_tsv -i Demo_input/demo_input_catalog.tsv -p Project -o Test_Output_catalog_tsv -g GRCh37 -t WGS

#005 catalog_csv<br>
python mSigPortal_Profiler_Extraction.py -f catalog_csv -i Demo_input/demo_input_catalog.csv -p Project -o Test_Output_catalog_csv -g GRCh37 -t WGS
