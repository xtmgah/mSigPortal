# mSigPortal
Integrative mutational signature portal (MsigPortal) for cancer genomic study

Here are some basic steps:

1. Download the script and two Packages from Github:
https://github.com/xtmgah/SigProfilerPlotting<br>
https://github.com/xtmgah/SigProfilerMatrixGenerator<br>
https://github.com/xtmgah/mSigPortal


2. Install the two Packages locally based on the following order:

pip install -e Path/SigProfilerPlotting-master/
pip install -e Path/SigProfilerMatrixGenerator-master/


3. Download Reference Genome, e.g.,'GRCh37’  

$ python 
>> from SigProfilerMatrixGenerator import install as genInstall 
>> genInstall.install('GRCh37', rsync=False, bash=True)


4. Run the Script mSigPortal_Profiler_Extraction.py

python mSigPortal_Profiler_Extraction.py -f vcf -F PASS@alt_allele_in_normal@- -i Demo_input/demo_input_multi.vcf -p Project -o Test_Output -g GRCh37 -t WGS

 
