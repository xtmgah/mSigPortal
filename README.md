# mSigPortal
Integrative mutational signature portal (MsigPortal) for cancer genomic study

Here are some basic steps:

1. Download the script and three Packages from Github: <br>

https://github.com/xtmgah/SigProfilerPlotting<br>
https://github.com/xtmgah/SigProfilerMatrixGenerator<br>
https://github.com/xtmgah/SigProfilerClusters<br>

https://github.com/xtmgah/mSigPortal


2. Install them locally based on the following order: <br>

pip install -e Path/SigProfilerClusters-master/ <br>
pip install -e Path/SigProfilerPlotting-master/ <br>
pip install -e Path/SigProfilerMatrixGenerator-master/ <br>


3. Download Reference Genome, e.g.,'GRCh37’  

$ python <br>
from SigProfilerMatrixGenerator import install as genInstall <br>
genInstall.install('GRCh37', rsync=False, bash=True)<br>
genInstall.install('GRCh38', rsync=False, bash=True)<br>

 

4. Run Script:
python3 mSigPortal_Profiler_Extraction.py -f vcf -i Demo_input_Test/demo_input_multi.vcf.gz -p Project -o test-15 -g GRCh37 -t WGS -C True

 
