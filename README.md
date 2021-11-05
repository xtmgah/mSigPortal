# mSigPortal
Integrative mutational signature portal (MsigPortal) for cancer genomic study

Here are some basic steps:

1. Download the script and two Packages from Github:

https://github.com/xtmgah/SigProfilerPlotting<br>
https://github.com/xtmgah/SigProfilerMatrixGenerator<br>
https://github.com/xtmgah/mSigPortal

2-1 Download the updated script: SigProfilerMatrixGeneratorFunc.py
2-2 Download the updated script: sigProfilerPlotting.py<br>

3. Install the two Packages locally based on the following order:

pip install -e Path/SigProfilerPlotting-master/ <br>
pip install -e Path/SigProfilerMatrixGenerator-master/


4. Download Reference Genome, e.g.,'GRCh37’  

$ python <br>
from SigProfilerMatrixGenerator import install as genInstall <br>
genInstall.install('GRCh37', rsync=False, bash=True)<br>

 

5. Run Script:
python mSigPortal_Profiler_Extraction.py -f vcf -i Demo_input/demo_input_multi.vcf -p Project -o Test_Output -g GRCh37 -t WGS -c True

 
