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


5-1 Replace the script [SigProfilerMatrixGeneratorFunc.py] in https://github.com/AlexandrovLab/SigProfilerMatrixGenerator/tree/master/SigProfilerMatrixGenerator/scripts/SigProfilerMatrixGeneratorFunc.py, with the updated one from our github portal.<br>
5-2 Replace the script [sigProfilerPlotting.py] in https://github.com/AlexandrovLab/SigProfilerPlotting/tree/master/sigProfilerPlotting/sigProfilerPlotting.py , with the updated one from our github portal<br>
 

6. Run Script:
python mSigPortal_Profiler_Extraction.py -f vcf -i Demo_input/demo_input_multi.vcf -p Project -o Test_Output -g GRCh37 -t WGS -c True

 
