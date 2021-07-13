# GEO_syn_parser
RNAseq data parser from GEO

**Clone the repository and do as follows:**
1. install a miniconda3, then create a environment by doing:
```
conda create -n py3_geo python=3
conda activate py3_geo
```
2. install requirements
```
pip install -r requirements.txt
```
3. test your environment deploy, it is ok, if no error report when you run:
```
python env.py
```
4. run parser:
```
python scrna_parser_runner.py -h
```
Examples can be found in run.sh.
An example bash script that was used to run the parser for TISMO on the DFCI cluster can be found in run_parser2.sh file