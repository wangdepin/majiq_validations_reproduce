# Validation Tools

This package implements different scripts that are used for validating/testing MAJIQ. Some of those evaluations include the evaluation of competing methods under the same conditions.

This package has two parts
+ Root directory where different scripts are located. They implement the main pipeline for generating the evaluation. Each script is one figure or result.
+ Tools directory. This directory include one file per tool that can be tested. Currently there is full evaluation for rMATS, SUPPA2, leafcutter and all MAJIQ differential splicing methods. There is partial implementation for Whippet. The `skeleton.py` file includes an explanation of a basic tool py file structure. 

## Importing Tools in script

The structure is created to automatically add a tool if a py file exist in the tools folder. For that reason there are two snippets.
+ importing the objects from __init__.py:
```
from tools import operator, all_stats
```
+ importing dinamically the tools python files if they are required by the input.
  ```
  for xx in args.tool:
        if xx[0] not in all_stats:
            print ('ERROR tool %s is not available' % xx[0])
            exit(-1)

        module_ = __import__('tools.' + xx[0].lower(), fromlist=xx[0].title())
        class_ = getattr(module_, xx[0].title())
        operator[xx[0]] = class_()
        tools[(xx[0], xx[1])] = xx[2:]
  ```
  This case tool is a parsed parameter that includes type of tool, id, and input files ( see run_iir.py file for example)

## Existing scripts
This is the list of existing scripts, We will organize them by function, so putting together the ones that are parts of the same process.

### Simulation validation
This pipeline implements a validation of different tools calculating FDR/FPR metrics. This is done using simulated data done by Greg Grant. How this simulation is done, check Norton et al 2017, Bioinformatics.

Let´s assume we have RNASeq simulation and the TPR or Transcript counts matching that RNA simulation. Now we face a major problem, how to transform transcript counts from the ground truth to something is comparable to all the tools. Each tool represent splicing events in different ways so we need to be able to transform the information to all those _languages_.

For that reason we do this in two steps.
1. We parse the transcript counts to junction/exon inclusion based on the counts. Per junction we sum all the transcript counts that include that junction, normalized by the transcript length. Creating a ground truth dictionary. This is run by subcommand `run_validation.py gen_truth`
2. We parse each tool output file and given their definition, we check which would be the proper dpsi quantification. Note that this uses the tools event definition to check their accuracy. We do not have a realistic unified ground truth for everyone since each tool defines things differently. Majiq defines LSVs, rMATS classical events, leafcutter intronic cluster. Some of them detect complexity, others do not. This is executed using subcommand `run_validation.py validate`

*Note: We are working on a variation of this based on differentially spliced genes, in order to be able to compare the tools in the same space*

#### Tools generated for the simulation validation
+ `add_transcipt_len.py`: Validation step one requires the transcript lenght to be able to normilize, we include that length in the transcript gene translation file using this script
+ `merge_column_files.py`: BEERS simulation includes the ground truth separated by sample in different files. This script merges all of them in a multicolumn file, that will be used in step 1.


### Reproducibility evaluation
Running reproducibility analysis is done in two scripts
1. `run_repro.py`: calculates reproducibility and generates a log file and a pickle files
2. `plot_rr.py`: Creates plot using the pickle files created in the previous step.

### IIR evaluation

### RTPCR evaluations

### Other tools
+ `extract_ground_truth.py`: Using a voila file as input creates a new voila file with PSI generated by the ground truth. This ground truth values are extracted from files creating during the simulation validation pipelines. 
+ `dpsi_to_settings.py`: This script uses a list of majiq deltapsi runs to create a tool independent list of runs that can be used _a posteriori_. We will call that file _run\_settings_ file.
+ `gen_tool_runs.py`: Using a _run\_settings_ file as input it generates the pipeline to run an specific tool. Currently this is configured to create pipelines for MAJIQ, rMATS, SUPPA, Whippet and leafcutter.
+ `split_file.py`: This file is used by suppa pipeline to split temporary files, as it is required on suppa usage manual.
+ `extract_introns.py`: This script reads a gff3 file and extracts all the annotated introns including gene id transcript Id and intron coordinates.
+ `plot_l1_distances.py`: This script generates l1 distances plot between a group replica EPSI and the median of EPSI for the whole group.
+ `
