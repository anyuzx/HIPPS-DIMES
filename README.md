# HIPPS-DIMES
Maximum Entropy Based HI-C/Distance Map - Polymer Physics - Structures Method

# Description

This python script can be used to generate ensemble of genome structures from Hi-C contact map or mean spatial distance map. The method is based on Maximum Entropy principle and the relation between the contact probability and the mean spatial distance from polymer theory. The application of the method can be found in this work https://www.biorxiv.org/content/10.1101/2020.05.21.109421v1.

# Documentation

## Install

```bash
pip install --editable .
```

This command will install the required packages, and install the script as a python module. Once installed, you can call `HippsDimes` directly in the terminal to run the script

## How to use

To display help information, use

```bash
HippsDimes --help
```

### Output files

This script will generate several files:

- A text file for the final simulated mean distance map
- A text file for the connectivity matrix
- A `.xyz` formatted file for the ensemble of genome structures generated (can be turned off)
- A csv formatted file for cost versus iteration data (can be turned off)

### Examples



### Explanantion of the arguments and options

#### Argument

- `INPUT`: File path for the input file. The input file can be a Hi-C contact map or a mean spatial distance map as measured in Multiplexed FISH experiment.
- `OUTPUT_PREFIX`: Prefix for outputfiles. For instance, if one specify it to be `TEST`, then all the output files will start with `TEST_`.

#### Options

- `-e` or `--ensemble`: Number of individual conformations to be generated. This script will generate an ensemble of structures consistent with the input Hi-C contact map or the mean spatial distance map. Each individual conformations are different from each other. You can specify how many such individual conformations you want to generate.
- `-a` or `--alpha`: Value of the contact map to distance map conversion exponent. If the input file is Hi-C contact map, the method first convert the contact map to a mean spatial distance map. The equation of the conversion is d_{ij} ~ c_{ij}^{1/\alpha}. The default value of \alpha is 4.0, estimated in this work 10.1126/science.aaf8084
- `-s` or `--selection`: Specify chromosome or region. This option only works when the input file has [`cooler`](https://github.com/open2c/cooler) format. The value of this option is passed to the `cooler.Cooler.matrix().fetch()` method. For details, please refer their [documentation](https://cooler.readthedocs.io/en/latest/concepts.html#matrix-selector)
- `-i` or `--iteration`: The method relies on iterative scaling to find the optimal parameters. This option specifies the number of iterations. Generally, the more iterations the model runs, the better results are. However, the convergence of the model slow down when iteration increases. For larger size of contact map and the mean distance map, the number of iterations needed to good convergence is larger.
- `-r` or `--learning-rate`: The learning rate for iterative scaling. Higher learning rate achieves faster convergence. However, the model can crash if learning rate is too large. The default value is 10. One should play around this option to see what works best.
- `--input-type`: The type of the input file. To use the script, the type must be specified. The method can work on both the contact map (`cmap`) or distance map (`dmap`).
- `--input-format`: The format of the input file. If the type of input file is Hi-C contact map, then the script support `cooler` format Hi-C contact map file or a pure text based file. In the text based file, each line corresponds to the row of the contact map. If the type of input file is mean distance map, then the script only support the text based file in which each line represents the row of the mean distance map.
- `--log`: A log file will be written if this option is specified. The log file contains the data of cost versus iteration
- `--no-xyzs`: Turn off writing x,y,z coordinates of genome structures to files.
