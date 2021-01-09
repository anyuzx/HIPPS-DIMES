
# HIPPS-DIMES
Maximum Entropy Based HI-C/Distance Map - Polymer Physics - Structures Method

# Description

This python script can be used to generate ensemble of genome structures from Hi-C contact map or mean spatial distance map. The method is based on Maximum Entropy principle and the relation between the contact probability and the mean spatial distance from polymer theory. The application of the method can be found in this work https://www.biorxiv.org/content/10.1101/2020.05.21.109421v1.

In general, this script accepts the input file of a Hi-C contact or a mean spatial distance map (can be measured in Multiplex FISH experiment), and generate an ensemble of individual conformations consistent with the inputted contact map or distance map. The output conformations are in `.xyz` format, which users can use to calculate quantities of interest and can be viewed using `VMD` or other compatible softwares.

# Documentation

## Requirements

- Basic understanding of `bash`, `Python` and how to use terminal
- Python 3.8+
- Python Package Installer [`pip`](https://pip.pypa.io/en/stable/). Follow the instruction https://pip.pypa.io/en/stable/installing/ to install it if it is not installed already. Make sure the version of `pip3` is up-to-date.
- Linux or Mac OS system

## Install

First, open a terminal window and download this repository by put the following command in the terminal and hit Enter,

```bash
git clone https://github.com/anyuzx/HIPPS-DIMES
```

Next, install required packages using the command below,

```bash
cd HIPPS-DIMES
pip3 install --editable .
```

This command will install the required packages, and install the script as a python module. Note that you need to install `pip3` if it is not installed already (Follow the instruction on the official document here https://pip.pypa.io/en/stable/installing/). Once installed, you can call `HippsDimes` directly in the terminal to run the script. The packages installed are:

* `Click`
* `Numpy`
* `Scipy`
* `Pandas`
* `Tqdm`
* `Cooler`

## How to use

To get started, please go through the jupyter notebook **`walkthrough.ipynb`** in this repository.

In addition, it is helpful to view the help information for each arguments and options. To display help information, use

```bash
HippsDimes --help
```

### Input files

This script accept input files in two formats. If the input file is a Hi-C contact map, it can be in `cooler` format or pure text format. If the input file is a mean spatial distance map, the script only accepts a pure text formatted file. The text format for a matrix is the following: each row of the file corresponds to the row of the matrix. Values are space-separated. An example of such file is,

``` text
1  2  3
2  1  2
3  2  1
```

### Output files

This script will generate several files:

- A text file for the final simulated mean distance map
- A text file for the final simulated contact map
- A text file for the connectivity matrix
- A `.xyz` formatted file for the ensemble of genome structures generated (can be turned off)
- A csv formatted file for cost versus iteration data (can be turned off)

### Examples

#### Example 1

First, download a cooler format Hi-C contact map from [here](https://drive.google.com/file/d/1eIxGv1JbIrEAVoUSQK_O_ebIjWo6toTJ/view?usp=sharing) (**The file size is about 116 Mb**). This Hi-C contact map is for Chicken cell mitotic chromosome, originally retrieved from [GEO repository](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102740). Rename it to `hic_example.cool`. Then execute the following command,

```bash
HippsDimes hic_example.cool test --input-type cmap --input-format cooler -s chr7:10M-15M -i 10 -e 10
```

This command tells the script to load the Hi-C contact map `hic_example.cool` and perform the iterative scaling algorithm. The argument `test` instructs the files names of output files start with `test_`. Option `--input-type cmap` specifies that the input file is a contact map. Option `--input-format cooler` specifies that the input file is a `cooler` file. Option `-s chr7:10M-15M` specifies that the algorithm is performed on the region 10 Mbps - 15 Mbps on Chromosome 7. Note that these three options are required and cannot be neglected. **Some option arguments are optional, some are required. Please refer to the section below and use `HippsDimes --help` for details**

When the program finishes, the script will generate several output files: `test.xyz`, `test_connectivity_matrix.txt`, and `test_dmap_final.txt`. `test.xyz` contains 10 sets of individual conformations of x, y, z coordinates and can be viewed using `VMD` or other compatible visualization softwares.

#### Example 2

In this example, we use Hi-C contact map for HeLa cell line Chromosome 14 at time point of 12 hours after the release from prometaphase. For the purpase of demonstration, you can download the Hi-C `.cool` file from [here](https://drive.google.com/file/d/1j-zfDUP6LOZGCxz9uA3LaMI372ct1cU_/view?usp=sharing) which is origannly retreived from [GEO repository](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102740) under accession number GSE102740. **Before you download, note that the file has size of about 655 Mb**. Once downloaded, execute the following command,

```bash
HippsDimes GSM3909682_TB-HiC-Dpn-R2-T12_hg19.1000.multires.cool::6 test --input-type cmap --input-format cooler -s chr14:20M-107M -i 10000 -e 10
```

Similar to the first example, this command tells the script to load the Hi-C cooler file `GSM3909682_TB-HiC-Dpn-R2-T12_hg19.1000.multires.cool` and its group 6 (data for several different resolution is stored in different groups) and perform the HIPPS/DIMES algorithm. In this example, we change the number of iterations to be 10000 by using the option `-i 10000`. On a AMD Ryzen 5 3600 CPU machine, it takes about 3-4 mins to finish the program. Once it is finished, several ouput files are generated.

### Explanantion of the arguments and options

#### Argument

- `INPUT`: File path for the input file. The input file can be a Hi-C contact map or a mean spatial distance map as measured in Multiplexed FISH experiment.
- `OUTPUT_PREFIX`: Prefix for outputfiles. For instance, if one specify it to be `TEST`, then all the output files will start with `TEST_`.

#### Options

- `-e` or `--ensemble`: Number of individual conformations to be generated. This script will generate an ensemble of structures consistent with the input Hi-C contact map or the mean spatial distance map. Each individual conformations are different from each other. You can specify how many such individual conformations you want to generate. If not specified, its value would be 1000.
- `-a` or `--alpha`: Value of the contact map to distance map conversion exponent. If the input file is Hi-C contact map, the method first convert the contact map to a mean spatial distance map. The equation of the conversion is d_{ij} ~ c_{ij}^{1/\alpha}. The default value of \alpha is 4.0, estimated in this work 10.1126/science.aaf8084. If not specified, its value is 4.0
- `-s` or `--selection`: Specify chromosome or region. This option is only required and works when the input file has [`cooler`](https://github.com/open2c/cooler) format. The value of this option is passed to the `cooler.Cooler.matrix().fetch()` method. For details, please refer their [documentation](https://cooler.readthedocs.io/en/latest/concepts.html#matrix-selector).
- `-i` or `--iteration`: The method relies on iterative scaling to find the optimal parameters. This option specifies the number of iterations. Generally, the more iterations the model runs, the better results are. However, the convergence of the model slow down when iteration increases. For larger size of contact map and the mean distance map, the number of iterations needed to good convergence is larger. If not specified, its default value is 10000.
- `-r` or `--learning-rate`: The learning rate for iterative scaling. Higher learning rate achieves faster convergence. However, the model can crash if learning rate is too large. One should play around this option to see what works best. If not specified, its default value is 10.0
- `--input-type`: The type of the input file. To use the script, the type must be specified. The method can work on both the contact map (`cmap`) or distance map (`dmap`). This option is required.
- `--input-format`: The format of the input file. If the type of input file is Hi-C contact map, then the script support `cooler` format Hi-C contact map file or a pure text based file. In the text based file, each line corresponds to the row of the contact map. If the type of input file is mean distance map, then the script only support the text based file in which each line represents the row of the mean distance map. This option is required.
- `--log`: A log file will be written if this option is specified. The log file contains the data of cost versus iteration.
- `--no-xyzs`: Turn off writing x,y,z coordinates of genome structures to files.
