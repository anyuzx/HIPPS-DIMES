![Imgur](https://i.imgur.com/IyTmFyE.jpg)

This python program is the implementation of the HIPPS-DIMES method[^1][^2]. HIPPS-DIMES is a computational method based on the maximum entropy principle, with experimental measured contact map or pair-wise distances as constraints, to generate a unique ensemble of <ins>3D chromatin structures</ins>. In a nutshell, this program accepts the input file of a mean spatial distance map (which can be measured in Multiplexed FISH experiment) or a Hi-C contact map (which is converted to distance map internally), and generate an ensemble of individual chromatin conformations that are consistent with the input. The output conformations are stored as `.xyz` format files, and can be used to calculate quantities of interest and can be visualized using `VMD` or other compatible softwares. 

_The theory and applications of this method can be found in our work published in Physical Review X ([link](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.11.011051)) and in Nature Communications ([link](https://www.nature.com/articles/s41467-023-36412-4))_

# Documentation

## Requirements

- Python 3.8+

## Install

First, download this repository using,

```bash
git clone https://github.com/anyuzx/HIPPS-DIMES
```

Next, go into the repository folder and install required packages using the command below,

```bash
cd HIPPS-DIMES
pip3 install --editable .
```

This command will install the required packages, and install the script as a
python module. Note that you need to install `pip3` if it is not installed
already (Follow the instruction on the official document here
https://pip.pypa.io/en/stable/installing/). Once installed, you can call
`HippsDimes` directly in the terminal to run the script. The packages installed
are:

- `Click`
- `Numpy`
- `Scipy`
- `Pandas`
- `Tqdm`
- `Cooler`
- `rich`

## How to use

To get started, please go through the jupyter notebook **`walkthrough.ipynb`**
in this repository.

In addition, it is helpful to view the help information for each arguments and
options. To display help information, use

```bash
HippsDimes --help
```

### Input files

This script accept input files in two formats. If the input file is a Hi-C
contact map, it can be in either `.cool` format or pure text format. If the
input file is a mean spatial distance map, the script only accepts a pure text
formatted file. The text format for a matrix is the following: each row of the
file corresponds to the row of the matrix. Values are space-separated. The
content of the file should look like this,

```text
1  2  3
2  1  2
3  2  1
```

### Output files

This script will generate several files:

- A text file for the final simulated mean distance map
- A text file for the final simulated contact map
- A text file for the connectivity matrix
- A `.xyz` formatted file for the ensemble of genome structures generated (can
  be turned off)
- A csv formatted file for cost versus iteration data (can be turned off)

### Examples

#### Example 1

First, download a cooler format Hi-C contact map from
[here](https://drive.google.com/file/d/1eIxGv1JbIrEAVoUSQK_O_ebIjWo6toTJ/view?usp=sharing)
(**The file size is about 116 Mb**). This Hi-C contact map is for Chicken cell
mitotic chromosome, originally retrieved from
[GEO repository](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102740).
Rename it to `hic_example.cool`. Then execute the following command,

```bash
HippsDimes hic_example.cool test --input-type cmap --input-format cooler -s chr7:10M-15M -i 10 -e 10
```

This command tells the script to load the Hi-C contact map `hic_example.cool`
and perform the iterative scaling algorithm. The argument `test` instructs the
files names of output files start with `test_`. Option `--input-type cmap`
specifies that the input file is a contact map. Option `--input-format cooler`
specifies that the input file is a `cooler` file. Option `-s chr7:10M-15M`
specifies that the algorithm is performed on the region 10 Mbps - 15 Mbps on
Chromosome 7. Note that these three options are required and cannot be
neglected. **Some option arguments are optional, some are required. Please refer
to the section below and use `HippsDimes --help` for details**

When the program finishes, the script will generate several output files:
`test.xyz`, `test_connectivity_matrix.txt`, and `test_dmap_final.txt`.
`test.xyz` contains 10 sets of individual conformations of x, y, z coordinates
and can be viewed using `VMD` or other compatible visualization softwares.

#### Example 2

In this example, we use Hi-C contact map for HeLa cell line Chromosome 14 at
time point of 12 hours after the release from prometaphase. For the purpase of
demonstration, you can download the Hi-C `.cool` file from
[here](https://drive.google.com/file/d/1j-zfDUP6LOZGCxz9uA3LaMI372ct1cU_/view?usp=sharing)
which is origannly retreived from
[GEO repository](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102740)
under accession number GSE102740. **Before you download, note that the file has
size of about 655 Mb**. Once downloaded, execute the following command,

```bash
HippsDimes GSM3909682_TB-HiC-Dpn-R2-T12_hg19.1000.multires.cool::6 test --input-type cmap --input-format cooler -s chr14:20M-107M -i 10000 -e 10
```

Similar to the first example, this command tells the script to load the Hi-C
cooler file `GSM3909682_TB-HiC-Dpn-R2-T12_hg19.1000.multires.cool` and its group
6 (data for several different resolution is stored in different groups) and
perform the HIPPS/DIMES algorithm. In this example, we change the number of
iterations to be 10000 by using the option `-i 10000`. On a AMD Ryzen 5 3600 CPU
machine, it takes about 3-4 mins to finish the program. Once it is finished,
several ouput files are generated.

#### Additional examples

The jupyter notebook `walkthrough.ipynb` in this repository contains additional examples. In particular, if you would like see an example of direct application of HIPPS-DIMES on imaging data, please go through the notebook.

### Explanantion of the arguments and options

#### Argument

- `INPUT`: File path for the input file. The input file can be a Hi-C contact
  map or a mean spatial distance map as measured in Multiplexed FISH experiment.
- `OUTPUT_PREFIX`: Prefix for outputfiles. For instance, if one specify it to be
  `TEST`, then all the output files will start with `TEST_`.

#### Options

- `-k` or `--connectivity-matrix`: Provide the path to the existing connectivity
  matrix one would like to use as initialization. Useful if restarting using the
  result from the previous run.
- `-e` or `--ensemble`: Number of individual conformations to be generated. This
  script will generate an ensemble of structures consistent with the input Hi-C
  contact map or the mean spatial distance map. Each individual conformations
  are different from each other. You can specify how many such individual
  conformations you want to generate. If not specified, its value would be 1000.
- `-a` or `--alpha`: Value of the contact map to distance map conversion
  exponent. If the input file is Hi-C contact map, the method first convert the
  contact map to a mean spatial distance map. The equation of the conversion is
  d*{ij} ~ c*{ij}^{1/\alpha}. The default value of \alpha is 4.0, estimated in
  this work 10.1126/science.aaf8084. If not specified, its value is 4.0
- `-s` or `--selection`: Specify chromosome or region. This option is only
  required and works when the input file has
  [`cooler`](https://github.com/open2c/cooler) format. The value of this option
  is passed to the `cooler.Cooler.matrix().fetch()` method. For details, please
  refer their
  [documentation](https://cooler.readthedocs.io/en/latest/concepts.html#matrix-selector).
- `-m` or `--method`: Specify the method used for optimization. The default
  method is Iterative Scaling (IS). Currently, Iterative scaling (IS), gradient
  descent (GD) and direct inversion (DI) are supported.
- `-l` or `--lamd`: Specify the weight for L1 or L2 regularization. Default
  value is zero, meaning no regularization. Regularization is typically used to
  avoid over-fitting.
- `-r` or `--reg`: Specify the type of regularization. Default is L2
  regularization. L1 and L2 are supported.
- `-i` or `--iteration`: The method relies on iterative scaling to find the
  optimal parameters. This option specifies the number of iterations. Generally,
  the more iterations the model runs, the better results are. However, the
  convergence of the model slow down when iteration increases. For larger size
  of contact map and the mean distance map, the number of iterations needed to
  good convergence is larger. If not specified, its default value is 10000.
- `-r` or `--learning-rate`: Learning rate. This hyperparameter controls the
  speed of convergence. If its value is too small, then convergence is very
  slow. If its value is too large, the program may never converge. Typically,
  learning rate can be set to be 1-30 if use Iterative scaling method. It should
  be a very small value (such as 1e-8) when using gradient descent optimization.
  The default value is 10.0.
- `--input-type`: The type of the input file. To use the script, the type must
  be specified. The method can work on both the contact map (`cmap`) or distance
  map (`dmap`). This option is required.
- `--input-format`: The format of the input file. If the type of input file is
  Hi-C contact map, then the script support `cooler` format Hi-C contact map
  file or a pure text based file. In the text based file, each line corresponds
  to the row of the contact map. If the type of input file is mean distance map,
  then the script only support the text based file in which each line represents
  the row of the mean distance map. This option is required.
- `--log`: A log file will be written if this option is specified. The log file
  contains the data of cost versus iteration.
- `--no-xyzs`: Turn off writing x,y,z coordinates of genome structures to files.
- `--ignore-missing-data`: Turn on this argument will let the program ignore the
  missing elements or infinite number in the contact map or distance map
- `--balance`: Turn on the matrix balance for contact map. Only effective when
  `input_type == cmap` and `input_format == cooler`
- `--not-normalize`: Turn off the auto normalization of the contact map. Only
  effective when `input_type == cmap`
- `--enforce-nonnegative-connectivity-matrix`: Constrain all the "spring
  constants" to be nonnegative

### Tips for using this program

- In practice, contact map or distance map larger than 5000x5000 is too large
  for the method to converge. If your matrix is larger than 5000x5000, I suggest
  that you can either perform a coarse-graining on the original matrix to get a
  smaller one or you can use the model on a subregion of the contact
  map/distance map.
- When using Iterative scaling for optimization, the learning rate typically can
  be set between 1 and 50. You should try different values to see what is the
  optimal learning rate to use. For gradient descent, the learning rate
  typically needed to be set very small, such as 1e-7.
- If your contact map/distance map has a lot of missing or zero entries. You can
  try to turn on the option `--ignore-missing-data`. This will tell the code not
  considering these missing entries. Thus giving you a less biased result
- Whenever the contact map is feeded, the programe will normalize the contact
  map by dividing it by its maximum value entry. If you don't want this, you can
  set option `--not-normalize`. This will tell the code not normalize the
  contact map at all
- Note that when feeding the contact map, there is no physical length scale
  associated with it. Thus we cannot set a unit to the resulting distance matrix
  or the structures. In this sense, the structures generated are dimensionless.
  But one can use additional information to set the length scale of the problem.
  For instance, if you have a reasonable estimate of average distance between
  two nearest loci, then you can use this distance as the measure to rescale the
  structure to be consistent with it.

## How to cite

If you used this program in your publication, please cite the following
reference:

* _Shi, Guang, and D. Thirumalai. "From Hi-C Contact Map to Three-dimensional
Organization of Interphase Human Chromosomes." Physical Review X 11.1
(2021): 011051._

* _Shi, Guang, and D. Thirumalai. "A maximum-entropy model to predict 3D structural ensembles of chromatins from pairwise distances: Applications to Interphase Chromosomes and Structural Variants." bioRxiv (2022)._

[^1]: _Shi, Guang, and D. Thirumalai. "From Hi-C Contact Map to Three-dimensional
Organization of Interphase Human Chromosomes." Physical Review X 11.1
(2021): 011051._
[^2]: _Shi, G., Thirumalai, D. A maximum-entropy model to predict 3D structural ensembles of chromatin from pairwise distances with applications to interphase chromosomes and structural variants. Nat Commun 14, 1150 (2023)._
