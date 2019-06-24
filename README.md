# Anisotropic osmosis filter for shadow removal in images

**Author**: Simone Parisotto

**Other authors**: L. Calatroni, M. Caliari, C.B. Schoenlieb, J. Weickert

**Version 1.0**

**Date: 19/09/2018**

This is a companion software for the paper [Anisotropic osmosis filter for shadow removal in images](https://iopscience.iop.org/article/10.1088/1361-6420/ab08d2),
also available at [arXiv](https://arxiv.org/pdf/1809.06298.pdf) repository.
Please cite as:

```
@article{ParCalCalSchWei18,
 author        = {Parisotto, S. and Calatroni, L. and Caliari, M. and Sch\"{o}nlieb, C.B. and Weickert, J.},
 title         = {{Anisotropic osmosis filtering for shadow removal in images}},
 journal       = {Inverse Problems},
 year          = {2019},
 doi           = {10.1088/1361-6420/ab08d2}
}
``` 

#### Example
For a shadowed image (left), identify the shadow boundary (middle-left) and compare the results between the isotropic osmosis filter (middle-right) in [[Vogel](https://link.springer.com/chapter/10.1007/978-3-642-38267-3_31), [Weickert](https://link.springer.com/chapter/10.1007/978-3-642-40395-8_3)] and the proposed anisotropic osmosis model (right) in [[ParCalCalSchWei18](https://arxiv.org/abs/1809.06298)].

<img src="https://raw.githubusercontent.com/simoneparisotto/Anisotropic-osmosis-filter/master/runme_syntethic/results/case11/u_start_65.png" width="150px"> <img src="https://raw.githubusercontent.com/simoneparisotto/Anisotropic-osmosis-filter/master/runme_syntethic/results/case11/theta65.png" width="150px"> <img src="https://raw.githubusercontent.com/simoneparisotto/Anisotropic-osmosis-filter/master/runme_syntethic/results/case11/u_CLA_65.png" width="150px"> <img src="https://raw.githubusercontent.com/simoneparisotto/Anisotropic-osmosis-filter/master/runme_syntethic/results/case11/u_MIR_65.png" width="150px">

######  Software acknowledgements
* [Anisotropic Diffusion](https://uk.mathworks.com/matlabcentral/fileexchange/47244-anisotropic-diffusion-stable-scheme)
* The action of the matrix exponential function: [expleja](https://bitbucket.org/expleja/expleja)
* [Tensor Voting with steerable filters](https://uk.mathworks.com/matlabcentral/fileexchange/47398-tensor-voting-with-steerable-filters)
* Export figures: [export_fig](https://github.com/altmany/export\_fig)
* Cyclic color map: [phasemap](https://uk.mathworks.com/matlabcentral/fileexchange/57020-cyclic-color-map)

### License
[BSD 3-Clause License](https://opensource.org/licenses/BSD-3-Clause)