# Anisotropic osmosis filter for shadow removal in images


**Version 1.0**

**Date: 19/09/2018**

**Author**: Simone Parisotto

**Other authors**: L. Calatroni, M. Caliari, C.B. Schönlieb, J. Weickert

This is a companion software for the [submission](https://arxiv.org/pdf/1809.06298.pdf):

```
@article{ParCalCalSchWei18,
 author        = {Parisotto, S. and Calatroni, L. and Caliari, M. and Sch\"{o}nlieb, C.B. and Weickert, J.},
 title         = {{Anisotropic osmosis filtering for shadow removal in images}},
 year          = {2018},
 month         = {sep}, 
 journal       = {ArXiv e-prints},
 archivePrefix = {arXiv},
 eprint        = {1809.06298},
}
```

Please use the following entry to cite this code:

```
@Misc{ParCalCalSchWei18code,
 author       = {Parisotto, S. and Calatroni, L. and Caliari, M. and Sch\"{o}nlieb, C.B. and Weickert, J.},
 title        = {{Anisotropic osmosis filter for shadow removal in images}},
 howpublished = {GitHub repository},
 month        = {September},
 url          = {https://github.com/simoneparisotto/Anisotropic-osmosis-filter/},
 year         = {2018}
}
```

#### Example
For a shadowed image (left), identify the shadow boundary (middle-left) and compare the results between the isotropic osmosis filter in [[Vogel](https://link.springer.com/chapter/10.1007/978-3-642-38267-3_31), [Weickert](https://link.springer.com/chapter/10.1007/978-3-642-40395-8_3)] and the proposed anisotropic osmosis model [[ParCalCalSchWei18](https://arxiv.org/abs/1809.06298)].

<img src="https://raw.githubusercontent.com/simoneparisotto/Anisotropic-osmosis-filter/master/runme_syntethic/results/case11/u_start_65.png" width="150px"> <img src="https://raw.githubusercontent.com/simoneparisotto/Anisotropic-osmosis-filter/master/runme_syntethic/results/case11/theta65.png" width="150px"> <img src="https://raw.githubusercontent.com/simoneparisotto/Anisotropic-osmosis-filter/master/runme_syntethic/results/case11/u_CLA_65.png" width="150px"> <img src="https://raw.githubusercontent.com/simoneparisotto/Anisotropic-osmosis-filter/master/runme_syntethic/results/case11/u_MIR_65.png" width="150px">

######  Software acknowledgements
* [Anisotropic Diffusion](https://uk.mathworks.com/matlabcentral/fileexchange/47244-anisotropic-diffusion-stable-scheme)
* The action of the matrix exponential function: [expleja](https://bitbucket.org/expleja/expleja)
* [Tensor Voting with steerable filters](https://uk.mathworks.com/matlabcentral/fileexchange/47398-tensor-voting-with-steerable-filters)
* Export figures: [export_fig](https://github.com/altmany/export\_fig)
* Cyclic color map: [phasemap](https://uk.mathworks.com/matlabcentral/fileexchange/57020-cyclic-color-map)