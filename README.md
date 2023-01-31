# Diffractive Voronoi Lensless Imaging

### [Project Page](https://vccimaging.org/Publications/Fu2022Diffractive/) | [Paper](https://doi.org/10.1364/OE.475004)

[Qiang Fu](https://fuqiangx.github.io/), [Dong-Ming Yan](https://sites.google.com/site/yandongming/dong-ming-yans-homepage), [Wolfgang Heidrich](https://vccimaging.org/People/heidriw/)

We demonstrate the optimization of Voronoi-Fresnel phase using our proposed two-step quasi-Centroidal Voronoi Tessellation optimization scheme. In the first step, we fix the number of Voronoi sites K, and adopt a modified Lloyd iteration routine to maximize the panchromatic Modulation Transfer Function volume (MTFv). In the second step, we run the same optimization for different K, and find the best number of Voronoi sites. Please try the demo code ```demo_optimization.m``` for a small scale optimization (this is part of the implementation of Fig. 3 in our paper).

If you find our work useful in your research, please cite the following paper.

```
@article{fu2022diffractive,
  title={Diffractive lensless imaging with optimized {Voronoi-Fresnel} phase},
  author={Fu, Qiang and Yan, Dong-Ming and Heidrich, Wolfgang},
  journal={Optics Express},
  volume={30},
  number={25},
  pages={45807--45823},
  year={2022},
  publisher={Optica Publishing Group}
}
```

## Requirements
This code is developed with Matlab R2020b. We use external third-party codes for the Voronoi diagram, which are included in the sub-directory ```3rdparty```. Terms of related licenses apply.

## License
Our code is licensed under BSL-1. By downloading the software, you agree to the terms of this License. 

## Questions
Should you have any questions, please feel free to contact Qiang Fu at qiang[dot]fu[at]kaust[dot]edu[dot]sa.
