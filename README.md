![Forest compostional dynamics](https://github.com/bioatmosphere/UVAFME-VOC/blob/master/figures/Fig_Forest_Succession_Dynamics.jpg)

# UVAFME-VOC

![GitHub repo size](https://img.shields.io/github/repo-size/bioatmosphere/UVAFME-VOC)
![GitHub stars](https://img.shields.io/github/stars/bioatmosphere/UVAFME-VOC?style=social)
![GitHub forks](https://img.shields.io/github/forks/bioatmosphere/UVAFME-VOC?style=social)
![Twitter Follow](https://img.shields.io/twitter/follow/bioatmo_sphere?style=social)

UVAFME-VOC is an Individual-based forest volatile organic compounds emission model that is developed from the state-of-the-art forest gap model--University of Virginia Forest Model Enhanced (UVAFME). **Though with an VOC component, this model can TOTALLY be treated and applied as an individual-based forest gap model**, as the VOC component is literally added on top of UVAFME without changing any part of it. 

Please check out the  manuscript by Wang, Shugart, and Lerdau in [**Ecological Modelling**](https://doi.org/10.1016/j.ecolmodel.2017.02.006) for details.

UVAFME-VOC (v1.0) is written in Fortran 2003 (a Python-based version that would make it more widely accessible is under conceiving). Any questions related to the explicit VOC simulation in particular and the individual-based forest model in general can  be directed to B. Wang @ wbwenwen@gmail.com or bw8my@virginia.edu

### Run UVAFME-VOC

**Get the code/repo**

```shell
git clone https://github.com/bioatmosphere/UVAFME-VOC
```
**Understand the directory structure**:

1. **src/**: 24 in total including all source code (in .f90/.F90) and a Makefile.

2. **input_data/**: all input files needed

3. **output_data/**: where outputs reside

4. **file_list.txt**: a .txt file lising directories of all files needed

**Compile and Run the Program**

Navigate to the src/ folder and run the following commands from a shell:

```shell
# compile and move the successfully compiled .exe file out of /src above to its parent folder
# NOTE: a compiler (e.g. intel) is required
make UVAFME.exe
mv UVAFME.exe ..
# run the program
./UVAFME.exe file_list.txt
```

### Publications arising from this model

**Wang, B.**, Shugart, H. H., & Lerdau, M. T. (2017). [An individual-based model of forest volatile organic compound emissions—UVAFME-VOC v1.0](https://doi.org/10.1016/j.ecolmodel.2017.02.006). **Ecological Modelling**, 350, 69-78.

**Wang, B.**, Shugart, H. H., Shuman, J. K., & Lerdau, M. T. (2016). [Forests and ozone: productivity, carbon storage, and feedbacks](https://www.nature.com/articles/srep22133). **Scientific Reports**, 6, 22133.

**Wang, B.**, Shuman, J., Shugart, H. H., & Lerdau, M. T. (2018). [Biodiversity matters in feedbacks between climate change and air quality: a study using an individual‐based model](https://doi.org/10.1002/eap.1721). **Ecological Applications**.

Shugart, H. H., **Wang, B.**, Fischer, R., Ma, J., Fang, J., Yan, X., ... & Armstrong, A. H. (2018). [Gap models and their individual-based relatives in the assessment of the consequences of global change](https://doi.org/10.1088/1748-9326/aaaacc). **Environmental Research Letters**, 13, 033001.

**Wang, B.**, Shugart, H.H. & Lerdau, M.T. (2019). [Complexities between plants and the atmosphere](https://doi.org/10.1038/s41561-019-0413-8). **Nature Geoscience** 12, 693–694 
