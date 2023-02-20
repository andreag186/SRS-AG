#
# Enabling Provenance of Scientific Workflows via Jupyter:

#
# _A Genomic Workflow Case Study_

![Shape1](RackMultipart20230220-1-kopjg7_html_cb55ddb5edd60516.gif)

Author: Andrea Maria Grecu

(Student ID: 293287544, UPI: agre945, Current Programme: Postgraduate Diploma in Science)

Supervisor(s): Dr Allen Rodrigo and Dr Teng Li

Thomas Building, School of Biological Sciences, University of Auckland, Private Bag 92019, Auckland 1142 New Zealand.

## Career Development Statement

Within the three months I have spent working within Dr Allen Rodrigo's computational laboratory, I have learnt many fundamental skills required to work within the bioinformatics field. I have gained familiarity with up to 5 programming languages, the Nesi platform, bioinformatic methods/pipelines and a general understanding of what it is to conduct research. This summer research scholarship has afforded me the opportunity to realise my career.

## Abstract

Computational methods within bioinformatics often fail in reproducibility due to insufficient or inaccessible documentation. It is imperative to the progression of bioinformatics research that an effort is made to produce interactive methods accessible to a wider audience with limited computer literacy. Implementing provenance tracking, which involves meticulously tracking all data transformations within a workflow, enables anyone to reuse a published computational method.

Orr et al. (2020) described a phylogenomic method to measure somatic mutations within a phenotypically mosaic plant individual and further predict the somatic mutation rate. The bioinformatic workflow (aka. pipeline) and input data used are open-access and found within the original journal publication. We attempted to replicate the results using the pipeline published and failed due to workflow decay of the pipeline since publication.

This study aims to partially replicate the original pipeline created by Orr et al. (2020) on the interactive web-based computing platform; Jupyter Notebook. Jupyter Notebook supports various coding languages, combines text and code, and saves outputs and visualisations. Our notebook provides a computational narrative with clear explanations of each step, assumptions, logic, input data, and expected output. Our notebook facilitates the application of the pipeline to new data with detailed annotations and code shortcuts, requiring minimal computer literacy. This project showcases provenance tracking of bioinformatic pipelines on Jupyter to increase reproducibility and prevent workflow decay, serving as a model for future applications.

## Introduction

Reproducibility is a cornerstone of the scientific method ensuring the validity of scientific findings and building confidence in the knowledge generated by science. The publication of a 2016 survey revealing that 70% of researchers had tried and failed to replicate the results of another scientist's experiment sparked widespread attention in the scientific community, with many declaring a 'reproducibility crisis' (Baker, 2016).

This crisis affects numerous scientific disciplines; however, bioinformatics and wider scientific computing have received particularly negative notoriety (Ivie & Thain, 2019). Defining reproducibility in computational science can be challenging. For the purposes of this study, we will consider a computational method as reproducible if the results can be replicated using the same data **and** code, as defined by Whitaker (2017).

Provenance tracking in a computational workflow involves creating a detailed record of the data's history and tracking its changes as it moves through the workflow while providing detailed annotations on data transformations at each step. Incorporating provenance tracking in computational methods is required to prevent workflow breakdown (Kanwal et al., 2017).

In this study, we explore the process of incorporating provenance tracking in the preexisting genomic workflow or pipeline published by Orr et al. (2020) to improve reproducibility.

## Methods

A Linux system with Ubuntu was used in conjunction with the Conda package manager, Jupyter and Github Desktop to create the Jupyter Notebook and Github Repository linked [here](https://github.com/andreag186/SRS-AG). The notes below describe the procedure required to set up a Jupyter Notebook with a Conda Environment as a kernel and enabled version control via Github in sufficient detail for replication.

**Procedure** :

a. Download the necessary repository files from a given repository using the clone feature with [GitHub Desktop](https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/adding-and-cloning-repositories/cloning-and-forking-repositories-from-github-desktop)(Gilroy & Kaplan, 2019).

b. Install Conda Package Manager on the Linux system

i. Download the Anaconda Installer [here](https://www.anaconda.com/products/distribution)

ii. Run the following command in a terminal window

bash Anaconda-latest-Linux-x86\_64.sh

c. Install [Jupyter](https://jupyter.org/install) using Conda by running the command below in the terminal

conda install jupyter

d. Create a new Conda environment using the command below in the terminal, replacing "myenv" with your specific environment.

conda create --name myenv

e. Create a kernel for the new environment using the commands below in the terminal, replacing "myenv" with your environment name and "My Environment" with the desired kernel name to be displayed on Jupyter. See more [here](https://arshren.medium.com/how-to-setup-conda-environments-and-add-kernels-for-jupyter-notebook-f2ebf968a409).

conda install ipykernel

python -m ipykernel install --user --name myenv --display-name "My Environment"

f. Activate the conda environment by running the command below in the terminal (replace myenv with your environment name)

_Remember to do this before opening Jupyter Notebook_

conda activate myenv

g. Launch Jupyter Lab using the command below in the terminal.

jupyter lab

h. Create a new Jupyter notebook using the "Python 3" option and select the kernel created in step e.

i. Use GitHub for version control by creating a new repository on GitHub and linking it to the local repository containing the files and notebooks (Gilroy & Kaplan, 2019). See more [here](https://docs.qubole.com/en/latest/user-guide/notebooks-and-dashboards/notebooks/jupyter-notebooks/managing-jupy-notebook-versions/link-jupy-notebook-github.html).

j. Commit changes to the local repository and push them to the remote repository using the command line or GitHub Desktop. See more [here](https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/making-changes-in-a-branch/pushing-changes-to-github).

## Results and Discussion

**Orr et al. (2020). - Somatic Variation Detection Pipeline**

The study published by Orr et al. (2020) aimed to develop a method for measuring somatic mutations occurring across the branches of an individual plant, further using these measurements to estimate the somatic mutation rate of the individual. The somatic mutation rate estimated can be extrapolated to the species and higher taxa of the individual.

The individual studied by Orr et al. (2020) was a long-lived _Eucalyptus melliodora_ tree, described as phenotypically mosaic due to only one of the individual branches displaying resistance to defoliation by the Australian Christmas Beetle, the genus _Anoplognathus_.

There was no available published genome of _E. melliodora_; thus, Orr et al. (2020) developed a pseudo-reference genome using the genome of the closely related _Eucalyptus grandis_ as a reference. The computational method for the creation of a pseudo-reference genome is included in the pipeline published to an open-access repository.

The results produced by Orr et al. (2020) suggested that _E. melliodora_ evolved mechanisms to lower their mutation rate as a species. By publishing the computational methods used to obtain these results, further study on the evolution of somatic mutation rates in different plant taxa may be facilitated.

**Factors Preventing Reproducibility**

We attempted to run the computational method or pipeline designed by Orr et al. (2020) in the published repository ([found here](https://github.com/adamjorr/somatic-variation)) using the instructions and input [raw reads](https://www.ncbi.nlm.nih.gov/bioproject/553104) provided and the [reference genome](https://www.ncbi.nlm.nih.gov/nuccore/AUSX00000000) provided by Myburg et al. (2014). The repository claims to contain all the necessary scripts required to run the pipeline, organised into several folders according to their function in the pipeline.

The repository suggests running the pipeline 'using the makefiles provided within the analysis folder(s)' provided. A makefile is a text file containing a set of rules and dependencies for building and compiling programs, frequently employed to create bioinformatics workflows (Spjuth et al., 2015). To utilise this method, the user must download the relevant folders to their workspace, provide the required data in the appropriate format within the data/ folder, and navigate to the corresponding directory to call the makefile Once in the correct directory, a makefile is called using the command below;

make

This study focuses on the 'Analysis Makefile', the first makefile that the repository suggests running. There are 28 rules within this makefile, representing 28 different steps producing associated output file(s) (Figure 1). The steps of the makefile must be run sequentially, with outputs from previous steps used as inputs for subsequent steps in addition to makefiles not supporting multi-threaded jobs (Spjuth et al., 2015).

The repository does not provide a method to run the steps individually, with the makefile running until failure and repeating previously completed steps when called again without script modification. Furthermore, each step takes approximately two days to complete due to the large data size, with errors often occurring hours or days after running the makefile.

This pipeline structure is not only frustrating to run, but greatly hinders reproducibility.

The Analysis makefile directly calls upon five scripts within the scripts folder provided by the repository. Across the pipeline, scripts provided are written in 5+ programming languages, each with its assumptions and parameters, most of which are not detailed. These scripts also call upon various third-party tools (Figure 2), introducing a high risk of workflow decay due to unpredictable accessibility or software updates of these tools (Hettne et al., 2012). We encountered a problem with the 'load-into-counting.py' script from the Khmer software, which was incompatible with Python 3.8.16. We resolved the issue by reverting to an earlier version of Khmer.

The instructions provided in the repository do not fully detail the computational capacity required to run the pipeline successfully, only detailing thread usage for a few scripts. For instance, we found that the create\_consensus script required 500GB of RAM to run successfully and took over five days to run using the '48 threads' suggested by the instructions.

While the repository provides some initialising instructions to install tools and set a path to the script folder, the computational environment in which they conducted the pipeline is not specified. When reproducing computational methods, the most common error occurs when running the workflow in a different environment, thus making it crucial to specify the original environment in which a pipeline is designed (Vitek & Kalibera, 2011).

**Jupyter Notebook; A Platform Enabling Provenance Tracking**

To overcome issues preventing reproducibility in the original pipeline published by Orr et al. 2020, we created a notebook using the Jupyter platform. The notebook provides a computational narrative with detailed explanations of each step, including inputs, outputs, code, assumptions/parameters, and the role of the step with respect to the pipeline and greater aim. Our notebook should enable the pipeline to be easily applied to new data sets by users with differing computational skills.

Code cells within the pipeline enable users to check computational capacity, install third-party tools and replicate the work environment, input new data, change script parameters, execute individual steps, and inspect outputs directly from the notebook (Figure 3).

We suggest Jupyter Notebook as an effective platform for provenance tracking and structuring complex genomic workflows. Using one platform to edit and run scripts in many languages with version control and store outputs in a user-friendly text-based interface sets Jupyter notebooks apart from other workflow management tools.

Issues previously identified in literature with Jupyter notebooks, such as a limited capacity for running tasks in parallel, incompatibility with air-gapped secure networks and automatic printing of large outputs (i.e. sequences), are rapidly being solved due to the platform's growing popularity (Colonnelli et al., 2022).

**Reproducibility Crisis**

The scientific method is often attributed to rapid research progress in Western Europe during the 1600s, from Copernicus to Newton (Betz, 2011). However, the scientific method's basic principles have existed in indigenous knowledge systems since the beginning of human life (Liebenberg, 1990). Observing, testing, and sharing results to confirm them is the essence of the scientific method and the foundation of human knowledge.

Reproducibility failures can motivate scientific inquiry, as seen in the case of Orr et al. 's . (2020) non-reproducible pipeline,which has inspired the "meta-research" described in this report, to improve its computational method and subsequently similar studies (Redish et al., 2018; Ioannidis, 2010). However, replicating studies face challenges in getting published due to publisher bias towards 'novel research' (Fidler & Wilcox, 2018). In the context of computational methods, inadequate details can hinder the replication of final results and impede others' comprehension of the research.

Although time-consuming, presenting comprehensive details and provenance tracking in a method is critical, particularly for complex bioinformatic pipelines like this one. It guarantees reproducibility, simplifies modifying parameters to explore diverse aspects, and enables scientific advancements based on findings.

By creating a user-friendly descriptive Jupyter notebook, we have made the pipeline more accessible, overcoming the limitations of the original repository. However, additional work is necessary to establish the reproducibility of the broader pipeline with a range of users and input data types.

## Figures
[Analysis Makefile Input-Output Figure.pdf](https://github.com/andreag186/SRS-AG/files/10788113/Analysis.Makefile.Input-Output.Figure.pdf)

**Figure 1. Input Files and Output Files associated with the Analysis Makefile of the Original Method.** Graphical representation of the input files required and output files produced by running the command "make" as the suggested method for executing the first makefile (Analysis Makefile) in the Orr et al. (2020) pipeline. The folders depicted in coloured rectangles are either provided in the GitHub repository or produced via the makefile (as depicted with the star symbol). The files are depicted in white rectangles and linked by arrows to their respective folder and other files they depend on.

![](RackMultipart20230220-1-kopjg7_html_932843848f692892.png)

**Figure 2. Scripts and Third Party Tools associated with the Analysis Makefile of the Original Method.** Graphical representation of all the scripts (provided by the GitHub repository or within the package) and third-party tool packages called upon the command "make" as the suggested method for executing the first makefile (Analysis Makefile) in the Orr et al. (2020) pipeline. Black arrows indicate which variables call one another.

![Shape 2](RackMultipart20230220-1-kopjg7_html_32fc62a78f152592.gif)

**A**

 ![](RackMultipart20230220-1-kopjg7_html_f0b478ca532ab52a.png)

![Shape 2](RackMultipart20230220-1-kopjg7_html_32fc62a78f152592.gif)

**B**

 ![](RackMultipart20230220-1-kopjg7_html_cc3c82c3f4de2b5e.png)

![Shape 2](RackMultipart20230220-1-kopjg7_html_32fc62a78f152592.gif)

**C**

 ![](RackMultipart20230220-1-kopjg7_html_527af384dff8a75f.png)

![Shape 2](RackMultipart20230220-1-kopjg7_html_6e969b2b62a0ae61.gif) ![Shape 2](RackMultipart20230220-1-kopjg7_html_4d3840cdd448d04b.gif)

**D**

**E**

 ![](RackMultipart20230220-1-kopjg7_html_d159c16eff979a02.png) ![](RackMultipart20230220-1-kopjg7_html_37d6e45e88eb916.png)

![Shape 2](RackMultipart20230220-1-kopjg7_html_32fc62a78f152592.gif)

**F**

 ![](RackMultipart20230220-1-kopjg7_html_db3e6832c869e8d2.png)

**Figure 3. Features of Jupyter Notebook Pipeline**

Screenshots from code cells of the Jupyter Notebook showing (A) computational capacity checking, (B) installation of a third-party tool via the package manager Conda, (C) altering input files and script parameters, (D) execution of a single step integrating Boolean expressions and an if/elif loop (E) for the clean\_reads scripts which run steps 1-3 and (F) printing and inspection of output(s).

## References

Baker, M. (2016). 1,500 scientists lift the lid on reproducibility. _Nature_, _533_(7604), 452–454. https://doi.org/10.1038/533452a

Bertrand, B., & Harrisson, A. (2019). Building and Packaging EPICS Modules With Conda. _International Conference on Accelerator and Large Experimental Physics Control Systems_, _17_, 223–227. https://doi.org/10.18429/JACoW-ICALEPCS2019-MOPHA014

Betz, F. (2011). _Origin of Scientific Method_. ResearchGate; unknown. https://www.researchgate.net/publication/251204905\_Origin\_of\_Scientific\_Method

Colonnelli, I., Aldinucci, M., Cantalupo, B., Padovani, L., Rabellino, S., Spampinato, C., Morelli, R., Di Carlo, R., Magini, N., & Cavazzoni, C. (2022). Distributed workflows with Jupyter. _Future Generation Computer Systems_, _128_, 282–298. https://doi.org/10.1016/j.future.2021.10.007

Fidler, F., & Wilcox, J. (2018). _Reproducibility of Scientific Results (Stanford Encyclopedia of Philosophy)_. Stanford.edu. https://plato.stanford.edu/entries/scientific-reproducibility/#ReplRepeReprScieResu

Hettne, K., Wolstencroft, K., Khalid Belhajjame, Goble, C., Mina, E., Harish Dharuri, Verdes-Montenegro, L., Garrido, J., David de Roure, & Roos, M. (2012). Best practices for workflow design: how to prevent workflow decay. _Research Explorer the University of Manchester_, 23. https://research.manchester.ac.uk/en/publications/best-practices-for-workflow-design-how-to-prevent-workflow-decay

Ioannidis, J. P. A. (2010). Meta-research: The art of getting it wrong. _Research Synthesis Methods_, _1_(3-4), 169–184. https://doi.org/10.1002/jrsm.19

Ivie, P., & Thain, D. (2019). _Reproducibility in Scientific Computing | ACM Computing Surveys_. ACM Computing Surveys (CSUR). https://dl.acm.org/doi/abs/10.1145/3186266?casa\_token=6q40fcniYeYAAAAA:1Q7pC8bnm6h1J9lt7yp6uKz7mnoho5Yn4qVzKoS9yGaa-AeeHnwyUqvfPhKqHW-zsvO8b7u\_JLSPNA

Kanwal, S., Khan, F. Z., Lonie, A., & Sinnott, R. O. (2017). Investigating reproducibility and tracking provenance – A genomic workflow case study. _BMC Bioinformatics_, _18_(1). https://doi.org/10.1186/s12859-017-1747-0

Liebenberg, L. (1990). _The Art of Tracking_. David Phillip Pub. https://www.academia.edu/download/32658743/Liebenberg-1990-The-Art-of-Tracking.pdf

Madrid-Márquez, L., Rubio-Escudero, C., Pontes, B., González-Pérez, A., Riquelme, J. C., & Sáez, M. E. (2022). MOMIC: A Multi-Omics Pipeline for Data Analysis, Integration and Interpretation. _Applied Sciences_, _12_(8), 3987. https://doi.org/10.3390/app12083987

Myburg, A. A., Grattapaglia, D., Tuskan, G. A., Hellsten, U., Hayes, R. D., Grimwood, J., Jenkins, J., Lindquist, E., Tice, H., Bauer, D., Goodstein, D. M., Dubchak, I., Poliakov, A., Mizrachi, E., Kullan, A. R. K., Hussey, S. G., Pinard, D., van der Merwe, K., Singh, P., & van Jaarsveld, I. (2014). The genome of Eucalyptus grandis. _Nature_, _510_(7505), 356–362. https://doi.org/10.1038/nature13308

Redish, A. D., Kummerfeld, E., Morris, R. L., & Love, A. C. (2018). Reproducibility failures are essential to scientific inquiry. _Proceedings of the National Academy of Sciences_, _115_(20), 5042–5046. https://doi.org/10.1073/pnas.1806370115

Spjuth, O., Bongcam-Rudloff, E., Hernández, G. C., Forer, L., Giovacchini, M., Guimera, R. V., Kallio, A., Korpelainen, E., Kańduła, M. M., Krachunov, M., Kreil, D. P., Kulev, O., Łabaj, P. P., Lampa, S., Pireddu, L., Schönherr, S., Siretskiy, A., & Vassilev, D. (2015). Experiences with workflows for automating data-intensive bioinformatics. _Biology Direct_, _10_(1). https://doi.org/10.1186/s13062-015-0071-8

Vitek, J., & Kalibera, T. (2011). _Repeatability, reproducibility, and rigor in systems research | Proceedings of the ninth ACM international conference on Embedded software_. ACM Conferences. https://dl.acm.org/doi/abs/10.1145/2038642.2038650?casa\_token=4VOGgOl41RcAAAAA:DepYJ-KNOVPvhoK\_\_DewU4Po-BetD26q1-ulRK23vsp6oLXQ9BYnU0e5bJDULnKAQxiBMxlvc7cxMg

Whitaker, K. (2017). Showing your working: a how to guide to reproducible research. _Figshare_. https://doi.org/10.6084/m9.figshare.5443201.v1