# CIDR metagenomics_workflow

This release has been designed to work from inside of a singularity container. 

# Metagenomics singularity/apptainer setup and execution
A guide for running CIDR metagenomics from singularity/apptainer images. The tools here can be cloned and built independent of containers from their GitHub repos, [organism_query](https://github.com/GSTT-CIDR/organism_query)_and [CIDR_metagenomics (not up yet)](https://github.com/).

## Rationale
Using Singularity/Aptainer for the distribution of CIDR metagenomics software across the NHS metagenomics NoE enables standardised and repeatable use of the informatics workflows across trial sites. It reduces intrusion on partner's hardware, hands-on set up time/support and streamlines updating. The container operates locally and without root privilege, **so no data will leave the sequencing device.**

## Apptainer installation
Install singularity on the host machine through either:
1. OS package manager: yum/apt [Apptainer guide](https://apptainer.org/docs/admin/main/installation.html)
2. Mamba: ```mamba install apptainer```
3. From source: [Apptainer guide](https://apptainer.org/docs/admin/main/installation.html)

## Building apptainer image from definition
**If you have already received the _.sif_ image, proceed to the next step _Running apptainer images_.** 
1. Clone the container definition repo: ```git clone``` 
2. Navigate to the directory containing the container definition
3. Build the image : ```sudo apptainer build --bind /:/mnt organism_query.sif metag_v1.def```
3. Build the image : ```sudo apptainer build --bind /:/mnt cidr_metagenomics.sif cidr_metagenomics.def```
     - _You may need to increase your shell's ulimit cap if building fails_: ```ulimit -n 20000```
## Running apptainer images
After execution, the luncher GUI should appear. Refer to the WEBPAGE for documentation on how to run.
```
apptainer run --bind /tmp/.X11-unix:/tmp/.X11-unix --bind /:/mnt --env DISPLAY=$DISPLAY organism_query.sif 
```

##For the metagenomics
apptainer run --bind /tmp/.X11-unix:/tmp/.X11-unix --bind ../:/mnt --bind /data:/data --env DISPLAY=$DISPLAY cidr_metagenomics.sif