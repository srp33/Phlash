# Phlash - server

## Project setup

- For the Docker image to build properly with GeneMarkS, you need to download it from the [GeneMark website](http://topaz.gatech.edu/GeneMark/license_download.cgi).
    - Download "GeneMarkS v.4.30" for "LINUX 64".
    - After registering, you should be given the opportunity to download files called "genemark_suite_linux_64.tar.gz" and "gm_key_64.gz".
    - Store the .gz files in this directory.
- Execute the `build_docker` script to build the Docker image, then execute the `run_docker` script to start the container for the back-end.
    - The server is configured to run on port 5000 although this may be changed in in the `run_docker` script.
