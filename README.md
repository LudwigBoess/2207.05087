# CRESCENDO

In this repository you will find all scripts and dependencies to reproduce the figures presented in [Böss et al (submitted)](https://arxiv.org/abs/2207.05087).

To initialize the dependencies run `julia build.jl`.

This will install all packages as well as a custom registry and will download and unpack the test data.
Please note that the download is 3.3GB, and the unpacked test data is 12GB large!

For the simulation snapshots send an email to `lboess@usm.lmu.de` and I will provide them to you.

You can get all plots by running from the main repository directory: `julia src/fig1.jl`.
If you use an interactive session make sure to run it from the main folder, as this is where the environment is loaded.