[![pipeline status](https://gitlab.tudelft.nl/idema-group/flippy/badges/master/pipeline.svg)](https://gitlab.tudelft.nl/idema-group/flippy/-/commits/master)
[![coverage report](https://gitlab.tudelft.nl/idema-group/flippy/badges/master/coverage.svg)](https://gitlab.tudelft.nl/idema-group/flippy/-/commits/master)
[![release version](https://img.shields.io/badge/dynamic/json?url=https://gitlab.tudelft.nl/idema-group/flippy/-/raw/master/VERSION.json&query=$.*&color=blue&label=version)](https://gitlab.tudelft.nl/idema-group/flippy/-/releases)
[![licence](https://img.shields.io/badge/licence-MIT-green)](https://gitlab.tudelft.nl/idema-group/flippy/-/blob/master/LICENSE)
[![EMail](https://img.shields.io/badge/EMail-D14836?logo=Mail.ru&logoColor=white&logoWidth=20)](mailto:flippy@mailbox.org)
# flippy
c++ package for dynamically triangulated membrane simulations.

<img src="assets/flippy.png" alt="flippy" width="300"/>

## Gallery

<img src="https://surfdrive.surf.nl/files/index.php/s/6HCtX7B4NwE6w48/download" alt="rbc" width="300">
<img src="https://surfdrive.surf.nl/files/index.php/s/1O3oJBNMT0vLxIv/download" alt="bubble_collision" width="300">
<img src="https://surfdrive.surf.nl/files/index.php/s/mua39VAvdxmobD4/download" alt="cell division" width="300">

## Current Version

current release [![release version](https://img.shields.io/badge/dynamic/json?url=https://gitlab.tudelft.nl/idema-group/flippy/-/raw/master/VERSION.json&query=$.*&color=blue&label=version)](https://gitlab.tudelft.nl/idema-group/flippy/-/releases) is experimental. This means that the API may change significantly. 

## Support
Flippy is still in active development and the documentation is almost non-existent, but I try my best to reply to e-mails and provide support for scientists who want to use flippy.
### for questions about general usage
please use the support email [![EMail](https://img.shields.io/badge/EMail-D14836?logo=Mail.ru&logoColor=white&logoWidth=20)](mailto:flippy@mailbox.org).
### for bugfixes 
please create an [issue](https://gitlab.tudelft.nl/idema-group/flippy/-/issues).
### for feature requests
again the [issues](https://gitlab.tudelft.nl/idema-group/flippy/-/issues) page can be used but be aware that new features will be slow to come.
### want to join the team?
write me an email!


## Release Version Nomenclature

This repository follows [Semantic Versioning](https://semver.org/) guidelines.

Given a version number *MAJOR*.*MINOR*.*PATCH*, the flippy project will increment the:

- *MAJOR* version when you make incompatible API changes,
- *MINOR* version when you add functionality in a backwards compatible manner, and
- *PATCH* version when you make backwards compatible bug fixes.

Additional labels for pre-release and build metadata are available as extensions to the *MAJOR*.*MINOR*.*PATCH* format.

## licence 

*flippy* is under MIT License, which means you can do pretty much anything with it. For more information read the `LICENCE` file.

## How to get it

*flippy* is a headers only library, so all you need to is to download the `flippy` subfolder and copy it into your project.

## Documentation
You can find flippy's [user manual](https://gitlab.tudelft.nl/idema-group/flippy/-/wikis/User_Manual) and automatically generated code [documentation](https://gitlab.tudelft.nl/idema-group/flippy/-/wikis/Documentation) over on the [wiki](https://gitlab.tudelft.nl/idema-group/flippy/-/raw/master/docs/mainpage.md).

## Examples of usage

This is a simple example of a fluctuating sphere in under 100 lines. The code is not as elegant as it could be, but it is general enough to be used as a starting pint in a real project. This code assumes that the `flippy` project is in the same folder as the main file.
This code saves a data file in the folder where the binary file will be executed, under the name `test_run.json` which is a full snapshot of the final configuration.

```c++:demo/simplest_MC/main.cpp
```
This code can be found in the subfolder `demos/simplest_MC` together with a python script that visualizes the data, and a simple `CMake` file that can be used to build the code.
