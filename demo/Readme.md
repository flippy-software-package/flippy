# flippy demos

each sub-folder contains a simulation based on flippy.
In each folder you will (at least) find:
 - `main.cpp` file containing the `c++` code of the simulation
 - `CMakeLists.txt` file that contains build instructions for the `c++` code
 - `Readme.md` file that explains what the simulation is about
 
In order to compile the simulations you can use the cmake file. From shell this can be done in few steps:
1. open command line in the demonstration folder or navigate to it.
2. create a build directory
```shell
mkdir cmake-build
```
3. change to the build directory
```shell
cd cmake-build
```
4. create build files with cmake
```shell
cmake ..
```
5. build the executable 
```shell
cmake --build . --config Release
```
The executable file will be now located  in the `cmake-build` folder and have the same name as specified in `CMakeLists.txt` in the `project(...)` statement. One caveat for `MSVC` users, the executable will be located in the `Release` sub-folder.


## how to use shell
### MacOS
Open the `terminal` app which comes preinstalled. If you have never used this app before then you will need to first install the command line tools by executing
```shell
xcode-select —install
```
a dialog window will appear, which you should confirm by clicking install.
This step will take a bit and will install to command line tools and a `c++` compiler on your mac.

### Windows
Easiest way to get shell on windows (that I am aware of) is to install the [`PowerShell`](https://www.microsoft.com/en-us/p/powershell/9mz1snwt0n5d#activetab=pivot:overviewtab).
### linux
If you are a linux user you know this.