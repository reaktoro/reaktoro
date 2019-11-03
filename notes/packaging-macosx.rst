Packaging Reaktoro in MacOSX
============================

Firstly, make sure you have Xcode installed. Then, open Xcode, and go to
`Xcode → Open Developer Tool → More Developer Tools...` and search for
*Auxiliary Tools for Xcode*. Download `Auxiliary Tools for Xcode - Late July 2012`,
since the later versions do not ship with `PackageMaker`.

Open the downloaded file, search for `PackageMaker` and copy it to Applications.
After this, CPack will be able to pack Reaktoro executables and libraries into a MacOSX installer.

The guidelines for setting up CPack variables can be found [here](https://cmake.org/Wiki/CMake:Component_Install_With_CPack).
Read also [here](https://developer.apple.com/library/mac/documentation/Darwin/Conceptual/KEXTConcept/KEXTConceptPackaging/packaging_tutorial.html)
for more information about PackageMaker.
