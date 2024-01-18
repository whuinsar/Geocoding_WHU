# Geocoding_WHU
The tailored sources, compiled binaries and datasets for testing, coupled with our work "Fast and accurate SAR geocoding with a plane approximation."

# Surfing the sources & compiling
Development environment: Visual Studio 2019 + Windows 10 + CMake + vcpkg
With vcpkg toolchain selected, the vcpkg.json shall be parsed to install all necessary packages.

Three implementations, the plain forwarding method, Eineder's polynomial approximation method, and ours are provided in the subdirectory 'geocoding.' 

# Testing
1. Download and unzip the last release with Bandzip; the default Windows zip tools fail to uncompress this file.
2. In the Powershell terminal, navigate to the subdirectory `ALOS-1` (or `Sentinel-1`, `TSX`), and execute the bash script `geocoding.bat.`
3. Use the MATLAB script in the top-level folder to check the mosaicked DEMS and the geocoding results. As we cannot upload the entire binary SLC file, one can evaluate its reliability by exploring the spatial distributions of the cites when multiple solutions are found, investigating the layover artifacts and holes resulting from the shadow within the DEM of the slant-range coordinate system, or scrutinizing the layover range intervals.
