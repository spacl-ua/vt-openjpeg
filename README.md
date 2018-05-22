VT-OPENJPEG Package
================================
VT-OpenJPEG is an open-source implementation of visually-optimized JPEG 2000 encoding derived from OpenJPEG (https://github.com/uclouvain/openjpeg).
It also contains a Matlab toolbox to measure visibiility threshold (VT) of each wavelet coefficients subbands.
The VT-based JPEG2000 encoder is described in references \[1,4] and the steps to measure VTs  are described in references \[1,2,3].

Installation
------------
There are two ways to download the latest version of the vt-openjpeg project
* (Option 1) Download the zip file at https://github.com/yzlinaz/vt-openjpeg/archive/master.zip
* or (option 2) Git clone with a git client:
  ```
  git clone  https://github.com/yzlinaz/vt-openjpeg.git
  ```

The easiest way of building the project is to use the cmake & make system provided by the original OpenJPEG project, i.e. at the current folder type:

```
cmake .
make
```



Optionally, you can make a system-wide installation with

```
sudo make install
```

The compressor and decompressors are created at:
* `bin/opj_compress` (compressor)
* `bin/opj_decompress` (decompressor)

More build and installation options (e.g. Using third-party library for TIFF/PNG format) can be achieved by following the instruction provided in INSTALL

Quick Manual
---------------------------------

The simplest way to run the encoder is by specifiying the desired visual quality. Visual quality is determined by Just Noticeable Distortion (JND). For example, to encode an image at JND=1, the following commandline can be used:

```
opj_compress -i test.bmp -o test.jp2 -I -JND 1 -monitor ../data/asus-pa328q.txt    
```    

Since JND is dependent on display characteristics, the software allows specification of the monitor parameters. In the above example, the monitor parameters are provided in ../data/asus-pa328q.txt
For more advanced usage please refer to the wiki page: https://github.com/yzlinaz/vt-openjpeg/wiki/Usage


Details on folders hierarchy:
-----------------------------
_NOTE_: The symbol "+" denotes the particular folder has been modified from the original OpenJPEG implementation to enable VT-optimized encoding.

* `src/` (+)
  * `lib/`
    * `openjp2/` (+): contains the sources of the openjp2 library (Part 1 & 2), which contains modification for VT-optimized encoder.
    * `openjpwl/`: contains the additional sources if you want to build a JPWL-flavoured library.
    * `openjpip/`: complete client-server architecture for remote browsing of JPEG 2000 images.
    * `openjp3d/`: JP3D implementation
    * `openmj2/`: MJ2 implementation

  * `bin/`: contains all applications that use the openjpeg library
    * `hvs/`: contains C code for calculation of quantization stepsize given JND level.
    * `common/`: common files to all applications
    * `jp2/` (+): a basic codec, which is able to encode/decode  VT-optimized codestream besides the original OpenJPEG codec's functionality
    * `mj2/`: motion JPEG 2000 executables
    * `jpip/`: OpenJPIP applications (server and dec server)
      * `java/`: a Java client viewer for JPIP
    * `jp3d/`: JP3D applications
      * `tcltk/`: a test tool for JP3D
    * `wx/`
      * `OPJViewer/`: gui for displaying j2k files (based on wxWidget)
* `matlab/` (+): contains the Matlab toolbox for VT measurement
* `vt_tests/` (+): contains shell script example that illustrates interface of the encoder.
* `data/` (+): contains example VT measurement data.
* `wrapping/`
  * `java/`: java jni to use openjpeg in a java program
* `thirdparty/`: thirdparty libraries used by some applications. These libraries will be built only if there are not found on the system. Note that libopenjpeg itself does not have any dependency.
* `doc/`: doxygen documentation setup file and man pages
* `tests/`: configuration files and utilities for the openjpeg test suite. All test images are located in 'http://openjpeg.googlecode.com/svn/data' folder.
* `cmake/`: cmake related files


Other informative files in this project:
----------------------------------------

* See `LICENSE` for license and copyright information.
* See `INSTALL` for installation procedures.
* See `NEWS` for user visible changes in successive releases.
* See `CHANGES`  for per-revision changes.


References:
----------------------------------------

\[1] H. Oh, A. Bilgin and M. W. Marcellin, "Visually Lossless Encoding for JPEG2000," in IEEE Transactions on Image Processing, vol. 22, no. 1, pp. 189-201, January, 2013.

\[2] A. B. Watson and D. G. Pelli, "QUEST: A Bayesian adaptive psychometirc method," Perception and psychophysics, vol. 33, no. 2, February, 1983.

\[3] "Psychtoolbox-3 software. Available online: [http://psychtoolbox.org/]."

\[4] F. Liu, Y. Lin, E. L. Ahanonu, M. W. Marcellin, A. Ashok, E. A. Krupinski and A. Bilgin, "Visibility Thresholds for Visually Lossy JPEG2000," Proc. SPIE 9971, Applications of Digital Image Processing XXXIX, 99711P, September, 2016.
