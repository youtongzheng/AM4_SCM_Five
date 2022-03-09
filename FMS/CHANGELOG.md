# ChangeLog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0),
and this project uses `yyyy.rr[.pp]`, where `yyyy` is the year a patch is released,
`rr` is a sequential release number (starting from `01`), and an optional two-digit
sequential patch number (starting from `01`).

## [2020.03] - 2020-10-08
### Added
- FMS2_IO: Adds header_buffer_val to the fms2io namelist which sets the netcdf header size in bytes. The default value is 16kb
- FMS2_IO: Adds netcdf_default_format to the fms2io namelist which allows the user to change the netcdf file type. The default value is 64bit.
- FMS2_IO: Adds support to read netcdf string global attributes
- FMS2_IO: Adds an optional argument to open_file, `dont_add_res_to_filename`, which indicates that the filename should not be modified (default adds .res to restart file name)
- FMS2_IO: Modifies the `register_variable_attribute` and `register_global_attribute` interfaces by adding str_len as an argument. This is a workaround to get fms2io to work with PGI because they don't support class (*) with len=*.
- FMS2_IO: Adds unit test that tests `write_data` and `read_data` when using a domain with a mask table 
- FMS2_IO: Adds fms2io’s version of get_mosaic_tile_grid
- MPP_IO: Adds `-Duse_mpp_io` compile option for data_override, interpolator, amip_interp, diag_manager, topography, and xgrid to select using mpp_io instead of fms2_io
- MPP_INIT: Adds unit tests for routines/functions that are called in mpp_init
- CMAKE: A cmake build system has been added with a CI build using cmake
### Changed
- FMS2_IO: Improves performance of previous release by gathering the domain decomposed data into one global buffer and doing one write rather than doing multiple reads
- DATA_OVERRIDE: Changes line in time_interp_external2 to enable 3D overrides

### Fixed
- DATA_OVERRIDE[2]: Fixes a crash when doing ongrid data_override calls with a domain with halos
- DIAG_MANAGER[2]: Fixes an issue where time_bnds were written incorrectly for the last time stamp 
- DIAG_MANAGER: Regional diagnostics with a mask table now work
- FMS2_IO: Unit test includes fms2io_init call to improve functionality 
- MPP: BOZ literals that are used in variable declaration are converted to integers using the int() function.
- MPP_DOMAINS2: Fixed unit test
- FMS_IO: Changes the logic in get_tile_string to fix bug where tile numbers 9 and 99 produce an inappropriate error

### Tag Commit Hashes
- 2020.03-beta4 (4d38679c1e18e920feb03d69f8a9762eb6a047aa)
- 2020.03-beta3 (521a15135a99d1f2da7d82f238353945f82ce1dd)
- 2020.03-beta2 (3dae0dfa405d555ecc09bbd2d60a1be24461f69e)
- 2020.03-beta1 (f7f1c1c73c1f478a53e84caee6aff2fa840ad086)
- 2020.03-alpha1 (2dd30b7ca0ac75a4a38b969e4a6d446eb395b4dd)


## [2020.02] - 2020-05-01
### Added
- FMS2_IO:  An fms2_io_nml namelist has been created.  It includes the variable ncchksz. This is the replacement for the environment variable NC_BLKSZ set in model run scripts and used by mpp_io.  The default value is 64 KB.  Any time a file is opened in fms2_io (nf90_open or nf90_create), the optional argument `chunksize=ncchksz` is passed to the NetCDF library.  NetCDF attempts to use this value to control the blocksize utilized for reads and writes of data from the filesystem.
- FMS2_IO: Adds support to `compute_global_cheksum.inc` for `real32`, assuming the flag `-DOVERLOAD_R4` is used when compiling.

### Changed
- MPP_DOMAINS - nesting:  The logic supporting nested domains for mosaic grids has been overhauled and extended.  FMS now supports multiple nests and telescoping nests (nest embedded within a nest).  The requirement for a nest to lie wholly within a single tile has been relaxed and a first-level nest may cross tile boundaries, but may not contain a tile corner.  Communications for two-way nesting have also been improved.  The 2019 December Public Release and 2020.02 GFDL Release within the  [GFDL_atmos_cubed_sphere] (https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere) have been updated and are compatible with this release of FMS.
- FMS2_IO:  The intent of fileobj is changed from (in) to (inout) in netcdf_restore_state_wrap, restore_domain_state, and netcdf_restore_state because the file object type has a pointer that is being reassigned in one of the routines lower in the call stack.

### Deprecated
- MPP_DOMAINS - nesting:  The initial nesting implementation is no longer supported.  Please see the Changed::MPP_DOMAINS sub-entry under.

### Removed
- GENERAL:  References to the macro _ALLOCATABLE have been replaced with “allocatable”, _ALLOCATED has been replaced with “allocated”, and _NULL has been removed.  It is now assumed that all compilers support the Fortran 2003 standard.  The macros still exist in fms_platforms.h for compatibility within other components.
- DIAG_MANAGER:  “fms_platform.h” is no longer included in any of the diag_manager routines.  Instead, fms_platform_mod is now being use-associated where necessary.  This fixes an issue for debuggers not providing correct line numbers. 

### Tag Commit Hashes
- 2020.02-beta1 (bbc6f8d33cfb75a411bbcd3f8423fa74b8b7cdfd)
- 2020.02-beta2 (6242941a632f6e261234f3a575e59efa1bfb1b36)
- 2020.02-beta3 (e6fd03b070eb1ba8ffbca740fd69681ee16d6be7)


## [2020.01] - 2020-03-13
### Added
- Adds the modules `axis_utils2`, `mosaic2`, and `time_interp_external2` that use `fms2_io`
- Adds unit tests for thread affinity in test_fms/affinity
- Autotools unit tests are now run with `srun, mpirun`, or `aprun` (whichever is present on your system)
- `fms2_io` provides three new derived types, which target the different I/O paradigms used in GFDL GCMs:
  - type(FmsNetcdfFile_t) - This type provides a thin wrapper over the netCDF4 library, but allows the user to assign a “pelist” to the file.  If a pelist is assigned, only the first rank on the list directly interacts with the netCDF library, and performs broadcasts to relay the information to the rest of the ranks on the list.  This derived type also allows the user to perform “compressed” reads/writes and restart variable “registration” (i.e. the ability to store pointers to user-defined/allocated buffers) and reads/writes in a manner similar to the existing calls in fms_io.
  - type(FmsNetcdfUnstructuredDomainFile_t) - This type extends the FmsNetcdfFile_t type to support “domain-decomposed” reads/writes on a user defined mpp_domains unstructured grid.  Users are required to inform the FMS2_io module which dimensions correspond to the unstructured grid using the appropriate `register_axis` call before any domain-decomposed reads/writes can be performed.
  - type(FmsNetcdfDomainFile_t) - This type extends the FmsNetcdfFile_t type to support “domain-decomposed” reads/writes on a user-defined `mpp_domains` two-dimensional lon-lat or cubed-sphere grid.  Users are required to inform the `fms2_io` module which dimensions (‘x’ or ‘y’) correspond to the two-dimensional grid using the appropriate `register_axis` calls before any domain-decomposed reads/writes can be performed.
- `fms2_io` requires the user to manage objects of the types described above.  Calls to `open_file` (`close_file`) act as constructors (destructors) for each of the objects, and must be explicitly made by the user.  Each object must be constructed before any I/O can be performed.  Examples describing how to use these new types are available in [test_fms/fms2_io](https://github.com/NOAA-GFDL/FMS/tree/master/test_fms/fms2_io)
- `fms2_io` treats *_FillValue* attributes as *valid_max* or *valid_min* range specifiers if none of *valid_range*, *valid_min*, and *valid_max* are specified as described in these [netcdf conventions](https://www.unidata.ucar.edu/software/netcdf/docs/attribute_conventions.html)

### Changed
- The `diag_manager` IO is handled by fms2_io instead of mpp_io. Default behavior assumes that the mpp_io namelist variable is set to *cf_compliant = .true.*
- The user must specify the diagnostic attributes that they want to write to the output files.  Example: If there is no *units* attribute, then the variable metadata will not include *units*, and it will not automatically add *units = “none”*.
- Calls to `register_diag_axis` for an X or Y axis that is shifted from the *CENTER* position need to include the optional argument *domain_position* and should be equal to *EAST* or *NORTH* based on the position relative to the domain. EAST and *NORTH* are exposed through `diag_manager_mod`.
- Changed the handling of *average_T* and *time_bnds* variables so that they are set to values that are sent in and are not manipulated as was the case in mpp_io.
- `interpolator`, `xgrid`, `data_override`, and `amip_interp` now call  `fms2_io` routines
- Support for enabling/disabling quad-precision (used in certain calculations) has been changed to a hard on/off switch.  Default behavior is quad-precision disabled. To enable, add the following CPP macro -DENABLE_QUAD_PRECISION.  This change was necessary to remove guessing at the proper setting via a mix of compiler vendor and system-defined environment variables which resulted in different behaviors on machines unbeknownst to the user.

### Deprecated
- fms2_io does **NOT** use the *scale_factor*, *add_offset*, or other attributes to manipulate the data.  The variable/data is returned to the caller as it appears in the file.  All post-read data manipulations should be handled by the caller.
### Removed
- Removes the use of bats when running unit tests

### Tag Commit Hashes
- 2020.01-alpha1 (09dc8e9e0f1c852e9e9190834176d16943cd3729)
- 2020.01-beta1 (e1c0d9d01d844938adc0d18afa09532f336bcdfe)
- 2020.01-beta2 (b68de5382a5ce631ddd6167de8d85f7c9ae54351)

## [2019.01] - 2019-11-26
### Added
- switch from "city" versioning style to `yyyy.<2_digit_version_number>[.<2-digit-patch number>]` style
- main development branch is `master` instead of `dev/master`
- affinity handling moved to the affinity directory
- support for building with autotools
- [fms_platform.h](include/fms_platform.h) contains directives that support building on macOS
- unit and build tests are available in the [test_fms](test_fms) directory
- updated [fv3gfs/makefile](fv3gfs/makefile) for use with current EMC build system
### Fixed
- Fixed `time_interp_missing` parameter in [time_interp/time_interp_external.F90](time_interp/time_interp_external.F90) to be within range when compiled in mixed-mode.
- reverted `QUAD_PRECISION` cpp macro behavior to pre-Xanadu behavior in [include/fms_platform.h](include/fms_platform.h)
- Fixed a GNU compiler issue with the logical check to set the netCDF fill value in `mpp_io_write` by separating the logical `.AND.` into nested `if` statements.
