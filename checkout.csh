#!/bin/tcsh -f
# Checkout Script for Experiment 'SCM_am4_xanadu_2020.03'
# ------------------------------------------------------------------------------
# The script created at 2021-12-06T09:35:18 via:
# /ncrc/home2/fms/local/opt/fre-commands/bronx-18/bin/fremake --link --ncores=8 --platform=ncrc4.intel18 --target=prod-openmp --walltime=120 --xmlfile=/ncrc/home2/Youtong.Zheng/awg/xanadu_2020.03/scm_xanadu_2020.03_xml/awg_xanadu_scm.xml SCM_am4_xanadu_2020.03
# ------------------------------------------------------------------------------

source $MODULESHOME/init/csh
echo Using source directory = /ncrc/home2/Youtong.Zheng/awg/xanadu_2020.03/SCM_am4_xanadu_2020.03/src...
cd /ncrc/home2/Youtong.Zheng/awg/xanadu_2020.03/SCM_am4_xanadu_2020.03/src

module avail git >& .git_avail
if (! -z .git_avail) then
    module load git
endif

unalias *

# ---------------- component 'fms'
echo "Cloning https://github.com/NOAA-GFDL/FMS.git on branch/tag main"
set git_output=`git clone -q --recursive -b main https://github.com/NOAA-GFDL/FMS.git >& /dev/stdout`
if ( $? != 0 ) then
     echo "$git_output" | sed 's/^/**GIT ERROR** /' > /dev/stderr
     exit 1
endif
# Additional checkout commands from XML file

            ( cd FMS  && git checkout  2020.03 )
                       
          

# ---------------- component 'mom6'
echo "Cloning http://gitlab.gfdl.noaa.gov/fms/ocean_shared.git on branch/tag master"
set git_output=`git clone -q --recursive -b master http://gitlab.gfdl.noaa.gov/fms/ocean_shared.git >& /dev/stdout`
if ( $? != 0 ) then
     echo "$git_output" | sed 's/^/**GIT ERROR** /' > /dev/stderr
     exit 1
endif
# Additional checkout commands from XML file

          git clone -b dev/gfdl/2018.04.06 https://github.com/NOAA-GFDL/MOM6-examples.git mom6
          pushd mom6
          git checkout dev/gfdl/2018.04.06  #needed for older git on zeus
          git submodule init src/MOM6 src/SIS2 src/icebergs tools/python/MIDAS
          git clone --recursive https://github.com/NOAA-GFDL/MOM6.git src/MOM6 
          git clone             https://github.com/NOAA-GFDL/SIS2.git src/SIS2
          git clone             https://github.com/NOAA-GFDL/icebergs.git src/icebergs
          git submodule update #This gets the right version of submodules
          popd

          pushd mom6
          set platform_domain = `perl -T -e "use Net::Domain(hostdomain) ; print hostdomain"`
          if ("${platform_domain}" =~ *"fairmont.rdhpcs.noaa.gov"* ) then
            ln -s /scratch4/GFDL/gfdlscr/pdata/gfdl_O/datasets/ .datasets
          else if ("${platform_domain}" =~ *"ccs.ornl.gov"* ) then
            ln -s /lustre/atlas/proj-shared/cli061/pdata/gfdl_O/datasets/ .datasets
          else
            ln -s /lustre/f2/pdata/gfdl/gfdl_O/datasets/ .datasets
          endif
          popd   

          test -e mom6/.datasets
          if ($status != 0) then
            echo ""; echo "" ; echo "   WARNING:  .datasets link in MOM6 examples directory is invalid"; echo ""; echo ""
          endif

        

# ---------------- component 'ice_sis'
echo "Cloning http://gitlab.gfdl.noaa.gov/fms/ice_sis.git on branch/tag master"
set git_output=`git clone -q --recursive -b master http://gitlab.gfdl.noaa.gov/fms/ice_sis.git >& /dev/stdout`
if ( $? != 0 ) then
     echo "$git_output" | sed 's/^/**GIT ERROR** /' > /dev/stderr
     exit 1
endif
echo "Cloning http://gitlab.gfdl.noaa.gov/fms/ice_param.git on branch/tag master"
set git_output=`git clone -q --recursive -b master http://gitlab.gfdl.noaa.gov/fms/ice_param.git >& /dev/stdout`
if ( $? != 0 ) then
     echo "$git_output" | sed 's/^/**GIT ERROR** /' > /dev/stderr
     exit 1
endif
# Additional checkout commands from XML file

          ( cd ice_sis    && git checkout xanadu )
          ( cd ice_param  && git checkout xanadu )
         
        

# ---------------- component 'land_null'
echo "Cloning http://gitlab.gfdl.noaa.gov/fms/land_null.git on branch/tag master"
set git_output=`git clone -q --recursive -b master http://gitlab.gfdl.noaa.gov/fms/land_null.git >& /dev/stdout`
if ( $? != 0 ) then
     echo "$git_output" | sed 's/^/**GIT ERROR** /' > /dev/stderr
     exit 1
endif
# Additional checkout commands from XML file

          ( cd land_null && git checkout xanadu )
         
        

# ---------------- component 'FMScoupler'
echo "Cloning https://github.com/NOAA-GFDL/FMScoupler.git on branch/tag 2020.03"
set git_output=`git clone -q --recursive -b 2020.03 https://github.com/NOAA-GFDL/FMScoupler.git >& /dev/stdout`
if ( $? != 0 ) then
     echo "$git_output" | sed 's/^/**GIT ERROR** /' > /dev/stderr
     exit 1
endif
# Additional checkout commands from XML file

           git clone --recursive -b 2020.01 https://github.com/NOAA-GFDL/atmos_drivers.git
           git clone --recursive -b xanadu http://gitlab.gfdl.noaa.gov/fms/atmos_phys.git
           git clone --recursive -b user/znt/scm_20201026 http://gitlab.gfdl.noaa.gov/fms/atmos_fv_dynamics.git
           git clone -q --recursive -b user/znt/scm_20210616 http://gitlab.gfdl.noaa.gov/fms/atmos_scm.git
          
          

exit 0
