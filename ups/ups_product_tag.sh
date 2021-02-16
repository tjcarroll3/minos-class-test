#!/bin/bash
#

# ups_product_tag.sh <ups version>

usage()
{
   echo "USAGE: `basename ${0}` <ups version> <ddt_ups_version>"
   echo "Please run from the ups directory of the repository"
}


# -------------------------------------------------------------------
# start processing
# -------------------------------------------------------------------

ups_version=${1}
ddt_ups_version=${2}

if [ -z ${ups_version} ] || [ -z ${ddt_ups_version} ]
then 
    usage
    echo "Provided novasoft version: $ups_version"
    echo "Provided novaddt  version: $ddt_ups_version"
    exit 1
fi

# -------------------------------------------------------------------
# need to edit product_deps, build_novasoft.sh and bootstrap.sh
# to have the correct product version
# also need to tag the repository
# -------------------------------------------------------------------

# first edit the product_deps, bootstrap.sh, and build_novasoft.sh files
sed -i~ 's/parent novasoft develop/parent novasoft '${ups_version}'/'       product_deps
sed -i~ 's/pkgver=devel/pkgver='${ups_version}'/'                         bootstrap.sh	
sed -i~ 's/trunk/tags\/\${pkgver}/'              			  bootstrap.sh	
sed -i~ 's/pkgver=devel/pkgver='${ups_version}'/' 			  build_novasoft.sh

# now commit these changes to the repository
svn commit -m"preparing for ups version tag $ups_version"

# copy the trunk to the tag
svn cp svn+ssh://p-novaart@cdcvs.fnal.gov/cvs/projects/novaart/pkgs.svn/trunk \
svn+ssh://p-novaart@cdcvs.fnal.gov/cvs/projects/novaart/pkgs.svn/tags/${ups_version} \
-m'create tag for ups version ${ups_version}'

# the tag is good to go, so now revert the trunk to use devel instead of ${ups_version}
sed -i~ 's/parent novasoft '${ups_version}'/parent novasoft devel/'       product_deps
sed -i~ 's/novaddt          '${ddt_ups_version}'/novaddt          devel/' product_deps
sed -i~ 's/tags\/\${pkgver}/trunk/'                                       bootstrap.sh
sed -i~ 's/pkgver='${ups_version}'/pkgver=devel/'                         bootstrap.sh
sed -i~ 's/pkgver='${ups_version}'/pkgver=devel/'                         build_novasoft.sh
sed -i~ 's/novaddt        '${ddt_ups_version}'/novaddt        devel/'     ../CMakeLists.txt

svn commit -m'reset the scripts in the trunk to setup devel as the ups version'
svn commit -m'reset the CMakeLists.txt in the trunk to setup devel as the ups version for the ddt' ../CMakeLists.txt

exit 0

