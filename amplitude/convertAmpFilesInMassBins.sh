#!/bin/bash
##########################################################################
#
#    Copyright 2010
#
#    This file is part of rootpwa
#
#    rootpwa is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    rootpwa is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
#
##########################################################################
#-------------------------------------------------------------------------
# File and Version Information:
# $Rev::                             $: revision of last commit
# $Author::                          $: author of last commit
# $Date::                            $: date of last commit
#
# Description:
#      converts all .amp files in mass bin directories into ROOT files
#
#      uses PWA environment variable(s)
#      PWA_ENV_SET
#      ROOTPWA
#      PWA_DATA_DIR
#
#
# Author List:
#      Boris Grube          TUM            (original author)
#
#
#-------------------------------------------------------------------------


echo ">>> info: ${0} started on $(date)"
echo ">>> info: converting .amp files in ${PWA_DATA_DIR}"
echo
if [[ -z "${PWA_ENV_SET}" ]]
then
    echo "!!! error: PWA environment is not setup. source setupThis.sh!"
    exit 1
fi


RECREATE_LINKS="false"


# find mass bin directories
MASS_BINS=$(find ${PWA_DATA_DIR} -type d -regex '.*/[0-9]+.[0-9]+' -printf '%f\n' | sort -n)
if [[ -z "${MASS_BINS}" ]]
then
    echo "!!! error: cannot find any mass bins in ${PWA_DATA_DIR}"
    exit 1
fi


for MASS_BIN in ${MASS_BINS}
do
    DIR=${PWA_DATA_DIR}/${MASS_BIN}
		for AMP_FILE in ${DIR}/{,PSP}AMPS/*.amp
		do
				ROOT_FILE=${AMP_FILE%.amp}.root
				echo ">>> info: converting '${AMP_FILE}' to '${ROOT_FILE}'"
				CMD="root -b -q rootlogon.C convertAmpToTree.C+\(\\\"${AMP_FILE}\\\",\\\"${ROOT_FILE}\\\"\)"
				echo ${CMD}
				time eval ${CMD}
				echo
				if [[ "${RECREATE_LINKS}" == "true" ]]
				then
						ln --symbolic --force --verbose ../$(basename ${ROOT_FILE}) $(dirname ${ROOT_FILE})/SYM
				fi
		done

done


echo ">>> info: ${0} finished on $(date)"
exit 0