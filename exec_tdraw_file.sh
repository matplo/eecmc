#!/bin/bash

function usage() {
	echo "Usage: $0 <template_config.yaml> <file_input> [file_output]"
	exit 1
}

template_config=$1
if [ -z "$template_config" ]; then
		usage
		exit 1
fi

file_input=$2
if [ -z "$file_input" ]; then
		usage
		exit 1
fi

file_output=$3
if [ -z "$file_output" ]; then
		file_output=${file_input%.*}_h.root
fi

echo
echo "template_config: $template_config"
echo "file_input: $file_input"
echo "file_output: $file_output"
echo

# Run the script

function thisdir()
{
	SOURCE="${BASH_SOURCE[0]}"
	while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
	  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
	  SOURCE="$(readlink "$SOURCE")"
	  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
	done
	DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
	echo ${DIR}
}
THISD=$(thisdir)

cd ${THISD}

yaspenv_shell=$(which yaspenv.sh)
if [ -z ${yaspenv_shell} ]; then
  echo "Error: yaspenv.sh not found"
  exit 1
fi

SYS_YASP_DIR=$(${yaspenv_shell} yasp -q feature yasp_dir 2>&1 | tail -n 1)
if [ -z ${SYS_YASP_DIR} ]; then
  echo "Error: SYS_YASP_DIR not found"
  exit 1
fi

echo "[i] using ${SYS_YASP_DIR}"
source ${SYS_YASP_DIR}/venvyasp/bin/activate
module use ${SYS_YASP_DIR}/software/modules
module load yasp
module load bundle/hepbase

if [ ! -e ${THISD}/eecmc.module ]; then
  echo "Error: ${THISD}/eecmc.module not found"
  exit 1
fi
module use ${THISD}
module load eecmc.module
module list

#python tdraw_file.py $template_config $file_input $file_output
tmp_file=$(mktemp)
lre $template_config --define input=$file_input output=$file_output > $tmp_file

draw_from_yaml.py -c $tmp_file
