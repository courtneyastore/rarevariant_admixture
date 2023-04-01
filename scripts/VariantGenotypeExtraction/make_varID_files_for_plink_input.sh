#!/usr/bin/env bash

set -e -u -o pipefail

bindir="$(dirname "$0")"


export var_file

get_input() {

    while getopts f:h option
    do
        case "${option}"
        in
			f) var_file=${OPTARG};;
            h) echo "Use -f flag to input your file";;
            *) echo "test";;
        esac
    done
}

create_varID_file(){

	while IFS=$'\t' read -r v a2;
	do
		out_file="${v}.txt"

		echo -e "${v}\t${a2}" > "${out_file}"
		
	done < "$var_file"
	
}

main () {
    get_input "$@" || exit 2
    create_varID_file
    echo "All done :)"
}

main "$@"
