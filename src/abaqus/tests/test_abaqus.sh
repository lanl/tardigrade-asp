# Make bash script more like high-level languages.
set -Eeuo pipefail

# Get this script's file name
script=`basename "$0"`

# Parse arguments
if [ "$#" -ne 3 ]; then
    echo "${script} USAGE:"
    echo "./${script} ABAQUS_PROGRAM USER_SUBROUTINE INPUT_FILE"
    echo "  REQUIRED POSITIONAL ARGUMENTS:"
    echo "    ABAQUS_PROGRAM: absolute path to abaqus executable"
    echo "    USER_SUBROUTINE: absolute path to user subroutine file"
    echo "    INPUT_FILE: input file base name without extension"
    exit 1
fi
abaqus_program=$1
input=$2
user=$3

# Verify abaqus path 
if [ ! "$(command -v ${abaqus_program})" ]; then
    echo "command ${abaqus_program} not found"
    exit 2
fi

# Verify expected input file path
if [ ! -f "../${input}.inp" ]; then
    echo "input file ../${input}.inp not found"
    exit 3
fi

# Verify user subroutine file
if [ ! -f "${user}" ]; then
    echo "user file ${user} not found"
    exit 4
fi

# Clean up build directory and copy input file
rm ${input}.{com,dat,log,msg,odb,prt,sim,sta} || true
cp ../${input}.inp .

# Run abaqus against the input file and redirect interactive output to log file
# Use interactive to avoid sleep and wait statements
${abaqus_program} -job ${input} -user ${user} -interactive >> ${input}.log 2>&1

# Test the output
grep "Hello Abaqus" ${input}.log 
grep "COMPLETED SUCCESSFULLY" ${input}.sta
