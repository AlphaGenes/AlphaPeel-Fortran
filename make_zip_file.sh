mkdir AlphaPeel

# Assumes that the program and manual have both been built.

# To build the program run:
#cmake . ; make

# to build the manual:
#( cd docs ; make latexpdf)

cp example_multi_locus_spec.txt AlphaPeel
cp example_single_locus_spec.txt AlphaPeel

cp -r example AlphaPeel

# Copy in the documentation.
cp docs/build/latex/AlphaPeel.pdf AlphaPeel

if [ $? != 0 ]; then                   # last command: echo
    echo "The manual needs to be built." # last command: [
    exit 1
fi

# Copy in the binaries
cp binaries/* AlphaPeel

# Create a version file

version=`git describe --tags --abbrev=0`
commit=`git rev-parse --short HEAD`

echo Version: $version > AlphaPeel/version.txt
echo Commit: $commit >> AlphaPeel/version.txt

zip -r AlphaPeel.zip AlphaPeel
