echo "Creating project with cmake."

mkdir -p Release.Linux
cd Release.Linux

# Build types are Debug or Release
cmake ../src -Dall=ON -DCMAKE_BUILD_TYPE=Release
