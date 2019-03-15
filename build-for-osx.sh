echo "Creating project with cmake."
mkdir -p Release.OSX
cd Release.OSX
cmake ../src -Dall=ON -DCMAKE_BUILD_TYPE=Release
