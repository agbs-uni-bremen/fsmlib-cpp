echo "Creating project with cmake."
mkdir -p Debug.OSX
cd Debug.OSX
cmake ../src -Dall=ON -DCMAKE_BUILD_TYPE=Debug
