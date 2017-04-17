# build Xcode environment
echo "creating Xcode project."
mkdir -p xcode
cd xcode
touch rttmbtlib_dummy.cpp rttwcet_dummy.cpp
echo "configure build using 'cmake -DGENERATOR=Xcode ..."
cmake ../src -DGENERATOR=Xcode -Dall=ON -DCMAKE_BUILD_TYPE=Debug -G Xcode
