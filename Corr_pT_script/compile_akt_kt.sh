g++ *.cc \
 -O2 -ansi -W -Wall -Wshadow -m64 -Wno-shadow \
 -o ${1}.exe \
 -I$PYTHIA8/include -L$PYTHIA8/lib/ -L$PYTHIA8/lib/archive -lpythia8 \
 `fastjet-config --cxxflags --libs` \
 `root-config --cflags --ldflags --glibs `

