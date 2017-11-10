CXXFLAGS=-Wno-write-strings
CXX=g++
H5CXX=HDF5_CC=/usr/bin/g++ h5cc $(CXXFLAGS)

#all: testRepro testRepro2 testRepro3 testReproHDF5
all: af_run test_aster test_aster_allOrbit test_read_area

#testRepro.o: testRepro.cpp
#	$(CXX) -o $@ -c $<

#testRepro2.o: testRepro2.cpp
#	$(CXX) -o $@ -c $<

#testRepro3.o: testRepro3.cpp
#	$(CXX) -o $@ -c $<

#testReproHDF5.o: testReproHDF5.cpp
#	$(H5CXX) -c $< -o $@ 

test_read_area.o: test_read_area.cpp
	$(H5CXX) -c $< -o $@

af_run.o: af_run.cpp
	$(H5CXX) -c $< -o $@

test_aster.o: test_aster.cpp
	$(H5CXX) -c $< -o $@

test_aster_allOrbit.o: test_aster_allOrbit.cpp
	$(H5CXX) -c $< -o $@

reproject.o: reproject.cpp
	$(CXX) -o $@ -c $<

io.o: io.cpp
	$(H5CXX) -c $< -o $@

# testRepro: testRepro.o reproject.o
#	$(CXX) -o ../$@ $+ -lm

# testRepro2: testRepro2.o reproject.o
#	$(CXX) -o ../$@ $+ -lm

# testRepro3: testRepro3.o reproject.o
#	$(CXX) -o ../$@ $+ -lm

# testReproHDF5: testReproHDF5.o reproject.o io.o
#	$(H5CXX) -o ../$@ $+ -lm

test_read_area: test_read_area.o reproject.o io.o
	$(H5CXX) -o ../$@ $+ -lm

af_run: af_run.o reproject.o io.o
	$(H5CXX) -o ../$@ $+ -lm

test_aster: test_aster.o reproject.o io.o
	$(H5CXX) -o ../$@ $+ -lm

test_aster_allOrbit: test_aster_allOrbit.o reproject.o io.o
	$(H5CXX) -o ../$@ $+ -lm

clean:
	rm *.o ../test_read_area ../test_aster ../test_aster_allOrbit ../af_run
#	rm *.o ../testRepro ../testRepro2 ../testRepro3 ../testReproHDF5
