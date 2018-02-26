
CXXFLAGS=-std=c++11 -Wno-write-strings -fopenmp

# Use this on BW
CXX=CC ${CXXFLAGS}
H5CXX=CC ${CXXFLAGS}
GDALDIR = /projects/sciteam/jq0/TerraFusion/gdal

# Use this set for local build
#CXX=g++
#H5CXX=HDF5_CC=/usr/bin/g++ h5cc $(CXXFLAGS)

#all: testRepro testRepro2 testRepro3 testReproHDF5
all: AFtool test_aster test_aster_allOrbit test_read_area test_MISR_MODIS test_userdefinedgrids test_MISR_offset

AFtool.o: AFtool.cpp
	$(H5CXX) -c $< -o $@

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

#af_run.o: af_run.cpp
#	$(H5CXX) -c $< -o $@

test_aster.o: test_aster.cpp
	$(H5CXX) -c $< -o $@

test_MISR_MODIS.o: test_MISR_MODIS.cpp
	$(H5CXX) -c $< -o $@

test_aster_allOrbit.o: test_aster_allOrbit.cpp
	$(H5CXX) -c $< -o $@

test_userdefinedgrids.o: test_userdefinedgrids.cpp
	$(H5CXX) -c $< -o $@

test_MISR_offset.o: test_MISR_offset.cpp
	$(H5CXX) -c $< -o $@

reproject.o: reproject.cpp
	$(CXX) -o $@ -c $<

io.o: io.cpp
	$(H5CXX) -c $< -o $@

gdalio.o: gdalio.cpp
	$(CXX) -I. -I$(GDALDIR)/include -o $@ -c $< 	

misrutil.o: misrutil.cpp
	$(CXX) -c $< -o $@

# testRepro: testRepro.o reproject.o
#	$(CXX) -o ../$@ $+ -lm

# testRepro2: testRepro2.o reproject.o
#	$(CXX) -o ../$@ $+ -lm

# testRepro3: testRepro3.o reproject.o
#	$(CXX) -o ../$@ $+ -lm

# testReproHDF5: testReproHDF5.o reproject.o io.o
#	$(H5CXX) -o ../$@ $+ -lm

AFtool: AFtool.o AF_InputParmeterFile.o AF_debug.o reproject.o io.o
	$(H5CXX) -o ../$@ $+ -lm

test_read_area: test_read_area.o reproject.o io.o
	$(H5CXX) -o ../$@ $+ -lm

#af_run: af_run.o reproject.o io.o
#	$(H5CXX) -o ../$@ $+ -lm

test_aster: test_aster.o reproject.o io.o
	$(H5CXX) -o ../$@ $+ -lm

test_MISR_MODIS: test_MISR_MODIS.o reproject.o io.o
	$(H5CXX) -o ../$@ $+ -lm

test_aster_allOrbit: test_aster_allOrbit.o reproject.o io.o
	$(H5CXX) -o ../$@ $+ -lm

test_userdefinedgrids: test_userdefinedgrids.o reproject.o io.o gdalio.o
	$(H5CXX) -o ../$@ $+ -lm -L$(GDALDIR)/lib -lgdal

test_MISR_offset: test_MISR_offset.o io.o misrutil.o gdalio.o
	$(H5CXX) -o ../$@ $+ -lm -L$(GDALDIR)/lib -lgdal

clean:
	rm *.o ../AFtool ../test_read_area ../test_aster ../test_aster_allOrbit ../test_MISR_MODIS ../test_userdefinedgrids ../test_MISR_offset
#	rm *.o ../testRepro ../testRepro2 ../testRepro3 ../testReproHDF5
