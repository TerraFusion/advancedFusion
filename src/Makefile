#This Makefile template can be used to compile the AF tool on Linux
#To make this work on your own Linux, you need to install HDF5(either 1.8 or 1.10) and GDAL. 

#Replace /mnt/wrk/myang6/terra-fusion/af/hdf5-1.8.16 with your own HDF5 path
#Replace /mnt/wrk/myang6/terra-fusion/af/gdalbin with your own GDAL path

CXXFLAGS=-std=c++11 -Wno-write-strings -I/mnt/wrk/myang6/terra-fusion/af/hdf5-1.8.16/hdf5/include -L/mnt/wrk/myang6/terra-fusion/af/hdf5-1.8.16/hdf5/lib -lhdf5 -lhdf5_hl -lz

# Use this on BW
#CXX=CC ${CXXFLAGS}
#H5CXX=CC ${CXXFLAGS}

# Use this set for local build
CXX=g++ -g -std=c++11 -Wno-write-strings
GDALDIR = /mnt/wrk/myang6/terra-fusion/af/gdalbin
H5CXX=g++ -g $(CXXFLAGS)


all: AFtool test_aster test_aster_allOrbit test_read_area test_MISR_MODIS test_modis2aster test_Clipping_MISR_MODIS test_userdefinedgrids test_MISR_offset

AFtool.o: AFtool.cpp
	$(H5CXX) -c $< -o $@

test_read_area.o: test_read_area.cpp
	$(H5CXX) -c $< -o $@

test_aster.o: test_aster.cpp
	$(H5CXX) -c $< -o $@

test_MISR_MODIS.o: test_MISR_MODIS.cpp
	$(H5CXX) -c $< -o $@

test_modis2aster.o: test_modis2aster.cpp
	$(H5CXX) -c $< -o $@

test_Clipping_MISR_MODIS.o: test_Clipping_MISR_MODIS.cpp
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


AFtool: AFtool.o reproject.o io.o  misrutil.o gdalio.o AF_InputParmeterFile.o AF_debug.o AF_output_util.o AF_output_MODIS.o AF_output_MISR.o AF_output_ASTER.o
	$(H5CXX) -o ../$@ $+ -lm -L$(GDALDIR)/lib -lgdal -fopenmp

test_read_area: test_read_area.o reproject.o io.o
	$(H5CXX) -o ../$@ $+ -lm

test_aster: test_aster.o reproject.o io.o
	$(H5CXX) -o ../$@ $+ -lm

test_MISR_MODIS: test_MISR_MODIS.o reproject.o io.o
	$(H5CXX) -o ../$@ $+ -lm

test_modis2aster: test_modis2aster.o reproject.o io.o
	$(H5CXX) -o ../$@ $+ -lm

test_Clipping_MISR_MODIS: test_Clipping_MISR_MODIS.o reproject.o io.o
	$(H5CXX) -o ../$@ $+ -lm

test_aster_allOrbit: test_aster_allOrbit.o reproject.o io.o
	$(H5CXX) -o ../$@ $+ -lm

test_userdefinedgrids: test_userdefinedgrids.o reproject.o io.o gdalio.o
	$(H5CXX) -o ../$@ $+ -lm -L$(GDALDIR)/lib -lgdal -fopenmp

test_MISR_offset: test_MISR_offset.o io.o misrutil.o gdalio.o
	$(H5CXX) -o ../$@ $+ -lm -L$(GDALDIR)/lib -lgdal -fopenmp

clean:
	rm *.o ../AFtool ../test_read_area ../test_aster ../test_aster_allOrbit ../test_MISR_MODIS ../test_modis2aster ../test_Clipping_MISR_MODIS ../test_userdefinedgrids ../test_MISR_offset
#	rm *.o ../testRepro ../testRepro2 ../testRepro3 ../testReproHDF5
