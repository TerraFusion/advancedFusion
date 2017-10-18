CC=gcc
H5CC=h5cc

#all: testRepro testRepro2 testRepro3 testReproHDF5
all: af_run test_aster test_aster_allOrbit test_read_area

#testRepro.o: testRepro.c
#	$(CC) -o $@ -c $<
#testRepro2.o: testRepro2.c
#	$(CC) -o $@ -c $<
#testRepro3.o: testRepro3.c
#	$(CC) -o $@ -c $<
#testReproHDF5.o: testReproHDF5.c
#	$(H5CC) -c $< -o $@ 
test_read_area.o: test_read_area.c
	$(H5CC) -c $< -o $@
af_run.o: af_run.c
	$(H5CC) -c $< -o $@
test_aster.o: test_aster.c
	$(H5CC) -c $< -o $@
test_aster_allOrbit.o: test_aster_allOrbit.c
	$(H5CC) -c $< -o $@
reproject.o: reproject.c
	$(CC) -o $@ -c $<
io.o: io.c
	$(H5CC) -c $< -o $@
# testRepro: testRepro.o reproject.o
#	$(CC) -o ../$@ $+ -lm
# testRepro2: testRepro2.o reproject.o
#	$(CC) -o ../$@ $+ -lm
# testRepro3: testRepro3.o reproject.o
#	$(CC) -o ../$@ $+ -lm
# testReproHDF5: testReproHDF5.o reproject.o io.o
#	$(H5CC) -o ../$@ $+ -lm
test_read_area: test_read_area.o reproject.o io.o
	$(H5CC) -o ../$@ $+ -lm
af_run: af_run.o reproject.o io.o
	$(H5CC) -o ../$@ $+ -lm
test_aster: test_aster.o reproject.o io.o
	$(H5CC) -o ../$@ $+ -lm
test_aster_allOrbit: test_aster_allOrbit.o reproject.o io.o
	$(H5CC) -o ../$@ $+ -lm
clean:
	rm *.o ../testRepro ../testRepro2 ../testRepro3 ../testReproHDF5 ../test_aster ../test_aster_allOrbit ../af_run
