include ../../config.mk

iflags=-I../../tgmg -I../../shared
lflags=-L../../shared -L.
objs=common.o fluid_2d.o fluid_2d_sio.o mgs_common.o mgs_mac.o mgs_fem.o \
	 bi_interp.o visco_impl.o tri_solve.o gas_layer.o gas_layer_data.o \
	 gas_layer_solve.o fileinfo.o
src=$(patsubst %.o,%.cc,$(objs))
execs=mg_test fluid_test edges tri_test h_analysis tt_analysis p_strip

all:
	$(MAKE) -C ../../shared
	$(MAKE) -C ../../tgmg lib
	$(MAKE) executables
	
executables: $(execs)

depend: $(src)
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

-include Makefile.dep

libf2d.a: $(objs)
	rm -f libf2d.a
	ar rs libf2d.a $^

fluid_test: fluid_test.cc libf2d.a
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lf2d $(lp_lflags)

mg_test: mg_test.cc mgs_common.o mgs_mac.o mgs_fem.o common.o
	$(cxx) $(cflags) $(iflags) $(lflags) -o $@ $^

tri_test: tri_test.cc tri_solve.o
	$(cxx) $(cflags) $(iflags) $(lflags) -o $@ $^ $(lp_lflags)

edges: edges.cc common.o
	$(cxx) $(cflags) $(iflags) $(lflags) -o $@ $^ -lgpmtx $(png_lflags)

h_analysis: h_analysis.cc common.o
	$(cxx) $(cflags) $(iflags) $(lflags) -o $@ $^ -lgpmtx $(png_lflags)

tt_analysis: tt_analysis.cc libf2d.a
	$(cxx) $(cflags) $(iflags) $(lflags) -o $@ $^ -lgpmtx -lf2d $(png_lflags) $(lp_lflags)

p_strip: p_strip.cc common.o
	$(cxx) $(cflags) $(iflags) $(lflags) -o $@ $^ -lgpmtx $(png_lflags)

data:
	rsync -rvz olympus.seas.harvard.edu:/data/files/splash/\* .

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

clean:
	rm -f $(execs) $(objs) libf2d.a

.PHONY: clean data all executables depend
