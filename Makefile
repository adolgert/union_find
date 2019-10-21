# This just exists because fingers know to type "make".
.PHONY: all, clean

all:
	scons

clean:
	scons -c
	rm -rf .scons*

debug:
	gdb --args ./cluster_tbb_test --log_level=all --run_test=known_single_thread_tbb0
