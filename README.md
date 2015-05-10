## bench
a collection of a few simple fpu-intensive benchmarks

http://kluge.in-chemnitz.de/docs/notes/benchmark.php

### fpu
This benchmark measures the speed of individual arithmetic
operations such as add, mul, div, sqrt, sin, etc. To compile:

```
$ cd fpu
$ make
```

This will also generate assembly output for inspection.
To run the programme the cpu clock speed needs to be provided as well
as the number of operations:

```
user@pentium_mmx:~ $ ./fpu 200 100	# cpu clock is 200MHz, do 100 mill ops
 ...
 calculating: x<k> + a + a + ... + a - (x<k>+n*a), with a=1.521700
 addx87: 1.507s, res=-1.206e-04,  66.349 Mflops, 0.332 flops/cycle
 add1:   1.632s, res=-1.207e-04,  61.268 Mflops, 0.306 flops/cycle
 add2:   0.848s, res= 7.831e-06,   0.118 Gflops, 0.590 flops/cycle
 add<1>: 1.570s, res=-1.207e-04,  63.709 Mflops, 0.319 flops/cycle
 add<2>: 0.847s, res= 7.816e-06,   0.118 Gflops, 0.590 flops/cycle
 add<3>: 0.586s, res= 1.878e-05,   0.171 Gflops, 0.853 flops/cycle
 add<4>: 0.596s, res=-1.319e-05,   0.168 Gflops, 0.839 flops/cycle
 add<5>: 0.552s, res=-2.651e-05,   0.181 Gflops, 0.905 flops/cycle
 ...

 user@nehalem:~ $ ./fpu 2667 100	# cpu clock is 2.67GHz, do 100 mill ops
 ...
 calculating: x<k> + a + a + ... + a - (x<k>+n*a), with a=1.521700
 addx87: 0.118s, res=-1.206e-04,   0.845 Gflops, 0.317 flops/cycle
 add1:   0.115s, res=-2.540e-01,   0.871 Gflops, 0.327 flops/cycle
 add2:   0.057s, res= 1.175e-02,   1.743 Gflops, 0.654 flops/cycle
 add<1>: 0.114s, res=-2.540e-01,   0.874 Gflops, 0.328 flops/cycle
 add<2>: 0.057s, res= 1.175e-02,   1.748 Gflops, 0.656 flops/cycle
 add<3>: 0.038s, res= 3.682e-02,   2.622 Gflops, 0.983 flops/cycle
 add<4>: 0.038s, res=-2.601e-02,   2.623 Gflops, 0.984 flops/cycle
 add<5>: 0.038s, res=-5.064e-02,   2.627 Gflops, 0.985 flops/cycle
```

- `add<2>` means two independent additions are executed consecutively,
	to make use of pipelining of instructions
- `res=...` can be ignored, it is only used so the compiler
	does not optimise loops away.
- `... Gflops` is the number of floating point operations per second
- `... flops/cycle` is the number of floating point operations per CPU	
	cycle


For a list of timing of assembly instructions see:
- http://kluge.in-chemnitz.de/docs/notes/asm_timing.php
- http://agner.org/optimize/instruction_tables.pdf


### linalg

This benchmark solves a dense linear equation system, similar
to the [Linpack Benchmark](http://www.top500.org/project/linpack)
which is used to compile the
[Top 500 Supercomputer list](http://www.top500.org/).

This benchmark requires optimised versions of the BLAS routines,
like [ATLAS](http://math-atlas.sourceforge.net/),
[GotoBLAS2] (https://www.tacc.utexas.edu/research-development/tacc-software/gotoblas2)
and [Eigen](http://eigen.tuxfamily.org/).


To compile:
```
$ cd linalg
$ make
```
