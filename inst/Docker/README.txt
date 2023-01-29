nvidpt_bioc.dfile
-----------------
The nvidpt_bioc.dfile is a Dockerfile that built vjcitn/nvidpt_bioc:0.0.1 in dockerhub

It is capable of running scviR immediately.

Use

docker run --gpus all --ipc=host --ulimit memlock=-1 --ulimit stack=67108864 -ti vjcitn/nvidpt_bioc:0.0.1 bash

on a CUDA-capable machine
