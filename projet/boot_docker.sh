#!/bin/bash -
xhost +local:docker
sudo docker run --rm -it \
    -e DISPLAY="$DISPLAY" \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    my_app:latest
