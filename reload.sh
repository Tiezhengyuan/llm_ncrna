#/bin/bash

lsmod | grep nvidia

sudo rmmod nvidia_uvm
sudo rmmod nvidia_drm
sudo rmmod nvidia_modeset
sudo rmmod nvidia

# Check processes using the kernel module
lsof /dev/nvidia*

modprobe nvidia_drm
modprobe nvidia_modeset
modprobe nvidia_uvm
modprobe nvidia
lsmod | grep nvidia

nvidia-smi
