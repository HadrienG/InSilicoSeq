#!/bin/bash

# set image registry and image name
REGISTRY=""
IMAGE_NAME="insilicoseq"

# Version-build: sync version form setup.py
VERSION=$(python3 setup.py --version)

# create custom builder for multi-arch image. Default docker-builder can not create multi-arch builds.
docker buildx \
create \
--name inSilicoSeq-builder \
--use

# build + push multi-arch image for AMD and ARM architectures
docker build . \
--builder inSilicoSeq-builder \
--platform linux/amd64,linux/arm64 \
-t "${REGISTRY}/${IMAGE_NAME}:${VERSION}" \
--provenance=false \
--push

# remove custom inSilicoSeq-builder
docker buildx \
rm inSilicoSeq-builder