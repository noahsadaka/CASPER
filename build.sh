#!/usr/bin/env bash
cmake -B build -S .
cmake --build build
cp build/spicetest.so .

