#!/bin/bash
[[ ! -z "$(ls ?filter/*.fits 2>/dev/null)" ]] && mv ?filter/*.fits ./ && rm -rf ?filter/
rm -f *_after.fits
