#!/usr/bin/bash

find ../plots -type f -regex ".*\.png" | parallel mogrify -trim {}
find ../plots -type f -regex ".*\.png" | parallel img2pdf -o {}.pdf {}
find ../plots -depth -name "*.png.pdf" -exec sh -c 'mv "$1" "${1%.png.pdf}.pdf"' _ {} \;
find ../plots -type f -regex ".*\.png" -exec rm {} \;
