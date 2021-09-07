#!/bin/sh
python -m pdoc --pdf pystablemotifs > Manual.md
pandoc --metadata=title:"pystablemotifs Documentation"  --pdf-engine=xelatex --variable=mainfont:"DejaVu Sans" --toc --toc-depth=4 --output=Manual.pdf  Manual.md
