#!/bin/sh
python -m pdoc --pdf PyStableMotifs > Manual.md
pandoc --metadata=title:"PyStableMotifs Documentation"  --pdf-engine=xelatex --variable=mainfont:"DejaVu Sans" --toc --toc-depth=4 --output=Manual.pdf  Manual.md
