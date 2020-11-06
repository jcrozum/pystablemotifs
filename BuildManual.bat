python -m pdoc --pdf StableMotifs > Manual.md
pandoc --metadata=title:"StableMotifs Documentation"  --pdf-engine=xelatex --variable=mainfont:"DejaVu Sans" --toc --toc-depth=4 --output=Manual.pdf  Manual.md
