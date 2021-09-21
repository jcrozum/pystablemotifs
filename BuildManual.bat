:: uses pdoc3, not pdoc, so must run `pip install pdoc3`
python -m pdoc --pdf pystablemotifs > Manual.md
pandoc --metadata=title:"pystablemotifs Documentation"  --pdf-engine=xelatex --variable=mainfont:"DejaVu Sans" --toc --toc-depth=4 --output=Manual.pdf  Manual.md
python -m pdoc --html pystablemotifs --force
