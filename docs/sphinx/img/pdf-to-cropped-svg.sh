#/usr/bin

# References:
#  - Inkscape: https://inkscape.org/sk/doc/inkscape-man.html
#  - pdfcrop: http://alexsleat.co.uk/2011/01/25/using-pdfcrop-to-remove-white-margins-ubuntu/

for f in "$@"; do pdfcrop $f $f; inkscape $f --export-text-to-path --export-plain-svg=${f%.pdf}.svg; done
