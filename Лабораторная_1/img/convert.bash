
for f in device.pdf; do
    convert -density 600 $f -alpha remove png/$f.png
done
