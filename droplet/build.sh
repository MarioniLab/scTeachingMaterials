cat workshop_droplet.Rmd | grep -v "^#A#" | sed "s/^#H# //" > practical.Rmd
cat workshop_droplet.Rmd | grep -v "^#H#" | sed "s/^#A# //" > answers.Rmd
rm workshop_droplet.Rmd

echo "rmarkdown::render('prepare.Rmd')" | R --slave --no-save
rm prepare.Rmd

rm answers.md
rm test.sh
