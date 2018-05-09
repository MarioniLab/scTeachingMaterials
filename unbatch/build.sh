cat workshop_unbatch.Rmd | grep -v "^#A#" | sed "s/^#H# //" > practical.Rmd
cat workshop_unbatch.Rmd | grep -v "^#H#" | sed "s/^#A# //" > answers.Rmd
rm workshop_unbatch.Rmd

echo "rmarkdown::render('prepare.Rmd')" | R --slave --no-save
rm prepare.Rmd
