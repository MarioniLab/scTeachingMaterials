cat workshop_kolod.Rmd | grep -v "^#A#" | sed "s/^#H# //" > practical.Rmd
cat workshop_kolod.Rmd | grep -v "^#H#" | sed "s/^#A# //" > answers.Rmd
rm workshop_kolod.Rmd

echo "rmarkdown::render('prepare.Rmd')" | R --slave --no-save
rm prepare.Rmd

rm answers.md
rm test.sh
