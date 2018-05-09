echo "rmarkdown::render('prepare.Rmd', clean=FALSE)" | R --slave --no-save
mv prepare.knit.md prepare.md

cat workshop_unbatch.Rmd | grep -v "^#H#" | sed "s/^#A# //" > answers.Rmd
echo "rmarkdown::render('answers.Rmd', clean=FALSE)" | R --slave --no-save
mv answers.knit.md answers.md
