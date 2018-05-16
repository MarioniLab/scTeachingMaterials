cat workshop_unbatch.Rmd | grep -v "^#A#" | sed "s/^#H# //" > practical.Rmd
cat workshop_unbatch.Rmd | grep -v "^#H#" | sed "s/^#A# //" > answers.Rmd
rm workshop_unbatch.Rmd

rm answers.md
rm test.sh
