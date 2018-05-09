rm LICENSE
rm README.md
rm -rf lectures
rm -rf .git

for x in basic droplet unbatch
do
    cd $x
    bash build.sh
    rm build.sh
    cd -
done

rm -- $0
