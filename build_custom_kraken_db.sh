
count=0

while read p;

do 

((count=count+1))
echo $count

sed -i "/>/ c\>REFSEQ_$i|kraken:taxid|1$i" $p

done<refseq_fna

count=0

while read p;

do 

((count=count+1))
echo $count

sed -i "/>/ c\>MOUSE_$i|kraken:taxid|1$i" $p

done<mouse_locs

# add to library - set this directory
dbdir=~/almeida_genomes/kraken2_db
find -name 'UMGS*.fa' -print0 | xargs -0 -I{} -P 8 -n1 kraken2-build --add-to-library {} --db $dbdir

mkdir ~/almeida_genomes/kraken2_db/taxonomy/
# construct nodes.dmp and names.dmp taxonomy files
echo -e "1\t|\t1\t|\tno rank\t|\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|" > ~/almeida_genomes/kraken2_db/taxonomy/nodes.dmp
while read line; do
    UMGSID=1"$line"
    echo -e "$UMGSID\t|\t1\t|\tspecies\t|\t|\t0\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|" >> ~/almeida_genomes/kraken2_db/taxonomy/nodes.dmp
done < ~/almeida_genomes/umgsids.txt

echo -e "1\t|\tall\t|\t\t|\tsynonym\t|" > ~/almeida_genomes/kraken2_db/taxonomy/names.dmp
echo -e "1\t|\troot\t|\t\t|\tscientific name\t|" >> ~/almeida_genomes/kraken2_db/taxonomy/names.dmp
while read line; do
    taxid=1"$line"
    UMGSID="$line"
    echo -e "$taxid\t|\tAlmeida_UMGSID_$UMGSID\t|\t\t|\tAlmeida name\t|" >> ~/almeida_genomes/kraken2_db/taxonomy/names.dmp
done < ~/almeida_genomes/umgsids.txt

# then build in parallel
kraken2-build --build --threads 16 --db ~/almeida_genomes/kraken2_db/
     