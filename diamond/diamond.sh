diamond makedb --in $1 --db $1
diamond blastx --db $1 --out $3 --outfmt 8 --query $2
