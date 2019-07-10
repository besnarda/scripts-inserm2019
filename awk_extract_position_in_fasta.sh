# this awk function need a double input ids.txt with name of sequence start end and a fasta file where to extract it
awk '
function char(){
                 n=split(X[p],Z)
                 for(i=1;i<=n;i+=2){print p":"Z[i]"-"Z[i+1] RS substr(s,Z[i],Z[i+1]-Z[i]+1)}
               }

        FNR==NR{
                 X[">"$1] = X[">"$1] ? X[">"$1] FS $2 FS $3 : $2 FS $3
                 next
               }

           !/>/{s = s $0}

            />/{if(s){ char();s = ""}p=$1}

            END{char()}
     ' 
